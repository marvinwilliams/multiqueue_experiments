#!/bin/bash

set -euo pipefail

function usage {
  echo "Usage:  $0 [<options>] <command> [<args>]

  Options:
    -h|help         Display this message

  Commands:
    create          Create a new experiment from definition
    run             Run an experiment or steps from it
    fetch           Aggregate experiment results" >&2
}

function create_usage {
  echo "Usage:  $0 create [<options>] [<definition>...]

  <definition>      Path to the experiment definition in JSON format

  Options:
    -h|help         Display this message
    -d|directory    Directory to create the experiments in
    --              Signal end of options" >&2
}

function run_usage {
  echo "Usage:  $0 run [<options>] <experiment> [<step>...]

  <experiment>      Path to the experiment
  <step>            Steps to run. Possible steps are
                      p[repare]: Prepare the experiment
                      e[xecute]: Execute the runs of this experiment
                      c[leanup]: Perform experiment cleanup
                    Defaults to 'p e c'

  Options:
    -h|help         Display this message
    --              Signal end of options" >&2
}

function fetch_usage {
  echo "Usage:  $0 fetch [<options>] [<experiment>...]

  <experiment>      Path to the experiment to fetch the results from

  Options:
    -h|help         Display this message
    -f|--filter     Filter passed to jq applied to each fetched run. Defaults to '.'
    --              Signal end of options" >&2
}

exp_base_dir=$(pwd)
positional_args=()
filter='.'
outfile="fetched_results.json"
run_dir_format="%04d"

if command -v hostname; then
  hostname=$(hostname -s)
elif [[ -r "/proc/sys/kernel/hostname" ]]; then
  hostname=$(cat /proc/sys/kernel/hostname)
else
  echo "Warning: Could not determine hostname, defaulting to 'unknown'" >&2
  hostname="unknown"
fi

function create_experiment {
  local exp_file exp_json exp_name exp_dir runs_json num_runs runs_dir i run_json run_dir resource_json res_name res_path
  for exp_file in "${positional_args[@]}"; do
    if [[ ! -f ${exp_file} ]]; then
      echo "Warning: Experiment file '${exp_file}' not found, skipping" >&2
      continue
    fi

    exp_json=$(jq -c '.' ${exp_file})

    exp_name=$(jq -r '.properties.name // ""' <<<${exp_json})

    if [[ -z $exp_name ]]; then
      echo "Error: Key 'name' nof found in 'properties' or \"name\" is empty" >&2
      return 1
    fi

    exp_dir="${exp_base_dir}/${exp_name}"

    if [[ -d ${exp_dir} ]]; then
      echo "Warning: Experiment directory '${exp_dir}' already exists, skipping" >&2
      return 1
    fi

    echo ":: Creating experiment '${exp_name}' in '${exp_dir}'"

    mkdir -p "${exp_dir}"

    runs_json="$(jq -c '.runs // []' <<<${exp_json})"
    num_runs=$(jq 'length' <<<${runs_json})
    echo "-> Writing static properties"
    jq --argjson runs $num_runs '{properties: .properties, prepare: (.prepare // []), cleanup: (.cleanup // []), runs: $runs}' <<<$exp_json > "${exp_dir}/static.json"

    echo "-> Adding resources"
    for resource_json in $(jq -c '.resources // [] | .[]' <<<${exp_json}); do
      res_name=$(jq -r '.name' <<<${resource_json})
      res_path=$(jq -r '.path' <<<${resource_json})
      if [[ $(jq '(.symlink // false) == true' <<<${resource_json}) == "true" ]]; then
        ln -s -r ${res_path} ${exp_dir}/${res_name}
      else
        cp ${res_path} "${exp_dir}/${res_name}"
      fi
    done

    runs_dir="${exp_dir}/runs"
    mkdir ${runs_dir}

    i=0
    while IFS='' read -r run_json; do
      run_dir=${runs_dir}/$(printf "${run_dir_format}" $i)
      i=$((i+1))
      echo -ne "-> Create run directories ($i/$num_runs) \r"
      mkdir ${run_dir}
      jq '{properties: (.properties // {}), commands: (.commands // [])}' <<<$run_json > "${run_dir}/static.json"

      while IFS='' read -r resource_json; do
        res_name=$(jq -r '.name' <<<${resource_json})
        res_path=$(jq -r '.path' <<<${resource_json})
        if jq '.symlink' > /dev/null <<<${resource_json}; then
          ln -s -r ${res_path} ${run_dir}/${res_name}
        else
          cp ${res_path} "${run_dir}/${res_name}"
        fi
      done < <(jq -c '.resources // [] | .[]' <<<${run_json})
    done < <(jq -c '.[]' <<<${runs_json})
    echo ""
  done
  echo ":: Done"
}

function prepare_experiment {
  local exp_dir="$1"
  local cmd exit_status

  echo ":: Preparing experiment"

  if [[ ! -f "${exp_dir}/static.json" ]]; then
    echo "Error: 'static.json' not found in '${exp_dir}'" >&2
    return 1
  fi

  while IFS='' read -r cmd_json; do
    cmd=()
    while IFS='' read -r arg; do
      cmd+=("$arg")
    done < <(jq -r '.[]' <<<${cmd_json})
    echo ""
    echo "-> Executing '${cmd[@]}'"
    set +e
    (cd ${exp_dir}; "${cmd[@]}")
    exit_status=$?
    set -e
    if [[ ${exit_status} -ne 0 ]]; then
      echo "Error: Preparation failed" >&2
      return 1
    fi
  done < <(jq -c '.prepare[]' < "${exp_dir}/static.json")
}

function execute_runs {
  local exp_dir="$1"

  if [[ ! -f "${exp_dir}/static.json" ]]; then
    echo "Error: 'static.json' not found in '${exp_dir}'" >&2
    return 1
  fi

  local num_runs=$(jq '.runs' < "${exp_dir}/static.json")
  local runs_dir="${exp_dir}/runs"
  local i run_dir exit_status cmd

  for ((i=0;i<num_runs;i++)); do
    run_dir=${runs_dir}/$(printf "${run_dir_format}" $i)
    if [[ ! -d "${run_dir}" ]]; then
      echo "Warning: Run directory '${run_dir}' not found, skipping" >&2
      continue
    fi
    if [[ ! -f "${run_dir}/static.json" ]]; then
      echo "Warning: Run directory '${run_dir}' contains no 'static.json', skipping" >&2
      continue
    fi
    num_cmds=$(jq -r '.commands | length' < "${run_dir}/static.json")
    echo ":: Start run $((i+1))/${num_runs} (${num_cmds} commands)"
    exit_status="null"
    while IFS='' read -r cmd_json; do
      cmd=()
      while IFS='' read -r arg; do
        cmd+=("$arg")
      done < <(jq -r '.[]' <<<${cmd_json})
      echo "-> Execute '${cmd[@]}'"
      set +e
      (cd "${run_dir}"; "${cmd[@]}" >> "stdout.log" 2>> "stderr.log")
      exit_status=$?
      set -e
      if [[ ${exit_status} -ne 0 ]]; then
        echo "Warning: Command failed, skipping remaining commands" >&2
        break
      fi
    done < <(jq -c '.commands[]' < "${run_dir}/static.json")
    i=$((i+1))
    jq -n --arg hostname ${hostname} --argjson exit_status ${exit_status} '{host: $hostname, exit_status: $exit_status}' > "${run_dir}/result.json"
    echo ""
  done
  echo ":: Done"
}

function cleanup_experiment {
  local exp_dir="$1"

  if [[ ! -f "${exp_dir}/static.json" ]]; then
    echo "Error: 'static.json' not found in '${exp_dir}'" >&2
    return 1
  fi

  local cmd

  echo ":: Cleaning up experiment"

  while IFS='' read -r cmd_json; do
    cmd=()
    while IFS='' read -r arg; do
      cmd+=("$arg")
    done < <(jq -r '.[]' <<<${cmd_json})
    echo ""
    echo "-> Executing '${cmd[@]}'"
    set +e
    (cd ${exp_dir}; "${cmd[@]}")
    set -e
  done < <(jq -c '.prepare[]' < "${exp_dir}/static.json")
}

function run_experiment {
  if [[ "${#positional_args[@]}" -eq 0 ]]; then
    echo "Error: No experiment given" >&2
    return 129
  fi
  local exp_dir="${positional_args[@]:0}"
  if [[ ! -d "${exp_dir}" ]]; then
    echo "Error: '${exp_dir}' not found" >&2
    return 1
  fi
  local steps
  if [[ "${#positional_args[@]}" -eq 1 ]]; then
    steps=("p" "e" "c")
  else
    steps=("${positional_args[@]:1}")
  fi
  for step in "${steps[@]}"; do
    case ${step} in
      p|prepare)
        prepare_experiment "${exp_dir}"
        if [[ $? -ne 0 ]]; then
          return $?
        fi
        ;;
      e|execute)
        execute_runs "${exp_dir}"
        ;;
      c|cleanup)
        cleanup_experiment "${exp_dir}"
        ;;
      *)
        echo "Error: Unknown step '${step}'" >&2
        break
    esac
  done
}

function fetch_results {
  local num_exps i exp_dir num_runs runs_dir j run_dir
  num_exps="${#positional_args[@]}"
  exec 4>&1
  i=1
  for exp_dir in "${positional_args[@]}"; do
    echo -ne ":: Fetching results ($i/${num_exps}) \r" >&4
    if [[ ! -d "${exp_dir}" ]]; then
      echo "Warning: '${exp_dir}' not found, skipping" >&2
      continue
    fi
    if [[ ! -f "${exp_dir}/static.json" ]]; then
      echo "Warning: 'static.json' not found in '${exp_dir}', skipping" >&2
      continue
    fi

    num_runs=$(jq '.runs' < "${exp_dir}/static.json")
    runs_dir="${exp_dir}/runs"

    for ((j=0;j<num_runs;j++)); do
      run_dir=${runs_dir}/$(printf "${run_dir_format}" $j)
      if [[ ! -d "${run_dir}" ]]; then
        echo "Warning: Run directory '${run_dir}' not found, skipping" >&2
        continue
      fi
      if [[ ! -f "${run_dir}/result.json" ]]; then
        echo "Warning: Run directory '${run_dir}' contains no 'result.json', skipping" >&2
        continue
      fi
      jq -c "${filter}" "${run_dir}/result.json"
    done
    i=$((i+1))
  done | jq --slurp '.' > "${outfile}"
  echo "" >&4
  echo ""
  echo ":: Done"
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      exit 129
      ;;
    *)
      command_verb="$1"
      shift
      break
      ;;
  esac
  shift
done

if [[ -z ${command_verb} ]]; then
  usage
  exit 1
fi

case ${command_verb} in
  create)
    while [[ $# -gt 0 ]]; do
      case $1 in
        -d|--exp-dir)
          shift
          if [[ $# -eq 0 ]]; then
            echo "Missing argument for option: -d|--exp-dir" >&2
            create_usage
            exit 129
          fi
          exp_base_dir="$1"
          ;;
        --)
          shift
          break
          ;;
        -*)
          echo "Unknown option: $1" >&2
          create_usage
          exit 129
          ;;
        *)
          positional_args+=("$1")
          ;;
      esac
      shift
    done
    while [[ $# -gt 0 ]]; do
      positional_args+=("$1")
    done
    create_experiment
    ;;
  run)
    while [[ $# -gt 0 ]]; do
      case $1 in
        --)
          shift
          break
          ;;
        -*)
          echo "Unknown option: $1" >&2
          run_usage
          exit 129
          ;;
        *)
          positional_args+=("$1")
          ;;
      esac
      shift
    done
    while [[ $# -gt 0 ]]; do
      positional_args+=("$1")
    done
    run_experiment
    ;;
  fetch)
    while [[ $# -gt 0 ]]; do
      case $1 in
        -f|--filter)
          shift
          if [[ $# -eq 0 ]]; then
            echo "Error: Missing argument for option: -f|--filter" >&2
            fetch_usage
            exit 129
          fi
          filter="$1"
          ;;
        -o|--output)
          shift
          if [[ $# -eq 0 ]]; then
            echo "Error: Missing argument for option: -o|--output" >&2
            fetch_usage
            exit 129
          fi
          outfile="$1"
          ;;
        --)
          shift
          break
          ;;
        -*)
          echo "Error: Unknown option: $1" >&2
          fetch_usage
          exit 129
          ;;
        *)
          positional_args+=("$1")
          ;;
      esac
      shift
    done
    while [[ $# -gt 0 ]]; do
      positional_args+=("$1")
    done
    fetch_results
    ;;
  *)
    echo "Unknown command: ${command_verb}" >&2
    usage
    exit 1
    ;;
esac
