#!/bin/bash

set -euo pipefail

function usage {
  echo "Usage:  $0 [<options>] <command> [<args>]

  Options:
    -h|--help         Display this message

  Commands:
    create            Create a new experiment from definition
    run               Run an experiment
    fetch             Aggregate experiment results" >&2
}

function create_usage {
  echo "Usage:  $0 create [<options>] <definition>

  <definition>        Path to the experiment definition in JSON format

  Options:
    -d|directory      Directory to create the experiment in
    -P|--no-prepare   Don't execute the prepare commands
    -R|--no-resources Don't add resources
    -h|--help         Display this message
    --                Signal end of options" >&2
}

function run_usage {
  echo "Usage:  $0 run [<options>] <experiment>

  <experiment>      Path to the experiment

  Options:
    -p|--prepare      (Re-) run the prepare commands
    -r|--resources    (Re-) add resources
    -E|--no-execute   Don't execute the runs itself
    -C|--no-cleanup   Don't perform cleanup
    -h|--help         Display this message
    --                Signal end of options" >&2
}

function fetch_usage {
  echo "Usage:  $0 fetch [<options>] [<experiment>...]

  <experiment>      Path to the experiment to fetch the results from

  Options:
    -f|--filter     Filter passed to jq applied to each fetched run. Defaults to '.'
    -h|--help       Display this message
    --              Signal end of options" >&2
}

run_dir_format="%04d"

exp_base_dir=$(pwd)
positional_args=()
filter='.'
outfile="fetch_results.json"
with_prepare="true"
with_resources="true"
with_execute="true"
with_cleanup="true"

if command -v hostname; then
  hostname=$(hostname -s)
elif [[ -r "/proc/sys/kernel/hostname" ]]; then
  hostname=$(cat /proc/sys/kernel/hostname)
else
  echo "Warning: Could not determine hostname, defaulting to 'unknown'" >&2
  hostname="unknown"
fi

function prepare_experiment {
  local exp_dir="$1"

  if [[ ! -d "${exp_dir}" ]]; then
    echo "Error: Experiment '${exp_dir}' not found" >&2
    return 1
  fi

  if [[ ! -f "${exp_dir}/static.json" ]]; then
    echo "Error: Experiment file '${exp_dir}/static.json' not found" >&2
    return 1
  fi

  echo ":: Preparing experiment..."

  local cmd exit_status
  while IFS='' read -r cmd_json; do
    cmd=()
    while IFS='' read -r arg; do
      cmd+=("$arg")
    done < <(jq -r '.[]' <<<${cmd_json})
    set +e
    (cd "${exp_dir}"; "${cmd[@]}")
    exit_status=$?
    set -e
    if [[ ${exit_status} -ne 0 ]]; then
      echo "Error: Preparation failed" >&2
      return 1
    fi
  done < <(jq -c '.prepare_commands[]?' < "${exp_dir}/static.json")
}

function add_resources {
  local exp_dir="$1"

  if [[ ! -d "${exp_dir}" ]]; then
    echo "Error: Experiment '${exp_dir}' not found" >&2
    return 1
  fi

  if [[ ! -f "${exp_dir}/static.json" ]]; then
    echo "Error: Experiment file '${exp_dir}/static.json' not found" >&2
    return 1
  fi

  echo ":: Adding resources..."

  local res_name res_path 
  while IFS='' read -r resource_json; do
    res_name=$(jq -r '.name' <<<${resource_json})
    res_path=$(jq -r '.path' <<<${resource_json})
    if jq -e '(.symlink // false)' <<<${resource_json} > /dev/null; then
      ln -s -r "${res_path}" ${exp_dir}/${res_name}
    else
      cp "${res_path}" "${exp_dir}/${res_name}"
    fi
  done < <(jq -c '.resources[]?' < "${exp_dir}/static.json")

  local runs_dir="${exp_dir}/runs"

  local num_runs=$(jq -r '.runs // 0' < "${exp_dir}/static.json")

  local i run_dir
  for ((i=0;i<num_runs;i++)); do
    run_dir="${runs_dir}/$(printf ${run_dir_format} $i)"
    if [[ ! -d "${run_dir}" ]]; then
      echo "Warning: Run directory '${run_dir}' not found, skipping" >&2
      continue
    fi
    if [[ ! -f "${run_dir}/static.json" ]]; then
      echo "Warning: Run directory '${run_dir}' contains no 'static.json', skipping" >&2
      continue
    fi
    while IFS='' read -r resource_json; do
      res_name=$(jq -r '.name' <<<${resource_json})
      res_path=$(jq -r '.path' <<<${resource_json})
      if jq -e '(.symlink // false)' <<<${resource_json} > /dev/null; then
        ln -s -r "${res_path}" "${run_dir}/${res_name}"
      else
        cp "${res_path}" "${run_dir}/${res_name}"
      fi
    done < <(jq -c '.resources[]?' < "${run_dir}/static.json")
  done
}

function create_experiment {
  if [[ ${#positional_args[@]} -eq 0 ]]; then
    echo "Error: Specify at least one experiment definition"
    return 129
  fi
  if [[ ${#positional_args[@]} -gt 1 ]]; then
    echo "Error: Excessive argument: ${positional_args[@]:1}"
    return 129
  fi
  local exp_file="${positional_args[@]:0}"
  if [[ ! -f ${exp_file} ]]; then
    echo "Error: Experiment file '${exp_file}' not found" >&2
    return 1
  fi

  local exp_json=$(jq -c '.' ${exp_file})

  local exp_name=$(jq -r '.properties.name // ""' <<<${exp_json})

  if [[ -z $exp_name ]]; then
    echo "Error: Key 'name' nof found in 'properties' or \"name\" is empty" >&2
    return 1
  fi

  local exp_dir="${exp_base_dir}/${exp_name}"

  if [[ -d ${exp_dir} ]]; then
    echo "Error: Experiment directory '${exp_dir}' already exists" >&2
    return 1
  fi

  echo ":: Creating experiment '${exp_name}' (${exp_dir})..."

  mkdir -p "${exp_dir}"

  local runs_json="$(jq -c '.runs // []' <<<${exp_json})"
  local num_runs=$(jq 'length' <<<${runs_json})
  jq --argjson runs $num_runs '{properties: .properties, prepare_commands: (.prepare_commands // []), resources: (.resources // []), cleanup_commands: (.cleanup_commands // []), runs: $runs}' <<<$exp_json > "${exp_dir}/static.json"

  local runs_dir="${exp_dir}/runs"
  mkdir "${runs_dir}"

  local i=0 run_dir
  while IFS='' read -r run_json; do
    run_dir=${runs_dir}/$(printf "${run_dir_format}" $i)
    mkdir ${run_dir}
    jq '{properties: (.properties // {}), resources: (.resources // []), commands: (.commands // [])}' <<<$run_json > "${run_dir}/static.json"
    i=$((i+1))
  done < <(jq -c '.[]' <<<${runs_json})

  if "${with_prepare}"; then
    prepare_experiment "${exp_dir}"
  fi
  if "${with_resources}"; then
    add_resources "${exp_dir}"
  fi
  echo ":: Done"
}

function execute_runs {
  local exp_dir="$1"

  if [[ ! -d "${exp_dir}" ]]; then
    echo "Error: Experiment '${exp_dir}' not found" >&2
    return 1
  fi

  if [[ ! -f "${exp_dir}/static.json" ]]; then
    echo "Error: 'static.json' not found in '${exp_dir}'" >&2
    return 1
  fi

  local num_runs=$(jq '.runs' < "${exp_dir}/static.json")
  local runs_dir="${exp_dir}/runs"
  local i run_dir num_cmds cmd_json exit_status j cmd 

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
    echo ":: Start run $((i+1))/${num_runs}..."
    exit_status="null"
    j=1
    while IFS='' read -r cmd_json; do
      cmd=()
      while IFS='' read -r arg; do
        cmd+=("$arg")
      done < <(jq -r '.[]' <<<${cmd_json})
      echo "($j/${num_cmds}) ${cmd[@]}"
      set +e
      (cd "${run_dir}"; "${cmd[@]}" >> "stdout.log" 2>> "stderr.log")
      exit_status=$?
      set -e
      if [[ ${exit_status} -ne 0 ]]; then
        echo "Warning: Command failed, skipping remaining commands" >&2
        break
      fi
      j=$((j+1))
    done < <(jq -c '.commands[]' < "${run_dir}/static.json")
    jq -n --arg hostname ${hostname} --argjson exit_status ${exit_status} '{host: $hostname, exit_status: $exit_status}' > "${run_dir}/result.json"
  done
}

function cleanup_experiment {
  local exp_dir="$1"

  if [[ ! -d "${exp_dir}" ]]; then
    echo "Error: Experiment '${exp_dir}' not found" >&2
    return 1
  fi

  if [[ ! -f "${exp_dir}/static.json" ]]; then
    echo "Error: 'static.json' not found in '${exp_dir}'" >&2
    return 1
  fi

  echo ":: Cleaning up experiment"

  local cmd_json cmd
  while IFS='' read -r cmd_json; do
    cmd=()
    while IFS='' read -r arg; do
      cmd+=("$arg")
    done < <(jq -r '.[]' <<<${cmd_json})
    set +e
    (cd "${exp_dir}"; "${cmd[@]}")
    set -e
  done < <(jq -c '.cleanup_commands[]?' < "${exp_dir}/static.json")
}

function run_experiment {
  if [[ "${#positional_args[@]}" -eq 0 ]]; then
    echo "Error: No experiment given" >&2
    return 129
  fi
  if [[ ${#positional_args[@]} -gt 1 ]]; then
    echo "Error: Excessive argument: ${positional_args[@]:1}"
    return 129
  fi

  local exp_dir="${positional_args[@]:0}"

  if "${with_prepare}"; then
    prepare_experiment "${exp_dir}"
  fi
  if "${with_resources}"; then
    add_resources "${exp_dir}"
  fi
  if "${with_execute}"; then
    execute_runs "${exp_dir}"
  fi
  if "${with_cleanup}"; then
    cleanup_experiment "${exp_dir}"
  fi
  echo ":: Done"
}

function fetch_results {
  local num_exps i exp_dir num_runs runs_dir exp_properties j run_dir
  num_exps="${#positional_args[@]}"

  echo ":: Fetching results"
  for exp_dir in "${positional_args[@]}"; do
    if [[ ! -d "${exp_dir}" ]]; then
      echo "Warning: '${exp_dir}' not found, skipping" >&2
      continue
    fi
    if [[ ! -f "${exp_dir}/static.json" ]]; then
      echo "Warning: 'static.json' not found in '${exp_dir}', skipping" >&2
      continue
    fi

    num_runs=$(jq '.runs' < "${exp_dir}/static.json")
    exp_properties=$(jq '.properties' < "${exp_dir}/static.json")
    runs_dir="${exp_dir}/runs"

    for ((j=0;j<num_runs;j++)); do
      run_dir=${runs_dir}/$(printf "${run_dir_format}" $j)
      if [[ ! -d "${run_dir}" ]]; then
        echo "Warning: Run directory '${run_dir}' not found, skipping" >&2
        continue
      fi
      if [[ ! -f "${run_dir}/static.json" ]]; then
        echo "Warning: Run directory '${run_dir}' contains no 'static.json', skipping" >&2
        continue
      fi
      if [[ ! -f "${run_dir}/result.json" ]]; then
        echo "Warning: Run directory '${run_dir}' contains no 'result.json', skipping" >&2
        continue
      fi
      jq -s -c --argjson exp_properties "${exp_properties}" "(\$exp_properties + .[0].properties + .[1]) | ${filter}" "${run_dir}/static.json" "${run_dir}/result.json"
    done
  done | jq --slurp > "${outfile}"
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
        -P|--no-prepare)
          with_prepare="false"
          ;;
        -R|--no-resources)
          with_resources="false"
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
    with_prepare="false"
    with_resources="false"
    while [[ $# -gt 0 ]]; do
      case $1 in
        --)
          shift
          break
          ;;
        -p|--prepare)
          with_prepare="true"
          ;;
        -r|--resources)
          with_resources="true"
          ;;
        -E|--no-execute)
          with_execute="false"
          ;;
        -C|--no-cleanup)
          with_cleanup="false"
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
