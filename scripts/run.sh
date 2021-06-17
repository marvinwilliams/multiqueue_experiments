#!/bin/bash

set -euo pipefail

cmake -B build -DCMAKE_BUILD_TYPE=Release -Dmultiqueue_BUILD_TESTS=OFF
cmake --build build --target all -j $(nproc)

mkdir -p experiments

./scripts/benchmark_all.sh
