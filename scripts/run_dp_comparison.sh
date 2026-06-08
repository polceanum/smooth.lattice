#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p bin results/local
CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:-"-O3 -std=c++17"}
./build.sh
$CXX $CXXFLAGS benchmarks/dp_pointer_baseline.cpp -o bin/dp_pointer_baseline
{
  echo "# DP pointer baseline versus smooth.lattice kernels"
  echo "# Generated $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  for N in 1000000 10000000; do
    echo "## P=2,3,5,7,11 N=$N"
    ./bin/dp_pointer_baseline 2,3,5,7,11 "$N"
    python3 smooth_toolkit.py unrank 2,3,5,7,11 "$N"
  done
} | tee results/local/dp_comparison.txt
