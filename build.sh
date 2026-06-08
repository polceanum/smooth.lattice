#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p bin
CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:-"-O3 -std=c++17"}

$CXX $CXXFLAGS src/smooth_3prime_beatty_ranker.cpp -o bin/smooth_3prime_beatty_ranker
$CXX $CXXFLAGS src/smooth_layer_compressed_general.cpp -o bin/smooth_layer_compressed_general
$CXX $CXXFLAGS src/smooth_sums_only_scalable.cpp -o bin/smooth_sums_only_scalable
$CXX $CXXFLAGS src/smooth_interval_audit_exps_k6.cpp -o bin/smooth_interval_audit_exps_k6

echo "Built core kernels in ./bin"
