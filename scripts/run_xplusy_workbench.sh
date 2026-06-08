#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p bin results/local
CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:-"-O3 -std=c++17"}
$CXX $CXXFLAGS benchmarks/smooth_xplusy_baseline.cpp -o bin/smooth_xplusy_baseline
$CXX $CXXFLAGS benchmarks/smooth_xplusy_fj_loh_workbench.cpp -o bin/smooth_xplusy_fj_loh_workbench
./bin/smooth_xplusy_baseline | tee results/local/xplusy_baseline.txt
./bin/smooth_xplusy_fj_loh_workbench | tee results/local/xplusy_fj_loh.txt
