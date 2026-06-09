#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p bin
CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:-"-O3 -std=c++17"}

if ! printf '#include <boost/multiprecision/cpp_int.hpp>\nint main(){return 0;}\n' | $CXX $CXXFLAGS -x c++ - -o bin/.boost_multiprecision_check >/dev/null 2>&1; then
  cat >&2 <<'EOF'
ERROR: Boost.Multiprecision headers were not found by the active C++ compiler.

This project uses the header-only Boost.Multiprecision cpp_int type.
Install/provide Boost headers, then rerun ./build.sh. Examples:

  Debian/Ubuntu: sudo apt-get install libboost-dev
  macOS Homebrew: brew install boost

You may also pass an explicit include path, for example:

  CXXFLAGS="-O3 -std=c++17 -I/path/to/boost/include" ./build.sh
EOF
  exit 1
fi
rm -f bin/.boost_multiprecision_check

$CXX $CXXFLAGS src/smooth_3prime_beatty_ranker.cpp -o bin/smooth_3prime_beatty_ranker
$CXX $CXXFLAGS src/smooth_layer_compressed_general.cpp -o bin/smooth_layer_compressed_general
$CXX $CXXFLAGS src/smooth_sums_only_scalable.cpp -o bin/smooth_sums_only_scalable
$CXX $CXXFLAGS src/smooth_interval_audit_exps_k6.cpp -o bin/smooth_interval_audit_exps
$CXX $CXXFLAGS src/smooth_interval_audit_exps_k6.cpp -o bin/smooth_interval_audit_exps_k6

echo "Built core kernels in ./bin"
