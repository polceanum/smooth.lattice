#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/local
./build.sh
python3 smooth_toolkit.py audit 2,3,5,7,11 1000000000000 --timeout 300 | tee results/local/audit_5prime_1e12.json
python3 smooth_toolkit.py audit 2,3,5,7,11,13 1000000000000 --timeout 300 | tee results/local/audit_6prime_1e12.json
