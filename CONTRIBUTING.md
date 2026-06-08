# Contributing

This repository is currently a research prototype.

Please keep changes narrow and reproducible:

1. Run `./build.sh` and `python3 tests/test_smoke.py`.
2. If changing a fast kernel, compare against exact heap/DP generation for small ranks.
3. If changing an audited path, include the audit output.
4. Do not claim broad state-of-the-art status without adding an appropriate comparator.
5. Keep dependencies minimal.
