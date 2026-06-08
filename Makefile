.PHONY: build test clean audit5 audit6 bench-dp

build:
	./build.sh

test: build
	python3 tests/test_smoke.py

audit5: build
	python3 smooth_toolkit.py audit 2,3,5,7,11 1000000000000 --timeout 300

audit6: build
	python3 smooth_toolkit.py audit 2,3,5,7,11,13 1000000000000 --timeout 300

bench-dp: build
	bash scripts/run_dp_comparison.sh

clean:
	rm -rf bin __pycache__ tests/__pycache__ results/local
