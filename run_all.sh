#!/usr/bin/env bash
# Full reproduction: install, test, run, verify.
set -euo pipefail
cd "$(dirname "$0")"
python3 -m pip install -q -r requirements.txt
python3 -m pytest tests/ -q
python3 run_all.py "$@"
python3 verify.py
