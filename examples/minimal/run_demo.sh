#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
PYTHON=${PYTHON:-python3}
RSCRIPT=${RSCRIPT:-Rscript}
"$PYTHON" "$ROOT/examples/minimal/make_demo_h5ad.py"
"$RSCRIPT" "$ROOT/examples/minimal/run_demo.R" "$ROOT"
printf 'Demo completed: %s\n' "$ROOT/examples/minimal/output/demo_records.tsv"
