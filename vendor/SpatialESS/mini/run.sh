#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
Rscript "$ROOT/run_minimal.R"
printf "SpatialESS mini reproduction completed.\n"
