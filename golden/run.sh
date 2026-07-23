#!/usr/bin/env bash
# Golden-run regression gate — env wrapper around gate.py.
#
# Sets the lean toolchain (Rust render bins on TIMSIM_BIN; the imspy-free NECRO venv that carries
# necroflow + timsim-predict + timsim-eval + DiaNN reachability) then runs both closed loops and diffs
# the pinned baseline. Override any path via the env vars below; defaults suit the dev box.
set -euo pipefail

# The venv with necroflow + timsim-predict console scripts (timsim-rt/ccs/fragments) + timsim-eval.
NECRO="${NECRO:-/scratch/timsim-demo/timsim-2/timsim-necro/NECRO}"
# Rust render binaries (timsim-render, timsim-render-thermo, timsim-spectra, …).
export TIMSIM_BIN="${TIMSIM_BIN:-/scratch/timsim-demo/timsim-cli/target/release}"
# DiaNN + .NET (the Thermo .raw leg reads natively only via .NET 8; the Bruker .d leg needs neither).
export TIMSIM_DIANN="${TIMSIM_DIANN:-/home/administrator/dia-nn/diann-2.5.0/diann-linux}"
export DOTNET_ROOT="${DOTNET_ROOT:-$HOME/.dotnet}"

# NECRO's bin first so the DAG's bare `timsim-rt`/`python -m timsim_eval` resolve to the lean venv.
export PATH="$NECRO/bin:$DOTNET_ROOT:$PATH"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$NECRO/bin/python" "$HERE/gate.py" "$@"
