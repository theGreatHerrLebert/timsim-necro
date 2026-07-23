# timsim-necro — set up the lean simulator. Everything is ingested from published/git sources; no rustims.
PREFIX ?= $(HOME)/.local

.PHONY: predict-deps rust-bins setup
predict-deps:   ## the lean prediction + orchestration stack (necroflow + timsim-predict -> pepdl -> mscorepy)
	pip install -r requirements.txt

rust-bins:      ## the timsim-cli protocol/render binaries — installed from git, crates.io deps only
	cargo install --git https://github.com/theGreatHerrLebert/timsim-cli --features tdf,thermo --root $(PREFIX)

setup: predict-deps rust-bins   ## everything the simulator needs — no rustims, no imspy
