help: ## Display this help screen
	@grep -h \
		-E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

clippy: ## Run clippy checks over all workspace members
	@cargo check --all-features
	@cargo clippy --all-features --all-targets -- -D warnings

doc: ## Generate and tests docs including private items
	@cargo doc --no-deps --all --document-private-items

fmt: ## Check whether the code is formated correctly
	@cargo check --all-features
	@cargo fmt --all -- --check

test: ## Run tests for all the workspace members
	# Run light tests
	@cargo test --release --all --all-features --exclude integration-tests --exclude circuit-benchmarks
	# Run heavy tests serially to avoid OOM
	@cargo test --release --all --all-features --exclude integration-tests --exclude circuit-benchmarks serial_ -- --ignored --test-threads 1

test_doc: ## Test the docs
	@cargo test --release --all --all-features --doc

test_benches: ## Compiles the benchmarks
	@cargo test --verbose --release --all-features -p circuit-benchmarks --no-run

test-all: fmt doc clippy test_doc test_benches test ## Run all the CI checks locally (in your actual toolchain)

evm_bench: ## Run Evm Circuit benchmarks
	@cargo test --profile bench bench_evm_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

state_bench: ## Run State Circuit benchmarks
	@cargo test --profile bench bench_state_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

bit_keccak_bench: ## Run Bit Keccak Circuit benchmarks
	@cargo test --profile bench bench_bit_keccak_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

packed_keccak_bench: ## Run Packed Keccak Circuit benchmarks
	@cargo test --profile bench bench_packed_keccak_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

packed_multi_keccak_bench: ## Run Packed Multi Keccak Circuit benchmarks
	@cargo test --profile bench bench_packed_multi_keccak_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

bytecode_bench: ## Run Bytecode Circuit benchmarks
	@cargo test --profile bench bench_bytecode_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

circuit_benches: evm_bench state_bench ## Run All Circuit benchmarks


.PHONY: clippy doc fmt test test_benches test-all evm_bench state_bench circuit_benches help
