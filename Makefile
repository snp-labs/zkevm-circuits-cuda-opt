


DEGREE		:=18

test-all :
	cd zkevm-circuits ; make test-all

test_benches: ## Compiles the benchmarks
	export DEGREE=${DEGREE} ; cd zkevm-circuits ; cargo test --verbose --release --all-features -p circuit-benchmarks --no-run

evm_bench: ## Run Evm Circuit benchmarks
	export DEGREE=${DEGREE} ; cd zkevm-circuits ; cargo test --profile bench bench_evm_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

state_bench: ## Run State Circuit benchmarks
	export DEGREE=${DEGREE} ; cd zkevm-circuits ; cargo test --profile bench bench_state_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

circuit_benches: evm_bench state_bench ## Run All Circuit benchmarks

