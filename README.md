# Circuits for zkEVM


## Running benchmarks

-   EVM Circuit prover benches. -> `make evm_bench`.
-   EVM Circuit prover benches (using cuda gpu). -> `make evm_bench_cuda`.

-   State Circuit prover benches. -> `make state_bench`
-   State Circuit prover benches (using cuda gpu). -> `make state_bench_cuda`

You can also run all benchmarks by running: `make circuit_benches` or `make circuit_benches_cuda`.


There are currently several benchmarks to run in the workspace in regards to the circuits. 
All use the DEGREE env var to specify the degree of the K parameter that you want to use for your circuit in the bench process.

You can change the DEGREE env value in the Makefile  