

CU_KERNEL_DIR       :=${PWD}/cuda-kernel-src
CU_KERNEL           :=${CU_KERNEL_DIR}/kernel.ptx

DEGREE		:=12

test-all :
	cd zkevm-circuits ; make test-all

## Run Evm Circuit benchmarks
evm_bench_cuda: ${CU_KERNEL} ; 
	export DEGREE=${DEGREE}  \
	CUDA=1 \
	CU_KERNEL="${CU_KERNEL}" ; \
	cd zkevm-circuits ; cargo test --profile bench bench_evm_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

evm_bench:
	export DEGREE=${DEGREE} ; cd zkevm-circuits ; cargo test --profile bench bench_evm_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture


## Run State Circuit benchmarks
state_bench_cuda: ${CU_KERNEL} ;
	export DEGREE=${DEGREE}  \
	CUDA=1 \
	CU_KERNEL="${CU_KERNEL}" ; \
	cd zkevm-circuits ; cargo test --profile bench bench_state_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

state_bench: 
	export DEGREE=${DEGREE} ; cd zkevm-circuits ; cargo test --profile bench bench_state_circuit_prover -p circuit-benchmarks --features benches  -- --nocapture

## Run All Circuit benchmarks
circuit_benches: evm_bench state_bench 

circuit_benches_cuda: evm_bench_cuda state_bench_cuda 



${CU_KERNEL} : ${CU_KERNEL_DIR}/field.h  ${CU_KERNEL_DIR}/evaluate_h.cu ;
	nvcc -ptx ${CU_KERNEL_DIR}/evaluate_h.cu -o $@ 


clean :
	cd zkevm-circuits ; cargo clean 
	rm ${CU_KERNEL}