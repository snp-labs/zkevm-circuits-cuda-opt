//! State circuit benchmarks

#[cfg(test)]
mod tests {
    use ark_std::{end_timer, start_timer};
    use bus_mapping::evm::OpcodeId;
    use eth_types::Field;
    use halo2_proofs::plonk::{create_proof, keygen_pk, keygen_vk, verify_proof};
    use halo2_proofs::poly::kzg::commitment::{KZGCommitmentScheme, ParamsKZG, ParamsVerifierKZG};
    use halo2_proofs::poly::kzg::multiopen::{ProverSHPLONK, VerifierSHPLONK};
    use halo2_proofs::poly::kzg::strategy::SingleStrategy;
    use halo2_proofs::{
        halo2curves::bn256::{Bn256, Fr, G1Affine},
        poly::commitment::ParamsProver,
        transcript::{
            Blake2bRead, Blake2bWrite, Challenge255, TranscriptReadBuffer, TranscriptWriterBuffer,
        },
    };
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use std::env::var;
    use zkevm_circuits::bytecode_circuit::bytecode_unroller::{unroll, UnrolledBytecode};
    use zkevm_circuits::bytecode_circuit::dev::BytecodeCircuitTester;

    #[cfg_attr(not(feature = "benches"), ignore)]
    #[test]
    fn bench_bytecode_circuit_prover() {
        let degree: u32 = var("DEGREE")
            .unwrap_or_else(|_| "15".to_string())
            .parse()
            .expect("Cannot parse DEGREE env var as u32");

        // Contract code size exceeds 24576 bytes may not be deployable on Mainnet.
        const MAX_BYTECODE_LEN: usize = 24576;

        let num_rows = 1 << degree;
        const NUM_BLINDING_ROWS: usize = 7 - 1;
        let max_bytecode_row_num = num_rows - NUM_BLINDING_ROWS;
        let bytecode_len = std::cmp::min(MAX_BYTECODE_LEN, max_bytecode_row_num);
        let bytecodes_num: usize = max_bytecode_row_num / bytecode_len;

        // Create the circuit
        let bytecode_circuit = BytecodeCircuitTester::<Fr>::new(
            fillup_codebytes(bytecodes_num, bytecode_len),
            2usize.pow(degree),
        );

        // Initialize the polynomial commitment parameters
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        // Bench setup generation
        let setup_message = format!("Setup generation with degree = {}", degree);
        let start1 = start_timer!(|| setup_message);
        let general_params = ParamsKZG::<Bn256>::setup(degree, &mut rng);
        let verifier_params: ParamsVerifierKZG<Bn256> = general_params.verifier_params().clone();
        end_timer!(start1);

        // Initialize the proving key
        let vk = keygen_vk(&general_params, &bytecode_circuit).expect("keygen_vk should not fail");
        let pk =
            keygen_pk(&general_params, vk, &bytecode_circuit).expect("keygen_pk should not fail");
        // Create a proof
        let mut transcript = Blake2bWrite::<_, G1Affine, Challenge255<_>>::init(vec![]);

        // Bench proof generation time
        let proof_message = format!("Bytecode Proof generation with {} rows", degree);
        let start2 = start_timer!(|| proof_message);
        create_proof::<
            KZGCommitmentScheme<Bn256>,
            ProverSHPLONK<'_, Bn256>,
            Challenge255<G1Affine>,
            XorShiftRng,
            Blake2bWrite<Vec<u8>, G1Affine, Challenge255<G1Affine>>,
            BytecodeCircuitTester<Fr>,
        >(
            &general_params,
            &pk,
            &[bytecode_circuit],
            &[&[]],
            rng,
            &mut transcript,
        )
        .expect("proof generation should not fail");
        let proof = transcript.finalize();
        end_timer!(start2);

        // Bench verification time
        let start3 = start_timer!(|| "Bytecode Proof verification");
        let mut verifier_transcript = Blake2bRead::<_, G1Affine, Challenge255<_>>::init(&proof[..]);
        let strategy = SingleStrategy::new(&general_params);

        verify_proof::<
            KZGCommitmentScheme<Bn256>,
            VerifierSHPLONK<'_, Bn256>,
            Challenge255<G1Affine>,
            Blake2bRead<&[u8], G1Affine, Challenge255<G1Affine>>,
            SingleStrategy<'_, Bn256>,
        >(
            &verifier_params,
            pk.get_vk(),
            strategy,
            &[&[]],
            &mut verifier_transcript,
        )
        .expect("failed to verify bench circuit");
        end_timer!(start3);
    }

    /// fill bytecodes_num * bytecode_len bytes to the witness table
    fn fillup_codebytes<F: Field>(
        bytecodes_num: usize,
        bytecode_len: usize,
    ) -> Vec<UnrolledBytecode<F>> {
        let mut codebytes = vec![];
        (0..bytecodes_num).for_each(|_| {
            let bytecodes = (0..bytecode_len)
                .map(|v| {
                    OpcodeId::try_from(v as u8)
                        .unwrap_or(OpcodeId::STOP)
                        .as_u8()
                })
                .collect::<Vec<u8>>();
            let unrolled_bytes = unroll::<F>(bytecodes);
            codebytes.push(unrolled_bytes);
        });
        codebytes
    }
}
