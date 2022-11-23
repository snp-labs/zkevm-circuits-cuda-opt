//! The Super Circuit is a circuit that contains all the circuits of the
//! zkEVM in order to achieve two things:
//! - Check the correct integration between circuits via the shared lookup
//!   tables, to verify that the table layouts match.
//! - Allow having a single circuit setup for which a proof can be generated
//!   that would be verified under a single aggregation circuit for the first
//!   milestone.
//!
//! The current implementation contains the following circuits:
//!
//! - [x] EVM Circuit
//! - [ ] State Circuit
//! - [x] Tx Circuit
//! - [x] Bytecode Circuit
//! - [x] Copy Circuit
//! - [ ] Keccak Circuit
//! - [ ] MPT Circuit
//! - [ ] PublicInputs Circuit
//!
//! And the following shared tables, with the circuits that use them:
//!
//! - [x] Copy Table
//!   - [x] Copy Circuit
//!   - [x] EVM Circuit
//! - [ ] Rw Table
//!   - [ ] State Circuit
//!   - [ ] EVM Circuit
//!   - [ ] Copy Circuit
//! - [x] Tx Table
//!   - [x] Tx Circuit
//!   - [x] EVM Circuit
//!   - [x] Copy Circuit
//!   - [ ] PublicInputs Circuit
//! - [x] Bytecode Table
//!   - [x] Bytecode Circuit
//!   - [x] EVM Circuit
//!   - [x] Copy Circuit
//! - [ ] Block Table
//!   - [ ] EVM Circuit
//!   - [ ] PublicInputs Circuit
//! - [ ] MPT Table
//!   - [ ] MPT Circuit
//!   - [ ] State Circuit
//! - [x] Keccak Table
//!   - [ ] Keccak Circuit
//!   - [ ] EVM Circuit
//!   - [x] Bytecode Circuit
//!   - [x] Tx Circuit
//!   - [ ] MPT Circuit

use crate::bytecode_circuit::bytecode_unroller::{
    unroll, Config as BytecodeConfig, UnrolledBytecode,
};
use crate::copy_circuit::CopyCircuit;
use crate::evm_circuit::{table::FixedTableTag, EvmCircuit};
use crate::keccak_circuit::keccak_packed_multi::KeccakPackedConfig as KeccakConfig;
use crate::state_circuit::StateCircuitConfig;
use crate::table::{BlockTable, BytecodeTable, CopyTable, MptTable, RwTable, TxTable};
use crate::tx_circuit::{TxCircuit, TxCircuitConfig};
use crate::util::Challenges;
use crate::witness::{block_convert, Block, MptUpdates};

use bus_mapping::mock::BlockData;
use eth_types::geth_types::{self, GethData};
use eth_types::Field;

use halo2_proofs::arithmetic::CurveAffine;
use halo2_proofs::halo2curves::{
    bn256::Fr,
    group::{Curve, Group},
    secp256k1::Secp256k1Affine,
};
use halo2_proofs::{
    circuit::{Layouter, SimpleFloorPlanner, Value},
    plonk::{Circuit, ConstraintSystem, Error, Expression},
};

use rand::RngCore;
use std::array;
use strum::IntoEnumIterator;

/// Mock randomness used for `SuperCircuit`.
pub const MOCK_RANDOMNESS: u64 = 0x100;

/// Configuration of the Super Circuit
#[derive(Clone)]
pub struct SuperCircuitConfig<F: Field, const MAX_TXS: usize, const MAX_CALLDATA: usize> {
    tx_table: TxTable,
    rw_table: RwTable,
    mpt_table: MptTable,
    bytecode_table: BytecodeTable,
    block_table: BlockTable,
    copy_table: CopyTable,
    evm_circuit: EvmCircuit<F>,
    state_circuit: StateCircuitConfig<F>,
    tx_circuit: TxCircuitConfig<F>,
    bytecode_circuit: BytecodeConfig<F>,
    copy_circuit: CopyCircuit<F>,
    keccak_circuit: KeccakConfig<F>,
}

/// The Super Circuit contains all the zkEVM circuits
#[derive(Default, Debug)]
pub struct SuperCircuit<F: Field, const MAX_TXS: usize, const MAX_CALLDATA: usize> {
    // EVM Circuit
    /// Block witness. Usually derived via
    /// `evm_circuit::witness::block_convert`.
    pub block: Block<F>,
    /// Inputs for the keccak circuit
    pub keccak_inputs: Vec<Vec<u8>>,
    /// Passed down to the evm_circuit. Usually that will be
    /// `FixedTableTag::iter().collect()`.
    pub fixed_table_tags: Vec<FixedTableTag>,
    // Tx Circuit
    /// The transaction circuit that will be used in the `synthesize` step.
    pub tx_circuit: TxCircuit<F, MAX_TXS, MAX_CALLDATA>,
    // Bytecode Circuit
    // bytecodes: Vec<UnrolledBytecode<F>>,
    /// The maximium size for the underlying bytecode circuit.
    pub bytecode_size: usize,
}

impl<F: Field, const MAX_TXS: usize, const MAX_CALLDATA: usize>
    SuperCircuit<F, MAX_TXS, MAX_CALLDATA>
{
    /// Return the number of rows required to verify a given block
    pub fn get_num_rows_required(block: &Block<F>) -> usize {
        let num_rows_evm_circuit = {
            let mut cs = ConstraintSystem::default();
            let config = Self::configure(&mut cs);
            config.evm_circuit.get_num_rows_required(block)
        };
        let num_rows_tx_circuit = TxCircuitConfig::<F>::get_num_rows_required(MAX_TXS);
        num_rows_evm_circuit.max(num_rows_tx_circuit)
    }
}

impl<F: Field, const MAX_TXS: usize, const MAX_CALLDATA: usize> Circuit<F>
    for SuperCircuit<F, MAX_TXS, MAX_CALLDATA>
{
    type Config = SuperCircuitConfig<F, MAX_TXS, MAX_CALLDATA>;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<F>) -> Self::Config {
        let tx_table = TxTable::construct(meta);
        let rw_table = RwTable::construct(meta);
        let mpt_table = MptTable::construct(meta);
        let bytecode_table = BytecodeTable::construct(meta);
        let block_table = BlockTable::construct(meta);
        let q_copy_table = meta.fixed_column();
        let copy_table = CopyTable::construct(meta, q_copy_table);

        let power_of_randomness = array::from_fn(|i| {
            Expression::Constant(F::from(MOCK_RANDOMNESS).pow(&[1 + i as u64, 0, 0, 0]))
        });

        let keccak_circuit = KeccakConfig::configure(meta, power_of_randomness[0].clone());
        let keccak_table = keccak_circuit.keccak_table.clone();

        let evm_circuit = EvmCircuit::configure(
            meta,
            power_of_randomness.clone(),
            &tx_table,
            &rw_table,
            &bytecode_table,
            &block_table,
            &copy_table,
            &keccak_table,
        );
        let state_circuit =
            StateCircuitConfig::configure(meta, power_of_randomness.clone(), &rw_table, &mpt_table);
        let challenges = Challenges::mock(power_of_randomness[0].clone());

        Self::Config {
            tx_table: tx_table.clone(),
            rw_table,
            mpt_table,
            bytecode_table: bytecode_table.clone(),
            block_table,
            copy_table,
            evm_circuit,
            state_circuit,
            copy_circuit: CopyCircuit::configure(
                meta,
                &tx_table,
                &rw_table,
                &bytecode_table,
                copy_table,
                q_copy_table,
                power_of_randomness[0].clone(),
            ),
            tx_circuit: TxCircuitConfig::new(
                meta,
                tx_table,
                keccak_table.clone(),
                challenges.clone(),
            ),
            bytecode_circuit: BytecodeConfig::configure(
                meta,
                bytecode_table,
                keccak_table,
                challenges,
            ),
            keccak_circuit,
        }
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<F>,
    ) -> Result<(), Error> {
        let challenges = Challenges::mock(Value::known(self.block.randomness));

        // --- EVM Circuit ---
        let rws = self.block.rws.table_assignments();
        config
            .evm_circuit
            .load_fixed_table(&mut layouter, self.fixed_table_tags.clone())?;
        config.evm_circuit.load_byte_table(&mut layouter)?;
        config.rw_table.load(
            &mut layouter,
            &rws,
            self.block.state_circuit_pad_to,
            self.block.randomness,
        )?;
        config.state_circuit.load(&mut layouter)?;
        config
            .block_table
            .load(&mut layouter, &self.block.context, self.block.randomness)?;
        config
            .copy_table
            .load(&mut layouter, &self.block, self.block.randomness)?;
        config
            .evm_circuit
            .assign_block(&mut layouter, &self.block)?;
        // --- State Circuit ---
        config.mpt_table.load(
            &mut layouter,
            &MptUpdates::mock_from(&rws),
            self.block.randomness,
        )?;
        config.state_circuit.assign(
            &mut layouter,
            &rws,
            self.block.state_circuit_pad_to,
            self.block.randomness,
        )?;
        // --- Tx Circuit ---
        config.tx_circuit.load(&mut layouter)?;
        self.tx_circuit
            .assign(&config.tx_circuit, &mut layouter, &challenges)?;
        // --- Bytecode Circuit ---
        let bytecodes: Vec<UnrolledBytecode<F>> = self
            .block
            .bytecodes
            .iter()
            .map(|(_, b)| unroll(b.bytes.clone()))
            .collect();
        config.bytecode_circuit.load(&mut layouter)?;
        config.bytecode_circuit.assign(
            &mut layouter,
            self.bytecode_size,
            &bytecodes,
            &challenges,
        )?;
        // --- Keccak Table ---
        config.keccak_circuit.load(&mut layouter)?;
        config.keccak_circuit.assign_from_witness(
            &mut layouter,
            &self.keccak_inputs,
            self.block.randomness,
        )?;
        // --- Copy Circuit ---
        config
            .copy_circuit
            .assign_block(&mut layouter, &self.block, self.block.randomness)?;
        Ok(())
    }
}

impl<const MAX_TXS: usize, const MAX_CALLDATA: usize> SuperCircuit<Fr, MAX_TXS, MAX_CALLDATA> {
    /// From the witness data, generate a SuperCircuit instance with all of the
    /// sub-circuits filled with their corresponding witnesses.
    ///
    /// Also, return with it the minimum required SRS degree for the circuit and
    /// the Public Inputs needed.
    pub fn build(
        geth_data: GethData,
        rng: &mut (impl RngCore + Clone),
    ) -> Result<(u32, Self, Vec<Vec<Fr>>), bus_mapping::Error> {
        let txs = geth_data
            .eth_block
            .transactions
            .iter()
            .map(geth_types::Transaction::from)
            .collect();

        let mut builder =
            BlockData::new_from_geth_data(geth_data.clone()).new_circuit_input_builder();

        builder
            .handle_block(&geth_data.eth_block, &geth_data.geth_traces)
            .expect("could not handle block tx");
        let keccak_inputs = builder.keccak_inputs()?;
        let mut block = block_convert(&builder.block, &builder.code_db);
        block.randomness = Fr::from(MOCK_RANDOMNESS);
        block.state_circuit_pad_to = 1;

        let fixed_table_tags: Vec<FixedTableTag> = FixedTableTag::iter().collect();
        let log2_ceil = |n| u32::BITS - (n as u32).leading_zeros() - (n & (n - 1) == 0) as u32;

        let num_rows_required = Self::get_num_rows_required(&block);

        let k = log2_ceil(
            64 + fixed_table_tags
                .iter()
                .map(|tag| tag.build::<Fr>().count())
                .sum::<usize>(),
        );
        let bytecodes_len = block
            .bytecodes
            .iter()
            .map(|(_, bytecode)| bytecode.bytes.len())
            .sum::<usize>();
        let k = k.max(log2_ceil(64 + bytecodes_len));
        let k = k.max(log2_ceil(64 + num_rows_required));
        log::debug!("super circuit uses k = {}", k);

        let aux_generator = <Secp256k1Affine as CurveAffine>::CurveExt::random(rng).to_affine();
        let chain_id = block.context.chain_id;
        let tx_circuit = TxCircuit::new(aux_generator, chain_id.as_u64(), txs);

        let circuit = Self {
            block,
            fixed_table_tags,
            tx_circuit,
            keccak_inputs,
            // Instead of using 1 << k - NUM_BLINDING_ROWS, we use a much smaller number of enabled
            // rows for the Bytecode Circuit because otherwise it penalizes significantly the
            // MockProver verification time.
            bytecode_size: bytecodes_len + 64,
        };

        // SignVerifyChip -> ECDSAChip -> MainGate instance column
        let instances = vec![vec![]];
        Ok((k, circuit, instances))
    }
}

#[cfg(test)]
mod super_circuit_tests {
    use super::*;
    use ethers_signers::{LocalWallet, Signer};
    use halo2_proofs::dev::MockProver;
    use mock::{TestContext, MOCK_CHAIN_ID};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use std::collections::HashMap;

    use eth_types::{address, bytecode, geth_types::GethData, Word};

    // High memory usage test.  Run in serial with:
    // `cargo test [...] serial_ -- --ignored --test-threads 1`
    #[ignore]
    #[test]
    fn serial_test_super_circuit() {
        let mut rng = ChaCha20Rng::seed_from_u64(2);

        let chain_id = (*MOCK_CHAIN_ID).as_u64();

        let bytecode = bytecode! {
            GAS
            STOP
        };

        let wallet_a = LocalWallet::new(&mut rng).with_chain_id(chain_id);

        let addr_a = wallet_a.address();
        let addr_b = address!("0x000000000000000000000000000000000000BBBB");

        let mut wallets = HashMap::new();
        wallets.insert(wallet_a.address(), wallet_a);

        let mut block: GethData = TestContext::<2, 1>::new(
            None,
            |accs| {
                accs[0]
                    .address(addr_b)
                    .balance(Word::from(1u64 << 20))
                    .code(bytecode);
                accs[1].address(addr_a).balance(Word::from(1u64 << 20));
            },
            |mut txs, accs| {
                txs[0]
                    .from(accs[1].address)
                    .to(accs[0].address)
                    .gas(Word::from(1_000_000u64));
            },
            |block, _tx| block.number(0xcafeu64),
        )
        .unwrap()
        .into();

        block.sign(&wallets);

        let (k, circuit, instance) =
            SuperCircuit::<_, 1, 32>::build(block, &mut ChaCha20Rng::seed_from_u64(2)).unwrap();
        let prover = MockProver::run(k, &circuit, instance).unwrap();
        let res = prover.verify_par();
        if let Err(err) = res {
            eprintln!("Verification failures:");
            eprintln!("{:#?}", err);
            panic!("Failed verification");
        }
    }
}
