use std::collections::HashMap;

use bus_mapping::circuit_input_builder::{self, CopyEvent};
use eth_types::{Address, Field, ToLittleEndian, ToScalar, Word};
use halo2_proofs::halo2curves::bn256::Fr;
use itertools::Itertools;

use crate::{evm_circuit::util::RandomLinearCombination, table::BlockContextFieldTag};

use super::{tx::tx_convert, Bytecode, RwMap, Transaction};

/// Block is the struct used by all circuits, which constains all the needed
/// data for witness generation.
#[derive(Debug, Default, Clone)]
pub struct Block<F> {
    /// The randomness for random linear combination
    pub randomness: F,
    /// Transactions in the block
    pub txs: Vec<Transaction>,
    /// Read write events in the RwTable
    pub rws: RwMap,
    /// Bytecode used in the block
    pub bytecodes: HashMap<Word, Bytecode>,
    /// The block context
    pub context: BlockContext,
    /// Copy events for the EVM circuit's copy table.
    pub copy_events: Vec<CopyEvent>,
    /// Pad evm circuit to make selectors fixed, so vk/pk can be universal.
    pub evm_circuit_pad_to: usize,
    /// Length to rw table rows in state circuit
    pub state_circuit_pad_to: usize,
    /// Inputs to the SHA3 opcode
    pub sha3_inputs: Vec<Vec<u8>>,
}

/// Block context for execution
#[derive(Debug, Default, Clone)]
pub struct BlockContext {
    /// The address of the miner for the block
    pub coinbase: Address,
    /// The gas limit of the block
    pub gas_limit: u64,
    /// The number of the block
    pub number: Word,
    /// The timestamp of the block
    pub timestamp: Word,
    /// The difficulty of the blcok
    pub difficulty: Word,
    /// The base fee, the minimum amount of gas fee for a transaction
    pub base_fee: Word,
    /// The hash of previous blocks
    pub history_hashes: Vec<Word>,
    /// The chain id
    pub chain_id: Word,
}

impl BlockContext {
    /// Assignments for block table
    pub fn table_assignments<F: Field>(&self, randomness: F) -> Vec<[F; 3]> {
        [
            vec![
                [
                    F::from(BlockContextFieldTag::Coinbase as u64),
                    F::zero(),
                    self.coinbase.to_scalar().unwrap(),
                ],
                [
                    F::from(BlockContextFieldTag::Timestamp as u64),
                    F::zero(),
                    self.timestamp.to_scalar().unwrap(),
                ],
                [
                    F::from(BlockContextFieldTag::Number as u64),
                    F::zero(),
                    self.number.to_scalar().unwrap(),
                ],
                [
                    F::from(BlockContextFieldTag::Difficulty as u64),
                    F::zero(),
                    RandomLinearCombination::random_linear_combine(
                        self.difficulty.to_le_bytes(),
                        randomness,
                    ),
                ],
                [
                    F::from(BlockContextFieldTag::GasLimit as u64),
                    F::zero(),
                    F::from(self.gas_limit),
                ],
                [
                    F::from(BlockContextFieldTag::BaseFee as u64),
                    F::zero(),
                    RandomLinearCombination::random_linear_combine(
                        self.base_fee.to_le_bytes(),
                        randomness,
                    ),
                ],
                [
                    F::from(BlockContextFieldTag::ChainId as u64),
                    F::zero(),
                    RandomLinearCombination::random_linear_combine(
                        self.chain_id.to_le_bytes(),
                        randomness,
                    ),
                ],
            ],
            {
                let len_history = self.history_hashes.len();
                self.history_hashes
                    .iter()
                    .enumerate()
                    .map(|(idx, hash)| {
                        [
                            F::from(BlockContextFieldTag::BlockHash as u64),
                            (self.number - len_history + idx).to_scalar().unwrap(),
                            RandomLinearCombination::random_linear_combine(
                                hash.to_le_bytes(),
                                randomness,
                            ),
                        ]
                    })
                    .collect()
            },
        ]
        .concat()
    }
}

impl From<&circuit_input_builder::Block> for BlockContext {
    fn from(block: &circuit_input_builder::Block) -> Self {
        Self {
            coinbase: block.coinbase,
            gas_limit: block.gas_limit,
            number: block.number,
            timestamp: block.timestamp,
            difficulty: block.difficulty,
            base_fee: block.base_fee,
            history_hashes: block.history_hashes.clone(),
            chain_id: block.chain_id,
        }
    }
}

/// Convert a block struct in bus-mapping to a witness block used in circuits
pub fn block_convert(
    block: &circuit_input_builder::Block,
    code_db: &bus_mapping::state_db::CodeDB,
) -> Block<Fr> {
    Block {
        randomness: Fr::from(0xcafeu64),
        context: block.into(),
        rws: RwMap::from(&block.container),
        txs: block
            .txs()
            .iter()
            .enumerate()
            .map(|(idx, tx)| tx_convert(tx, idx + 1))
            .collect(),
        bytecodes: block
            .txs()
            .iter()
            .flat_map(|tx| {
                tx.calls()
                    .iter()
                    .map(|call| call.code_hash)
                    .unique()
                    .into_iter()
                    .map(|code_hash| {
                        let bytecode =
                            Bytecode::new(code_db.0.get(&code_hash).cloned().unwrap_or_default());
                        (bytecode.hash, bytecode)
                    })
            })
            .collect(),
        copy_events: block.copy_events.clone(),
        sha3_inputs: block.sha3_inputs.clone(),
        ..Default::default()
    }
}
