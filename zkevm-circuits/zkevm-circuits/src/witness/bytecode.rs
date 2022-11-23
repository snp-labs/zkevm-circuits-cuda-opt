use bus_mapping::evm::OpcodeId;
use eth_types::{Field, ToLittleEndian, Word};
use halo2_proofs::circuit::Value;
use sha3::{Digest, Keccak256};

use crate::{
    evm_circuit::util::RandomLinearCombination, table::BytecodeFieldTag, util::Challenges,
};

/// Bytecode
#[derive(Clone, Debug)]
pub struct Bytecode {
    /// Hash of bytecode
    pub hash: Word,
    /// Raw bytes
    pub bytes: Vec<u8>,
}

impl Bytecode {
    /// Construct from bytecode bytes
    pub fn new(bytes: Vec<u8>) -> Self {
        let hash = Word::from_big_endian(Keccak256::digest(&bytes).as_slice());
        Self { hash, bytes }
    }

    /// Assignments for bytecode table
    pub fn table_assignments<F: Field>(
        &self,
        challenges: &Challenges<Value<F>>,
    ) -> Vec<[Value<F>; 5]> {
        let n = 1 + self.bytes.len();
        let mut rows = Vec::with_capacity(n);
        let hash = challenges.evm_word().map(|challenge| {
            RandomLinearCombination::random_linear_combine(self.hash.to_le_bytes(), challenge)
        });

        rows.push([
            hash,
            Value::known(F::from(BytecodeFieldTag::Length as u64)),
            Value::known(F::zero()),
            Value::known(F::zero()),
            Value::known(F::from(self.bytes.len() as u64)),
        ]);

        let mut push_data_left = 0;
        for (idx, byte) in self.bytes.iter().enumerate() {
            let mut is_code = true;
            if push_data_left > 0 {
                is_code = false;
                push_data_left -= 1;
            } else if (OpcodeId::PUSH1.as_u8()..=OpcodeId::PUSH32.as_u8()).contains(byte) {
                push_data_left = *byte as usize - (OpcodeId::PUSH1.as_u8() - 1) as usize;
            }
            rows.push([
                hash,
                Value::known(F::from(BytecodeFieldTag::Byte as u64)),
                Value::known(F::from(idx as u64)),
                Value::known(F::from(is_code as u64)),
                Value::known(F::from(*byte as u64)),
            ])
        }
        rows
    }
}

impl From<&eth_types::bytecode::Bytecode> for Bytecode {
    fn from(b: &eth_types::bytecode::Bytecode) -> Self {
        Bytecode::new(b.to_vec())
    }
}
