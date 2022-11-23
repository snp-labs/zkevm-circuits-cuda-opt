//! Witness for all circuits.
//! The `Block<F>` is the witness struct post-processed from geth traces and
//! used to generate witnesses for circuits.

mod block;
pub use block::{block_convert, Block, BlockContext};
mod bytecode;
pub use bytecode::Bytecode;
mod call;
pub use call::Call;
mod mpt;
pub use mpt::{MptUpdate, MptUpdateRow, MptUpdates};
mod rw;
pub use rw::{Rw, RwMap, RwRow};
mod step;
pub use step::ExecStep;
mod tx;
pub use tx::Transaction;
