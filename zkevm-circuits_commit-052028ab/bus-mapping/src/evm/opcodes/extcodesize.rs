use crate::circuit_input_builder::{CircuitInputStateRef, ExecStep};
use crate::evm::Opcode;
use crate::operation::{TxAccessListAccountOp, RW};
use crate::Error;
use eth_types::{GethExecStep, ToAddress, ToWord};

#[derive(Debug, Copy, Clone)]
pub(crate) struct Extcodesize;

impl Opcode for Extcodesize {
    fn gen_associated_ops(
        state: &mut CircuitInputStateRef,
        geth_steps: &[GethExecStep],
    ) -> Result<Vec<ExecStep>, Error> {
        // TODO: finish this, only access list part is done
        let geth_step = &geth_steps[0];
        let mut exec_step = state.new_step(geth_step)?;

        let external_address = geth_steps[0].stack.last()?.to_address();
        state.stack_read(
            &mut exec_step,
            geth_step.stack.last_filled(),
            external_address.to_word(),
        )?;
        let is_warm = state.sdb.check_account_in_access_list(&external_address);
        state.push_op_reversible(
            &mut exec_step,
            RW::WRITE,
            TxAccessListAccountOp {
                tx_id: state.tx_ctx.id(),
                address: external_address,
                is_warm: true,
                is_warm_prev: is_warm,
            },
        )?;

        state.stack_write(
            &mut exec_step,
            geth_steps[1].stack.nth_last_filled(0),
            geth_steps[1].stack.nth_last(0)?,
        )?;

        Ok(vec![exec_step])
    }
}
