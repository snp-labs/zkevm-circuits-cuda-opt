use crate::{
    evm_circuit::{
        execution::ExecutionGadget,
        param::{N_BYTES_MEMORY_ADDRESS, N_BYTES_MEMORY_WORD_SIZE, STACK_CAPACITY},
        step::ExecutionState,
        util::{
            common_gadget::RestoreContextGadget,
            constraint_builder::{
                ConstraintBuilder, StepStateTransition,
                Transition::{Delta, To},
            },
            math_gadget::{IsZeroGadget, MinMaxGadget},
            memory_gadget::{MemoryAddressGadget, MemoryExpansionGadget},
            not, CachedRegion, Cell,
        },
        witness::{Block, Call, ExecStep, Transaction},
    },
    table::CallContextFieldTag,
    util::Expr,
};
use bus_mapping::{circuit_input_builder::CopyDataType, evm::OpcodeId};
use eth_types::Field;
use halo2_proofs::{circuit::Value, plonk::Error};

#[derive(Clone, Debug)]
pub(crate) struct ReturnGadget<F> {
    opcode: Cell<F>,

    range: MemoryAddressGadget<F>,

    is_success: Cell<F>,
    restore_context: RestoreContextGadget<F>,

    copy_length: MinMaxGadget<F, N_BYTES_MEMORY_ADDRESS>,
    copy_rw_increase: Cell<F>,
    copy_rw_increase_is_zero: IsZeroGadget<F>,

    return_data_offset: Cell<F>,
    return_data_length: Cell<F>,

    memory_expansion: MemoryExpansionGadget<F, 1, N_BYTES_MEMORY_WORD_SIZE>,
}

// TODO: rename this is reflect the fact that is handles REVERT as well.
impl<F: Field> ExecutionGadget<F> for ReturnGadget<F> {
    const NAME: &'static str = "RETURN";

    const EXECUTION_STATE: ExecutionState = ExecutionState::RETURN;

    fn configure(cb: &mut ConstraintBuilder<F>) -> Self {
        let opcode = cb.query_cell();
        cb.opcode_lookup(opcode.expr(), 1.expr());

        let offset = cb.query_cell();
        let length = cb.query_rlc();
        cb.stack_pop(offset.expr());
        cb.stack_pop(length.expr());
        let range = MemoryAddressGadget::construct(cb, offset, length);

        let is_success = cb.call_context(None, CallContextFieldTag::IsSuccess);
        cb.require_boolean("is_success is boolean", is_success.expr());
        cb.require_equal(
            "if is_success, opcode is RETURN. if not, opcode is REVERT",
            opcode.expr(),
            is_success.expr() * OpcodeId::RETURN.expr()
                + not::expr(is_success.expr()) * OpcodeId::REVERT.expr(),
        );

        // There are 4 cases non-mutually exclusive, A to D, to handle, depending on if
        // the call is, or is not, a create, root, or successful. See the specs at
        // https://github.com/privacy-scaling-explorations/zkevm-specs/blob/master/specs/opcode/F3RETURN.md
        // for more details.
        let is_create = cb.curr.state.is_create.expr();
        let is_root = cb.curr.state.is_root.expr();

        // These are globally defined because they are used across multiple cases.
        let copy_rw_increase = cb.query_cell();
        let copy_rw_increase_is_zero = IsZeroGadget::construct(cb, copy_rw_increase.expr());

        let memory_expansion = MemoryExpansionGadget::construct(
            cb,
            cb.curr.state.memory_word_size.expr(),
            [range.address()],
        );

        // Case A in the specs.
        cb.condition(is_create.clone() * is_success.expr(), |cb| {
            cb.require_equal(
                "increase rw counter once for each memory to bytecode byte copied",
                copy_rw_increase.expr(),
                range.length(),
            );
        });
        cb.condition(
            is_create.clone() * is_success.expr() * not::expr(copy_rw_increase_is_zero.expr()),
            |_cb| {
                // TODO: copy_table_lookup for contract creation.
            },
        );

        // Case B in the specs.
        cb.condition(is_root.expr(), |cb| {
            cb.require_next_state(ExecutionState::EndTx);
            cb.call_context_lookup(
                false.expr(),
                None,
                CallContextFieldTag::IsPersistent,
                is_success.expr(),
            );
            cb.require_step_state_transition(StepStateTransition {
                program_counter: To(0.expr()),
                stack_pointer: To(STACK_CAPACITY.expr()),
                rw_counter: Delta(
                    cb.rw_counter_offset()
                        + not::expr(is_success.expr())
                            * cb.curr.state.reversible_write_counter.expr(),
                ),
                gas_left: Delta(-memory_expansion.gas_cost()),
                reversible_write_counter: To(0.expr()),
                memory_word_size: To(0.expr()),
                ..StepStateTransition::default()
            });
        });

        // Case C in the specs.
        // TODO: have copy_table_lookup update rw_counter expression so that this can go
        // at the end of the constraints.
        let restore_context = cb.condition(not::expr(is_root.expr()), |cb| {
            RestoreContextGadget::construct(
                cb,
                is_success.expr(),
                2.expr() * not::expr(is_create.clone()) + copy_rw_increase.expr(),
                range.offset(),
                range.length(),
                memory_expansion.gas_cost(),
            )
        });

        // Case D in the specs.
        let (return_data_offset, return_data_length, copy_length) = cb.condition(
            not::expr(is_create.clone()) * not::expr(is_root.clone()),
            |cb| {
                let [return_data_offset, return_data_length] = [
                    CallContextFieldTag::ReturnDataOffset,
                    CallContextFieldTag::ReturnDataLength,
                ]
                .map(|field_tag| cb.call_context(None, field_tag));
                let copy_length =
                    MinMaxGadget::construct(cb, return_data_length.expr(), range.length());
                cb.require_equal(
                    "increase rw counter twice for each memory to memory byte copied",
                    copy_length.min() + copy_length.min(),
                    copy_rw_increase.expr(),
                );
                (return_data_offset, return_data_length, copy_length)
            },
        );
        cb.condition(
            not::expr(is_create.clone())
                * not::expr(is_root.clone())
                * not::expr(copy_rw_increase_is_zero.expr()),
            |cb| {
                cb.copy_table_lookup(
                    cb.curr.state.call_id.expr(),
                    CopyDataType::Memory.expr(),
                    cb.next.state.call_id.expr(),
                    CopyDataType::Memory.expr(),
                    range.offset(),
                    range.address(),
                    return_data_offset.expr(),
                    copy_length.min(),
                    0.expr(),
                    copy_rw_increase.expr(),
                );
            },
        );

        // Without this, copy_rw_increase would be unconstrained for non-create root
        // calls.
        cb.condition(not::expr(is_create) * is_root, |cb| {
            cb.require_zero(
                "rw counter is 0 if there is no copy event",
                copy_rw_increase.expr(),
            );
        });

        Self {
            opcode,
            range,
            is_success,
            copy_length,
            copy_rw_increase,
            copy_rw_increase_is_zero,
            return_data_offset,
            return_data_length,
            restore_context,
            memory_expansion,
        }
    }

    fn assign_exec_step(
        &self,
        region: &mut CachedRegion<'_, '_, F>,
        offset: usize,
        block: &Block<F>,
        _: &Transaction,
        call: &Call,
        step: &ExecStep,
    ) -> Result<(), Error> {
        self.opcode.assign(
            region,
            offset,
            Value::known(F::from(step.opcode.unwrap().as_u64())),
        )?;

        let [memory_offset, length] = [0, 1].map(|i| block.rws[step.rw_indices[i]].stack_value());
        let range = self
            .range
            .assign(region, offset, memory_offset, length, block.randomness)?;
        self.memory_expansion
            .assign(region, offset, step.memory_word_size(), [range])?;

        self.is_success
            .assign(region, offset, Value::known(call.is_success.into()))?;

        if !call.is_root && !call.is_create {
            for (cell, value) in [
                (&self.return_data_length, call.return_data_length.into()),
                (&self.return_data_offset, call.return_data_offset.into()),
            ] {
                cell.assign(region, offset, Value::known(value))?;
            }

            self.copy_length.assign(
                region,
                offset,
                F::from(call.return_data_length),
                F::from(length.as_u64()),
            )?;
        }

        let copy_rw_increase = if call.is_create {
            length.as_u64()
        } else if !call.is_root {
            2 * std::cmp::min(call.return_data_length, length.as_u64())
        } else {
            0
        };
        self.copy_rw_increase
            .assign(region, offset, Value::known(F::from(copy_rw_increase)))?;
        self.copy_rw_increase_is_zero
            .assign(region, offset, F::from(copy_rw_increase))?;

        if !call.is_root {
            self.restore_context
                .assign(region, offset, block, call, step, 3)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::test_util::run_test_circuits;
    use eth_types::{
        address, bytecode, evm_types::OpcodeId, geth_types::Account, Address, Bytecode, ToWord,
        Word,
    };
    use itertools::Itertools;
    use mock::TestContext;

    const CALLEE_ADDRESS: Address = Address::repeat_byte(0xff);
    const CALLER_ADDRESS: Address = Address::repeat_byte(0x34);

    fn callee_bytecode(is_return: bool, offset: u64, length: u64) -> Bytecode {
        let memory_address = 2;
        let memory_value = Word::MAX;
        let mut code = bytecode! {
            PUSH32(memory_value)
            PUSH1(memory_address)
            MSTORE
            PUSH32(length)
            PUSH32(offset)
        };
        code.write_op(if is_return {
            OpcodeId::RETURN
        } else {
            OpcodeId::REVERT
        });
        code
    }

    fn caller_bytecode(return_data_offset: u64, return_data_length: u64) -> Bytecode {
        bytecode! {
            PUSH32(return_data_length)
            PUSH32(return_data_offset)
            PUSH32(0) // call data length
            PUSH32(0) // call data offset
            PUSH32(0) // value
            PUSH32(CALLEE_ADDRESS.to_word())
            PUSH32(4000) // gas
            CALL
            STOP
        }
    }

    #[test]
    fn test_return_root() {
        let test_parameters = [(0, 0), (0, 10), (300, 20), (1000, 0)];
        for ((offset, length), is_return) in
            test_parameters.iter().cartesian_product(&[true, false])
        {
            let code = callee_bytecode(*is_return, *offset, *length);
            assert_eq!(
                run_test_circuits(
                    TestContext::<2, 1>::simple_ctx_with_bytecode(code).unwrap(),
                    None
                ),
                Ok(()),
                "(offset, length, is_return) = {:?}",
                (*offset, *length, *is_return)
            );
        }
    }

    #[test]
    fn test_return_nonroot() {
        let test_parameters = [
            ((0, 0), (0, 0)),
            ((0, 10), (0, 10)),
            ((0, 10), (0, 20)),
            ((0, 20), (0, 10)),
            ((64, 1), (0, 10)), // Expands memory in RETURN/REVERT opcode
            ((0, 10), (1000, 0)),
            ((1000, 0), (0, 10)),
            ((1000, 0), (1000, 0)),
        ];
        for (((callee_offset, callee_length), (caller_offset, caller_length)), is_return) in
            test_parameters.iter().cartesian_product(&[true, false])
        {
            let callee = Account {
                address: CALLEE_ADDRESS,
                code: callee_bytecode(*is_return, *callee_offset, *callee_length).into(),
                nonce: Word::one(),
                ..Default::default()
            };
            let caller = Account {
                address: CALLER_ADDRESS,
                code: caller_bytecode(*caller_offset, *caller_length).into(),
                nonce: Word::one(),
                ..Default::default()
            };

            let test_context = TestContext::<3, 1>::new(
                None,
                |accs| {
                    accs[0]
                        .address(address!("0x000000000000000000000000000000000000cafe"))
                        .balance(Word::from(10u64.pow(19)));
                    accs[1].account(&caller);
                    accs[2].account(&callee);
                },
                |mut txs, accs| {
                    txs[0]
                        .from(accs[0].address)
                        .to(accs[1].address)
                        .gas(100000u64.into());
                },
                |block, _tx| block.number(0xcafeu64),
            )
            .unwrap();

            assert_eq!(
                run_test_circuits(test_context, None),
                Ok(()),
                "(callee_offset, callee_length, caller_offset, caller_length, is_return) = {:?}",
                (
                    *callee_offset,
                    *callee_length,
                    *caller_offset,
                    *caller_length,
                    *is_return
                )
            );
        }
    }
}
