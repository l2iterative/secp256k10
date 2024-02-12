use risc0_zkvm::guest::env;
use secp256k10_guest::{Evaluator, EvaluationResult, Hint};
use serde::{Serialize, Deserialize};
use l2r0_small_serde::from_slice_compact;

#[derive(Debug, Serialize, Deserialize)]
pub struct Task {
    pub r: Vec<u8>,
    pub s: Vec<u8>,
    pub v: u8,
    pub z: Vec<u8>,
}

fn main() {
    let task_slice: Vec<u32> = env::read();
    let task: Task = from_slice_compact(&task_slice).unwrap();
    let hint: Hint = env::read();

    let eval = Evaluator::new(
        &task.r,
        &task.s,
        &task.z,
        task.v,
        hint,
    );

    let cycle_count_before = env::cycle_count();

    let res = eval.evaluate();
    assert!(matches!(res, EvaluationResult::Ok(_)));

    let pk_recovered = match res {
        EvaluationResult::Ok(v) => v,
        EvaluationResult::Err(_) => {
            unreachable!()
        }
    };
    let cycle_count_after = env::cycle_count();

    env::commit_slice(&pk_recovered);
    eprintln!("before evaluation: {}, after evaluation: {}", cycle_count_before, cycle_count_after);
}