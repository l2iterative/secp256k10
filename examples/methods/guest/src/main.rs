use risc0_zkvm::guest::env;
use secp256k10_guest::{Evaluator, EvaluationResult, Hint};
use serde::{Serialize, Deserialize};
use l2r0_small_serde::from_slice_compact;
use l2r0_profiler_guest::*;

#[derive(Debug, Serialize, Deserialize)]
pub struct Task {
    pub r: Vec<u8>,
    pub s: Vec<u8>,
    pub v: u8,
    pub z: Vec<u8>,
}

fn main() {
    l2r0_profiler_guest::init_trace_logger();
    start_timer!("Total");

    start_timer!("Read the task slice");
    let task_slice: Vec<u32> = env::read();
    stop_start_timer!("Convert the task slice to task");
    let task: Task = from_slice_compact(&task_slice).unwrap();
    stop_start_timer!("Read the hint");
    let hint: Hint = env::read();

    stop_start_timer!("Evaluation");
    let eval = Evaluator::new(
        &task.r,
        &task.s,
        &task.z,
        task.v,
        hint,
    );

    let res = eval.evaluate();
    assert!(matches!(res, EvaluationResult::Ok(_)));

    let pk_recovered = match res {
        EvaluationResult::Ok(v) => v,
        EvaluationResult::Err(_) => {
            unreachable!()
        }
    };
    stop_timer!();

    env::commit_slice(&pk_recovered);
    stop_timer!();
}