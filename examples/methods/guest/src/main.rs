use risc0_zkvm::guest::env;
use secp256k10_guest::{Evaluator, ComputeHintBuffer, Hint};
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
    stop_start_timer!("Get the compute hint length");
    let compute_hint_length: u32 = env::read();
    stop_start_timer!("Initialize the compute hint");

    stop_start_timer!("Evaluation");
    let mut compute_hint_provider = ComputeHintBuffer::new(compute_hint_length as usize);
    let eval = Evaluator::new(
        &task.r,
        &task.s,
        &task.z,
        task.v,
        hint,
    );

    let res = eval.evaluate(&mut compute_hint_provider);
    assert!(matches!(res, Ok(_)));

    let pk_recovered = match res {
        Ok(v) => v,
        Err(_) => {
            unreachable!()
        }
    };
    stop_timer!();

    env::commit_slice(&pk_recovered);
    stop_timer!();
}