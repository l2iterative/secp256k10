use ark_ec::AffineRepr;
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_secp256k1::Fr;
use l2r0_small_serde::to_vec_compact;
use methods::METHOD_ELF;
use risc0_zkvm::{ExecutorEnv, ExecutorImpl};
use secp256k10_host::{Hint, HintBuilder, PrivateKey, PublicKey};
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::rc::Rc;

use l2r0_profiler_host::CycleTracer;

#[derive(Debug, Serialize, Deserialize)]
pub struct Task {
    pub r: Vec<u8>,
    pub s: Vec<u8>,
    pub v: u8,
    pub z: Vec<u8>,
}

fn main() {
    let mut prng = rand::thread_rng();

    let pk = PrivateKey::sample(&mut prng);
    let vk = PublicKey::from(&pk);

    let z = Fr::rand(&mut prng);

    let sig = pk.sign(&mut prng, &z);

    let (hint, Some(compute_hint)) = HintBuilder::build(
        &sig.r.into_bigint().to_bytes_le(),
        &sig.s.into_bigint().to_bytes_le(),
        &z.into_bigint().to_bytes_le(),
        sig.v,
    ) else {
        unreachable!()
    };
    assert!(matches!(hint, Hint::Ok));

    let compute_hint_vec = compute_hint.to_vec();

    let task = Task {
        r: sig.r.into_bigint().to_bytes_le(),
        s: sig.s.into_bigint().to_bytes_le(),
        v: sig.v,
        z: z.into_bigint().to_bytes_le(),
    };

    let task_to_slice = to_vec_compact(&task).unwrap();

    let cycle_tracer = Rc::new(RefCell::new(CycleTracer::default()));

    let env = ExecutorEnv::builder()
        .write(&task_to_slice)
        .unwrap()
        .write(&hint)
        .unwrap()
        .write(&(compute_hint_vec.len() as u32))
        .unwrap()
        .write_slice(&compute_hint_vec)
        .trace_callback(|e| {
            cycle_tracer.borrow_mut().handle_event(e);
            Ok(())
        })
        .build()
        .unwrap();

    let mut exec = ExecutorImpl::from_elf(env, METHOD_ELF).unwrap();
    let session = exec.run().unwrap();

    println!("number of cycles: {}", session.user_cycles);

    cycle_tracer.borrow().print();

    let journal = session.journal.unwrap().bytes;
    assert_eq!(
        journal[0..32],
        vk.0.x().unwrap().into_bigint().to_bytes_le()
    );
    assert_eq!(
        journal[32..64],
        vk.0.y().unwrap().into_bigint().to_bytes_le()
    );
}
