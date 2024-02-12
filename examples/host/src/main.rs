use ark_ec::AffineRepr;
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_secp256k1::Fr;
use l2r0_small_serde::to_vec_compact;
use methods::{METHOD_ELF, METHOD_ID};
use risc0_zkvm::{default_prover, ExecutorEnv};
use secp256k10_host::{Hint, HintBuilder, PrivateKey, PublicKey};
use serde::{Deserialize, Serialize};

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

    let hint = HintBuilder::build(
        &sig.r.into_bigint().to_bytes_le(),
        &sig.s.into_bigint().to_bytes_le(),
        &z.into_bigint().to_bytes_le(),
        sig.v,
    );
    assert!(matches!(hint, Hint::Ok(_)));

    let task = Task {
        r: sig.r.into_bigint().to_bytes_le(),
        s: sig.s.into_bigint().to_bytes_le(),
        v: sig.v,
        z: z.into_bigint().to_bytes_le(),
    };

    let task_to_slice = to_vec_compact(&task).unwrap();

    let env = ExecutorEnv::builder()
        .write(&task_to_slice)
        .unwrap()
        .write(&hint)
        .unwrap()
        .build()
        .unwrap();

    let prover = default_prover();

    let receipt = prover.prove(env, METHOD_ELF).unwrap();
    receipt.verify(METHOD_ID).unwrap();

    let journal = receipt.journal.bytes;
    assert_eq!(
        journal[0..32],
        vk.0.x().unwrap().into_bigint().to_bytes_le()
    );
    assert_eq!(
        journal[32..64],
        vk.0.y().unwrap().into_bigint().to_bytes_le()
    );
}
