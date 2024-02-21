use crate::bytes_to_u32_digits;
use crate::{HintBuilder, PrivateKey, PublicKey};
use ark_ec::AffineRepr;
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_secp256k1::Fr;
use rand::thread_rng;
use secp256k10_guest::{ComputeHintStore, Evaluator, Hint};

#[test]
fn test_signature() {
    let mut prng = thread_rng();

    let pk = PrivateKey::sample(&mut prng);
    let _vk = PublicKey::from(&pk);

    let z = Fr::rand(&mut prng);

    let _sig = pk.sign(&mut prng, &z);
}

#[test]
fn generate_hint() {
    let mut prng = thread_rng();

    let pk = PrivateKey::sample(&mut prng);
    let z = Fr::rand(&mut prng);

    let sig = pk.sign(&mut prng, &z);

    let (hint, _) = HintBuilder::build(
        &sig.r.into_bigint().to_bytes_le(),
        &sig.s.into_bigint().to_bytes_le(),
        &z.into_bigint().to_bytes_le(),
        sig.v,
    );
    assert!(matches!(hint, Hint::Ok))
}

#[test]
fn evaluate_hint() {
    let mut prng = thread_rng();

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
    let mut compute_hint_provider = ComputeHintStore::new(&compute_hint_vec);

    let eval = Evaluator::new(
        &sig.r.into_bigint().to_bytes_le(),
        &sig.s.into_bigint().to_bytes_le(),
        &z.into_bigint().to_bytes_le(),
        sig.v,
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

    assert_eq!(
        pk_recovered[0..8],
        bytes_to_u32_digits(&vk.0.x().unwrap().into_bigint().to_bytes_le())
    );
    assert_eq!(
        pk_recovered[8..16],
        bytes_to_u32_digits(&vk.0.y().unwrap().into_bigint().to_bytes_le())
    );
}
