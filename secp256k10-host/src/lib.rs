use num_bigint::{BigInt, BigUint};
use std::sync::OnceLock;
use ark_ff::BigInteger;

static MODULUS_Q: OnceLock<BigUint> = OnceLock::new();
static MODULUS_N: OnceLock<BigUint> = OnceLock::new();
static N11: OnceLock<BigInt> = OnceLock::new();
static N12: OnceLock<BigInt> = OnceLock::new();
static N21: OnceLock<BigInt> = OnceLock::new();
static N22: OnceLock<BigInt> = OnceLock::new();

mod table_generation;
pub use table_generation::TableGeneration;

mod hinter;
pub use hinter::HintBuilder;

mod signer;
pub use signer::*;

#[cfg(test)]
mod integration_test;

pub use secp256k10_guest::Hint;

#[inline]
pub fn bytes_to_u32_digits(fe: &[u8]) -> [u32; 8] {
    let mut bytes = [0u8; 32];
    bytes.copy_from_slice(fe);
    bytemuck::cast::<[u8; 32], [u32; 8]>(bytes)
}
