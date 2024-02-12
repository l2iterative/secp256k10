use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{BigInteger, Field, PrimeField, UniformRand, Zero};
use ark_secp256k1::{Affine, Fr};
use num_bigint::BigUint;
use rand::{CryptoRng, RngCore};
use std::ops::Mul;

pub struct PrivateKey(pub Fr);
pub struct PublicKey(pub Affine);

pub struct Signature {
    pub r: Fr,
    pub s: Fr,
    pub v: u8,
}

impl From<&PrivateKey> for PublicKey {
    fn from(value: &PrivateKey) -> Self {
        Self(Affine::generator().mul(value.0).into_affine())
    }
}

impl PrivateKey {
    pub fn sample<R: CryptoRng + RngCore>(prng: &mut R) -> Self {
        Self(Fr::rand(prng))
    }

    pub fn sign<R: CryptoRng + RngCore>(&self, prng: &mut R, z: &Fr) -> Signature {
        loop {
            let k = Fr::rand(prng);
            if k.is_zero() {
                continue;
            }

            let point_r = Affine::generator().mul(&k).into_affine();
            let mut r = BigUint::from(point_r.x().unwrap().into_bigint());

            let r_backup = r.clone();
            r = r % BigUint::from(Fr::MODULUS);

            if r.is_zero() {
                continue;
            }

            let is_overflow = r_backup != r;

            let r = Fr::from_le_bytes_mod_order(&r.to_bytes_le());

            let s = (r * &self.0 + z) * k.inverse().unwrap();
            if s.is_zero() {
                continue;
            }

            let is_r_y_odd = point_r.y().unwrap().into_bigint().is_odd();
            let recid = (is_overflow as u8) * 2 + (is_r_y_odd as u8);
            return Signature { r, s, v: recid };
        }
    }
}
