use ark_ec::AffineRepr;
use ark_ff::{BigInteger, Field, Fp256, LegendreSymbol, PrimeField, Zero};
use ark_secp256k1::{Fq, Fr};
use num_bigint::BigUint;
use std::str::FromStr;
use std::sync::OnceLock;

static MODULUS_Q: OnceLock<BigUint> = OnceLock::new();
static MODULUS_N: OnceLock<BigUint> = OnceLock::new();

mod table_generation;

pub struct SignatureOkHint {
    pub r_y: [u8; 32],
    pub r_inv: [u8; 32],
    pub hints: Vec<[u8; 32]>,
}

pub enum SignatureFailedHint {
    RIsZero,
    SIsZero,
    RGeN,
    SGeN,
    RPlusNGeQ,
    InvalidV,
    YIsImaginary([u8; 32]),
    RecoveredKeyIsPointOfInfinity,
}

pub enum SignatureHint {
    Ok(SignatureOkHint),
    Failed(SignatureFailedHint),
}

// The implementation here follows the code in
// https://github.com/rust-bitcoin/rust-secp256k1/blob/master/secp256k1-sys/depend/secp256k1/src/modules/recovery/main_impl.h#L87
// https://github.com/RustCrypto/signatures/blob/master/ecdsa/src/recovery.rs
// which is under the MIT license
impl SignatureHint {
    pub fn new(r: &[u8; 32], s: &[u8; 32], z: &[u8; 32], recid: u8) -> Self {
        let r = BigUint::from_bytes_le(r);
        let s = BigUint::from_bytes_le(s);

        if r.is_zero() {
            return Self::Failed(SignatureFailedHint::RIsZero);
        }

        if s.is_zero() {
            return Self::Failed(SignatureFailedHint::SIsZero);
        }

        let q = MODULUS_Q.get_or_init(|| {
            BigUint::from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671663",
            )
            .unwrap()
        });

        let n = MODULUS_N.get_or_init(|| {
            BigUint::from_str(
                "115792089237316195423570985008687907852837564279074904382605163141518161494337",
            )
            .unwrap()
        });

        if r.ge(&n) {
            return Self::Failed(SignatureFailedHint::RGeN);
        }

        if s.ge(&n) {
            return Self::Failed(SignatureFailedHint::SGeN);
        }

        if !(recid < 4) {
            return Self::Failed(SignatureFailedHint::InvalidV);
        }

        let mut r_mod_q = r.clone();

        // x is reduced
        if recid & 2 == 1 {
            r_mod_q = r_mod_q + n;
            if r_mod_q.ge(&q) {
                return Self::Failed(SignatureFailedHint::RPlusNGeQ);
            }
        }

        let r_x = Fq::from(r_mod_q);

        let x3_plus_ax_plus_b = r_x.square() * r_x + Fq::from(7u8);
        let legendre = x3_plus_ax_plus_b.legendre();
        if legendre == LegendreSymbol::QuadraticNonResidue {
            // One QNR is 3
            let mut witness = x3_plus_ax_plus_b.double() + x3_plus_ax_plus_b;
            witness = witness.sqrt().unwrap();

            let mut witness_bytes = [0u8; 32];
            witness_bytes.copy_from_slice(&witness.into_bigint().to_bytes_le());
            return Self::Failed(SignatureFailedHint::YIsImaginary(witness_bytes));
        }

        let mut r_y = x3_plus_ax_plus_b.sqrt().unwrap();
        if recid & 1 == 1 {
            // y is odd
            if r_y.into_bigint().is_even() {
                r_y.neg_in_place();
            }
        } else {
            // y is even
            if r_y.into_bigint().is_odd() {
                r_y.neg_in_place();
            }
        }

        // Note: r_inv is over Fr
        let r_inv = Fr::from(r).inverse().unwrap();
        let z = Fr::from_le_bytes_mod_order(z);
        let s = Fr::from(s);

        let u1 = -(r_inv * z);
        let u2 = r_inv * s;

        let mut u1_sum = None;
        let mut hints = Vec::new();

        let u1_bits = u1.into_bigint().to_bits_le();
        for (i, bit) in u1_bits.iter().enumerate() {
            if *bit {
                let (x2, y2) = secp256k10_guest::G_TABLES[i];

                let x2 = Fq::from_le_bytes_mod_order(&x2);
                let y2 = Fq::from_le_bytes_mod_order(&y2);

                if u1_sum.is_none() {
                    u1_sum = Some((x2, y2));
                } else {
                    let (x1, y1) = u1_sum.as_ref().unwrap();

                    let slope = (y1 - &y2) * (x1 - &x2).inverse().unwrap();
                    let mut slope_bytes = [0u8; 32];
                    slope_bytes.copy_from_slice(&slope.into_bigint().to_bytes_le());
                    hints.push(slope_bytes);

                    let x3 = slope.square() - x1 - x2;
                    let y3 = slope * &(x1 - &x3) - y1;
                    u1_sum = Some((x3, y3));
                }
            }
        }

        // The implementation as follows relies on an observation about secp256k1
        // The doubling algorithm would fail if y = 0.
        //
        // Fortunately, since the curve equation is y^2 = x^3 + 7.
        // And -7 is a non-cubic residue.
        //
        // In other words, no reasonable points would have y = 0.

        let mut u2_sum = None;
        let mut u2_cur = (r_x, r_y);
        let u2_bits = u2.into_bigint().to_bits_le();
        for (i, bit) in u2_bits.iter().enumerate() {
            if i != 0 {
                let (x1, y1) = &u2_cur;

                let x1_sqr = x1.square();
                let y1_dbl = y1.double();

                let slope = (x1_sqr.double() + x1_sqr) * y1_dbl.inverse().unwrap();
                let mut slope_bytes = [0u8; 32];
                slope_bytes.copy_from_slice(&slope.into_bigint().to_bytes_le());
                hints.push(slope_bytes);

                let x3 = slope.square() - x1 - x1;
                let y3 = slope * &(x1 - &x3) - y1;
                u2_cur = (x3, y3);
            }

            if *bit {
                let (x2, y2) = &u2_cur;

                if u2_sum.is_none() {
                    u2_sum = Some((*x2, *y2));
                } else {
                    let (x1, y1) = u2_sum.unwrap();

                    let slope = (y1 - y2) * (x1 - x2).inverse().unwrap();
                    let mut slope_bytes = [0u8; 32];
                    slope_bytes.copy_from_slice(&slope.into_bigint().to_bytes_le());
                    hints.push(slope_bytes);

                    let x3 = slope.square() - &x1 - x2;
                    let y3 = slope * &(x1 - &x3) - y1;
                    u2_sum = Some((x3, y3));
                }
            }
        }

        // Adding u1_sum and u2_sum together would have a corner case if their x coordinates are the same.
        // In that case, if the y coordinate is the same, then we double it
        // Otherwise, return an error
        match (u1_sum, u2_sum) {
            (Some((u1x, u1y)), Some((u2x, u2y))) => {
                if u1x == u2x {
                    if u1y != u2y {
                        return Self::Failed(SignatureFailedHint::RecoveredKeyIsPointOfInfinity);
                    } else {
                        let u1x_sqr = u1x.square();
                        let u1y_dbl = u1y.double();

                        let slope = (u1x_sqr.double() + u1x_sqr) * u1y_dbl.inverse().unwrap();

                        let mut slope_bytes = [0u8; 32];
                        slope_bytes.copy_from_slice(&slope.into_bigint().to_bytes_le());
                        hints.push(slope_bytes);
                    }
                } else {
                    let slope = (u1y - u2y) * (u1x - u2x).inverse().unwrap();
                    let mut slope_bytes = [0u8; 32];
                    slope_bytes.copy_from_slice(&slope.into_bigint().to_bytes_le());
                    hints.push(slope_bytes);
                }
            }
            (_, _) => {
                // only possible when z is zero, which would not happen with non-negligible probability
                unreachable!()
            }
        }

        let mut r_y_bytes = [0u8; 32];
        r_y_bytes.copy_from_slice(&r_y.into_bigint().to_bytes_le());

        let mut r_inv_bytes = [0u8; 32];
        r_inv_bytes.copy_from_slice(&r_inv.into_bigint().to_bytes_le());

        let res = SignatureOkHint {
            r_y: r_y_bytes,
            r_inv: r_inv_bytes,
            hints,
        };

        return Self::Ok(res);
    }
}
