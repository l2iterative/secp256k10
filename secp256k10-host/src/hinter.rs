use crate::{bytes_to_u32_digits, ENDO_COEFF, N11, N12, N21, N22};
use crate::{MODULUS_N, MODULUS_Q};
use ark_ff::{BigInteger, Field, LegendreSymbol, PrimeField, Zero};
use ark_secp256k1::{Fq, Fr};
use num_bigint::{BigInt, BigUint, Sign, ToBigInt};
use num_traits::sign::Signed;
use secp256k10_guest::{ComputeHint, Hint};
use std::ops::{MulAssign, Neg};
use std::str::FromStr;

pub struct HintBuilder {}

// The implementation here follows the code in
// https://github.com/rust-bitcoin/rust-secp256k1/blob/master/secp256k1-sys/depend/secp256k1/src/modules/recovery/main_impl.h#L87
// https://github.com/RustCrypto/signatures/blob/master/ecdsa/src/recovery.rs
// which is under the MIT license
impl HintBuilder {
    pub fn build(r: &[u8], s: &[u8], z: &[u8], recid: u8) -> (Hint, Option<ComputeHint>) {
        // Step 0: prepare data
        let n = MODULUS_N.get_or_init(|| {
            BigUint::from_str(
                "115792089237316195423570985008687907852837564279074904382605163141518161494337",
            )
            .unwrap()
        });
        let q = MODULUS_Q.get_or_init(|| {
            BigUint::from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671663",
            )
            .unwrap()
        });
        let n11 =
            N11.get_or_init(|| BigInt::from_str("64502973549206556628585045361533709077").unwrap());
        let n12 = N12
            .get_or_init(|| BigInt::from_str("367917413016453100223835821029139468248").unwrap());
        let n21 = N21
            .get_or_init(|| BigInt::from_str("303414439467246543595250775667605759171").unwrap());
        let n22 =
            N22.get_or_init(|| BigInt::from_str("64502973549206556628585045361533709077").unwrap());

        // Step 1: check if r is zero
        let r = BigUint::from_bytes_le(r);
        if r.is_zero() {
            return (Hint::FormatError, None);
        }

        // Step 2: check if s is zero
        let s = BigUint::from_bytes_le(s);
        if s.is_zero() {
            return (Hint::FormatError, None);
        }

        // Step 3: check if r and s are smaller than n
        if r.ge(&n) {
            return (Hint::FormatError, None);
        }

        if s.ge(&n) {
            return (Hint::FormatError, None);
        }

        // Step 4: check if the recid is valid
        if !(recid < 4) {
            return (Hint::FormatError, None);
        }

        // Step 5: compute r_mod_q, which would be added by n if x is reduced
        let mut r_mod_q = r.clone();

        if recid & 2 != 0 {
            // x is reduced
            r_mod_q = r_mod_q + n;
            if r_mod_q.ge(&q) {
                return (Hint::FormatError, None);
            }
        }

        // Step 6: check if r_mod_q is a valid x coordinate, and if so, compute and obtain y
        let r_x = <Fq as From<BigUint>>::from(r_mod_q);

        let x3_plus_ax_plus_b = r_x.square() * r_x + Fq::from(7u8);
        let legendre = x3_plus_ax_plus_b.legendre();
        if legendre == LegendreSymbol::QuadraticNonResidue {
            // One QNR is 3
            let mut witness = x3_plus_ax_plus_b.double() + x3_plus_ax_plus_b;
            witness = witness.sqrt().unwrap();
            return (
                Hint::YIsImaginary(bytes_to_u32_digits(&witness.into_bigint().to_bytes_le())),
                None,
            );
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

        // Step 7: compute u1 and u2

        // Note: r_inv is over Fr
        let r_inv = Fr::from(r).inverse().unwrap();
        let z = Fr::from_le_bytes_mod_order(z);
        let s = Fr::from(s);

        let u1 = -(&r_inv * &z);
        let u2 = r_inv * s;

        // Step 8: compute sum_u1_k1 and sum_u1_k2
        let ((is_u1_k1_negative, u1_k1_abs), (is_u1_k2_negative, u1_k2_abs)) =
            scalar_decomposition(n11, n12, n21, n22, n, &u1);

        let mut sum_u1_k1 = if is_u1_k1_negative {
            (
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&secp256k10_guest::G_BASE_1.0)),
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&secp256k10_guest::G_BASE_1.2)),
            )
        } else {
            (
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&secp256k10_guest::G_BASE_1.0)),
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&secp256k10_guest::G_BASE_1.1)),
            )
        };
        let mut sum_u1_k2 = if is_u1_k2_negative {
            (
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&secp256k10_guest::G_BASE_2.0)),
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&secp256k10_guest::G_BASE_2.2)),
            )
        } else {
            (
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&secp256k10_guest::G_BASE_2.0)),
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&secp256k10_guest::G_BASE_2.1)),
            )
        };
        let mut hints = Vec::<[u32; 8]>::new();

        for i in 0..4 {
            for j in 0..8 {
                let bits = (u1_k1_abs[i] >> (j * 4)) & 0xF;
                if bits == 0 {
                    continue;
                }
                let (x2, y2) = secp256k10_guest::G_TABLES[i * 8 + j][(bits - 1) as usize];
                let x2 = Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&x2));
                let y2 = Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&y2));
                sum_u1_k1 =
                    point_add_and_push_hint((&sum_u1_k1.0, &sum_u1_k1.1), (&x2, &y2), &mut hints);
            }
        }

        if u1_k1_abs[4] == 1 {
            let (x2, y2) = secp256k10_guest::G_LAST_ENTRY;
            let x2 = Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&x2));
            let y2 = Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&y2));
            sum_u1_k1 =
                point_add_and_push_hint((&sum_u1_k1.0, &sum_u1_k1.1), (&x2, &y2), &mut hints);
        }

        if is_u1_k1_negative {
            sum_u1_k1.1.neg_in_place();
        }

        for i in 0..4 {
            for j in 0..8 {
                let bits = (u1_k2_abs[i] >> (j * 4)) & 0xF;
                if bits == 0 {
                    continue;
                }

                let (x2, y2) = secp256k10_guest::G_TABLES[i * 8 + j][(bits - 1) as usize];
                let x2 = Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&x2));
                let y2 = Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&y2));
                sum_u1_k2 =
                    point_add_and_push_hint((&sum_u1_k2.0, &sum_u1_k2.1), (&x2, &y2), &mut hints);
            }
        }

        let endo_coeff = ENDO_COEFF.get_or_init(|| {
            Fq::from_str(
                "60197513588986302554485582024885075108884032450952339817679072026166228089408",
            )
            .unwrap()
        });

        if u1_k2_abs[4] == 1 {
            let (x2, y2) = secp256k10_guest::G_LAST_ENTRY;
            let x2 = Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&x2));
            let y2 = Fq::from_le_bytes_mod_order(&bytemuck::cast_slice(&y2));
            sum_u1_k2 =
                point_add_and_push_hint((&sum_u1_k2.0, &sum_u1_k2.1), (&x2, &y2), &mut hints);
        }

        // Step 9: compute sum_u1

        sum_u1_k2.0.mul_assign(endo_coeff);
        if is_u1_k2_negative {
            sum_u1_k2.1.neg_in_place();
        }

        let sum_u1 = match (sum_u1_k1, sum_u1_k2) {
            ((x1, y1), (x2, y2)) => {
                if x1 == x2 {
                    if y1 == y2 {
                        point_double_add_push_hint((&x1, &y1), &mut hints)
                    } else {
                        unreachable!()
                    }
                } else {
                    point_add_and_push_hint((&x1, &y1), (&x2, &y2), &mut hints)
                }
            }
        };

        // Step 10: build the table for u2

        // The implementation as follows relies on an observation about secp256k1
        // The doubling algorithm would fail if y = 0.
        //
        // Fortunately, since the curve equation is y^2 = x^3 + 7.
        // And -7 is a non-cubic residue.
        //
        // In other words, no reasonable points would have y = 0.

        // Generate the 4-bit table for u2 base and u2-lambda base
        let mut table_u2 = [(Fq::ZERO, Fq::ZERO, Fq::ZERO); 15];

        table_u2[0].0 = r_x.clone();
        table_u2[0].1 = r_y.clone();
        table_u2[0].2 = &r_x * endo_coeff;

        {
            let (x3, y3) = point_double_add_push_hint((&table_u2[0].0, &table_u2[0].1), &mut hints);
            let x3prime = &x3 * &endo_coeff;
            table_u2[1] = (x3, y3, x3prime);
        }

        for i in 2..15 {
            let (x3, y3) = point_add_and_push_hint(
                (&table_u2[i - 1].0, &table_u2[i - 1].1),
                (&r_x, &r_y),
                &mut hints,
            );
            let x3prime = &x3 * &endo_coeff;
            table_u2[i] = (x3, y3, x3prime);
        }

        // Step 11: compute sum_u2

        let ((is_u2_k1_negative, u2_k1_abs), (is_u2_k2_negative, u2_k2_abs)) =
            scalar_decomposition(n11, n12, n21, n22, n, &u2);

        let mut sum_u2: Option<(Fq, Fq)> = None;

        if u2_k1_abs[4] == 1 {
            // 15 avoids the need to double 4 times
            let (x3, y3) = point_add_and_push_hint(
                (&table_u2[14].0, &table_u2[14].1),
                (&r_x, &r_y),
                &mut hints,
            );
            if is_u2_k1_negative {
                sum_u2 = Some((x3, y3.neg()));
            } else {
                sum_u2 = Some((x3, y3));
            }
        }

        if u2_k2_abs[4] == 1 {
            // This is due to a bound, in that if u2_k1 is big, then it appears that u2_k2 cannot be big
            assert_eq!(u2_k1_abs[4], 0);

            // 15 avoids the need to double 4 times
            let (x3, y3) = point_add_and_push_hint(
                (&table_u2[14].2, &table_u2[14].1),
                (&table_u2[0].2, &r_y),
                &mut hints,
            );
            if is_u2_k2_negative {
                sum_u2 = Some((x3, y3.neg()));
            } else {
                sum_u2 = Some((x3, y3));
            }
        }

        for i in 0..4 {
            for j in 0..8 {
                if !(i == 0 && j == 0) {
                    if sum_u2.is_some() {
                        let mut cur = sum_u2.unwrap();
                        cur = point_double_add_push_hint((&cur.0, &cur.1), &mut hints);
                        cur = point_double_add_push_hint((&cur.0, &cur.1), &mut hints);
                        cur = point_double_add_push_hint((&cur.0, &cur.1), &mut hints);
                        cur = point_double_add_push_hint((&cur.0, &cur.1), &mut hints);
                        sum_u2 = Some(cur);
                    }
                }

                let bits_k1 = ((u2_k1_abs[3 - i] >> ((7 - j) * 4)) & 0xF) as usize;
                let bits_k2 = ((u2_k2_abs[3 - i] >> ((7 - j) * 4)) & 0xF) as usize;

                if bits_k1 != 0 {
                    if sum_u2.is_some() {
                        let res = {
                            let entry = table_u2[bits_k1 - 1].clone();
                            let cur = sum_u2.as_ref().unwrap();
                            if is_u2_k1_negative {
                                point_sub_and_push_hint(
                                    (&cur.0, &cur.1),
                                    (&entry.0, &entry.1),
                                    &mut hints,
                                )
                            } else {
                                point_add_and_push_hint(
                                    (&cur.0, &cur.1),
                                    (&entry.0, &entry.1),
                                    &mut hints,
                                )
                            }
                        };
                        sum_u2 = Some(res);
                    } else {
                        let entry = table_u2[bits_k1 - 1].clone();
                        if is_u2_k1_negative {
                            sum_u2 = Some((entry.0, entry.1.neg()));
                        } else {
                            sum_u2 = Some((entry.0, entry.1));
                        }
                    }
                }

                if bits_k2 != 0 {
                    if sum_u2.is_some() {
                        let res = {
                            let entry = table_u2[bits_k2 - 1].clone();
                            let cur = sum_u2.as_ref().unwrap();
                            if is_u2_k2_negative {
                                point_sub_and_push_hint(
                                    (&cur.0, &cur.1),
                                    (&entry.2, &entry.1),
                                    &mut hints,
                                )
                            } else {
                                point_add_and_push_hint(
                                    (&cur.0, &cur.1),
                                    (&entry.2, &entry.1),
                                    &mut hints,
                                )
                            }
                        };
                        sum_u2 = Some(res);
                    } else {
                        let entry = table_u2[bits_k2 - 1].clone();
                        if is_u2_k2_negative {
                            sum_u2 = Some((entry.2, entry.1.neg()));
                        } else {
                            sum_u2 = Some((entry.2, entry.1));
                        }
                    }
                }
            }
        }

        let r_y_digits = bytes_to_u32_digits(&r_y.into_bigint().to_bytes_le());
        let r_inv_digits = bytes_to_u32_digits(&r_inv.into_bigint().to_bytes_le());

        // Step 12: compute the final sum

        // Adding u1_sum and u2_sum together would have a corner case if their x coordinates are the same.
        // In that case, if the y coordinate is the same, then we double it
        // Otherwise, return an error
        match (sum_u1, sum_u2) {
            ((u1x, u1y), Some((u2x, u2y))) => {
                if u1x == u2x {
                    if u1y != u2y {
                        // impossible to happen, would imply the ability to sign the point of infinity
                        unreachable!()
                    } else {
                        let _ = point_double_add_push_hint((&u1x, &u1y), &mut hints);
                    }
                } else {
                    let _ = point_add_and_push_hint((&u1x, &u1y), (&u2x, &u2y), &mut hints);
                }
            }
            (_, _) => {
                // only possible when z is zero, which would not happen with non-negligible probability
                unreachable!()
            }
        }

        let compute_hint = ComputeHint {
            r_y: r_y_digits,
            r_inv: r_inv_digits,
            hints,
        };

        return (Hint::Ok, Some(compute_hint));
    }
}

pub fn scalar_decomposition(
    n11: &BigInt,
    n12: &BigInt,
    n21: &BigInt,
    n22: &BigInt,
    n: &BigUint,
    u: &Fr,
) -> ((bool, [u32; 8]), (bool, [u32; 8])) {
    let u_bigint = BigUint::from(u.into_bigint()).to_bigint().unwrap();

    let beta_1: BigInt = (&u_bigint) * n22 / &n.to_bigint().unwrap();
    let beta_2: BigInt = (&u_bigint) * n12 / &n.to_bigint().unwrap();

    let b11: BigInt = (&beta_1) * n11;
    let b12: BigInt = (&beta_2) * n21;
    let b1: BigInt = (&b11) + (&b12);

    let b21: BigInt = (&beta_1) * n12;
    let b22: BigInt = (&beta_2) * n22;
    let b2: BigInt = (&b21) - &b22;

    let k1 = &u_bigint - &b1;
    let k2 = -b2;

    let mut k1_abs = [0u32; 8];
    let mut k2_abs = [0u32; 8];

    {
        let v = k1.abs().to_biguint().unwrap().to_u32_digits();
        for i in 0..v.len() {
            k1_abs[i] = v[i];
        }

        let v = k2.abs().to_biguint().unwrap().to_u32_digits();
        for i in 0..v.len() {
            k2_abs[i] = v[i];
        }
    }

    assert!(k1_abs[4] == 1 || k1_abs[4] == 0);
    assert_eq!(k1_abs[5], 0);
    assert_eq!(k1_abs[6], 0);
    assert_eq!(k1_abs[7], 0);

    assert!(k2_abs[4] == 1 || k2_abs[4] == 0);
    assert_eq!(k2_abs[5], 0);
    assert_eq!(k2_abs[6], 0);
    assert_eq!(k2_abs[7], 0);

    (
        (k1.sign() == Sign::Minus, k1_abs),
        (k2.sign() == Sign::Minus, k2_abs),
    )
}

pub fn point_add_and_push_hint(
    point_1: (&Fq, &Fq),
    point_2: (&Fq, &Fq),
    hints: &mut Vec<[u32; 8]>,
) -> (Fq, Fq) {
    let (x1, y1) = point_1;
    let (x2, y2) = point_2;

    let slope = (y1 - y2) * (x1 - x2).inverse().unwrap();
    hints.push(bytes_to_u32_digits(&slope.into_bigint().to_bytes_le()));

    let x3 = slope.square() - x1 - x2;
    let y3 = slope * &(x1 - &x3) - y1;

    return (x3, y3);
}

pub fn point_sub_and_push_hint(
    point_1: (&Fq, &Fq),
    point_2: (&Fq, &Fq),
    hints: &mut Vec<[u32; 8]>,
) -> (Fq, Fq) {
    let (x1, y1) = point_1;
    let (x2, y2) = point_2;

    let slope = (y1 + y2) * (x1 - x2).inverse().unwrap();
    hints.push(bytes_to_u32_digits(&slope.into_bigint().to_bytes_le()));

    let x3 = slope.square() - x1 - x2;
    let y3 = slope * &(x1 - &x3) - y1;

    return (x3, y3);
}

pub fn point_double_add_push_hint(point: (&Fq, &Fq), hints: &mut Vec<[u32; 8]>) -> (Fq, Fq) {
    let (x1, y1) = point;

    let x1_sqr = x1.square();
    let y1_dbl = y1.double();

    let slope = (x1_sqr.double() + x1_sqr) * y1_dbl.inverse().unwrap();
    hints.push(bytes_to_u32_digits(&slope.into_bigint().to_bytes_le()));

    let x3 = slope.square() - x1 - x1;
    let y3 = slope * &(x1 - &x3) - y1;

    (x3, y3)
}
