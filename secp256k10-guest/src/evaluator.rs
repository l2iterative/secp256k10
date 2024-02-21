use crate::utils::{add, bytes_to_u32_digits, mul_mod, mul_quotient, sub_and_borrow};
use crate::{ComputeHintStreamer, Hint};
use l2r0_profiler_guest::*;

pub struct Evaluator {
    pub r: [u32; 8],
    pub s: [u32; 8],
    pub z: [u32; 8],
    pub recid: u8,
    pub hint: Hint,
}

#[derive(Debug)]
pub enum EvaluationError {
    RIsZero,
    SIsZero,
    RGeN,
    SGeN,
    RPlusNGeQ,
    InvalidV,
    YIsImaginary,
    RecoveredKeyIsPointOfInfinity,
    WrongHint,
}

impl Evaluator {
    pub fn new(r: &[u8], s: &[u8], z: &[u8], recid: u8, hint: Hint) -> Self {
        Self {
            r: bytes_to_u32_digits(r),
            s: bytes_to_u32_digits(s),
            z: bytes_to_u32_digits(z),
            recid,
            hint,
        }
    }

    pub fn evaluate<O: ComputeHintStreamer>(
        &self,
        compute_hint: &mut O,
    ) -> Result<[u32; 16], EvaluationError> {
        start_timer!("Step 0: prepare data");

        let n = [
            0xd0364141u32,
            0xbfd25e8cu32,
            0xaf48a03bu32,
            0xbaaedce6u32,
            0xfffffffeu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
        ];

        let q = [
            0xfffffc2fu32,
            0xfffffffeu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
        ];

        let n_minus_one = [
            0xd0364140u32,
            0xbfd25e8cu32,
            0xaf48a03bu32,
            0xbaaedce6u32,
            0xfffffffeu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
        ];

        let q_minus_one = [
            0xfffffc2eu32,
            0xfffffffeu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
        ];

        let q_minus_two = [
            0xfffffc2du32,
            0xfffffffeu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
        ];

        let endo_coeff = [
            0x8e6afa40u32,
            0x3ec693d6u32,
            0xed0a766au32,
            0x630fb68au32,
            0x53cbcb16u32,
            0x919bb861u32,
            0x9a83f8efu32,
            0x851695d4u32,
        ];

        let n11 = [
            0x9284eb15u32,
            0xe86c90e4u32,
            0xa7d46bcdu32,
            0x3086d221u32,
            0u32,
            0u32,
            0u32,
            0u32,
        ];

        let n12 = [
            0x9d44cfd8u32,
            0x57c1108du32,
            0xa8e2f3f6u32,
            0x14ca50f7u32,
            1u32,
            0u32,
            0u32,
            0u32,
        ];

        let n21 = [
            0x0abfe4c3u32,
            0x6f547fa9u32,
            0x010e8828u32,
            0xe4437ed6u32,
            0u32,
            0u32,
            0u32,
            0u32,
        ];

        let n22 = [
            0x9284eb15u32,
            0xe86c90e4u32,
            0xa7d46bcdu32,
            0x3086d221u32,
            0u32,
            0u32,
            0u32,
            0u32,
        ];

        let u256_bound = [
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
            0xffffffffu32,
        ];

        stop_start_timer!("Step 1: check if r is zero");

        let mut r_is_zero = true;
        for i in 0..8 {
            if self.r[i] != 0 {
                r_is_zero = false;
            }
        }

        if r_is_zero {
            if self.hint != Hint::FormatError {
                return Err(EvaluationError::WrongHint);
            }
            return Err(EvaluationError::RIsZero);
        }

        stop_start_timer!("Step 2: check if s is zero");

        let mut s_is_zero = true;
        for i in 0..8 {
            if self.s[i] != 0 {
                s_is_zero = false;
            }
        }

        if s_is_zero {
            if self.hint != Hint::FormatError {
                return Err(EvaluationError::WrongHint);
            }
            return Err(EvaluationError::SIsZero);
        }

        stop_start_timer!("Step 3: check if r and s are smaller than n");

        let mut r_less_than_n = false;
        for i in 0..8 {
            if self.r[7 - i] > n[7 - i] {
                if self.hint != Hint::FormatError {
                    return Err(EvaluationError::WrongHint);
                }
                return Err(EvaluationError::RGeN);
            } else if self.r[7 - i] < n[7 - i] {
                r_less_than_n = true;
                break;
            }
        }

        if !r_less_than_n {
            if self.hint != Hint::FormatError {
                return Err(EvaluationError::WrongHint);
            }
            return Err(EvaluationError::RGeN);
        }

        let mut s_less_than_n = false;
        for i in 0..8 {
            if self.s[7 - i] > n[7 - i] {
                if self.hint != Hint::FormatError {
                    return Err(EvaluationError::WrongHint);
                }
                return Err(EvaluationError::SGeN);
            } else if self.s[7 - i] < n[7 - i] {
                s_less_than_n = true;
                break;
            }
        }

        if !s_less_than_n {
            return Err(EvaluationError::SGeN);
        }

        stop_start_timer!("Step 4: check if the recid is valid");

        if !(self.recid < 4) {
            return Err(EvaluationError::InvalidV);
        }

        stop_start_timer!("Step 5: compute r_mod_q, which would be added by n if x is reduced");
        let mut r_mod_q: [u32; 8] = self.r.clone();

        // x is reduced
        if self.recid & 2 != 0 {
            add::<8, 8>(&mut r_mod_q, &n);

            let mut r_mod_q_less_than_q = false;
            for i in 0..8 {
                if r_mod_q[7 - i] > n[7 - i] {
                    if self.hint != Hint::FormatError {
                        return Err(EvaluationError::WrongHint);
                    }
                    return Err(EvaluationError::RPlusNGeQ);
                } else if r_mod_q[7 - i] < n[7 - i] {
                    r_mod_q_less_than_q = true;
                    break;
                }
            }

            if !r_mod_q_less_than_q {
                if self.hint != Hint::FormatError {
                    return Err(EvaluationError::WrongHint);
                }
                return Err(EvaluationError::RPlusNGeQ);
            }
        }

        // no more format error
        if self.hint == Hint::FormatError {
            return Err(EvaluationError::WrongHint);
        }

        stop_start_timer!(
            "Step 6: check if r_mod_q is a valid x coordinate, and if so, compute and obtain y"
        );

        let x_sqr = mul_mod(&self.r, &self.r, &q);
        let x_cubic = mul_mod(&x_sqr, &self.r, &q);

        let mut x3_plus_ax_plus_b = [7u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];
        add::<8, 8>(&mut x3_plus_ax_plus_b, &x_cubic);

        let one = [1u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];
        x3_plus_ax_plus_b = mul_mod(&x3_plus_ax_plus_b, &one, &q);

        let two = [2u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];
        let three = [3u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];

        if let Hint::YIsImaginary(w) = &self.hint {
            let w_sqr = mul_mod(w, w, &q);

            let w_sqr_three = mul_mod(&w_sqr, &three, &q);
            for i in 0..8 {
                if w_sqr_three[i] != x3_plus_ax_plus_b[i] {
                    return Err(EvaluationError::WrongHint);
                }
            }

            return Err(EvaluationError::YIsImaginary);
        }

        let mut r_y = compute_hint.next();
        let r_y_sqr = mul_mod(&r_y, &r_y, &q);
        for i in 0..8 {
            if r_y_sqr[i] != x3_plus_ax_plus_b[i] {
                return Err(EvaluationError::WrongHint);
            }
        }

        if self.recid & 1 == 1 {
            // y is odd
            if r_y[0] & 1 == 0 {
                r_y = mul_mod(&r_y, &q_minus_one, &q);
            }
        } else {
            // y is even
            if r_y[0] & 1 == 1 {
                r_y = mul_mod(&r_y, &q_minus_one, &q);
            }
        }

        stop_start_timer!("Step 7: compute u1 and u2");

        let r_inv = compute_hint.next();

        let should_be_one = mul_mod(&r_inv, &r_mod_q, &n);
        for i in 0..8 {
            if should_be_one[i] != one[i] {
                return Err(EvaluationError::WrongHint);
            }
        }

        let mut u1 = mul_mod(&self.z, &r_inv, &n);
        u1 = mul_mod(&u1, &n_minus_one, &n);
        let u2 = mul_mod(&self.s, &r_inv, &n);

        stop_start_timer!("Step 8: compute sum_u1_k1 and sum_u1_k2");
        let ((is_u1_k1_negative, u1_k1_abs), (is_u1_k2_negative, u1_k2_abs)) =
            scalar_decomposition(&n11, &n12, &n21, &n22, &n, &n_minus_one, &u256_bound, &u1);

        let overflow = [0x000003d1u32, 0x1u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];

        let mut sum_u1_k1 = if is_u1_k1_negative {
            (crate::G_BASE_1.0, crate::G_BASE_1.2)
        } else {
            (crate::G_BASE_1.0, crate::G_BASE_1.1)
        };
        let mut sum_u1_k2 = if is_u1_k2_negative {
            (crate::G_BASE_2.0, crate::G_BASE_2.2)
        } else {
            (crate::G_BASE_2.0, crate::G_BASE_2.1)
        };

        // handle the lower 128 bits first and consider the highest bit a special case
        for i in 0..4 {
            for j in 0..8 {
                let bits = (u1_k1_abs[i] >> (j * 4)) & 0xF;
                if bits == 0 {
                    continue;
                }
                let (x1, y1) = &sum_u1_k1;
                let (x2, y2) = crate::G_TABLES[i * 8 + j][(bits - 1) as usize];
                sum_u1_k1 = point_add_and_get_hint(
                    &q,
                    &q_minus_one,
                    &overflow,
                    &one,
                    x1,
                    y1,
                    &x2,
                    &y2,
                    compute_hint,
                )?;
            }
        }

        if u1_k1_abs[4] == 1 {
            let (x2, y2) = crate::G_LAST_ENTRY;

            let (x1, y1) = &sum_u1_k1;
            sum_u1_k1 = point_add_and_get_hint(
                &q,
                &q_minus_one,
                &overflow,
                &one,
                x1,
                y1,
                &x2,
                &y2,
                compute_hint,
            )?;
        }

        if is_u1_k1_negative {
            sum_u1_k1.1 = mul_mod(&sum_u1_k1.1, &q_minus_one, &q);
        }

        for i in 0..4 {
            for j in 0..8 {
                let bits = (u1_k2_abs[i] >> (j * 4)) & 0xF;
                if bits == 0 {
                    continue;
                }

                let (x2, y2) = crate::G_TABLES[i * 8 + j][(bits - 1) as usize];
                let (x1, y1) = &sum_u1_k2;
                sum_u1_k2 = point_add_and_get_hint(
                    &q,
                    &q_minus_one,
                    &overflow,
                    &one,
                    x1,
                    y1,
                    &x2,
                    &y2,
                    compute_hint,
                )?;
            }
        }

        if u1_k2_abs[4] == 1 {
            let (x2, y2) = crate::G_LAST_ENTRY;
            let (x1, y1) = &sum_u1_k2;
            sum_u1_k2 = point_add_and_get_hint(
                &q,
                &q_minus_one,
                &overflow,
                &one,
                x1,
                y1,
                &x2,
                &y2,
                compute_hint,
            )?;
        }

        sum_u1_k2.0 = mul_mod(&sum_u1_k2.0, &endo_coeff, &q);
        if is_u1_k2_negative {
            sum_u1_k2.1 = mul_mod(&sum_u1_k2.1, &q_minus_one, &q);
        }

        stop_start_timer!("Step 9: compute sum_u1");

        let sum_u1 = {
            let mut x_is_the_same = true;
            for i in 0..8 {
                if sum_u1_k1.0[i] != sum_u1_k2.0[i] {
                    x_is_the_same = false;
                    break;
                }
            }

            if x_is_the_same {
                let mut y_is_the_same = true;
                for i in 0..8 {
                    if sum_u1_k1.1[i] != sum_u1_k2.1[i] {
                        y_is_the_same = false;
                        break;
                    }
                }

                // this would not be very likely to happen
                if x_is_the_same && y_is_the_same {
                    point_double_and_get_hint(
                        &q,
                        &q_minus_one,
                        &q_minus_two,
                        &overflow,
                        &one,
                        &two,
                        &three,
                        &sum_u1_k1.0,
                        &sum_u1_k1.1,
                        compute_hint,
                    )?
                } else {
                    // only possible when z is zero, which would not happen with non-negligible probability
                    unreachable!()
                }
            } else {
                point_add_and_get_hint(
                    &q,
                    &q_minus_one,
                    &overflow,
                    &one,
                    &sum_u1_k1.0,
                    &sum_u1_k1.1,
                    &sum_u1_k2.0,
                    &sum_u1_k2.1,
                    compute_hint,
                )?
            }
        };

        stop_start_timer!("Step 10: build the table for u2");
        let mut table_u2 = [([0u32; 8], [0u32; 8], [0u32; 8]); 15];

        table_u2[0].0 = r_mod_q;
        table_u2[0].1 = r_y;
        table_u2[0].2 = mul_mod(&r_mod_q, &endo_coeff, &q);

        {
            let (x3, y3) = point_double_and_get_hint(
                &q,
                &q_minus_one,
                &q_minus_two,
                &overflow,
                &one,
                &two,
                &three,
                &table_u2[0].0,
                &table_u2[0].1,
                compute_hint,
            )?;
            let x3prime = mul_mod(&x3, &endo_coeff, &q);
            table_u2[1] = (x3, y3, x3prime);
        }

        for i in 2..15 {
            let (x3, y3) = point_add_and_get_hint(
                &q,
                &q_minus_one,
                &overflow,
                &one,
                &table_u2[i - 1].0,
                &table_u2[i - 1].1,
                &r_mod_q,
                &r_y,
                compute_hint,
            )?;
            let x3prime = mul_mod(&x3, &endo_coeff, &q);
            table_u2[i] = (x3, y3, x3prime);
        }

        stop_start_timer!("Step 11: compute sum_u2");
        let ((is_u2_k1_negative, u2_k1_abs), (is_u2_k2_negative, u2_k2_abs)) =
            scalar_decomposition(&n11, &n12, &n21, &n22, &n, &n_minus_one, &u256_bound, &u2);

        let mut sum_u2: Option<([u32; 8], [u32; 8])> = None;

        if u2_k1_abs[4] == 1 {
            // 15 avoids the need to double 4 times
            let (x3, y3) = point_add_and_get_hint(
                &q,
                &q_minus_one,
                &overflow,
                &one,
                &table_u2[14].0,
                &table_u2[14].1,
                &r_mod_q,
                &r_y,
                compute_hint,
            )?;
            sum_u2 = Some((x3, y3));
        }

        if u2_k2_abs[4] == 1 {
            // 15 avoids the need to double 4 times
            let (x3, y3) = point_add_and_get_hint(
                &q,
                &q_minus_one,
                &overflow,
                &one,
                &table_u2[14].0,
                &table_u2[14].1,
                &table_u2[0].2,
                &r_y,
                compute_hint,
            )?;
            sum_u2 = Some((x3, y3));
        }

        for i in 0..4 {
            for j in 0..8 {
                if !(i == 0 && j == 0) {
                    if sum_u2.is_some() {
                        let mut cur = sum_u2.unwrap();
                        cur = point_double_and_get_hint(
                            &q,
                            &q_minus_one,
                            &q_minus_two,
                            &overflow,
                            &one,
                            &two,
                            &three,
                            &cur.0,
                            &cur.1,
                            compute_hint,
                        )?;
                        cur = point_double_and_get_hint(
                            &q,
                            &q_minus_one,
                            &q_minus_two,
                            &overflow,
                            &one,
                            &two,
                            &three,
                            &cur.0,
                            &cur.1,
                            compute_hint,
                        )?;
                        cur = point_double_and_get_hint(
                            &q,
                            &q_minus_one,
                            &q_minus_two,
                            &overflow,
                            &one,
                            &two,
                            &three,
                            &cur.0,
                            &cur.1,
                            compute_hint,
                        )?;
                        cur = point_double_and_get_hint(
                            &q,
                            &q_minus_one,
                            &q_minus_two,
                            &overflow,
                            &one,
                            &two,
                            &three,
                            &cur.0,
                            &cur.1,
                            compute_hint,
                        )?;
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
                                point_sub_and_get_hint(
                                    &q,
                                    &q_minus_one,
                                    &overflow,
                                    &one,
                                    &cur.0,
                                    &cur.1,
                                    &entry.0,
                                    &entry.1,
                                    compute_hint,
                                )?
                            } else {
                                point_add_and_get_hint(
                                    &q,
                                    &q_minus_one,
                                    &overflow,
                                    &one,
                                    &cur.0,
                                    &cur.1,
                                    &entry.0,
                                    &entry.1,
                                    compute_hint,
                                )?
                            }
                        };
                        sum_u2 = Some(res);
                    } else {
                        let entry = table_u2[bits_k1 - 1].clone();
                        if is_u2_k1_negative {
                            sum_u2 = Some((entry.0, mul_mod(&entry.1, &q_minus_one, &q)));
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
                                point_sub_and_get_hint(
                                    &q,
                                    &q_minus_one,
                                    &overflow,
                                    &one,
                                    &cur.0,
                                    &cur.1,
                                    &entry.2,
                                    &entry.1,
                                    compute_hint,
                                )?
                            } else {
                                point_add_and_get_hint(
                                    &q,
                                    &q_minus_one,
                                    &overflow,
                                    &one,
                                    &cur.0,
                                    &cur.1,
                                    &entry.2,
                                    &entry.1,
                                    compute_hint,
                                )?
                            }
                        };
                        sum_u2 = Some(res);
                    } else {
                        let entry = table_u2[bits_k2 - 1].clone();
                        if is_u2_k2_negative {
                            sum_u2 = Some((entry.2, mul_mod(&entry.1, &q_minus_one, &q)));
                        } else {
                            sum_u2 = Some((entry.2, entry.1));
                        }
                    }
                }
            }
        }

        stop_start_timer!("Step 12: compute the final sum");

        match (sum_u1, sum_u2) {
            ((u1x, u1y), Some((u2x, u2y))) => {
                if u1x == u2x {
                    if u1y != u2y {
                        return Err(EvaluationError::RecoveredKeyIsPointOfInfinity);
                    } else {
                        let res = point_double_and_get_hint(
                            &q,
                            &q_minus_one,
                            &q_minus_two,
                            &overflow,
                            &one,
                            &two,
                            &three,
                            &u1x,
                            &u1y,
                            compute_hint,
                        )?;
                        let mut pk = [0u32; 16];
                        pk[0..8].copy_from_slice(&res.0);
                        pk[8..16].copy_from_slice(&res.1);
                        stop_timer!();
                        return Ok(pk);
                    }
                } else {
                    let res = point_add_and_get_hint(
                        &q,
                        &q_minus_one,
                        &overflow,
                        &one,
                        &u1x,
                        &u1y,
                        &u2x,
                        &u2y,
                        compute_hint,
                    )?;
                    let mut pk = [0u32; 16];
                    pk[0..8].copy_from_slice(&res.0);
                    pk[8..16].copy_from_slice(&res.1);
                    stop_timer!();
                    return Ok(pk);
                }
            }
            (_, _) => {
                // only possible when z is zero, which would not happen with non-negligible probability
                unreachable!()
            }
        }
    }
}

pub fn scalar_decomposition(
    n11: &[u32; 8],
    n12: &[u32; 8],
    n21: &[u32; 8],
    n22: &[u32; 8],
    n: &[u32; 8],
    n_minus_one: &[u32; 8],
    u256_bound: &[u32; 8],
    u: &[u32; 8],
) -> ((bool, [u32; 8]), (bool, [u32; 8])) {
    let beta_1 = mul_quotient(u, &n22, &n, &n_minus_one);
    let beta_2 = mul_quotient(u, &n12, &n, &n_minus_one);

    let b11 = mul_mod(&beta_1, &n11, u256_bound);
    let b12 = mul_mod(&beta_2, &n21, u256_bound);
    let mut b1 = b11.clone();
    let carry = add(&mut b1, &b12);
    assert_eq!(carry, 0);

    let b21 = mul_mod(&beta_1, &n12, u256_bound);
    let b22 = mul_mod(&beta_2, &n22, u256_bound);
    let mut b22_larger_than_b21 = false;
    for i in 0..8 {
        if b22[7 - i] > b21[7 - i] {
            b22_larger_than_b21 = true;
            break;
        } else if b22[7 - i] < b21[7 - i] {
            break;
        }
    }

    let mut b2;
    let is_b2_negative;
    if b22_larger_than_b21 {
        b2 = b22.clone();
        let borrow = sub_and_borrow(&mut b2, &b21);
        assert_eq!(borrow, 0);
        is_b2_negative = true;
    } else {
        b2 = b21.clone();
        let borrow = sub_and_borrow(&mut b2, &b22);
        assert_eq!(borrow, 0);
        is_b2_negative = false;
    }

    let is_k2_negative = !is_b2_negative;
    let k2_abs = &b2;

    // by ceiling instead of rounding, the error would be larger.
    // |v1| + |v2| can only guarantee that the maximal values would be at most 129 bits.

    let mut k1_abs = [0u32; 8];
    let is_k1_negative;
    let mut u_larger_than_b1 = false;
    for i in 0..8 {
        if u[7 - i] > b1[7 - i] {
            u_larger_than_b1 = true;
            break;
        } else if u[7 - i] < b1[7 - i] {
            break;
        }
    }

    if u_larger_than_b1 {
        k1_abs[0..8].copy_from_slice(&u[0..8]);
        let borrow = sub_and_borrow::<8, 8>(&mut k1_abs, &b1);
        assert_eq!(borrow, 0);
        is_k1_negative = false;
    } else {
        k1_abs[0..8].copy_from_slice(&b1[0..8]);
        let borrow = sub_and_borrow::<8, 8>(&mut k1_abs, &u);
        assert_eq!(borrow, 0);
        is_k1_negative = true;
    }

    ((is_k1_negative, k1_abs), (is_k2_negative, *k2_abs))
}

pub fn point_double_and_get_hint(
    q: &[u32; 8],
    q_minus_one: &[u32; 8],
    q_minus_two: &[u32; 8],
    overflow: &[u32; 8],
    one: &[u32; 8],
    two: &[u32; 8],
    three: &[u32; 8],
    x1: &[u32; 8],
    y1: &[u32; 8],
    compute_hint_provider: &mut impl ComputeHintStreamer,
) -> Result<([u32; 8], [u32; 8]), EvaluationError> {
    let x1_sqr = mul_mod(x1, x1, q);
    let x1_sqr_three = mul_mod(&x1_sqr, three, q);

    let y1_dbl = mul_mod(y1, two, q);

    let slope = compute_hint_provider.next();
    let should_be_x1_sqr_three = mul_mod(&y1_dbl, &slope, &q);
    for i in 0..8 {
        if should_be_x1_sqr_three[i] != x1_sqr_three[i] {
            return Err(EvaluationError::WrongHint);
        }
    }

    let slope_square = mul_mod(&slope, &slope, &q);
    let x1_neg = mul_mod(&x1, &q_minus_one, &q);
    let x1_neg_two = mul_mod(&x1, &q_minus_two, &q);

    let mut x3 = slope_square.clone();
    let carry = add::<8, 8>(&mut x3, &x1_neg_two);
    if carry == 1 {
        add::<8, 8>(&mut x3, &overflow);
    }
    x3 = mul_mod(&x3, &one, &q);

    let mut x3_minus_x1 = x3.clone();
    let carry = add::<8, 8>(&mut x3_minus_x1, &x1_neg);
    if carry == 1 {
        add::<8, 8>(&mut x3_minus_x1, &overflow);
    }

    let x1_minus_x3 = mul_mod(&x3_minus_x1, &q_minus_one, &q);

    let mut y3 = mul_mod(&slope, &x1_minus_x3, &q);
    let y1_neg = mul_mod(&y1, &q_minus_one, &q);
    let carry = add::<8, 8>(&mut y3, &y1_neg);
    if carry == 1 {
        add::<8, 8>(&mut y3, &overflow);
    }
    y3 = mul_mod(&y3, &one, &q);

    Ok((x3, y3))
}

pub fn point_add_and_get_hint(
    q: &[u32; 8],
    q_minus_one: &[u32; 8],
    overflow: &[u32; 8],
    one: &[u32; 8],
    x1: &[u32; 8],
    y1: &[u32; 8],
    x2: &[u32; 8],
    y2: &[u32; 8],
    compute_hint_provider: &mut impl ComputeHintStreamer,
) -> Result<([u32; 8], [u32; 8]), EvaluationError> {
    let x2_neg = mul_mod(&x2, &q_minus_one, &q);
    let y2_neg = mul_mod(&y2, &q_minus_one, &q);

    let mut x1_minus_x2 = x1.clone();
    let carry = add::<8, 8>(&mut x1_minus_x2, &x2_neg);
    if carry == 1 {
        add::<8, 8>(&mut x1_minus_x2, &overflow);
    }

    let slope = compute_hint_provider.next();

    let mut should_be_y1 = mul_mod(&x1_minus_x2, &slope, &q);
    let carry = add::<8, 8>(&mut should_be_y1, &y2);
    if carry == 1 {
        add::<8, 8>(&mut should_be_y1, &overflow);
    }
    should_be_y1 = mul_mod(&should_be_y1, &one, &q);

    for i in 0..8 {
        if should_be_y1[i] != y1[i] {
            return Err(EvaluationError::WrongHint);
        }
    }

    let slope_square = mul_mod(&slope, &slope, &q);
    let x1_neg = mul_mod(&x1, &q_minus_one, &q);

    let mut x3 = slope_square.clone();
    let carry = add::<8, 8>(&mut x3, &x1_neg);
    if carry == 1 {
        add::<8, 8>(&mut x3, &overflow);
    }
    x3 = mul_mod(&x3, &one, &q);

    let carry = add::<8, 8>(&mut x3, &x2_neg);
    if carry == 1 {
        add::<8, 8>(&mut x3, &overflow);
    }
    x3 = mul_mod(&x3, &one, &q);

    let mut x3_minus_x2 = x3.clone();
    let carry = add::<8, 8>(&mut x3_minus_x2, &x2_neg);
    if carry == 1 {
        add::<8, 8>(&mut x3_minus_x2, &overflow);
    }

    let x2_minus_x3 = mul_mod(&x3_minus_x2, &q_minus_one, &q);

    let mut y3 = mul_mod(&x2_minus_x3, &slope, &q);
    let carry = add::<8, 8>(&mut y3, &y2_neg);
    if carry == 1 {
        add::<8, 8>(&mut y3, &overflow);
    }
    y3 = mul_mod(&y3, &one, &q);

    Ok((x3, y3))
}

pub fn point_sub_and_get_hint(
    q: &[u32; 8],
    q_minus_one: &[u32; 8],
    overflow: &[u32; 8],
    one: &[u32; 8],
    x1: &[u32; 8],
    y1: &[u32; 8],
    x2: &[u32; 8],
    y2: &[u32; 8],
    compute_hint_provider: &mut impl ComputeHintStreamer,
) -> Result<([u32; 8], [u32; 8]), EvaluationError> {
    let x2_neg = mul_mod(&x2, &q_minus_one, &q);

    let mut x1_minus_x2 = x1.clone();
    let carry = add::<8, 8>(&mut x1_minus_x2, &x2_neg);
    if carry == 1 {
        add::<8, 8>(&mut x1_minus_x2, &overflow);
    }

    let slope = compute_hint_provider.next();

    let should_be_y1_plus_y2 = mul_mod(&x1_minus_x2, &slope, &q);
    let mut y1_plus_y2 = y1.clone();
    let carry = add::<8, 8>(&mut y1_plus_y2, &y2);
    if carry == 1 {
        add::<8, 8>(&mut y1_plus_y2, &overflow);
    }
    y1_plus_y2 = mul_mod(&y1_plus_y2, &one, &q);

    for i in 0..8 {
        if should_be_y1_plus_y2[i] != y1_plus_y2[i] {
            return Err(EvaluationError::WrongHint);
        }
    }

    let slope_square = mul_mod(&slope, &slope, &q);
    let x1_neg = mul_mod(&x1, &q_minus_one, &q);

    let mut x3 = slope_square.clone();
    let carry = add::<8, 8>(&mut x3, &x1_neg);
    if carry == 1 {
        add::<8, 8>(&mut x3, &overflow);
    }
    x3 = mul_mod(&x3, &one, &q);

    let carry = add::<8, 8>(&mut x3, &x2_neg);
    if carry == 1 {
        add::<8, 8>(&mut x3, &overflow);
    }
    x3 = mul_mod(&x3, &one, &q);

    let mut x3_minus_x2 = x3.clone();
    let carry = add::<8, 8>(&mut x3_minus_x2, &x2_neg);
    if carry == 1 {
        add::<8, 8>(&mut x3_minus_x2, &overflow);
    }

    let x2_minus_x3 = mul_mod(&x3_minus_x2, &q_minus_one, &q);

    let mut y3 = mul_mod(&x2_minus_x3, &slope, &q);
    let carry = add::<8, 8>(&mut y3, &y2);
    if carry == 1 {
        add::<8, 8>(&mut y3, &overflow);
    }
    y3 = mul_mod(&y3, &one, &q);

    Ok((x3, y3))
}
