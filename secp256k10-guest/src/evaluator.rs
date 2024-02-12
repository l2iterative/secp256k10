use crate::utils::bytes_to_u32_digits;
use crate::Hint;
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

#[derive(Debug)]
pub enum EvaluationResult {
    Ok([u32; 16]),
    Err(EvaluationError),
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

    pub fn evaluate(&self) -> EvaluationResult {
        start_timer!("r is zero");
        let mut r_is_zero = true;
        for i in 0..8 {
            if self.r[i] != 0 {
                r_is_zero = false;
            }
        }

        if r_is_zero {
            if self.hint != Hint::FormatError {
                return EvaluationResult::Err(EvaluationError::WrongHint);
            }
            return EvaluationResult::Err(EvaluationError::RIsZero);
        }

        stop_start_timer!("s is zero");
        let mut s_is_zero = true;
        for i in 0..8 {
            if self.s[i] != 0 {
                s_is_zero = false;
            }
        }

        if s_is_zero {
            if self.hint != Hint::FormatError {
                return EvaluationResult::Err(EvaluationError::WrongHint);
            }
            return EvaluationResult::Err(EvaluationError::SIsZero);
        }

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

        stop_start_timer!("r less than n");

        let mut r_less_than_n = false;
        for i in 0..8 {
            if !r_less_than_n && self.r[7 - i] > n[7 - i] {
                if self.hint != Hint::FormatError {
                    return EvaluationResult::Err(EvaluationError::WrongHint);
                }
                return EvaluationResult::Err(EvaluationError::RGeN);
            } else if self.r[7 - i] < n[7 - i] {
                r_less_than_n = true;
            }
        }

        if !r_less_than_n {
            if self.hint != Hint::FormatError {
                return EvaluationResult::Err(EvaluationError::WrongHint);
            }
            return EvaluationResult::Err(EvaluationError::RGeN);
        }

        stop_start_timer!("s less than n");

        let mut s_less_than_n = false;
        for i in 0..8 {
            if !s_less_than_n && self.s[7 - i] > n[7 - i] {
                if self.hint != Hint::FormatError {
                    return EvaluationResult::Err(EvaluationError::WrongHint);
                }
                return EvaluationResult::Err(EvaluationError::SGeN);
            } else if self.s[7 - i] < n[7 - i] {
                s_less_than_n = true;
            }
        }

        if !s_less_than_n {
            return EvaluationResult::Err(EvaluationError::SGeN);
        }

        stop_start_timer!("recid < 4");

        if !(self.recid < 4) {
            return EvaluationResult::Err(EvaluationError::InvalidV);
        }

        let mut r_mod_q: [u32; 8] = self.r.clone();

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

        stop_start_timer!("if recid & 2 != 0, x is reduced and needs to be recovered");

        // x is reduced
        if self.recid & 2 != 0 {
            crate::utils::add::<8, 8>(&mut r_mod_q, &n);

            let mut r_mod_q_less_than_q = false;
            for i in 0..8 {
                if !r_mod_q_less_than_q && r_mod_q[7 - i] > n[7 - i] {
                    if self.hint != Hint::FormatError {
                        return EvaluationResult::Err(EvaluationError::WrongHint);
                    }
                    return EvaluationResult::Err(EvaluationError::RPlusNGeQ);
                } else if r_mod_q[7 - i] < n[7 - i] {
                    r_mod_q_less_than_q = true;
                }
            }

            if !r_mod_q_less_than_q {
                if self.hint != Hint::FormatError {
                    return EvaluationResult::Err(EvaluationError::WrongHint);
                }
                return EvaluationResult::Err(EvaluationError::RPlusNGeQ);
            }
        }

        // no more format error
        if self.hint == Hint::FormatError {
            return EvaluationResult::Err(EvaluationError::WrongHint);
        }

        stop_start_timer!("compute x3_plus_ax_plus_b");

        let x_sqr = crate::utils::mul_mod(&self.r, &self.r, &q);
        let x_cubic = crate::utils::mul_mod(&x_sqr, &self.r, &q);

        let mut x3_plus_ax_plus_b = [7u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];
        crate::utils::add::<8, 8>(&mut x3_plus_ax_plus_b, &x_cubic);

        let one = [1u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];
        x3_plus_ax_plus_b = crate::utils::mul_mod(&x3_plus_ax_plus_b, &one, &q);

        let two = [2u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];
        let three = [3u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];

        stop_start_timer!("check if y is imaginary");

        if let Hint::YIsImaginary(w) = &self.hint {
            let w_sqr = crate::utils::mul_mod(w, w, &q);

            let w_sqr_three = crate::utils::mul_mod(&w_sqr, &three, &q);
            for i in 0..8 {
                if w_sqr_three[i] != x3_plus_ax_plus_b[i] {
                    return EvaluationResult::Err(EvaluationError::WrongHint);
                }
            }

            return EvaluationResult::Err(EvaluationError::YIsImaginary);
        }

        let compute_hint = match &self.hint {
            Hint::FormatError => {
                unreachable!()
            }
            Hint::YIsImaginary(_) => {
                unreachable!()
            }
            Hint::Ok(h) => h,
            Hint::RecoveredKeyIsPointOfInfinity(h) => h,
        };

        stop_start_timer!("check r_y");

        let mut r_y = compute_hint.r_y;
        let r_y_sqr = crate::utils::mul_mod(&r_y, &r_y, &q);
        for i in 0..8 {
            if r_y_sqr[i] != x3_plus_ax_plus_b[i] {
                return EvaluationResult::Err(EvaluationError::WrongHint);
            }
        }

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

        stop_start_timer!("revert r_y if needed");

        if self.recid & 1 == 1 {
            // y is odd
            if r_y[0] & 1 == 0 {
                r_y = crate::utils::mul_mod(&r_y, &q_minus_one, &q);
            }
        } else {
            // y is even
            if r_y[0] & 1 == 1 {
                r_y = crate::utils::mul_mod(&r_y, &q_minus_one, &q);
            }
        }

        stop_start_timer!("check r_inv");

        let should_be_one = crate::utils::mul_mod(&compute_hint.r_inv, &r_mod_q, &n);
        for i in 0..8 {
            if should_be_one[i] != one[i] {
                return EvaluationResult::Err(EvaluationError::WrongHint);
            }
        }

        stop_start_timer!("compute u1");

        let mut u1 = crate::utils::mul_mod(&self.z, &compute_hint.r_inv, &n);
        u1 = crate::utils::mul_mod(&u1, &n_minus_one, &n);

        let overflow = [0x000003d1u32, 0x1u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];

        let mut iter = compute_hint.hints.as_slice().into_iter();

        let point_double = |x1: &[u32; 8],
                            y1: &[u32; 8],
                            iter: &mut core::slice::Iter<[u32; 8]>|
         -> Result<([u32; 8], [u32; 8]), EvaluationError> {
            let x1_sqr = crate::utils::mul_mod(x1, x1, &q);
            let x1_sqr_three = crate::utils::mul_mod(&x1_sqr, &three, &q);

            let y1_dbl = crate::utils::mul_mod(y1, &two, &q);

            let slope_res = iter.next();
            if slope_res.is_none() {
                return Err(EvaluationError::WrongHint);
            }

            let slope = slope_res.unwrap();
            let should_be_x1_sqr_three = crate::utils::mul_mod(&y1_dbl, &slope, &q);
            for i in 0..8 {
                if should_be_x1_sqr_three[i] != x1_sqr_three[i] {
                    return Err(EvaluationError::WrongHint);
                }
            }

            let slope_square = crate::utils::mul_mod(&slope, &slope, &q);
            let x1_neg = crate::utils::mul_mod(&x1, &q_minus_one, &q);
            let x1_neg_two = crate::utils::mul_mod(&x1, &q_minus_two, &q);

            let mut x3 = slope_square.clone();
            let carry = crate::utils::add::<8, 8>(&mut x3, &x1_neg_two);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut x3, &overflow);
            }
            x3 = crate::utils::mul_mod(&x3, &one, &q);

            let mut x3_minus_x1 = x3.clone();
            let carry = crate::utils::add::<8, 8>(&mut x3_minus_x1, &x1_neg);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut x3_minus_x1, &overflow);
            }

            let x1_minus_x3 = crate::utils::mul_mod(&x3_minus_x1, &q_minus_one, &q);

            let mut y3 = crate::utils::mul_mod(&slope, &x1_minus_x3, &q);
            let y1_neg = crate::utils::mul_mod(&y1, &q_minus_one, &q);
            let carry = crate::utils::add::<8, 8>(&mut y3, &y1_neg);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut y3, &overflow);
            }
            y3 = crate::utils::mul_mod(&y3, &one, &q);

            Ok((x3, y3))
        };

        let point_add = |x1: &[u32; 8],
                         y1: &[u32; 8],
                         x2: &[u32; 8],
                         y2: &[u32; 8],
                         iter: &mut core::slice::Iter<[u32; 8]>|
         -> Result<([u32; 8], [u32; 8]), EvaluationError> {
            let x2_neg = crate::utils::mul_mod(&x2, &q_minus_one, &q);
            let y2_neg = crate::utils::mul_mod(&y2, &q_minus_one, &q);

            let mut x1_minus_x2 = x1.clone();
            let carry = crate::utils::add::<8, 8>(&mut x1_minus_x2, &x2_neg);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut x1_minus_x2, &overflow);
            }

            let slope_res = iter.next();
            if slope_res.is_none() {
                return Err(EvaluationError::WrongHint);
            }

            let slope = slope_res.unwrap();

            let mut should_be_y1 = crate::utils::mul_mod(&x1_minus_x2, &slope, &q);
            let carry = crate::utils::add::<8, 8>(&mut should_be_y1, &y2);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut should_be_y1, &overflow);
            }
            should_be_y1 = crate::utils::mul_mod(&should_be_y1, &one, &q);

            for i in 0..8 {
                if should_be_y1[i] != y1[i] {
                    return Err(EvaluationError::WrongHint);
                }
            }

            let slope_square = crate::utils::mul_mod(&slope, &slope, &q);
            let x1_neg = crate::utils::mul_mod(&x1, &q_minus_one, &q);

            let mut x3 = slope_square.clone();
            let carry = crate::utils::add::<8, 8>(&mut x3, &x1_neg);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut x3, &overflow);
            }
            x3 = crate::utils::mul_mod(&x3, &one, &q);

            let carry = crate::utils::add::<8, 8>(&mut x3, &x2_neg);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut x3, &overflow);
            }
            x3 = crate::utils::mul_mod(&x3, &one, &q);

            let mut x3_minus_x2 = x3.clone();
            let carry = crate::utils::add::<8, 8>(&mut x3_minus_x2, &x2_neg);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut x3_minus_x2, &overflow);
            }

            let x2_minus_x3 = crate::utils::mul_mod(&x3_minus_x2, &q_minus_one, &q);

            let mut y3 = crate::utils::mul_mod(&x2_minus_x3, &slope, &q);
            let carry = crate::utils::add::<8, 8>(&mut y3, &y2_neg);
            if carry == 1 {
                crate::utils::add::<8, 8>(&mut y3, &overflow);
            }
            y3 = crate::utils::mul_mod(&y3, &one, &q);

            Ok((x3, y3))
        };

        stop_start_timer!("compute u1 * G");
        let mut u1_sum = None;
        for i in 0..8 {
            for j in 0..32 {
                let bit = (u1[i] & (1 << j)) != 0;
                if bit {
                    let (x2, y2) = crate::G_TABLES[i * 32 + j];

                    if u1_sum.is_none() {
                        u1_sum = Some((x2, y2));
                    } else {
                        let (x1, y1) = u1_sum.as_ref().unwrap();
                        let res = point_add(x1, y1, &x2, &y2, &mut iter);
                        if res.is_ok() {
                            u1_sum = Some(res.unwrap());
                        } else {
                            return EvaluationResult::Err(res.unwrap_err());
                        }
                    }
                }
            }
        }

        stop_start_timer!("compute u2");

        let u2 = crate::utils::mul_mod(&self.s, &compute_hint.r_inv, &n);

        stop_start_timer!("compute u2 * R");

        let mut u2_sum = None;
        let mut u2_cur = (r_mod_q, r_y);
        for i in 0..8 {
            for j in 0..32 {
                if i != 0 || j != 0 {
                    let (x1, y1) = &u2_cur;
                    let res = point_double(x1, y1, &mut iter);
                    if res.is_ok() {
                        u2_cur = res.unwrap();
                    } else {
                        return EvaluationResult::Err(res.unwrap_err());
                    }
                }

                let bit = (u2[i] & (1 << j)) != 0;
                if bit {
                    let (x2, y2) = &u2_cur;
                    if u2_sum.is_none() {
                        u2_sum = Some((*x2, *y2));
                    } else {
                        let (x1, y1) = u2_sum.unwrap();
                        let res = point_add(&x1, &y1, x2, y2, &mut iter);
                        if res.is_ok() {
                            u2_sum = Some(res.unwrap());
                        } else {
                            return EvaluationResult::Err(res.unwrap_err());
                        }
                    }
                }
            }
        }

        stop_timer!();

        match (u1_sum, u2_sum) {
            (Some((u1x, u1y)), Some((u2x, u2y))) => {
                if u1x == u2x {
                    if u1y != u2y {
                        return EvaluationResult::Err(
                            EvaluationError::RecoveredKeyIsPointOfInfinity,
                        );
                    } else {
                        let res = point_double(&u1x, &u1y, &mut iter);
                        if let Ok(res) = res {
                            let mut pk = [0u32; 16];
                            pk[0..8].copy_from_slice(&res.0);
                            pk[8..16].copy_from_slice(&res.1);
                            return EvaluationResult::Ok(pk);
                        } else {
                            return EvaluationResult::Err(res.unwrap_err());
                        }
                    }
                } else {
                    let res = point_add(&u1x, &u1y, &u2x, &u2y, &mut iter);
                    if let Ok(res) = res {
                        let mut pk = [0u32; 16];
                        pk[0..8].copy_from_slice(&res.0);
                        pk[8..16].copy_from_slice(&res.1);
                        return EvaluationResult::Ok(pk);
                    } else {
                        return EvaluationResult::Err(res.unwrap_err());
                    }
                }
            }
            (_, _) => {
                // only possible when z is zero, which would not happen with non-negligible probability
                unreachable!()
            }
        }
    }
}
