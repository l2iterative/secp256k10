use std::ops::{Add, Div, Mul, Rem, Sub};

#[cfg(target_os = "zkvm")]
extern "C" {
    fn sys_bigint(
        result: *mut [u32; 8],
        op: u32,
        x: *const [u32; 8],
        y: *const [u32; 8],
        modulus: *const [u32; 8],
    );
}

#[inline(always)]
pub fn add32_and_overflow(a: u32, b: u32, carry: u32) -> (u32, u32) {
    let v = (a as u64).wrapping_add(b as u64).wrapping_add(carry as u64);
    ((v >> 32) as u32, (v & 0xffffffff) as u32)
}
#[inline]
pub fn add<const I: usize, const J: usize>(accm: &mut [u32; I], new: &[u32; J]) -> u32 {
    let mut carry = 0;
    (carry, accm[0]) = add32_and_overflow(accm[0], new[0], carry);
    for i in 1..I {
        (carry, accm[i]) = add32_and_overflow(accm[i], new[i], carry);
    }
    for i in I..J {
        (carry, accm[i]) = add32_and_overflow(accm[i], carry, 0);
    }
    carry
}

#[inline]
pub fn sub_with_borrow(a: u32, b: u32, carry: u32) -> (u32, u32) {
    let res = ((a as u64).wrapping_add(0x100000000))
        .wrapping_sub(b as u64)
        .wrapping_sub(carry as u64);
    (
        (res & 0xffffffff) as u32,
        1u32.wrapping_sub((res >> 32) as u32),
    )
}

#[inline]
pub fn sub_and_borrow<const I: usize, const J: usize>(accu: &mut [u32; I], new: &[u32; J]) -> u32 {
    let (cur, borrow) = accu[0].overflowing_sub(new[0]);
    accu[0] = cur;

    let mut borrow = borrow as u32;
    for i in 1..J {
        (accu[i], borrow) = sub_with_borrow(accu[i], new[i], borrow);
    }
    for i in J..I {
        (accu[i], borrow) = sub_with_borrow(accu[i], borrow, 0);
    }
    borrow
}

#[cfg(not(target_os = "zkvm"))]
pub fn mul_mod(a: &[u32; 8], b: &[u32; 8], n: &[u32; 8]) -> [u32; 8] {
    let a = num_bigint::BigUint::from_bytes_le(bytemuck::cast_slice::<_, u8>(a));
    let b = num_bigint::BigUint::from_bytes_le(bytemuck::cast_slice::<_, u8>(b));
    let n = num_bigint::BigUint::from_bytes_le(bytemuck::cast_slice::<_, u8>(n));

    let res_digits = (a * b % n).to_u32_digits();

    let mut res = [0u32; 8];
    for (i, digit) in res_digits.iter().enumerate() {
        res[i] = *digit;
    }
    res
}

#[cfg(target_os = "zkvm")]
#[inline(always)]
pub fn mul_mod(a: &[u32; 8], b: &[u32; 8], n: &[u32; 8]) -> [u32; 8] {
    let mut res = [0u32; 8];

    unsafe {
        sys_bigint(
            &mut res as *mut [u32; 8],
            0u32,
            a as *const [u32; 8],
            b as *const [u32; 8],
            n as *const [u32; 8],
        );
    }

    return res;
}

#[cfg(not(target_os = "zkvm"))]
pub fn mul_quotient(a: &[u32; 8], b: &[u32; 8], n: &[u32; 8], n_minus_one: &[u32; 8]) -> [u32; 8] {
    let a = num_bigint::BigUint::from_bytes_le(bytemuck::cast_slice::<_, u8>(a));
    let b = num_bigint::BigUint::from_bytes_le(bytemuck::cast_slice::<_, u8>(b));
    let n = num_bigint::BigUint::from_bytes_le(bytemuck::cast_slice::<_, u8>(n));
    let n_minus_one = num_bigint::BigUint::from_bytes_le(bytemuck::cast_slice::<_, u8>(n_minus_one));

    let r = (&a).mul(&b).rem(&n);
    let r_prime = (&a).mul(&b).rem(&n_minus_one);

    let q;
    if r > r_prime {
        q = r.sub(&r_prime).add(1u8).mul(&n_minus_one).rem(&n);
    } else {
        q = r_prime.sub(&r);
    }

    let mut res = [0u32; 8];
    for (i, digit) in q.to_u32_digits().iter().enumerate() {
        res[i] = *digit;
    }
    res
}

#[cfg(target_os = "zkvm")]
pub fn mul_quotient(a: &[u32; 8], b: &[u32; 8], n: &[u32; 8], n_minus_one: &[u32; 8]) -> [u32; 8] {
    let mut r = [0u32; 8];
    let mut r_prime = [0u32; 8];

    unsafe {
        sys_bigint(
            &mut r as *mut [u32; 8],
            0u32,
            a as *const [u32; 8],
            b as *const [u32; 8],
            n as *const [u32; 8],
        );
    }

    unsafe {
        sys_bigint(
            &mut r_prime as *mut [u32; 8],
            0u32,
            a as *const [u32; 8],
            b as *const [u32; 8],
            n_minus_one as *const [u32; 8],
        );
    }

    let r_greater_than_r_prime = false;
    for i in 0..8 {
        if r[7 - i] > r_prime[7 - i] {
            r_greater_than_r_prime = true;
            break;
        }
    }

    if r_greater_than_r_prime {
        // compute r - r_prime + 1
        let borrow = sub_and_borrow(&mut r, &r_prime);
        assert!(borrow == 0);

        let carry = add::<8, 1>(&mut r, &[1u32]);
        assert!(carry == 0);

        unsafe {
            sys_bigint(
                &mut r_prime as *mut [u32; 8],
                0u32,
                r as *const [u32; 8],
                n_minus_one as *const [u32; 8],
                n as *const [u32; 8],
            );
        }

        return r_prime;
    } else {
        let borrow = sub_and_borrow(&mut r_prime, &r);
        assert!(borrow == 0);

        return r;
    }
}


#[inline]
pub fn bytes_to_u32_digits(fe: &[u8]) -> [u32; 8] {
    let mut bytes = [0u8; 32];
    bytes.copy_from_slice(fe);
    bytemuck::cast::<[u8; 32], [u32; 8]>(bytes)
}

#[cfg(test)]
mod test{
    use crate::utils::mul_quotient;

    #[test]
    fn test_mul_quotient() {
        let a = [
            0x5127940du32,
            0x3d06144fu32,
            0xd7f7043du32,
            0xe4dac1e9u32,
            0x9f36e134u32,
            0x4244884cu32,
            0x760ae8bcu32,
            0x0c4c9f9au32,
        ];

        let b = [
            0x6a63a921u32,
            0x803913d7u32,
            0x7623a1efu32,
            0x57f2ae75u32,
            0xa07d6d0cu32,
            0x85bd3eb9u32,
            0x047026cfu32,
            0x7590db0du32,
        ];

        let n = [
            0xf5ab2dadu32,
            0x8999e588u32,
            0x37587cbfu32,
            0x6dbb942cu32,
            0xc1f649a9u32,
            0x9617bfaau32,
            0xd6e25d26u32,
            0xec194960u32,
        ];

        let n_minus_one = [
            0xf5ab2dacu32,
            0x8999e588u32,
            0x37587cbfu32,
            0x6dbb942cu32,
            0xc1f649a9u32,
            0x9617bfaau32,
            0xd6e25d26u32,
            0xec194960u32,
        ];

        let res = mul_quotient(&a, &b, &n, &n_minus_one);
        assert_eq!(res,
            [
                0xd13515aau32,
                0x41757786u32,
                0xb3f1c4bau32,
                0x35e4f06fu32,
                0xa7a3e3e1u32,
                0x388466a2u32,
                0xd93d7f96u32,
                0x061fdcf6u32,
            ]
        );
    }
}