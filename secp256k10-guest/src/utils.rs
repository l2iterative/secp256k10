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

pub fn add32_and_overflow(a: u32, b: u32, carry: u32) -> (u32, u32) {
    let v = (a as u64).wrapping_add(b as u64).wrapping_add(carry as u64);
    ((v >> 32) as u32, (v & 0xffffffff) as u32)
}
#[inline(always)]
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
#[inline]
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

#[inline]
pub fn bytes_to_u32_digits(fe: &[u8]) -> [u32; 8] {
    let mut bytes = [0u8; 32];
    bytes.copy_from_slice(fe);
    bytemuck::cast::<[u8; 32], [u32; 8]>(bytes)
}