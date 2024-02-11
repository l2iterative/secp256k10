use ark_ec::AffineRepr;
use ark_ff::{BigInteger, PrimeField};
use ark_secp256k1::Affine;

pub struct TableGeneration {
    pub g_series: Vec<([u8; 32], [u8; 32])>,
}

impl TableGeneration {
    pub fn new() -> Self {
        let mut res = Vec::<([u8; 32], [u8; 32])>::new();

        let mut cur = Affine::generator();
        for _ in 0..256 {
            let x = cur.x;
            let y = cur.y;

            let mut x_bytes = [0u8; 32];
            let mut y_bytes = [0u8; 32];

            x_bytes.copy_from_slice(&x.into_bigint().to_bytes_le());
            y_bytes.copy_from_slice(&y.into_bigint().to_bytes_le());

            res.push((x_bytes, y_bytes));

            cur = (cur + &cur).into();
        }

        Self { g_series: res }
    }
}

#[cfg(test)]
mod test {
    use crate::table_generation::TableGeneration;
    use ark_ec::AffineRepr;
    use ark_ff::PrimeField;
    use ark_secp256k1::{Affine, Fq};

    #[test]
    fn check_consistency() {
        let hint = TableGeneration::new();

        let mut cur = Affine::generator();
        for i in 0..256 {
            let x = cur.x;
            let y = cur.y;

            let x_reconstructed = Fq::from_le_bytes_mod_order(&hint.g_series[i].0);
            let y_reconstructed = Fq::from_le_bytes_mod_order(&hint.g_series[i].1);

            assert_eq!(x, x_reconstructed);
            assert_eq!(y, y_reconstructed);

            cur = (cur + &cur).into();
        }
    }
}
