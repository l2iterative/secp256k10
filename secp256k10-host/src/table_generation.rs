use ark_ec::AffineRepr;
use ark_ff::{BigInteger, PrimeField};
use ark_secp256k1::Affine;

pub struct TableGeneration {
    pub g_series: Vec<([u32; 8], [u32; 8])>,
}

impl TableGeneration {
    pub fn new() -> Self {
        let mut res = Vec::<([u32; 8], [u32; 8])>::new();

        let mut cur = Affine::generator();
        for _ in 0..256 {
            let x = cur.x;
            let y = cur.y;

            let mut x_bytes = [0u8; 32];
            let mut y_bytes = [0u8; 32];

            x_bytes.copy_from_slice(&x.into_bigint().to_bytes_le());
            y_bytes.copy_from_slice(&y.into_bigint().to_bytes_le());

            res.push((bytemuck::cast(x_bytes), bytemuck::cast(y_bytes)));

            cur = (cur + &cur).into();
        }

        Self { g_series: res }
    }

    pub fn print(&self) {
        for entry in self.g_series.iter() {
            print!("(");
            print!("[");
            for v in entry.0.iter() {
                print!("{}u32,", v);
            }
            print!("],");
            print!("[");
            for v in entry.1.iter() {
                print!("{}u32,", v);
            }
            print!("]");
            print!("),");
        }
    }
}

#[cfg(test)]
mod test {
    use crate::table_generation::TableGeneration;
    use ark_ec::AffineRepr;
    use ark_ff::PrimeField;
    use ark_secp256k1::{Affine, Fq};
    use std::ops::Neg;

    #[test]
    fn check_consistency() {
        let hint = TableGeneration::new();
        // hint.print();

        let mut cur = Affine::generator();
        for i in 0..256 {
            let x = cur.x;
            let y = cur.y;

            let x_reconstructed =
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_series[i].0));
            let y_reconstructed =
                Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_series[i].1));

            assert_eq!(x, x_reconstructed);
            assert_eq!(y, y_reconstructed);

            cur = (cur + &cur).into();
        }
    }
}
