use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::{BigInteger, PrimeField};
use ark_secp256k1::Affine;
use std::ops::Add;

pub struct TableGeneration {
    pub g_series: [[([u32; 8], [u32; 8]); 7]; 32],
}

impl TableGeneration {
    pub fn new() -> Self {
        let mut res = [[([0u32; 8], [0u32; 8]); 7]; 32];

        let mut cur = Affine::generator();
        for i in 0..32 {
            if i != 0 {
                cur = cur
                    .into_group()
                    .double()
                    .double()
                    .double()
                    .double()
                    .into_affine();
            }

            let mut table_cur = cur.clone();
            res[i][0].0.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.x().unwrap().into_bigint().to_bytes_le(),
            ));
            res[i][0].1.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.y().unwrap().into_bigint().to_bytes_le(),
            ));

            table_cur = table_cur.add(&cur).into_affine();
            res[i][1].0.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.x().unwrap().into_bigint().to_bytes_le(),
            ));
            res[i][1].1.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.y().unwrap().into_bigint().to_bytes_le(),
            ));

            table_cur = table_cur.add(&cur).into_affine();
            res[i][2].0.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.x().unwrap().into_bigint().to_bytes_le(),
            ));
            res[i][2].1.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.y().unwrap().into_bigint().to_bytes_le(),
            ));

            table_cur = table_cur.add(&cur).into_affine();
            res[i][3].0.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.x().unwrap().into_bigint().to_bytes_le(),
            ));
            res[i][3].1.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.y().unwrap().into_bigint().to_bytes_le(),
            ));

            table_cur = table_cur.add(&cur).into_affine();
            res[i][4].0.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.x().unwrap().into_bigint().to_bytes_le(),
            ));
            res[i][4].1.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.y().unwrap().into_bigint().to_bytes_le(),
            ));

            table_cur = table_cur.add(&cur).into_affine();
            res[i][5].0.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.x().unwrap().into_bigint().to_bytes_le(),
            ));
            res[i][5].1.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.y().unwrap().into_bigint().to_bytes_le(),
            ));

            table_cur = table_cur.add(&cur).into_affine();
            res[i][6].0.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.x().unwrap().into_bigint().to_bytes_le(),
            ));
            res[i][6].1.copy_from_slice(&bytemuck::cast_slice(
                &table_cur.y().unwrap().into_bigint().to_bytes_le(),
            ));
        }

        Self { g_series: res }
    }

    pub fn print(&self) {
        print!("[");
        for i in 0..32 {
            println!("[");
            for j in 0..7 {
                println!("([");
                for v in self.g_series[i][j].0.iter() {
                    print!("{}u32,", v);
                }
                print!("],");
                print!("[");
                for v in self.g_series[i][j].1.iter() {
                    print!("{}u32,", v);
                }
                println!("]),");
            }
            println!("],");
        }
        print!("]");
    }
}

#[cfg(test)]
mod test {
    use crate::table_generation::TableGeneration;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ff::{BigInteger, PrimeField};
    use ark_secp256k1::{Affine, Fq, Fr};
    use std::ops::Mul;

    #[test]
    fn check_consistency() {
        let hint = TableGeneration::new();
        hint.print();

        for i in 0..32 {
            for j in 0..7 {
                let mut r_bigint = Fr::from((j + 1) as u8).into_bigint();
                r_bigint.muln((i as u32) * 4);

                let r = Fr::from_bigint(r_bigint).unwrap();
                let point_r = Affine::generator().mul(&r).into_affine();

                let x = point_r.x;
                let y = point_r.y;

                let x_reconstructed = Fq::from_le_bytes_mod_order(
                    &bytemuck::cast_slice::<u32, u8>(&hint.g_series[i][j].0),
                );
                let y_reconstructed = Fq::from_le_bytes_mod_order(
                    &bytemuck::cast_slice::<u32, u8>(&hint.g_series[i][j].1),
                );

                assert_eq!(x, x_reconstructed);
                assert_eq!(y, y_reconstructed);
            }
        }
    }
}
