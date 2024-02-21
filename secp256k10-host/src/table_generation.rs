use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_secp256k1::{Affine, Fr};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::ops::{Add, Mul, Neg};
use std::str::FromStr;

pub struct TableGeneration {
    pub g_series: [[([u32; 8], [u32; 8]); 15]; 32],
    pub g_last_one: ([u32; 8], [u32; 8]),
    pub g_base_1: ([u32; 8], [u32; 8], [u32; 8]),
    pub g_base_2: ([u32; 8], [u32; 8], [u32; 8]),
}

impl TableGeneration {
    pub fn new() -> Self {
        let mut res = [[([0u32; 8], [0u32; 8]); 15]; 32];

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

            for j in 0..15 {
                res[i][j].0.copy_from_slice(&bytemuck::cast_slice(
                    &table_cur.x().unwrap().into_bigint().to_bytes_le(),
                ));
                res[i][j].1.copy_from_slice(&bytemuck::cast_slice(
                    &table_cur.y().unwrap().into_bigint().to_bytes_le(),
                ));

                if j != 14 {
                    table_cur = table_cur.add(&cur).into_affine();
                }
            }
        }

        cur = cur
            .into_group()
            .double()
            .double()
            .double()
            .double()
            .into_affine();

        let mut g_last_one = ([0u32; 8], [0u32; 8]);
        g_last_one.0.copy_from_slice(&bytemuck::cast_slice(
            &cur.x().unwrap().into_bigint().to_bytes_le(),
        ));
        g_last_one.1.copy_from_slice(&bytemuck::cast_slice(
            &cur.y().unwrap().into_bigint().to_bytes_le(),
        ));

        let mut prng = ChaCha20Rng::seed_from_u64(0u64);
        let base_2 = Affine::rand(&mut prng);
        let base_1 = base_2
            .mul(
                Fr::from_str(
                    "78074008874160198520644763525212887401909906723592317393988542598630163514318",
                )
                .unwrap(),
            )
            .neg()
            .into_affine();

        let mut g_base_1 = ([0u32; 8], [0u32; 8], [0u32; 8]);
        g_base_1.0.copy_from_slice(&bytemuck::cast_slice(
            &base_1.x().unwrap().into_bigint().to_bytes_le(),
        ));
        g_base_1.1.copy_from_slice(&bytemuck::cast_slice(
            &base_1.y().unwrap().into_bigint().to_bytes_le(),
        ));
        g_base_1.2.copy_from_slice(&bytemuck::cast_slice(
            &base_1.y().unwrap().neg().into_bigint().to_bytes_le(),
        ));

        let mut g_base_2 = ([0u32; 8], [0u32; 8], [0u32; 8]);
        g_base_2.0.copy_from_slice(&bytemuck::cast_slice(
            &base_2.x().unwrap().into_bigint().to_bytes_le(),
        ));
        g_base_2.1.copy_from_slice(&bytemuck::cast_slice(
            &base_2.y().unwrap().into_bigint().to_bytes_le(),
        ));
        g_base_2.2.copy_from_slice(&bytemuck::cast_slice(
            &base_2.y().unwrap().neg().into_bigint().to_bytes_le(),
        ));

        Self {
            g_series: res,
            g_last_one,
            g_base_1,
            g_base_2,
        }
    }

    pub fn print(&self) {
        println!("Main table:");

        print!("[");
        for i in 0..32 {
            println!("[");
            for j in 0..15 {
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

        println!();

        println!("Last entry:");
        println!("([");
        for v in self.g_last_one.0.iter() {
            print!("{}u32,", v);
        }
        print!("],");
        print!("[");
        for v in self.g_last_one.1.iter() {
            print!("{}u32,", v);
        }
        println!("])");

        println!();

        println!("Base entry 1:");
        println!("([");
        for v in self.g_base_1.0.iter() {
            print!("{}u32,", v);
        }
        print!("],");
        print!("[");
        for v in self.g_base_1.1.iter() {
            print!("{}u32,", v);
        }
        print!("],");
        print!("[");
        for v in self.g_base_1.2.iter() {
            print!("{}u32,", v);
        }
        println!("])");

        println!("Base entry 2:");
        println!("([");
        for v in self.g_base_2.0.iter() {
            print!("{}u32,", v);
        }
        print!("],");
        print!("[");
        for v in self.g_base_2.1.iter() {
            print!("{}u32,", v);
        }
        print!("],");
        print!("[");
        for v in self.g_base_2.2.iter() {
            print!("{}u32,", v);
        }
        println!("])");
    }
}

#[cfg(test)]
mod test {
    use crate::table_generation::TableGeneration;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ff::{BigInteger, PrimeField};
    use ark_secp256k1::{Affine, Fq, Fr};
    use std::ops::{Mul, Neg};
    use std::str::FromStr;

    #[test]
    fn check_consistency() {
        let hint = TableGeneration::new();
        // hint.print();

        for i in 0..32 {
            for j in 0..15 {
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

        let mut r_bigint = Fr::from(1 as u8).into_bigint();
        r_bigint.muln(128);
        let r = Fr::from_bigint(r_bigint).unwrap();
        let point_r = Affine::generator().mul(&r).into_affine();

        let x = point_r.x;
        let y = point_r.y;

        let x_reconstructed =
            Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_last_one.0));
        let y_reconstructed =
            Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_last_one.1));

        assert_eq!(x, x_reconstructed);
        assert_eq!(y, y_reconstructed);

        let lambda = Fr::from_str(
            "78074008874160198520644763525212887401909906723592317393988542598630163514318",
        )
        .unwrap();

        let base_1_x =
            Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_base_1.0));
        let base_1_y =
            Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_base_1.1));
        let base_1_y_neg =
            Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_base_1.2));

        let base_2_x =
            Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_base_2.0));
        let base_2_y =
            Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_base_2.1));
        let base_2_y_neg =
            Fq::from_le_bytes_mod_order(&bytemuck::cast_slice::<u32, u8>(&hint.g_base_2.2));

        assert_eq!(base_1_y.neg(), base_1_y_neg);
        assert_eq!(base_2_y.neg(), base_2_y_neg);

        let base_1 = Affine::new(base_1_x, base_1_y);
        let base_2 = Affine::new(base_2_x, base_2_y);

        let vanishing = (base_1 + base_2.mul(&lambda).into_affine()).into_affine();
        assert!(vanishing.is_zero())
    }
}
