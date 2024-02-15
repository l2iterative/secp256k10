use serde::{Deserialize, Serialize};
use std::array::TryFromSliceError;

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ComputeHint {
    pub r_y: [u32; 8],
    pub r_inv: [u32; 8],
    pub hints: Vec<[u32; 8]>,
}

impl ComputeHint {
    pub fn size(&self) -> usize {
        return 16 + self.hints.len() * 8;
    }

    pub fn to_vec(&self) -> Vec<u32> {
        let mut res = Vec::new();
        res.extend_from_slice(&self.r_y);
        res.extend_from_slice(&self.r_inv);
        for hint in self.hints.iter() {
            res.extend_from_slice(hint);
        }
        res
    }
}

pub struct ComputeHintProvider<'a> {
    store: &'a [u32],
}

impl<'a> ComputeHintProvider<'a> {
    pub fn new(store: &'a [u32]) -> Self {
        Self { store }
    }
    pub fn get_r_y(&self) -> &'a [u32; 8] {
        <&[u32; 8]>::try_from(&self.store[0..8]).unwrap()
    }

    pub fn get_r_inv(&self) -> &'a [u32; 8] {
        <&[u32; 8]>::try_from(&self.store[8..16]).unwrap()
    }
    pub fn get_hints(&self, index: usize) -> Result<&'a [u32; 8], TryFromSliceError> {
        <&[u32; 8]>::try_from(&self.store[16 + index * 8..24 + index * 8])
    }
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Hint {
    FormatError,
    YIsImaginary([u32; 8]),
    Ok,
    RecoveredKeyIsPointOfInfinity,
}
