use serde::{Deserialize, Serialize};

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ComputeHint {
    pub r_y: [u32; 8],
    pub r_inv: [u32; 8],
    pub hints: Vec<[u32; 8]>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Hint {
    FormatError,
    YIsImaginary([u32; 8]),
    Ok(ComputeHint),
    RecoveredKeyIsPointOfInfinity(ComputeHint),
}
