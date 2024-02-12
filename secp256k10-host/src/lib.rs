use num_bigint::BigUint;
use std::sync::OnceLock;

static MODULUS_Q: OnceLock<BigUint> = OnceLock::new();
static MODULUS_N: OnceLock<BigUint> = OnceLock::new();

mod table_generation;
pub use table_generation::TableGeneration;

mod hinter;
pub use hinter::{ComputeHint, Hint};

mod evaluator;
pub use evaluator::{EvaluationError, EvaluationResult, Evaluator};

mod signer;
pub use signer::*;

#[cfg(test)]
mod integration_test;

pub(crate) mod utils;
