extern crate serde;

mod table;
pub use table::{G_TABLES, G_LAST_ENTRY};

mod hinter;
pub use hinter::{ComputeHint, ComputeHintProvider, Hint};

mod evaluator;
pub use evaluator::{EvaluationResult, Evaluator};

pub(crate) mod utils;
