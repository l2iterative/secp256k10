extern crate serde;

mod table;
pub use table::{G_BASE_1, G_BASE_2, G_LAST_ENTRY, G_TABLES};

mod hinter;
pub use hinter::{ComputeHint, ComputeHintStore, ComputeHintStreamer, Hint};

#[cfg(target_os = "zkvm")]
pub use hinter::ComputeHintBuffer;

mod evaluator;
pub use evaluator::Evaluator;

pub(crate) mod utils;
