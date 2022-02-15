//! Data structures and algorithms implementations related to
//! quadratic numbers (roots of quadratic equation).

mod integer;
mod surd;

use integer::QuadraticInt; // TODO (v0.2): expose quadratic int
pub use surd::{QuadraticSurd, QuadraticSurdBase};

pub type Quadratic32 = QuadraticSurd<f32>;
pub type Quadratic64 = QuadraticSurd<f64>;

#[cfg(feature = "num-bigint")]
pub type BigQuadratic = QuadraticSurd<num_bigint::BigInt>;
