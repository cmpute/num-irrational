pub mod cont_frac;
pub mod quad_surd;
pub mod symbols;
pub mod traits;

pub use cont_frac::{ContinuedFraction, GeneralContinuedFraction};
pub use quad_surd::QuadraticSurd;
pub use traits::{Approximation, Computable};

pub type Quadratic32 = QuadraticSurd<f32>;
pub type Quadratic64 = QuadraticSurd<f64>;

#[cfg(feature = "num-bigint")]
pub type BigQuadratic = QuadraticSurd<num_bigint::BigInt>;
