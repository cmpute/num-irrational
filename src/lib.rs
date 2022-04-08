//! This crate provides representations of irrational numbers within following categories:
//! - Math constants (`pi`, `e`, etc.)
//!   - Values and continued fraction representations
//! - [Quadratic Numbers](https://en.wikipedia.org/wiki/Algebraic_number#Examples)
//!   - [Quadratic Irrational](https://en.wikipedia.org/wiki/Quadratic_irrational_number): [QuadraticSurd]
//!   - [Quadratic Integer](https://en.wikipedia.org/wiki/Quadratic_integer): [QuadraticInt]
//!   - [Gaussian Integer](https://en.wikipedia.org/wiki/Gaussian_integer): [GaussianInt]
//! - [Continued Fraction](https://en.wikipedia.org/wiki/Continued_fraction)
//!   - [Simple continued fraction](https://en.wikipedia.org/wiki/Continued_fraction): [ContinuedFraction], [InfiniteContinuedFraction]
//!   - [General continued fraction](https://en.wikipedia.org/wiki/Generalized_continued_fraction): [GeneralContinuedFraction]
//!   - Transcendental functions represented in continued fractions
//!
//! It's based on the `num` creates.
//!
//! # Examples
//!
//! ```rust
//! use num_irrational::{FromSqrt, QuadraticSurd};
//! let sq2 = QuadraticSurd::from_sqrt(2i32).unwrap();
//! println!("Square root of 2: {}", sq2); // √2
//!
//! use num_irrational::Computable;
//! let sq2_approx = sq2.approximated(&100).value();
//! println!("Rational approximation with denominator under 100: {}", sq2_approx); // 99/70
//!
//! use core::convert::TryFrom;
//! use num_irrational::ContinuedFraction;
//! // let sq2_fraction = ContinuedFraction::from(sq2); // only if feature `complex` is disabled
//! let sq2_fraction = ContinuedFraction::try_from(sq2).unwrap();
//! println!("Continued Fraction: {}", sq2_fraction); // [1; (2)]
//! ```
//!
//! # Optional Features
//! - `complex`: Enable negative square root base support for [QuadraticSurd]. Note that this flag might
//! change some behavior of the operators on [QuadraticSurd].
//! - `num-complex`: Enable converting [QuadraticSurd] to `num_complex::Complex`. You probably want to enable
//! the `complex` feature at the same time.
//! - `num-bigint`: Enable using big integers as the internal representation.
//!

pub mod cont_frac;
pub mod quadratic;
pub mod symbols;
mod traits;

pub use cont_frac::general::GeneralContinuedFraction;
pub use cont_frac::simple::{ContinuedFraction, InfiniteContinuedFraction};
pub use quadratic::surd::QuadraticSurd;
pub use quadratic::integer::QuadraticInt;
pub use traits::*;

#[cfg(feature = "complex")]
pub use quadratic::integer::GaussianInt;

/// [QuadraticSurd] with 32-bit integers
pub type Quadratic32 = QuadraticSurd<i32>;
/// [QuadraticSurd] with 64-bit integers
pub type Quadratic64 = QuadraticSurd<i64>;

#[cfg(feature = "num-bigint")]
/// [QuadraticSurd] with big integers
pub type BigQuadratic = QuadraticSurd<num_bigint::BigInt>;

// XXX: add support for general algebraic numbers in future
// REF: https://github.com/programmerjake/algebraics/
