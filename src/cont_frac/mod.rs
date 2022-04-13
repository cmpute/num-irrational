//! Data structures and algorithms implementations related to
//! regular and generalized continued fraction
//!
//! There are three abstractions of the continued fraction
//! 1. [ContinuedFraction][ContinuedFraction] represents a simple continued fraction with limited length or periodic
//! 2. [InfiniteContinuedFraction][InfiniteContinuedFraction] represents a infinite simple continued fraction on an iterator
//! 3. [GeneralContinuedFraction][GeneralContinuedFraction] is a trait that provides method to operate on a general continued fraction.
//!
//! # References:
//! - <https://pi.math.cornell.edu/~gautam/ContinuedFractions.pdf>
//! - <https://crypto.stanford.edu/pbc/notes/contfrac/>
//! - <http://www.numbertheory.org/continued_fractions.html>
//! - <http://www.numbertheory.org/php/cfrac.html>
//! - <https://github.com/blynn/frac>
//!
// TODO: support the hurwitz complex continued fraction
// XXX: support async version of InfiniteContinuedFraction and GeneralContinuedFraction, after https://github.com/rust-lang/rust/issues/79024
//      when the underlying iterator of InfiniteContinuedFraction is a async iterator, the return from all the methods should also be a async iterator
//      some other candidates: rayon, dpc-pariter, parallel-stream

mod block;
mod general;
mod infinite;
mod simple;

pub use general::*;
pub use infinite::*;
pub use simple::*;
