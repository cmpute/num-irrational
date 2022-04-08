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
// TODO: support async version of InfiniteContinuedFraction

mod block;
mod general;
mod simple;
mod infinite;

pub use general::*;
pub use simple::*;
pub use infinite::*;
