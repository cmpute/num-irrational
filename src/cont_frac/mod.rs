//! Data structures and algorithms implementations related to
//! regular and generalized continued fraction
//!
//! There are three abstractions of the continued fraction
//! 1. [ContinuedFraction][simple::ContinuedFraction] represents a simple continued fraction with limited length or periodic
//! 2. [InfiniteContinuedFraction][simple::InfiniteContinuedFraction] represents a infinite simple continued fraction on an iterator
//! 3. [GeneralContinuedFraction][general::GeneralContinuedFraction] is a trait that provides method to operate on a general continued fraction.
//!
//! # References:
//! - <https://pi.math.cornell.edu/~gautam/ContinuedFractions.pdf>
//! - <https://crypto.stanford.edu/pbc/notes/contfrac/>
//! - <http://www.numbertheory.org/continued_fractions.html>
//! - <http://www.numbertheory.org/php/cfrac.html>
//! - <https://github.com/blynn/frac>
//!
// TODO: support the hurwitz complex continued fraction

// TODO (v0.1): selective expose structs and document them
mod block;
pub mod general;
pub mod simple;
