//! Data structures and algorithms implementations related to
//! regular and generalized continued fraction
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
mod general;
mod simple;

pub use general::*;
pub use simple::*;
