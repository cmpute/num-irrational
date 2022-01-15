//! This module contains the implementations of algorithms related to
//! regular and generalized continued fraction
//! 
//! References:
//! - <https://pi.math.cornell.edu/~gautam/ContinuedFractions.pdf>
//! - <https://crypto.stanford.edu/pbc/notes/contfrac/>
//! - <http://www.numbertheory.org/continued_fractions.html>
//! - <http://www.numbertheory.org/php/cfrac.html>
//! - <https://github.com/blynn/frac>

mod block;
mod simple;
mod general;

pub use simple::{ContinuedFraction, InfiniteContinuedFraction};
pub use general::GeneralContinuedFraction;
