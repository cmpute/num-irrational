use num_traits::ToPrimitive;

#[cfg(feature = "num-rational")]
use num_rational::Ratio;

/// A marker trait to symbol represents a irrational number
pub trait Irrational : ToPrimitive {}

// TODO: add Num trait to Irrational
// TODO: impl Irrational for Root<T> wher T: Integer

/// This module contains several predefined irrational math constants
pub mod symbols {

pub struct E {}
pub struct Pi {}

}

pub trait FromSqrt<T> : Sized  {
    type Err;

    fn from_sqrt(t: T) -> Result<Self, Self::Err>;
}

// note: we could also have FloatApproximation,
//       but it's only useful when we have a BigFloat type
// TODO: wait for https://github.com/rust-num/num-rational/issues/100
#[cfg(feature = "num-rational")]
pub trait RationalApproximation {
    /// Return an approximated rational representation of the number
    /// Here `limit` determines the (absolute) maximum of the numerator and denumerator 
    /// First return tells whether the approximation is exact
    fn approx_rational(&self, limit: T) -> (bool, Ratio<T>);
}

// TODO: implement the rational approximation for all irrational types
