use num_traits::ToPrimitive;
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

/// In case there are multiple solution for square root,
/// only canonical result will be returned 
pub trait FromSqrt<T> : Sized  {
    type Error;

    fn from_sqrt(t: T) -> Result<Self, Self::Error>;
}

#[derive(PartialEq, Debug)]
pub enum Approximation<T> {
    Approximated(T),
    Exact(T)
}

// note: we could also have FloatApproximation,
//       but it's only useful when we have a BigFloat type
// TODO: wait for https://github.com/rust-num/num-rational/issues/100
pub trait RationalApproximation<T> {
    /// Return an approximated rational representation of the number
    fn approx_rational(&self, limit: T) -> Approximation<Ratio<T>>;
}

// TODO: implement the rational approximation for all irrational types
