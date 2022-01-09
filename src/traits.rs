use num_traits::ToPrimitive;

/// A marker trait to symbol represents a irrational number
pub trait Irrational : ToPrimitive {}

// TODO: impl Num for Irrational

// TODO: impl Irrational for Root<T> wher T: Integer

/// This module contains several predefined irrational math constants
pub mod symbols {

pub struct E {}
pub struct Pi {}

}

pub trait FromSqrt<T> {
    fn from_sqrt(t: T) -> Self;
}

pub trait TryFromSqrt<T> : Sized {
    type Error;
    fn try_from_sqrt(t: T) -> Result<Self, Self::Error>;
}
