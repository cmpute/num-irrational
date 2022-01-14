use num_traits::{ToPrimitive, Signed, Unsigned};
use num_rational::Ratio;

#[cfg(feature = "num-bigint")]
use num_bigint::{BigInt, BigUint};

/// A marker trait to symbol represents a irrational number
pub trait Irrational : ToPrimitive {}

// TODO: add Num trait to Irrational
// TODO: impl Irrational for QuadraticSurd and ContinuedFraction

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
// TODO: implement limit as the limit of denominator
// TODO: wait for https://github.com/rust-num/num-rational/issues/100
pub trait RationalApproximation<T> {
    /// Return an approximated rational representation of the number
    fn approx_rational(&self, limit: &T) -> Approximation<Ratio<T>>;
}

// TODO: implement the rational approximation for all irrational types

pub trait WithSigned: Unsigned{
    type Signed;
    fn to_signed(self) -> Self::Signed;
}

pub trait WithUnsigned: Signed{
    type Unsigned;
    fn to_unsigned(self) -> Self::Unsigned;
}


macro_rules! impl_primitive_sign {
    ($TSigned:ty, $TUnsigned:ty) => {
        impl WithSigned for $TUnsigned {
            type Signed = $TSigned;
            fn to_signed(self) -> Self::Signed { self as $TSigned }
        }
        impl WithUnsigned for $TSigned {
            type Unsigned = $TUnsigned;
            fn to_unsigned(self) -> Self::Unsigned { self as $TUnsigned }
        }
    };
}
impl_primitive_sign!(i8, u8);
impl_primitive_sign!(i16, u16);
impl_primitive_sign!(i32, u32);
impl_primitive_sign!(i64, u64);

#[cfg(feature = "num-bigint")]
impl WithSigned for BigUint {
    type Signed = BigInt;
    fn to_signed(self) -> Self::Signed { BigUint::from(self) }
}

#[cfg(feature = "num-bigint")]
impl WithUnsigned for BigInt {
    type Unsigned = BigUint;
    fn to_unsigned(self) -> Self::Unsigned { self.data() }
}
