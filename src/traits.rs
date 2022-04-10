use num_rational::Ratio;

#[cfg(feature = "num-bigint")]
use num_bigint::{BigInt, BigUint};

/// Represents a possible error occured calling [FromSqrt]
#[derive(Clone, Debug, PartialEq)]
pub struct FromSqrtError<T> {
    pub data: T,
    pub kind: SqrtErrorKind,
}

/// Types of [FromSqrt] errors
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SqrtErrorKind {
    /// A proper representation will requires a integer type with more capacity
    Overflow,
    /// The square root is a complex number
    Complex,
    /// The square root of the target cannot be represented as a quadratic surd
    Unrepresentable,
}

/// Create a number instance from the square root of another number.
///
/// In case there are multiple solution for square root,
/// only canonical result will be returned
pub trait FromSqrt<T>: Sized {
    fn from_sqrt(t: T) -> Result<Self, FromSqrtError<T>>;
}

/// Represents a number conversion with probably exact value.
#[derive(PartialEq, Debug)]
pub enum Approximation<T> {
    Approximated(T),
    Exact(T),
}

impl<T> Approximation<T> {
    /// Get the computed value regardless of whether it's exact
    pub fn value(self) -> T {
        match self {
            Approximation::Approximated(v) => v,
            Approximation::Exact(v) => v,
        }
    }
}

// XXX: backport into num-traits (https://github.com/rust-num/num-rational/issues/100)
/// This trait represents a real number that is computable.
/// See [Wiki](https://en.wikipedia.org/wiki/Computable_number)
pub trait Computable<T> {
    /// Return an approximated rational representation of the number
    /// The `limit` argument specify the maximum value of denominator. This will
    /// ensures the error of the approximation is less than 1/limit^2.
    fn approximated(&self, limit: &T) -> Approximation<Ratio<T>>;
}

// TODO: implement the rational approximation for all irrational types

/// Represents a number type with corresponding signed version
///
/// Signed number can implement this trait as well, with `type Signed = Self`.
pub trait WithSigned {
    type Signed;
    fn to_signed(self) -> Self::Signed;
}

/// Represents a number type with corresponding unsigned version
///
/// Unsigned number can implement this trait as well, with `type Unsigned = Self`.
pub trait WithUnsigned {
    type Unsigned;
    fn to_unsigned(self) -> Self::Unsigned;
}

macro_rules! impl_primitive_sign {
    ($TSigned:ty, $TUnsigned:ty) => {
        impl WithSigned for $TUnsigned {
            type Signed = $TSigned;
            #[inline]
            fn to_signed(self) -> Self::Signed {
                self as $TSigned
            }
        }
        impl WithSigned for $TSigned {
            type Signed = $TSigned;
            #[inline]
            fn to_signed(self) -> Self {
                self
            }
        }
        impl WithUnsigned for $TSigned {
            type Unsigned = $TUnsigned;
            #[inline]
            fn to_unsigned(self) -> Self::Unsigned {
                self as $TUnsigned
            }
        }
        impl WithUnsigned for $TUnsigned {
            type Unsigned = $TUnsigned;
            #[inline]
            fn to_unsigned(self) -> Self {
                self
            }
        }
    };
}
impl_primitive_sign!(i8, u8);
impl_primitive_sign!(i16, u16);
impl_primitive_sign!(i32, u32);
impl_primitive_sign!(i64, u64);
impl_primitive_sign!(i128, u128);

#[cfg(feature = "num-bigint")]
impl WithSigned for BigUint {
    type Signed = BigInt;
    #[inline]
    fn to_signed(self) -> Self::Signed {
        BigInt::from(self)
    }
}

#[cfg(feature = "num-bigint")]
impl WithUnsigned for BigUint {
    type Unsigned = BigUint;
    #[inline]
    fn to_unsigned(self) -> Self {
        self
    }
}

#[cfg(feature = "num-bigint")]
impl WithUnsigned for BigInt {
    type Unsigned = BigUint;
    #[inline]
    fn to_unsigned(self) -> Self::Unsigned {
        self.to_biguint().unwrap()
    }
}

#[cfg(feature = "num-bigint")]
impl WithSigned for BigInt {
    type Signed = BigInt;
    #[inline]
    fn to_signed(self) -> Self::Signed {
        self
    }
}
