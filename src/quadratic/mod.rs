//! Data structures and algorithms implementations related to
//! quadratic numbers (roots of quadratic equation).

use num_integer::Roots;
use num_traits::{NumRef, Signed};

/// A helper trait to define valid type that can be used for quadratic numbers
pub trait QuadraticBase: Roots + NumRef + Clone + Signed {}
impl<T: Roots + NumRef + Clone + Signed> QuadraticBase for T {}
// FIXME: change to trait alias after https://github.com/rust-lang/rust/issues/41517

pub mod integer;
pub mod surd;
use core::ops::{Add, Sub};

/// This trait describes operations on quadratic field with self as the discriminant
pub trait QuadraticOps<Rhs = Self, Discr = Self, Output = Self>:
    Add<Rhs, Output = Output> + Sub<Rhs, Output = Output>
{
    type Scalar;

    /// Multiplication in the quadratic field defined by the discriminant `discr`
    fn mul(self, rhs: Rhs, discr: Discr) -> Output;
    /// Division in the quadratic field defined by the discriminant `discr`
    fn div(self, rhs: Rhs, discr: Discr) -> Output;

    /// Get the conjugate of the quadratic integer.
    ///
    /// The conjugate of a quadratic number `x + y√D` is `x - y√D`
    fn conj(self, discr: Discr) -> Output;

    /// Get the norm of the quadratic integer.
    ///
    /// The norm of a quadratic number `x + y√D` is `x² - Dy²`
    fn norm(self, discr: Discr) -> Self::Scalar;
}

// TODO: create QuadraticAssignOps
