//! Data structures and algorithms implementations related to
//! quadratic numbers (roots of quadratic equation).

pub mod integer;
pub mod surd;

/// This trait describes operations on quadratic field with self as the discriminant
pub trait QuadraticNum<Rhs = Self, Discr = Self> {
    type Output;
    type Element;

    fn add(self, rhs: Rhs, discr: Discr) -> Self::Output; // TODO: use Add trait
    fn sub(self, rhs: Rhs, discr: Discr) -> Self::Output; //       use Sub trait
    fn mul(self, rhs: Rhs, discr: Discr) -> Self::Output;
    fn div(self, rhs: Rhs, discr: Discr) -> Self::Output;

    /// Get the conjugate of the quadratic integer.
    /// 
    /// The conjugate of a quadratic number `x + y√D` is `x - y√D`
    fn conj(self, discr: Discr) -> Self::Output;

    /// Get the norm of the quadratic integer.
    /// 
    /// The norm of a quadratic number `x + y√D` is `x² - Dy²`
    fn norm(self, discr: Discr) -> Self::Element;
}

// TODO: implement QuadraticNum for QuadraticSurdCoeffs and construct quadraticsurd based on that
