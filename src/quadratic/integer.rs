//! Implementation of quadratic integers

/// Quadratic integer `a + bω` where `ω` = `sqrt(d)` or `(1+sqrt(d))/2`.
///
/// The different between [QuadraticInt] and [QuadraticSurd][crate::QuadraticSurd] is that the operations for the
/// latter will be in normal fields of real numbers or complex numbers, while the operations
/// for the former will be in the Quadratic Field (specifically in the quadratic integer ring ℤ\[ω\])
///
/// This struct can be used to represent [Gaussian Integers](https://en.wikipedia.org/wiki/Gaussian_integer) and
/// [Einstein Integers](https://en.wikipedia.org/wiki/Eisenstein_integer)
pub struct QuadraticInt<T> {
    a: T,
    b: T, // zero if the quad int is a conventional integer
    d: T, // zero if the quad int is a conventional integer
}

/// Const generic version of [QuadraticInt], the operations can only be performed between
/// the integers with the same base
pub struct FixedQuadraticInt<T, const D: i32> {
    a: T,
    b: T,
}
pub type GaussianInt<T> = FixedQuadraticInt<T, -1>;
pub type EisensteinInt<T> = FixedQuadraticInt<T, -3>;
