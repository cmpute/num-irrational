/// Quadratic integer `a + bω` where `ω` = `sqrt(D)` or `(1+sqrt(D))/2`.
///
/// The different between [QuadraticInt] and [QuadraticSurd][crate::QuadraticSurd] is that the operations for the
/// latter will be in normal fields of real numbers or complex numbers, while the operations
/// for the former will be in the Quadratic Field (specifically in the quadratic integer ring ℤ\[ω\])
///
/// This struct can be used to represent [Gaussian Integers](https://en.wikipedia.org/wiki/Gaussian_integer) and
/// [Einstein Integers](https://en.wikipedia.org/wiki/Eisenstein_integer)
pub struct QuadraticInt<T, const D: i64> {
    a: T,
    b: T, // zero if the quad int is a conventional integer
}
