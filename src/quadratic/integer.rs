
/// Quadratic integer `a + b sqrt(r)`.
pub struct QuadraticInt<T> {
    a: T,
    b: T, // zero if the quad int is a conventional integer
    r: T, // zero if the quad int is a conventional integer
}
