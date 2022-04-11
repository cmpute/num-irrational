//! Implementation of quadratic integers

use super::{QuadraticBase, QuadraticOps};
use core::ops::*;
use num_integer::Integer;
use num_traits::{NumRef, One, RefNum, Signed, Zero};

#[inline]
fn four<T: Add<Output = T> + One>() -> T {
    T::one() + T::one() + T::one() + T::one()
}

/// return -1 if v ≡ 0 mod 4, 0 if v ≡ 1 mod 4, 1 if v ≡ 2 or 3 mod 4
#[inline]
fn mod4d2<T: Signed>(v: &T) -> i8
where
    for<'r> &'r T: RefNum<T>,
{
    let m: T = v % four::<T>();
    if m.is_zero() {
        -1
    } else if v.is_negative() {
        let m2: T = (four::<T>() + m) / (T::one() + T::one());
        if m2.is_one() {
            1
        } else {
            0
        }
    } else {
        let m2: T = m / (T::one() + T::one());
        if m2.is_one() {
            1
        } else {
            0
        }
    }
}

// x/y rounded to the nearest integer. If frac(x/y) = 1/2, then round away from zero.
#[inline]
fn div_rem_round<T: Integer + Signed + NumRef>(x: T, y: &T) -> (T, T)
where
    for<'r> &'r T: RefNum<T>,
{
    let two = T::one() + T::one();
    let offset = y.abs() / two * x.signum();
    let (q, r) = (x + &offset).div_rem(y);
    (q, r - offset)
}

/// Underlying representation of a quadratic integer `a + bω` as (a,b), where `ω` is `√D` or `(1+√D)/2`.
///
/// Specifically, `ω=√D` when `D ≡ 2,3 mod 4`, and `ω=(1+√D)/2` when `D ≡ 1 mod 4`. Note that when `ω=(1+√D)/2`,
/// `ω² = (D-1)/4 + w`.
/// 
/// The division of two quadratic ints is defined as: calculate a/b = x + yω in the quadratic field, then round x
/// and y to nearest integer (round half-way cases away from zero, to be consistent with `round()` functions of
/// primitive types).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadraticIntCoeffs<T>(pub T, pub T);

impl<T: Add<Output = T>> Add<T> for QuadraticIntCoeffs<T> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        Self(self.0 + rhs, self.1)
    }
}
impl<T: Add<Output = T>> Add for QuadraticIntCoeffs<T> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl<T: Sub<Output = T>> Sub<T> for QuadraticIntCoeffs<T> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        Self(self.0 - rhs, self.1)
    }
}
impl<T: Sub<Output = T>> Sub for QuadraticIntCoeffs<T> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl<T: Signed + NumRef> Mul<T> for QuadraticIntCoeffs<T> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Self(self.0 * &rhs, self.1 * rhs)
    }
}

impl<T: Integer + Signed + NumRef> Div<T> for QuadraticIntCoeffs<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self(div_rem_round(self.0, &rhs).0, div_rem_round(self.1, &rhs).0)
    }
}
impl<T: Integer + Signed + NumRef + Clone> Rem<T> for QuadraticIntCoeffs<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn rem(self, rhs: T) -> Self::Output {
        Self(div_rem_round(self.0, &rhs).1, div_rem_round(self.1, &rhs).1)
    }
}

impl<T: Integer + Signed + NumRef> QuadraticIntCoeffs<T>
where
    for<'r> &'r T: RefNum<T>,
{
    /// Get the norm of the quadratic integer (taking reference)
    fn norm_ref(&self, discr: &T) -> T {
        if discr.is_zero() {
            return &self.0 * &self.0;
        }

        match mod4d2(discr) {
            0 => {
                // |a+bw| = (a+b/2)^2 - (b/2*sq)^2 = a^2+ab+b^2/4 - D*b^2/4 = a^2 + ab + b^2(1-D)/4
                &self.0 * &self.0
                    + &self.0 * &self.1
                    + &self.1 * &self.1 * ((T::one() - discr) / four::<T>())
            }
            1 => &self.0 * &self.0 - &self.1 * &self.1 * discr,
            _ => unreachable!(),
        }
    }

    /// Get the quotient and the remainder of self / rhs
    fn div_rem(self, rhs: Self, discr: &T) -> (Self, Self) {
        // This division algorithm is compatible with Gaussian and Eisenstein integers,
        // but it's not necessarily a proper Euclidean division algorithm for any
        // quadratic integers.

        let (a, b) = match mod4d2(discr) {
            0 => (
                // (a+bw)/(c+dw) = (a+bw)*conj(c+dw)/|c+dw| = (a+bw)(c+d-dw) / |c+dw|
                // = (ac+ad-adw+bcw+bdw-bdw^2)/|c+dw|
                // = (ac+ad-bd*(D-1)/4)/|c+dw| + (-ad+bc)w/|c+dw|
                &self.0 * &rhs.0 + &self.0 * &rhs.1
                    - &self.1 * &rhs.1 * ((discr - T::one()) / four::<T>()),
                &self.1 * &rhs.0 - &self.0 * &rhs.1,
            ),
            1 => (
                // (a+bw)/(c+dw) = (a+bw)*conj(c+dw)/|c+dw| = (a+bw)(c-dw) / |c+dw|
                // = (ac-bd*D)/|c+dw| + (bc-ad)w/|c+dw|
                &self.0 * &rhs.0 - &self.1 * &rhs.1 * discr,
                &self.1 * &rhs.0 - &self.0 * &rhs.1,
            ),
            _ => unreachable!(),
        };
        let n = rhs.norm_ref(discr);
        let (aq, ar) = div_rem_round(a, &n);
        let (bq, br) = div_rem_round(b, &n);
        let (q, r) = (Self(aq, bq), Self(ar, br));
        (q, QuadraticOps::mul(rhs, r, discr) / n)
    }
}

impl<T: Integer + Signed + NumRef> QuadraticOps<Self, &T, Self> for QuadraticIntCoeffs<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Scalar = T;

    #[inline]
    fn mul(self, rhs: Self, discr: &T) -> Self {
        match mod4d2(discr) {
            0 => Self(
                // (a+bw)*(c+dw) = ac + bdw^2 + (bc+ad)w = ac+bd(D-1)/4 + (bc+ad+bd)w
                &self.0 * &rhs.0 + &self.1 * &rhs.1 * ((discr - T::one()) / four::<T>()),
                self.0 * &rhs.1 + self.1 * (rhs.0 + rhs.1),
            ),
            1 => Self(
                &self.0 * &rhs.0 + &self.1 * &rhs.1 * discr,
                self.0 * rhs.1 + self.1 * rhs.0,
            ),
            _ => unreachable!(),
        }
    }
    #[inline]
    fn div(self, rhs: Self, discr: &T) -> Self {
        self.div_rem(rhs, discr).0
    }

    #[inline]
    fn conj(self, discr: &T) -> Self {
        if discr.is_zero() {
            return self;
        }

        match mod4d2(discr) {
            0 => Self(
                // conj(a+bw) = conj(a+b/2 + b/2*sq) = a+b/2 - b/2*sq = a+b - bw
                self.0 + &self.1,
                -self.1,
            ),
            1 => Self(self.0, -self.1),
            _ => unreachable!(),
        }
    }
    #[inline]
    fn norm(self, discr: &T) -> T {
        self.norm_ref(discr)
    }
}

/// Quadratic integer `a + bω`, where `ω` is `√D` or `(1+√D)/2`, based on [QuadraticIntCoeffs]
///
/// The different between [QuadraticInt] and [QuadraticSurd][crate::QuadraticSurd] is that the operations for the
/// latter will be in normal fields of real numbers or complex numbers, while the operations
/// for the former will be in the Quadratic Field (specifically in the quadratic integer ring ℤ\[ω\])
/// The arithmetic operations can only be performed between the integers with the same base.
#[derive(Debug, Clone, Copy)]
pub struct QuadraticInt<T> {
    coeffs: QuadraticIntCoeffs<T>, // (a, b), b = 0 is ensured when D = 0
    discr: T,                      // D
}

impl<T: Integer + Signed + NumRef> QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    /// Create a quadratic integer `a + bω`, where `ω` is `√r` or `(1+√r)/2`.
    /// Note that r must be not divisible by 4 (to be square free), otherwise the factor 4
    /// will be extracted from r to b.
    pub fn new(a: T, b: T, r: T) -> Self {
        if b.is_zero() || r.is_zero() {
            return Self {
                coeffs: QuadraticIntCoeffs(a, T::zero()),
                discr: T::zero(),
            };
        }
        #[cfg(not(feature = "complex"))]
        if r.is_negative() {
            panic!("Negative root is not supported without the `complex` feature");
        }

        let mut b = b;
        let mut discr = r;
        while (&discr % four::<T>()).is_zero() {
            discr = discr / four::<T>();
            b = b * (T::one() + T::one());
        }
        Self {
            coeffs: QuadraticIntCoeffs(a, b),
            discr,
        }
    }
    #[inline]
    pub fn from_coeffs(coeffs: QuadraticIntCoeffs<T>, r: T) -> Self {
        Self::new(coeffs.0, coeffs.1, r)
    }

    #[inline]
    pub fn conj(self) -> Self {
        Self {
            coeffs: self.coeffs.conj(&self.discr),
            discr: self.discr,
        }
    }
    #[inline]
    pub fn norm(self) -> T {
        self.coeffs.norm(&self.discr)
    }
    #[inline]
    pub fn is_rational(&self) -> bool {
        self.coeffs.1.is_zero()
    }

    /// Get the fundamental unit of the quadratic field ℚ[√d]
    fn unit(d: T) -> Self {
        // REF: http://www.numbertheory.org/gnubc/unit
        // REF: https://people.reed.edu/~jerry/361/lectures/rqunits.pdf
        unimplemented!()
    }
}

impl<T: QuadraticBase> QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    // Try to reduce the root base with a possible factor.
    fn reduce_root_hinted(self, hint: T) -> Self {
        let hint = hint.abs();
        if hint.is_zero() || hint.is_one() || hint.is_negative() {
            return self;
        }

        let (quo, rem) = self.discr.div_rem(&hint);
        if rem.is_zero() {
            // if hint is actually a factor
            let root = hint.clone().sqrt();
            if &root * &root == hint {
                // if hint is a square number, then remove the hint factor
                // note that r ≡ r/a^2 mod 4 (given a is odd)
                return Self::new(self.coeffs.0, self.coeffs.1 * root, quo);
            }
        }

        self
    }
}

impl<T: Neg<Output = T>> Neg for QuadraticInt<T> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            coeffs: QuadraticIntCoeffs(-self.coeffs.0, -self.coeffs.1),
            discr: self.discr,
        }
    }
}

// Reduce root base for binary operands. Return None if the two bases
// cannot be matched. This function assumes the root bases on both sides
// are not zero
#[inline]
#[cfg(not(feature = "complex"))]
fn reduce_bin_op<T: QuadraticBase>(
    lhs: QuadraticInt<T>,
    rhs: QuadraticInt<T>,
) -> Option<(QuadraticInt<T>, QuadraticInt<T>)>
where
    for<'r> &'r T: RefNum<T>,
{
    debug_assert!(!lhs.discr.is_zero());
    debug_assert!(!rhs.discr.is_zero());

    let result = if lhs.discr > rhs.discr {
        let hint = &lhs.discr / &rhs.discr;
        (lhs.reduce_root_hinted(hint), rhs)
    } else if rhs.discr > lhs.discr {
        let hint = &rhs.discr / &lhs.discr;
        (lhs, rhs.reduce_root_hinted(hint))
    } else {
        (lhs, rhs)
    };

    if result.0.discr == result.1.discr {
        Some(result)
    } else {
        None
    }
}

#[inline]
#[cfg(feature = "complex")]
fn reduce_bin_op<T: QuadraticBase>(
    lhs: QuadraticInt<T>,
    rhs: QuadraticInt<T>,
) -> Option<(QuadraticInt<T>, QuadraticInt<T>)>
where
    for<'r> &'r T: RefNum<T>,
{
    if lhs.discr.is_negative() ^ rhs.discr.is_negative() {
        return None;
    }

    let result = if lhs.discr.abs() > rhs.discr.abs() {
        let hint = &lhs.discr / &rhs.discr;
        (lhs.reduce_root_hinted(hint), rhs)
    } else if rhs.discr.abs() > lhs.discr.abs() {
        let hint = &rhs.discr / &lhs.discr;
        (lhs, rhs.reduce_root_hinted(hint))
    } else {
        (lhs, rhs)
    };

    if result.0.discr == result.1.discr {
        Some(result)
    } else {
        None
    }
}

#[inline]
fn reduce_bin_op_unwrap<T: QuadraticBase>(
    lhs: QuadraticInt<T>,
    rhs: QuadraticInt<T>,
) -> (QuadraticInt<T>, QuadraticInt<T>)
where
    for<'r> &'r T: RefNum<T>,
{
    reduce_bin_op(lhs, rhs).expect("two root bases are not compatible!")
}

impl<T: QuadraticBase> PartialEq for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    fn eq(&self, other: &Self) -> bool {
        // shortcuts
        if self.discr == other.discr {
            return self.coeffs == other.coeffs;
        }
        if self.discr.is_zero() || other.discr.is_zero() {
            return false;
        }

        let QuadraticIntCoeffs(la, lb) = &self.coeffs;
        let QuadraticIntCoeffs(ra, rb) = &other.coeffs;
        la == ra && lb * lb * &self.discr == rb * rb * &other.discr
    }
}

impl<T: QuadraticBase> Eq for QuadraticInt<T> where for<'r> &'r T: RefNum<T> {}

impl<T: Add<Output = T>> Add<T> for QuadraticInt<T> {
    type Output = Self;
    fn add(self, rhs: T) -> Self::Output {
        Self {
            coeffs: self.coeffs + rhs,
            discr: self.discr,
        }
    }
}
impl<T: QuadraticBase> Add for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        if rhs.is_rational() {
            return self + rhs.coeffs.0;
        }
        if self.is_rational() {
            return rhs + self.coeffs.0;
        }

        let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
        Self::from_coeffs(lhs.coeffs + rhs.coeffs, rhs.discr)
    }
}

impl<T: Sub<Output = T>> Sub<T> for QuadraticInt<T> {
    type Output = Self;
    fn sub(self, rhs: T) -> Self::Output {
        Self {
            coeffs: self.coeffs - rhs,
            discr: self.discr,
        }
    }
}
impl<T: QuadraticBase> Sub for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        if rhs.is_rational() {
            return self - rhs.coeffs.0;
        }
        if self.is_rational() {
            return -(rhs - self.coeffs.0);
        }

        let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
        Self::from_coeffs(lhs.coeffs - rhs.coeffs, rhs.discr)
    }
}

impl<T: QuadraticBase> Mul<T> for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        if rhs.is_zero() {
            Self::zero()
        } else {
            Self {
                coeffs: self.coeffs * rhs,
                discr: self.discr,
            }
        }
    }
}
impl<T: QuadraticBase> Mul for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        if rhs.is_rational() {
            return self * rhs.coeffs.0;
        }
        if self.is_rational() {
            return rhs * self.coeffs.0;
        }
        // TODO(v0.3): add is_pure check

        let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
        Self::from_coeffs(
            QuadraticOps::mul(lhs.coeffs, rhs.coeffs, &lhs.discr),
            rhs.discr,
        )
    }
}

impl<T: QuadraticBase> Div<T> for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self {
            coeffs: self.coeffs / rhs,
            discr: self.discr,
        }
    }
}
impl<T: QuadraticBase> Div for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        if rhs.is_rational() {
            return self / rhs.coeffs.0;
        }
        // TODO(v0.3): add is_pure check

        if self.is_rational() {
            Self::from_coeffs(
                QuadraticOps::div(self.coeffs, rhs.coeffs, &rhs.discr),
                rhs.discr,
            )
        } else {
            let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
            Self::from_coeffs(
                QuadraticOps::div(lhs.coeffs, rhs.coeffs, &lhs.discr),
                rhs.discr,
            )
        }
    }
}

impl<T: QuadraticBase> Rem<T> for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn rem(self, rhs: T) -> Self::Output {
        Self {
            coeffs: self.coeffs % rhs,
            discr: self.discr,
        }
    }
}
impl<T: QuadraticBase> Rem for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn rem(self, rhs: Self) -> Self::Output {
        if rhs.is_rational() {
            return self / rhs.coeffs.0;
        }
        // TODO(v0.3): add is_pure check

        let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
        Self::from_coeffs(lhs.coeffs.div_rem(rhs.coeffs, &lhs.discr).1, rhs.discr)
    }
}

impl<T: QuadraticBase> Zero for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    fn zero() -> Self {
        Self {
            coeffs: QuadraticIntCoeffs(T::zero(), T::zero()),
            discr: T::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.0.is_zero() && self.coeffs.1.is_zero()
    }
}

impl<T: QuadraticBase> One for QuadraticInt<T>
where
    for<'r> &'r T: RefNum<T>,
{
    fn one() -> Self {
        Self {
            coeffs: QuadraticIntCoeffs(T::one(), T::zero()),
            discr: T::zero(),
        }
    }
}

#[cfg(feature = "complex")]
mod complex {
    use super::*;

    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GaussianInt<T>(QuadraticIntCoeffs<T>);
    impl<T: Eq> Eq for GaussianInt<T> {}

    impl<T> GaussianInt<T> {
        pub const fn new(re: T, im: T) -> Self {
            Self(QuadraticIntCoeffs(re, im))
        }
    }

    impl<T: Signed> GaussianInt<T> {
        pub fn conj(self) -> Self {
            let QuadraticIntCoeffs(a, b) = self.0;
            Self(QuadraticIntCoeffs(a, -b))
        }
    }

    impl<T: Add<Output = T>> GaussianInt<T>
    where
        for<'r> &'r T: RefNum<T>,
    {
        pub fn norm_ref(&self) -> T {
            let QuadraticIntCoeffs(a, b) = &self.0;
            a * a + b * b
        }
        pub fn norm(self) -> T {
            self.norm_ref()
        }
    }

    impl<T: Add<Output = T>> Add for GaussianInt<T>
    where
        for<'r> &'r T: RefNum<T>,
    {
        type Output = Self;
        #[inline]
        fn add(self, rhs: Self) -> Self::Output {
            Self(self.0 + rhs.0)
        }
    }

    impl<T: Sub<Output = T>> Sub for GaussianInt<T>
    where
        for<'r> &'r T: RefNum<T>,
    {
        type Output = Self;
        #[inline]
        fn sub(self, rhs: Self) -> Self::Output {
            Self(self.0 - rhs.0)
        }
    }

    impl<T: Integer + Signed + NumRef> Mul for GaussianInt<T>
    where
        for<'r> &'r T: RefNum<T>,
    {
        type Output = Self;
        #[inline]
        fn mul(self, rhs: Self) -> Self::Output {
            Self(QuadraticOps::mul(self.0, rhs.0, &-T::one()))
        }
    }

    impl<T: Integer + Signed + NumRef> Div for GaussianInt<T>
    where
        for<'r> &'r T: RefNum<T>,
    {
        type Output = Self;
        #[inline]
        fn div(self, rhs: Self) -> Self::Output {
            Self(QuadraticOps::div(self.0, rhs.0, &-T::one()))
        }
    }

    impl<T: Integer + Signed + NumRef> Rem for GaussianInt<T>
    where
        for<'r> &'r T: RefNum<T>,
    {
        type Output = Self;
        #[inline]
        fn rem(self, rhs: Self) -> Self::Output {
            Self(self.0.div_rem(rhs.0, &-T::one()).1)
        }
    }

    // TODO (v0.3.1): implement num_integer::Integer for GaussianInt, especially is_even and is_odd can be derived from
    // https://crates.io/crates/gaussiant
}

#[cfg(feature = "complex")]
pub use complex::GaussianInt;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn div_rem_round_test() {
        const CASES: [(i8, i8, i8, i8); 18] = [
            // x, y, (x/y), (x%y)
            (-4, 4, -1, 0),
            (-3, 4, -1, 1),
            (-2, 4, -1, 2),
            (-1, 4, 0, -1),
            (0, 4, 0, 0),
            (1, 4, 0, 1),
            (2, 4, 1, -2),
            (3, 4, 1, -1),
            (4, 4, 1, 0),
            (-4, -4, 1, 0),
            (-3, -4, 1, 1),
            (-2, -4, 1, 2),
            (-1, -4, 0, -1),
            (0, -4, 0, 0),
            (1, -4, 0, 1),
            (2, -4, -1, -2),
            (3, -4, -1, -1),
            (4, -4, -1, 0),
        ];
        for (x, y, q, r) in CASES {
            assert_eq!(
                div_rem_round(x, &y),
                (q, r),
                "{} / {} != ({}, {})",
                x,
                y,
                q,
                r
            );
        }
    }

    #[test]
    fn test_ops_unit() {
        let zero = QuadraticInt::new(0, 0, 2);
        let e1 = QuadraticInt::new(1, 0, 2);
        let e2 = QuadraticInt::new(0, 1, 2);
        let e3 = QuadraticInt::new(-1, 0, 2);
        let e4 = QuadraticInt::new(0, -1, 2);
        let m12 = QuadraticInt::new(1, 1, 2);
        let m34 = QuadraticInt::new(-1, -1, 2);

        assert_eq!(e1 + zero, e1);
        assert_eq!(e1 - zero, e1);
        assert_eq!(e1 * zero, zero);
        assert_eq!(zero / e1, zero);

        assert_eq!(e1 + e2, m12);
        assert_eq!(e3 + e4, m34);
        assert_eq!(m12 - e1, e2);
        assert_eq!(m34 - e3, e4);
        assert_eq!(e1 * e2, e2);
        assert_eq!(e2 / e1, e2);
        assert_eq!(e2 * e3, e4);
        assert_eq!(e4 / e2, e3);
        assert_eq!(e3 * e4, e2);
        assert_eq!(e2 / e3, e4);
        assert_eq!(e4 * e1, e4);
        assert_eq!(e4 / e1, e4);

        #[cfg(feature = "complex")]
        {
            // unit complex number tests
            let zero = QuadraticInt::new(0, 0, 2);
            let e1 = QuadraticInt::new(1, 0, -1);
            let e2 = QuadraticInt::new(0, 1, -1);
            let e3 = QuadraticInt::new(-1, 0, -1);
            let e4 = QuadraticInt::new(0, -1, -1);
            let m12 = QuadraticInt::new(1, 1, -1);
            let m34 = QuadraticInt::new(-1, -1, -1);

            assert_eq!(e1 + zero, e1);
            assert_eq!(e1 - zero, e1);
            assert_eq!(e1 * zero, zero);
            assert_eq!(zero / e1, zero);

            assert_eq!(e1 + e2, m12);
            assert_eq!(e3 + e4, m34);
            assert_eq!(m12 - e1, e2);
            assert_eq!(m34 - e3, e4);
            assert_eq!(e1 * e2, e2);
            assert_eq!(e2 / e1, e2);
            assert_eq!(e2 * e3, e4);
            assert_eq!(e4 / e2, e3);
            assert_eq!(e3 * e4, e2);
            assert_eq!(e2 / e3, e4);
            assert_eq!(e4 * e1, e4);
            assert_eq!(e4 / e1, e4);

            assert_eq!(e1.norm(), 1);
            assert_eq!(e2.norm(), 1);
            assert_eq!(e3.norm(), 1);
            assert_eq!(e4.norm(), 1);
            assert_eq!(e1.conj(), e1);
            assert_eq!(e2.conj(), e4);
            assert_eq!(e3.conj(), e3);
            assert_eq!(e4.conj(), e2);
        }
    }

    #[test]
    fn test_ops() {
        let a = QuadraticInt::new(1, 3, 2);
        let b = QuadraticInt::new(-2, -1, 2);
        let c = QuadraticInt::new(-1, 2, 2);

        assert_eq!(a + b, QuadraticInt::new(-1, 2, 2));
        assert_eq!(a - b, QuadraticInt::new(3, 4, 2));
        assert_eq!(a * b, QuadraticInt::new(-8, -7, 2));
        assert_eq!(a / b, QuadraticInt::new(2, -3, 2));
        assert_eq!(a % b, QuadraticInt::new(-1, -1, 2));

        assert_eq!(b + c, QuadraticInt::new(-3, 1, 2));
        assert_eq!(b - c, QuadraticInt::new(-1, -3, 2));
        assert_eq!(b * c, QuadraticInt::new(-2, -3, 2));
        assert_eq!(b / c, QuadraticInt::new(-1, -1, 2));
        assert_eq!(b % c, QuadraticInt::new(1, 0, 2));
    }

    #[test]
    fn test_different_bases() {
        let a = QuadraticInt::new(1, 6, 5);
        let b = QuadraticInt::new(1, 3, 20);
        let c = QuadraticInt::new(1, 2, 45);

        assert_eq!(a, b);
        assert_eq!(a, c);
        assert_eq!(a + b, a + a);
        assert_eq!(a + c, a + a);
        assert_eq!(a - b, QuadraticInt::zero());
        assert_eq!(a - c, QuadraticInt::zero());
        assert_eq!(a * b, a * a);
        assert_eq!(a * c, a * a);
        assert_eq!(a / b, QuadraticInt::one());
        assert_eq!(a / c, QuadraticInt::one());
        assert_eq!(a % b, QuadraticInt::zero());
        assert_eq!(a % c, QuadraticInt::zero());
    }

    #[test]
    #[cfg(feature = "complex")]
    fn test_gaussian() {
        // test Gaussian integers
        let q12 = GaussianInt::new(1, 2);
        let qm12 = GaussianInt::new(-1, 2);
        let q1m2 = GaussianInt::new(1, -2);
        let q23 = GaussianInt::new(2, 3);
        let qm23 = GaussianInt::new(-2, 3);
        let q2m3 = GaussianInt::new(2, -3);

        assert_eq!(q12.conj(), q1m2);
        assert_eq!(q23.conj(), q2m3);
        assert_eq!(q12.norm(), 5);
        assert_eq!(q23.norm(), 13);

        for &v1 in &[q12, qm12, q1m2, q23, qm23, q2m3] {
            for &v2 in &[q12, qm12, q1m2, q23, qm23, q2m3] {
                assert_eq!(v1 + v2 - v1, v2, "add/sub test failed between {:?} and {:?}", v1, v2);
                assert_eq!(v1 * v2 / v1, v2, "mul/div test failed between {:?} and {:?}", v1, v2);
            }
        }
    }

    #[test]
    #[cfg(feature = "complex")]
    fn test_eisenstein() {
        // test Eisenstein integers
        let q12 = QuadraticInt::new(1, 2, -3);
        let qm12 = QuadraticInt::new(-1, 2, -3);
        let q3m2 = QuadraticInt::new(3, -2, -3);
        let q23 = QuadraticInt::new(2, 3, -3);
        let qm23 = QuadraticInt::new(-2, 3, -3);
        let q5m3 = QuadraticInt::new(5, -3, -3);

        assert_eq!(q12.conj(), q3m2);
        assert_eq!(q23.conj(), q5m3);
        assert_eq!(q12.norm(), 7);
        assert_eq!(q23.norm(), 19);
        assert_eq!(q12 * q12, QuadraticInt::new(-3, 8, -3));

        for &v1 in &[q12, qm12, q3m2, q23, qm23, q5m3] {
            for &v2 in &[q12, qm12, q3m2, q23, qm23, q5m3] {
                assert_eq!(v1 + v2 - v1, v2, "add/sub test failed between {:?} and {:?}", v1, v2);
                assert_eq!(v1 * v2 / v1, v2, "mul/div test failed between {:?} and {:?}", v1, v2);
            }
        }
    }
}
