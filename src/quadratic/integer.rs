//! Implementation of quadratic integers

use super::QuadraticNum;
use core::ops::*;
use num_traits::{NumRef, RefNum, Signed, One};

#[inline]
fn four<T: Add<Output=T> + One>() -> T {
    T::one() + T::one() + T::one() + T::one()
}

/// return -1 if v ≡ 0 mod 4, 0 if v ≡ 1 mod 4, 1 if v ≡ 2 or 3 mod 4
#[inline]
fn mod4d2<T: Signed + One>(v: &T) -> i8 where for<'r> &'r T: RefNum<T> {
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

fn div_round<T: Signed + NumRef + One>(x: T, y: &T) -> T
where for<'r> &'r T: RefNum<T>, {
    if x.is_negative() ^ y.is_negative() {
        (x - y / (T::one() + T::one())) / y
    } else {
        (x + y / (T::one() + T::one())) / y
    }
}

/// Underlying representation of a quadratic integer `a + bω` as (a,b), where `ω` is `√D` or `(1+√D)/2`.
/// 
/// Specifically, `ω=√D` when `D ≡ 2,3 mod 4`, and `ω=(1+√D)/2` when `D ≡ 1 mod 4`. Note that when `ω=(1+√D)/2`,
/// `ω² = (D-1)/4 + w`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QuadraticIntCoeffs<T> (pub T, pub T);

// TODO: constraint T to be Signed

impl<T: Signed + NumRef + One + Clone> QuadraticNum<Self, &T> for QuadraticIntCoeffs<T> where for<'r> &'r T: RefNum<T> {
    type Output = Self;
    type Element = T;

    #[inline]
    fn add(self, rhs: Self, discr: &T) -> Self {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
    #[inline]
    fn sub(self, rhs: Self, discr: &T) -> Self {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
    #[inline]
    fn mul(self, rhs: Self, discr: &T) -> Self {
        match mod4d2(discr) {
            0 => Self(
                // (a+bw)*(c+dw) = ac + bdw^2 + (bc+ad)w = ac+bd(D-1)/4 + (bc+ad+bd)w
                &self.0 * &rhs.0 + &self.1 * &rhs.1 * ((discr - T::one())/four::<T>()),
                self.0 * &rhs.1 + self.1 * (rhs.0 + rhs.1)
            ),
            1 => Self(
                &self.0 * &rhs.0 + &self.1 * &rhs.1 * discr,
                self.0 * rhs.1 + self.1 * rhs.0
            ),
            _ => unreachable!()
        }
    }
    #[inline]
    fn div(self, rhs: Self, discr: &T) -> Self {
        let (a, b) = match mod4d2(discr) {
            0 => (
                // (a+bw)/(c+dw) = (a+bw)*conj(c+dw)/|c+dw| = (a+bw)(c+d-dw) / |c+dw|
                // = (ac+ad-adw+bcw+bdw-bdw^2)/|c+dw|
                // = (ac+ad-bd*(D-1)/4)/|c+dw| + (-ad+bc)w/|c+dw|
                &self.0 * &rhs.0 + &self.0 * &rhs.1 - &self.1 * &rhs.1 * ((discr - T::one())/four::<T>()),
                &self.1 * &rhs.0 - &self.0 * &rhs.1
            ),
            1 => (
                // (a+bw)/(c+dw) = (a+bw)*conj(c+dw)/|c+dw| = (a+bw)(c-dw) / |c+dw|
                // = (ac-bd*D)/|c+dw| + (bc-ad)w/|c+dw|
                &self.0 * &rhs.0 - &self.1 * &rhs.1 * discr,
                &self.1 * &rhs.0 - &self.0 * &rhs.1
            ),
            _ => unreachable!()
        };
        let n = Self::norm(rhs, discr);
        Self(div_round(a, &n), div_round(b, &n))
    }
    #[inline]
    fn conj(self, discr: &T) -> Self::Output {
        match mod4d2(discr) {
            0 => Self(
                // conj(a+bw) = conj(a+b/2 + b/2*sq) = a+b/2 - b/2*sq = a+b - bw
                self.0 + &self.1, -self.1
            ),
            1 => Self(self.0, -self.1),
            _ => unreachable!()
        }
    }
    #[inline]
    fn norm(self, discr: &T) -> Self::Element {
        match mod4d2(discr) {
            0 => {
                // |a+bw| = (a+b/2)^2 - (b/2*sq)^2 = a^2+ab+b^2/4 - D*b^2/4 = a^2 + ab + b^2(1-D)/4
                &self.0 * &self.0 + &self.0 * &self.1 + &self.1 * &self.1 * ((T::one() - discr)/four::<T>())
            },
            1 => &self.0 * &self.0 - &self.1 * &self.1 * discr,
            _ => unreachable!()
        }
    }
}

/// Quadratic integer `a + bω`, where `ω` is `√D` or `(1+√D)/2`, based on [QuadraticIntCoeffs]
/// 
/// The different between [QuadraticInt] and [QuadraticSurd][crate::QuadraticSurd] is that the operations for the
/// latter will be in normal fields of real numbers or complex numbers, while the operations
/// for the former will be in the Quadratic Field (specifically in the quadratic integer ring ℤ\[ω\])
/// The arithmetic operations can only be performed between the integers with the same base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QuadraticInt<T> {
    coeffs: QuadraticIntCoeffs<T>,
    discr: T,
}

// TODO: define GaussianInt only when complex feature is enabled.

impl<T: Signed + NumRef + One + Clone> QuadraticInt<T> where for<'r> &'r T: RefNum<T> {
    /// Create a quadratic integer `a + bω`, where `ω` is `√r` or `(1+√r)/2`.
    /// Note that r must be not divisible by 4 (to be square free), otherwise the factor 4
    /// will be extracted from r to b.
    pub fn new(a: T, b: T, r: T) -> Self {
        let mut b = b;
        let mut discr = r;
        while (&discr % four::<T>()).is_zero() {
            discr = discr / four::<T>();
            b = b * (T::one() + T::one());
        }
        Self { coeffs: QuadraticIntCoeffs(a, b), discr }
    }

    #[inline]
    pub fn conj(self) -> Self {
        Self { coeffs: self.coeffs.conj(&self.discr), discr: self.discr }
    }
    #[inline]
    pub fn norm(self) -> T {
        self.coeffs.norm(&self.discr)
    }
}

macro_rules! forward_binop {
    (impl $trait:ident, $func:ident) => {
        impl<T: Signed + NumRef + One + Clone> $trait for QuadraticInt<T> where for<'r> &'r T: RefNum<T> {
            type Output = Self;
            #[inline]
            fn $func(self, rhs: Self) -> Self::Output {
                Self { coeffs: self.coeffs.$func(rhs.coeffs, &self.discr), discr: rhs.discr }
            }
        }
    };
}
forward_binop!(impl Add, add);
forward_binop!(impl Sub, sub);
forward_binop!(impl Mul, mul);
forward_binop!(impl Div, div);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ops_unit() {
        // unit complex number tests
        let e1 = QuadraticInt::new(1, 0, -1);
        let e2 = QuadraticInt::new(0, 1, -1);
        let e3 = QuadraticInt::new(-1, 0, -1);
        let e4 = QuadraticInt::new(0, -1, -1);
        let m12 = QuadraticInt::new(1, 1, -1);
        let m34 = QuadraticInt::new(-1, -1, -1);

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

    #[test]
    fn test_gaussian() {
        // test Gaussian integers
        let q12 = QuadraticInt::new(1, 2, -1);
        let qm12 = QuadraticInt::new(-1, 2, -1);
        let q1m2 = QuadraticInt::new(1, -2, -1);
        let q23 = QuadraticInt::new(2, 3, -1);
        let qm23 = QuadraticInt::new(-2, 3, -1);
        let q2m3 = QuadraticInt::new(2, -3, -1);

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
