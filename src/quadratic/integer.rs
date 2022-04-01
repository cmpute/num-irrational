//! Implementation of quadratic integers

use core::ops::*;
use num_traits::{FromPrimitive, NumRef, RefNum, Zero};

/// Quadratic integer `a + b√D`, where `ω` is `√D` or `(1+√D)/2`
/// 
/// Specifically, `ω=√D` when `D ≡ 2,3 mod 4`, and `ω=(1+√D)/2` when `D ≡ 1 mod 4`. Note that when `ω=(1+√D)/2`,
/// `ω² = (D-1)/4 + w`
/// 
/// The different between [QuadraticInt] and [QuadraticSurd][crate::QuadraticSurd] is that the operations for the
/// latter will be in normal fields of real numbers or complex numbers, while the operations
/// for the former will be in the Quadratic Field (specifically in the quadratic integer ring ℤ\[ω\])
/// The arithmetic operations can only be performed between the integers with the same base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QuadraticInt<T, const D: i32> {
    pub a: T,
    pub b: T,
}
//      https://www.cut-the-knot.org/arithmetic/int_domain.shtml

/// [Gaussian Integer](https://en.wikipedia.org/wiki/Gaussian_integer)
pub type GaussianInt<T> = QuadraticInt<T, -1>;

/// [Eisenstein Integer](https://en.wikipedia.org/wiki/Eisenstein_integer)
/// 
/// Note that the base for Eisenstein Integer here is `(1+√-3)/2` instead of the commonly used `(-1+√-3)/2`
pub type EisensteinInt<T> = QuadraticInt<T, -3>;

// TODO: constraint T to be Signed
impl<T: for<'r> Add<&'r T, Output=T> + Neg<Output=T>, const D: i32> QuadraticInt<T, D> {
    /// Get the conjugate of the quadratic integer.
    /// 
    /// The conjugate of a quadratic number `x + y√D` is `x - y√D`
    #[inline]
    pub fn conj(self) -> Self {
        match D.rem_euclid(4) {
            1 => Self {
                // conj(a+bw) = conj(a+b/2 + b/2*sq) = a+b/2 - b/2*sq = a+b - bw
                a: self.a + &self.b, b: -self.b
            },
            2 | 3 => Self { a: self.a, b: -self.b },
            _ => unreachable!()
        }
    }
}

impl<T: Add<Output=T> + Sub<Output=T> + for<'r> Mul<&'r T, Output=T> + FromPrimitive, const D: i32> QuadraticInt<T, D> 
where for<'r> &'r T: Mul<&'r T, Output=T> {
    /// Get the norm of the quadratic integer.
    /// 
    /// The norm of a quadratic number `x + y√D` is `x² - Dy²`
    #[inline]
    pub fn norm(&self) -> T {
        match D.rem_euclid(4) {
            1 => {
                // |a+bw| = (a+b/2)^2 - (b/2*sq)^2 = a^2+ab+b^2/4 - D*b^2/4 = a^2 + ab + b^2(1-D)/4
                &self.a * &self.a + &self.a * &self.b + &self.b * &self.b * &T::from_i32((1-D)/4).unwrap()
            },
            2 | 3 => &self.a * &self.a - &(&self.b * &self.b) * &T::from_i32(D).unwrap(),
            _ => unreachable!()
        }
    }
}

impl<T: Add<Output=T>, const D: i32> Add for QuadraticInt<T, D> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self { a: self.a + rhs.a, b: self.b + rhs.b }
    }
}

impl<T: Sub<Output=T>, const D: i32> Sub for QuadraticInt<T, D> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self { a: self.a - rhs.a, b: self.b - rhs.b }
    }
}

impl<T: Add<Output=T> + Mul<Output=T> + FromPrimitive, const D: i32> Mul for QuadraticInt<T, D>
where for<'r> &'r T: Mul<&'r T, Output=T>, {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        match D.rem_euclid(4) {
            1 => Self {
                // (a+bw)*(c+dw) = ac + bdw^2 + (bc+ad)w = ac+bd(D-1)/4 + (bc+ad+bd)w
                a: &self.a * &rhs.a + &self.b * &rhs.b * T::from_i32((D-1)/4).unwrap(),
                b: &self.a * &rhs.b + self.b * (rhs.a + rhs.b)
            },
            2 | 3 => Self {
                a: &self.a * &rhs.a + &self.b * &rhs.b * T::from_i32(D).unwrap(),
                b: self.a * rhs.b + self.b * rhs.a
            },
            _ => unreachable!()
        }
    }
}

fn div_round<T: Add<Output=T> + Sub<Output=T> + FromPrimitive + PartialOrd + Zero>(x: T, y: &T) -> T
where for<'r> &'r T: Div<&'r T, Output=T>, {
    if (x < T::zero()) ^ (y < &T::zero()) {
        &(x - y / &T::from_u8(2).unwrap()) / y
    } else {
        &(x + y / &T::from_u8(2).unwrap()) / y
    }
}

impl<T: Add<Output=T> + Sub<Output=T> + for<'r> Mul<&'r T, Output=T> + FromPrimitive + PartialOrd + Zero, const D: i32> Div for QuadraticInt<T, D>
where for<'r> &'r T: Mul<&'r T, Output=T> + Div<&'r T, Output=T>, {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let n = rhs.norm();
        match D.rem_euclid(4) {
            1 => Self {
                // (a+bw)/(c+dw) = (a+bw)*conj(c+dw)/|c+dw| = (a+bw)(c+d-dw) / |c+dw|
                // = (ac+ad-adw+bcw+bdw-bdw^2)/|c+dw|
                // = (ac+ad-bd*(D-1)/4)/|c+dw| + (-ad+bc)w/|c+dw|
                a: div_round(&self.a * &rhs.a + &self.a * &rhs.b - &self.b * &rhs.b * &T::from_i32((D-1)/4).unwrap(), &n),
                b: div_round(&self.b * &rhs.a - &self.a * &rhs.b, &n)
            },
            2 | 3 => Self {
                // (a+bw)/(c+dw) = (a+bw)*conj(c+dw)/|c+dw| = (a+bw)(c-dw) / |c+dw|
                // = (ac-bd*D)/|c+dw| + (bc-ad)w/|c+dw|
                a: div_round(&self.a * &rhs.a - &self.b * &rhs.b * &T::from_i32(D).unwrap(), &n),
                b: div_round(&self.b * &rhs.a - &self.a * &rhs.b, &n)
            },
            _ => unreachable!()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ops_unit() {
        // unit complex number tests
        let e1 = GaussianInt::<i32> { a: 1, b: 0 };
        let e2 = GaussianInt::<i32> { a: 0, b: 1 };
        let e3 = GaussianInt::<i32> { a: -1, b: 0 };
        let e4 = GaussianInt::<i32> { a: 0, b: -1 };
        let m12 = GaussianInt::<i32> { a: 1, b: 1 };
        let m34 = GaussianInt::<i32> { a: -1, b: -1 };

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
        let q12 = GaussianInt::<i32> { a: 1, b: 2 };
        let qm12 = GaussianInt::<i32> { a: -1, b: 2 };
        let q1m2 = GaussianInt::<i32> { a: 1, b: -2 };
        let q23 = GaussianInt::<i32> { a: 2, b: 3 };
        let qm23 = GaussianInt::<i32> { a: -2, b: 3 };
        let q2m3 = GaussianInt::<i32> { a: 2, b: -3 };

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
        let q12 = EisensteinInt::<i32> { a: 1, b: 2 };
        let qm12 = EisensteinInt::<i32> { a: -1, b: 2 };
        let q3m2 = EisensteinInt::<i32> { a: 3, b: -2 };
        let q23 = EisensteinInt::<i32> { a: 2, b: 3 };
        let qm23 = EisensteinInt::<i32> { a: -2, b: 3 };
        let q5m3 = EisensteinInt::<i32> { a: 5, b: -3 };

        assert_eq!(q12.conj(), q3m2);
        assert_eq!(q23.conj(), q5m3);
        assert_eq!(q12.norm(), 7);
        assert_eq!(q23.norm(), 19);
        assert_eq!(q12 * q12, EisensteinInt { a: -3, b: 8 });

        for &v1 in &[q12, qm12, q3m2, q23, qm23, q5m3] {
            for &v2 in &[q12, qm12, q3m2, q23, qm23, q5m3] {
                assert_eq!(v1 + v2 - v1, v2, "add/sub test failed between {:?} and {:?}", v1, v2);
                assert_eq!(v1 * v2 / v1, v2, "mul/div test failed between {:?} and {:?}", v1, v2);
            }
        }
    }
}
