use num_traits::{FromPrimitive, RefNum, NumRef};
use num_integer::{Integer, Roots, sqrt};
use num_bigint::BigInt;
use std::fmt;

/// A type representation quadratic surd number (a + b*sqrt(r)) / c
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct QuadraticSurd<T> {
    a: T,
    b: T, // zero when reduced if the surd is a rational number
    c: T, // positive when reduced
    r: T  // zero when reduced if the surd is a rational number
}

impl<T> QuadraticSurd<T> {
    #[inline]
    pub const fn new_raw(a: T, b: T, c: T, r: T) -> QuadraticSurd<T> {
        QuadraticSurd { a, b, c, r }
    }
}

impl<T: Integer + NumRef + Clone + Roots> QuadraticSurd<T> 
where for<'r> &'r T: RefNum<T>
{
    // TODO: it's possible to support r < 0, but special handling is needed
    // TODO: add method to reduce root r

    fn reduce(&mut self) {
        if self.c.is_zero() {
            panic!("denominator == 0");
        }

        // test if the surd is rational
        let root = sqrt(self.r.clone());
        if &root * &root == self.r {
            self.a = &self.a + &self.b * root;
            self.b = T::zero();
            self.r = T::zero();
        }

        // reduce common divisor
        let g = self.a.gcd(&self.b).gcd(&self.c);
        self.a = &self.a / &g;
        self.b = &self.b / &g;
        self.c = &self.c / g;

        // keep denom positive
        if self.c < T::zero() {
            self.a = T::zero() - &self.a;
            self.b = T::zero() - &self.b;
            self.c = T::zero() - &self.c;
        }
    }

    #[inline]
    pub fn new(a: T, b: T, c: T, r: T) -> Self {
        assert!(r > T::zero());

        let mut ret = QuadraticSurd::new_raw(a, b, c, r);
        ret.reduce();
        ret
    }

    #[inline]
    pub fn from_integer(target: T) -> Self {
        QuadraticSurd {
            a: target, b: T::zero(), c: T::one(), r: T::zero()
        }
    }

    // TODO: from_complex
    // TODO: from_rational

    #[inline]
    pub fn from_sqrt(target: T) -> Self {
        let root = sqrt(target.clone());
        if &root * &root == target {
            Self::from_integer(target)
        } else {
            QuadraticSurd {
                a: T::zero(), b: T::one(), c: T::one(), r: target
            }
        }
    }

    pub fn recip(&self) -> Self {
        let bb = &self.b * &self.b;
        QuadraticSurd::new(
            T::zero() - &self.c * &self.a,
            &self.c * &self.b,
            &bb * &self.r - &self.a * &self.a,
            self.r.clone()
        )
    }

    pub fn is_integer(&self) -> bool {
        self.c.is_one() && self.r.is_zero()
    }

    pub fn to_integer(&self) -> T {
        if self.is_integer() { self.a.clone() } else { self.trunc().a }
    }

    // TODO: to_rational, see https://github.com/rust-num/num-rational/issues/35

    pub fn trunc(&self) -> Self {
        let br = sqrt(&self.b * &self.b * &self.r);
        let num = if self.b >= T::zero() { &self.a + br } else { &self.a - br };
        return Self::from_integer(num / &self.c);
    }

    pub fn floor(&self) -> Self {
        let br = sqrt(&self.b * &self.b * &self.r);
        let num = if self.b >= T::zero() { &self.a + br } else { &self.a - br - T::one() };
        let num = if num >= T::zero() { num } else {num - &self.c + T::one()};
        return Self::from_integer(num / &self.c);
    }
}

impl<T: fmt::Display> fmt::Display for QuadraticSurd<T>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}√{}) / {}", self.a, self.b, self.r, self.c)
    }
}

// impl QuadraticSurd<i64> {
//     pub fn value(&self) -> f64 {
//         let fa = self.a as f64;
//         let fb = self.b as f64;
//         let fc = self.c as f64;
//         ((self.r as f64).sqrt() * fa + fb) / fc
//     }
// }

// TODO: impl ToPrimitive
// TODO: impl Roots for QuadraticSurd and Rational

#[cfg(test)]
mod tests {
    use super::*;

    pub const PHI: QuadraticSurd<i32> = QuadraticSurd::new_raw(1, 1, 2, 5); // 1.618
    pub const PHI_R: QuadraticSurd<i32> = QuadraticSurd::new_raw(-1, 1, 2, 5); // 0.618
    pub const N_PHI: QuadraticSurd<i32> = QuadraticSurd::new_raw(-1, -1, 2, 5); // -1.618
    pub const N_PHI_R: QuadraticSurd<i32> = QuadraticSurd::new_raw(1, -1, 2, 5); // -0.618

    #[test]
    fn conversion_test() {
        assert_eq!(PHI.floor().to_integer(), 1);
        assert_eq!(PHI_R.floor().to_integer(), 0);
        assert_eq!(PHI.to_integer(), 1);
        assert_eq!(PHI_R.to_integer(), 0);
        assert_eq!(N_PHI.floor().to_integer(), -2);
        assert_eq!(N_PHI_R.floor().to_integer(), -1);
        assert_eq!(N_PHI.to_integer(), -1);
        assert_eq!(N_PHI_R.to_integer(), 0);
    }

    #[test]
    fn arithmic_test() {
        
    }
}
