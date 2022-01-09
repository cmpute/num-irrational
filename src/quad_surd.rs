use core::ops::{Add, Sub, Mul, Div, Neg};
use num_traits::{FromPrimitive, ToPrimitive, PrimInt, RefNum, NumRef};
use num_integer::{Integer, Roots, sqrt};
use std::fmt;
use crate::traits::{Irrational, FromSqrt};

#[cfg(feature = "num-rational")]
use num_rational::Ratio;
#[cfg(feature = "num-bigint")]
use num_bigint::BigInt;
#[cfg(feature = "num-prime")]
use num_prime::PrimeBuffer;

/// A helper trait to define valid type that can be used for QuadraticSurd
pub trait QuadraticSurdBase: Integer + NumRef + Clone + Roots {}
impl<T: Integer + NumRef + Clone + Roots> QuadraticSurdBase for T {}

/// A type representation quadratic surd number (a + b*sqrt(r)) / c
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct QuadraticSurd<T>{
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

impl<T: Integer> QuadraticSurd<T> {
    /// Create a `QuadraticSurd` representation of an integer.
    /// The square root base will be zero.
    #[inline]
    pub fn from_integer(target: T) -> Self {
        QuadraticSurd {
            a: target, b: T::zero(), c: T::one(), r: T::zero()
        }
    }

    /// Create a `QuadraticSurd` representation of a rational number.
    /// The square root base will be zero.
    #[inline]
    #[cfg(feature = "num-rational")]
    pub fn from_rational(target: Ratio<T>) -> Self {
        QuadraticSurd {
            a: target.numer(), b: T::zero(), c: target.denom(), r: T::zero()
        }
    }

    #[inline]
    pub fn is_integer(&self) -> bool {
        self.c.is_one() && self.r.is_zero()
    }
}

impl<T: QuadraticSurdBase> QuadraticSurd<T>
where for<'r> &'r T: RefNum<T>
{
    // TODO: it's possible to support r < 0, but special handling is needed. And then we can support from_complex and to_complex

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

    /// Try to eliminate factors of root that is a squared number
    /// parameter `factor` is a hint of the factor
    #[cfg(not(feature = "num-prime"))]
    fn reduce_root(&mut self, factor: Option<T>) {
        if let Some(d) = factor {
            let (quo, rem) = self.r.div_rem(&d);
            if rem.is_zero() { // if d is actually a factor
                let root = sqrt(d.clone());
                if &root * &root == d { // if d is a square number
                    let g = root.gcd(&self.a).gcd(&self.c);
                    self.a = &self.a / &g;
                    self.b = &self.b * root / &g;
                    self.c = &self.c / g;
                    self.r = quo;
                }
            }
            debug_assert!(&self.r % &d == T::zero())
        }
    }

    #[cfg(feature = "num-prime")]
    fn reduce_root(&mut self, factor: Option<T>) {
        let mut new_r = r.clone();
        let mut leftout = T::one();

        if let Some(d) = factor {
            let (quo, rem) = self.r.div_rem(&d);
            if rem.is_zero() { // if d is actually a factor
                let root = sqrt(d.clone());
                if &root * &root == d { // if d is a square number
                    new_r = quo;
                    leftout = root;
                }
            }
        }

        panic!("not implemented") // TODO: implement with unified factorization interface
    }

    #[inline]
    pub fn new(a: T, b: T, c: T, r: T) -> Self {
        assert!(r > T::zero());

        let mut ret = QuadraticSurd::new_raw(a, b, c, r);
        ret.reduce();
        ret
    }

    /// Returns a reduced copy of self.
    ///
    /// This method will try to eliminate common divisors between a, b, c and
    /// also try to eliminate square factors of r.
    ///
    #[inline]
    pub fn reduced(&self) -> Self {
        let mut ret = self.clone();
        ret.reduce();
        ret.reduce_root(None);
        ret
    }

    /// Returns the reciprocal of the surd
    #[inline]
    pub fn recip(&self) -> Self {
        self.clone().into_recip()
    }

    pub fn into_recip(self) -> Self {
        let aa = &self.a * &self.a;
        let bb = &self.b * &self.b;
        QuadraticSurd::new(
            T::zero() - &self.c * self.a,
            self.c * self.b,
            &bb * &self.r - aa,
            self.r
        )
    }

    pub fn trunc(&self) -> Self {
        let br = sqrt(&self.b * &self.b * &self.r);
        let num = if self.b >= T::zero() { &self.a + br } else { &self.a - br };
        return Self::from_integer(num / &self.c);
    }

    pub fn fract(&self) -> Self {
        return self.clone() - self.trunc()
    }

    /// Converts to an integer, rounding towards zero
    #[inline]
    pub fn to_integer(&self) -> T {
        if self.is_integer() { self.a.clone() } else { self.trunc().a }
    }

    /// Converts to a rational, rounding square root towards zero
    #[inline]
    #[cfg(feature = "num-rational")]
    pub fn to_rational(&self) -> Ratio<T> {
        if self.r.is_zero() {
            Ratio::new_raw(self.a, self.c)
        } else {
            Ratio::new_raw(self.a + sqrt(&self.b * &self.b * &self.r), self.c)
        }
    }

    /// Rounds towards minus infinity
    #[inline]
    pub fn floor(&self) -> Self {
        let br = sqrt(&self.b * &self.b * &self.r);
        let num = if self.b >= T::zero() { &self.a + br } else { &self.a - br - T::one() };
        let num = if num >= T::zero() { num } else {num - &self.c + T::one()};
        return Self::from_integer(num / &self.c);
    }
}

fn quadsurd_to_f64<T: Integer + ToPrimitive> (a: T, b: T, c: T, r: T) -> Option<f64> {
    Some((a.to_f64()? + b.to_f64()? * r.to_f64()?.sqrt()) / c.to_f64()?)
}

#[cfg(feature = "num-rational")]
impl<T: PrimInt> QuadraticSurd<T> {
    pub fn to_accurate_rational(&self) -> Ratio<T> {
        panic!("not implemented") // TODO: implement
    }
}

#[cfg(all(feature = "num-rational", feature = "num-bigint"))]
impl QuadraticSurd<BigInt> {
    pub fn to_accurate_rational(&self, iterations: u32) -> Ratio<BigInt> {
        panic!("not implemented") // TODO: implement
    }
}

impl<T: fmt::Display> fmt::Display for QuadraticSurd<T>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}√{}) / {}", self.a, self.b, self.r, self.c)
    }
}

macro_rules! arith_impl {
    (impl $imp:ident, $method:ident) => {
        // Abstracts a/b `op` c/d = (a*lcm/b `op` c*lcm/d)/lcm where lcm = lcm(b,d)
        impl<T: QuadraticSurdBase> $imp<QuadraticSurd<T>> for QuadraticSurd<T>
        where for<'r> &'r T: RefNum<T> {
            type Output = QuadraticSurd<T>;
            fn $method(self, rhs: QuadraticSurd<T>) -> QuadraticSurd<T> {
                if self.r != rhs.r {
                    panic!("The root base should be the same!");
                }
                if self.c == rhs.c {
                    return QuadraticSurd::new(self.a.$method(rhs.a), self.b.$method(rhs.b), rhs.c, rhs.r);
                }

                let lcm = self.c.lcm(&rhs.c);
                let lhs_r = &lcm / self.c;
                let rhs_r = &lcm / rhs.c;
                QuadraticSurd::new(
                    (self.a / &lhs_r).$method(rhs.a / &rhs_r),
                    (self.b / &lhs_r).$method(rhs.b / &rhs_r),
                    lcm, rhs.r)
            }
        }
        // Abstracts the a/b `op` c/1 = (a*1 `op` b*c) / (b*1) = (a `op` b*c) / b pattern
        impl<T: QuadraticSurdBase> $imp<T> for QuadraticSurd<T>
        where for<'r> &'r T: RefNum<T> {
            type Output = QuadraticSurd<T>;
            #[inline]
            fn $method(self, rhs: T) -> QuadraticSurd<T> {
                QuadraticSurd::new(
                    self.a.$method(&self.c * rhs), self.b, self.c, self.r)
            }
        }
    };
}

arith_impl!(impl Add, add);
arith_impl!(impl Sub, sub);

impl<T: QuadraticSurdBase> Mul<T> for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    type Output = QuadraticSurd<T>;
    #[inline]
    fn mul(self, rhs: T) -> QuadraticSurd<T> {
        let gcd = self.c.gcd(&rhs);
        let rem = rhs / &gcd;
        QuadraticSurd::new(self.a * &rem, self.b * rem, self.c / gcd, self.r)
    }
}

impl<T: QuadraticSurdBase> Mul<QuadraticSurd<T>> for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    type Output = QuadraticSurd<T>;
    #[inline]
    fn mul(self, rhs: QuadraticSurd<T>) -> QuadraticSurd<T> {
        if self.r != rhs.r {
            panic!("The root base should be the same!");
        }

        let gcd_lr = (&self.a).gcd(&self.b).gcd(&rhs.c); // gcd between lhs numerator and rhs denominator
        let gcd_rl = (&rhs.a).gcd(&rhs.b).gcd(&self.c); // gcd between rhs numerator and lhs denominator

        let (la, lb, lc) = (&self.a / &gcd_lr, &self.b / &gcd_lr, &rhs.c / &gcd_rl);
        let (ra, rb, rc) = (&rhs.a / &gcd_rl, &rhs.b / &gcd_rl, &self.c / &gcd_lr);
        QuadraticSurd::new(
            &la * &ra + &lb * &rb * self.r,
            &la * &rb + &lb * &ra,
            lc * rc, rhs.r)
    }
}

impl<T: QuadraticSurdBase> Div<T> for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    type Output = QuadraticSurd<T>;
    #[inline]
    fn div(self, rhs: T) -> QuadraticSurd<T> {
        let gcd = self.a.gcd(&self.b).gcd(&rhs);
        let rem = rhs / &gcd;
        QuadraticSurd::new(self.a / &gcd, self.b / &gcd, self.c * rem, self.r)
    }
}

impl<T: QuadraticSurdBase> Div<QuadraticSurd<T>> for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    type Output = QuadraticSurd<T>;
    #[inline]
    fn div(self, rhs: QuadraticSurd<T>) -> QuadraticSurd<T> {
        if self.r != rhs.r {
            panic!("The root base should be the same!");
        }

        self.mul(rhs.into_recip())
    }
}

impl<T: Integer> From<T> for QuadraticSurd<T>
{
    fn from(x: T) -> QuadraticSurd<T> {
        QuadraticSurd::from_integer(x)
    }
}

impl<T: Integer + FromPrimitive> FromPrimitive for QuadraticSurd<T> {
    #[inline]
    fn from_i64(n: i64) -> Option<Self> {
        T::from_i64(n).map(Self::from_integer)
    }

    #[inline]
    fn from_u64(n: u64) -> Option<Self> {
        T::from_u64(n).map(Self::from_integer)
    }

    #[inline]
    fn from_f64(_: f64) -> Option<Self> {
        None // it makes no sense to convert a float number to a quadratic surd
    }
}

impl<T: QuadraticSurdBase + ToPrimitive> ToPrimitive for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    #[inline]
    fn to_i64(&self) -> Option<i64> {
        self.to_integer().to_i64()
    }

    #[inline]
    fn to_u64(&self) -> Option<u64> {
        self.to_integer().to_u64()
    }

    #[inline]
    fn to_f64(&self) -> Option<f64> {
        quadsurd_to_f64(self.a.to_i64()?, self.b.to_i64()?, self.c.to_i64()?, self.r.to_i64()?)
    }
}

impl<T: QuadraticSurdBase + ToPrimitive> Irrational for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {}

impl<T: Integer + Roots + Clone> FromSqrt<T> for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    #[inline]
    fn from_sqrt(target: T) -> Self {
        let root = sqrt(target.clone());
        if &root * &root == target {
            Self::from_integer(root)
        } else {
            QuadraticSurd {
                a: T::zero(), b: T::one(), c: T::one(), r: target
            }
        }
    }
}

#[cfg(feature = "num-rational")]
impl<T> FromSqrt<Ratio<T>> for QuadraticSurd<T> {
    #[inline]
    fn from_sqrt(target: Ratio<T>) -> Self {
        QuadraticSurd::from_integer(target.numer) / QuadraticSurd::from_integer(target.denom)
    }
}

// TODO: impl TryFromSqrt for Integer, QuadraticSurd and Rational.
// Note that square root of a QuadraticSurd can be (but not always) a QuadraticSurd as well

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

        assert!(matches!(PHI.to_f64(),     Some(v) if (v - 1.61803398874989f64).abs() < 1e-10));
        assert!(matches!(PHI_R.to_f64(),   Some(v) if (v - 0.61803398874989f64).abs() < 1e-10));
        assert!(matches!(N_PHI.to_f64(),   Some(v) if (v + 1.61803398874989f64).abs() < 1e-10));
        assert!(matches!(N_PHI_R.to_f64(), Some(v) if (v + 0.61803398874989f64).abs() < 1e-10));
    }

    #[test]
    fn arithmic_test() {
        assert_eq!((PHI + PHI_R).to_f64().unwrap(), 5f64.sqrt());
        assert_eq!((N_PHI + N_PHI_R).to_f64().unwrap(), -(5f64.sqrt()));
        assert_eq!((PHI + N_PHI).to_f64().unwrap(), 0f64);
        assert_eq!((PHI + N_PHI_R).to_f64().unwrap(), 1f64);
    }
}
