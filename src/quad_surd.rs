use core::ops::{Add, Sub, Mul, Div, Neg};
use num_traits::{FromPrimitive, ToPrimitive, RefNum, NumRef, Signed, Zero, One, CheckedAdd, CheckedMul};
use num_integer::{Integer, Roots, sqrt};
use std::fmt;
use crate::traits::{Irrational, FromSqrt, Approximation, RationalApproximation};
use crate::cont_frac::ContinuedFraction;

use num_rational::Ratio;
#[cfg(feature = "num-bigint")]
use num_bigint::BigInt;
#[cfg(feature = "num-prime")]
use num_prime::PrimeBuffer;

/// A helper trait to define valid type that can be used for QuadraticSurd
pub trait QuadraticSurdBase: Integer + NumRef + Clone + Roots + Signed {}
impl<T: Integer + NumRef + Clone + Roots + Signed> QuadraticSurdBase for T {}

/// A type representation quadratic surd number (a + b*sqrt(r)) / c
#[derive(PartialEq, Eq, Hash, Clone, Debug, Copy)]
pub struct QuadraticSurd<T>{
    a: T,
    b: T, // zero when reduced if the surd is a rational number
    c: T, // positive when reduced
    r: T  // zero when reduced if the surd is a rational number
}

impl<T> QuadraticSurd<T> {
    #[inline]
    pub(crate) const fn new_raw(a: T, b: T, c: T, r: T) -> Self {
        QuadraticSurd { a, b, c, r }
    }
}

impl<T: Integer> QuadraticSurd<T> {
    #[inline]
    pub fn is_integer(&self) -> bool {
        self.c.is_one() && self.b.is_zero()
    }

    #[inline]
    pub fn is_rational(&self) -> bool {
        self.b.is_zero() || self.r.is_zero()
    }

    #[inline]
    // TODO: disable this if complex is not supported. rather, panic in the creation
    // of the surd
    fn panic_if_complex(&self) {
        if self.r < T::zero() {
            panic!("it's a complex number!");
        }
    }
}

impl<T: Integer> From<T> for QuadraticSurd<T> {
    /// Create a `QuadraticSurd` representation of an integer.
    /// The square root base will be zero.
    #[inline]
    fn from(t: T) -> Self {
        QuadraticSurd {
            a: t, b: T::zero(), c: T::one(), r: T::zero()
        }
    }
}

impl<T: Integer> From<Ratio<T>> for QuadraticSurd<T> {
    /// Create a `QuadraticSurd` representation of a rational number.
    /// The square root base will be zero.
    #[inline]
    fn from(t: Ratio<T>) -> Self {
        let (a, c) = t.into();
        QuadraticSurd {
            a, b: T::zero(), c, r: T::zero()
        }
    }
}

impl<T: QuadraticSurdBase> QuadraticSurd<T>
where for<'r> &'r T: RefNum<T>
{
    // TODO: it's possible to support r < 0, but special handling is needed.
    //       And then we can support from_complex, to_complex, is_complex
    //       maybe enable this only if num-complex is enabled

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

        // TODO: reduce by iteratively check square numbers
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

        unimplemented!() // TODO: implement with unified factorization interface
    }

    #[inline]
    pub fn new(a: T, b: T, c: T, r: T) -> Self {
        assert!(r >= T::zero());

        let mut ret = QuadraticSurd::new_raw(a, b, c, r);
        ret.reduce();
        ret
    }

    /// Create a surd with format `a + b sqrt(r)` where a, b, r are rationals
    #[inline]
    pub fn new_split(a: Ratio<T>, b: Ratio<T>, r: Ratio<T>) -> Self {
        assert!(!r.is_negative());

        let surd_r = r.numer() * r.denom();
        let new_b_denom = b.denom() * r.denom();
        let surd_a = a.numer() * &new_b_denom;
        let surd_b = b.numer() * a.denom();
        let surd_c = a.denom() * new_b_denom;
        let mut ret = QuadraticSurd::new_raw(surd_a, surd_b, surd_c, surd_r);
        ret.reduce();
        ret
    }

    /// Get the root of a quadratic equation `ax^2 + bx + c` with rational coefficients
    /// This method only returns one of the root `(b + sqrt(b^2 - 4ac)) / 2a`, use `conjugate()` for the other root
    #[inline]
    pub fn from_root(a: Ratio<T>, b: Ratio<T>, c: Ratio<T>) -> Option<Self> {
        // degraded cases
        if a.is_zero() {
            if b.is_zero() { return None }
            return Some(Self::from(-c / b))
        }

        let two = T::one() + T::one();
        let four = &two * &two;
        let delta = &b * &b - Ratio::from(four) * &a * c;
        if delta.is_negative() { return None } // TODO: skip this with complex number support

        let aa = Ratio::from(two) * a;
        Some(Self::new_split(-b * aa.recip(), aa.recip(), delta))
    }

    /// Returns a reduced version of self.
    ///
    /// This method will try to eliminate common divisors between a, b, c and
    /// also try to eliminate square factors of r.
    ///
    #[inline]
    pub fn reduced(mut self) -> Self {
        self.reduce();
        self.reduce_root(None);
        self
    }

    /// Returns the reciprocal of the surd
    #[inline]
    pub fn recip(self) -> Self {
        let aa = &self.a * &self.a;
        let bb = &self.b * &self.b;
        QuadraticSurd::new(
            -(&self.c * self.a),
            self.c * self.b,
            &bb * &self.r - aa,
            self.r
        )
    }

    /// Return the conjugate of the surd, i.e. (a - b*sqrt(r)) / c
    #[inline]
    pub fn conjugate(self) -> Self {
        QuadraticSurd {
            a: self.a, b: -self.b, c: self.c, r: self.r
        }
    }

    pub fn trunc(&self) -> Self {
        self.panic_if_complex();

        let br = sqrt(&self.b * &self.b * &self.r);
        let num = if self.b >= T::zero() { &self.a + br } else { &self.a - br };
        return Self::from(num / &self.c);
    }

    pub fn fract(&self) -> Self {
        self.panic_if_complex();

        return self.clone() - self.trunc()
    }

    /// Converts to an integer, rounding towards zero
    #[inline]
    pub fn to_integer(&self) -> Approximation<T> {
        self.panic_if_complex();

        if self.is_integer() {
            Approximation::Exact(self.a.clone())
        } else {
            Approximation::Approximated(self.trunc().a)
        }
    }

    /// Converts to a rational, rounding square root towards zero
    #[inline]
    pub fn to_rational(&self) -> Approximation<Ratio<T>> {
        self.panic_if_complex();

        if self.r.is_zero() {
            Approximation::Exact(Ratio::new_raw(self.a.clone(), self.c.clone()))
        } else {
            Approximation::Approximated(Ratio::new_raw(
                &self.a + sqrt(&self.b * &self.b * &self.r), self.c.clone()))
        }
    }

    /// Rounds towards minus infinity
    #[inline]
    pub fn floor(&self) -> Self {
        let br = sqrt(&self.b * &self.b * &self.r);
        let num = if self.b >= T::zero() { &self.a + br } else { &self.a - br - T::one() };
        let num = if num >= T::zero() { num } else {num - &self.c + T::one()};
        return Self::from(num / &self.c);
    }
}

fn quadsurd_to_f64<T: Integer + ToPrimitive> (a: T, b: T, c: T, r: T) -> Option<f64> {
    Some((a.to_f64()? + b.to_f64()? * r.to_f64()?.sqrt()) / c.to_f64()?)
}

#[cfg(not(feature="num-complex"))]
impl<T: Integer> From<ContinuedFraction<T>> for &QuadraticSurd<T> {
    fn from(s: ContinuedFraction<T>) -> Self {
        unimplemented!()
    }
}

#[cfg(not(feature="num-complex"))]
impl<T: Integer> From<QuadraticSurd<T>> for &ContinuedFraction<T> {
    fn from(s: QuadraticSurd<T>) -> Self {
        // REF: http://www.numbertheory.org/courses/MP313/lectures/lecture17/page5.html
        unimplemented!()
    }
}

#[cfg(feature="num-complex")]
impl<T: Integer> TryFrom<QuadraticSurd<T>> for &ContinuedFraction<T> {
    fn try_from(s: QuadraticSurd<T>) -> Self {
        unimplemented!()
    }
}

// impl<T: Integer + Clone + CheckedAdd + CheckedMul> RationalApproximation<T> for &QuadraticSurd<T> {
//     fn approx_rational(self, limit: &T) -> Approximation<Ratio<T>> {
//         ContinuedFraction::from(self).approx_rational(limit)
//     }
// }

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

        self.mul(rhs.recip())
    }
}

impl<T: QuadraticSurdBase> Neg for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    type Output = QuadraticSurd<T>;
    #[inline]
    fn neg(self) -> QuadraticSurd<T> {
        QuadraticSurd::new_raw(-self.a, -self.b, self.c, self.r)
    }
}

impl<T: QuadraticSurdBase> Zero for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    #[inline]
    fn zero() -> Self {
        QuadraticSurd { a: T::zero(), b: T::zero(), c: T::one(), r: T::zero() }
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }
}

impl<T: QuadraticSurdBase> One for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    #[inline]
    fn one() -> Self {
        QuadraticSurd { a: T::one(), b: T::zero(), c: T::one(), r: T::zero() }
    }
    #[inline]
    fn is_one(&self) -> bool {
        self.a.is_one() && self.b.is_zero() && !self.c.is_zero()
    }
}

impl<T: Integer + FromPrimitive> FromPrimitive for QuadraticSurd<T> {
    #[inline]
    fn from_i64(n: i64) -> Option<Self> {
        T::from_i64(n).map(Self::from)
    }

    #[inline]
    fn from_u64(n: u64) -> Option<Self> {
        T::from_u64(n).map(Self::from)
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
        if self.r < T::zero() { return None; }
        match self.to_integer() {
            Approximation::Exact(v) => v.to_i64(),
            Approximation::Approximated(_) => None
        }
    }

    #[inline]
    fn to_u64(&self) -> Option<u64> {
        if self.r < T::zero() { return None; }
        match self.to_integer() {
            Approximation::Exact(v) => v.to_u64(),
            Approximation::Approximated(_) => None
        }
    }

    #[inline]
    fn to_f64(&self) -> Option<f64> {
        if self.r < T::zero() { return None; }
        quadsurd_to_f64(self.a.to_i64()?, self.b.to_i64()?, self.c.to_i64()?, self.r.to_i64()?)
    }
}

impl<T: QuadraticSurdBase + ToPrimitive> Irrational for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FromSqrtError {
    kind: SqrtErrorKind
}

#[derive(Copy, Clone, Debug, PartialEq)]
enum SqrtErrorKind {
    Overflow,
    Unrepresentable
}

impl<T: QuadraticSurdBase> FromSqrt<T> for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    type Error = FromSqrtError;

    #[inline]
    fn from_sqrt(target: T) -> Result<Self, Self::Error> {
        let root = sqrt(target.clone());
        if &root * &root == target {
            Ok(Self::from(root))
        } else {
            Ok(QuadraticSurd {
                a: T::zero(), b: T::one(), c: T::one(), r: target
            })
        }
    }
}

impl<T: QuadraticSurdBase + CheckedMul> FromSqrt<Ratio<T>> for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    type Error = FromSqrtError;

    #[inline]
    fn from_sqrt(target: Ratio<T>) -> Result<Self, Self::Error> {
        if target.is_integer() {
            return Self::from_sqrt(target.to_integer());
        }

        let dd = target.denom().checked_mul(target.denom());
        let nd = target.numer().checked_mul(target.denom());
        match (dd, nd) {
            (Some(new_c), Some(new_r)) => Ok(
                QuadraticSurd::new(T::zero(), T::one(), new_c, new_r)),
            (_, _) => Err(FromSqrtError { kind: SqrtErrorKind::Overflow })
        }
    }
}

impl<T: QuadraticSurdBase + CheckedMul + std::fmt::Debug> FromSqrt<QuadraticSurd<T>> for QuadraticSurd<T>
where for<'r> &'r T: RefNum<T> {
    type Error = FromSqrtError;

    #[inline] 
    fn from_sqrt(target: QuadraticSurd<T>) -> Result<Self, Self::Error> {
        if target.is_rational() {
            match target.to_rational() {
                Approximation::Exact(v) => return Self::from_sqrt(v), _ => panic!()
            };
        }

        // denote a_c = a/c, b_c = b/c
        // suppose (a_c + b_c * sqrt(r))^2 = target = x + y * sqrt(r)
        // let g = a_c/b_c, then g^2 - 2g(x/y) + r = 0
        // result is available only if g is rational

        let x = Ratio::new(target.a, target.c.clone());
        let y = Ratio::new(target.b, target.c);
        let x_y = x / &y;
        let delta2 = &x_y * &x_y - &target.r;
        if delta2.is_negative() {
            return Err(FromSqrtError { kind: SqrtErrorKind::Unrepresentable });
        } // TODO: if complex number is enabled, don't report error here
        let sqrt_delta = Self::from_sqrt(delta2)?;
        if !sqrt_delta.is_integer() {
            return Err(FromSqrtError { kind: SqrtErrorKind::Unrepresentable });
        }

        let delta2 = match sqrt_delta.to_integer() {
            Approximation::Exact(v) => v, _ => panic!()
        };
        let g = x_y - delta2; // XXX: shall we select another root x_y + delta2?
        let two = T::one() + T::one();
        let c2 = Ratio::from(two * g.numer() * g.denom()) / y;
        debug_assert!(c2.is_integer());
        let c = sqrt(c2.to_integer());
        debug_assert!(&c * &c == c2.to_integer());

        let (a, b) = g.into();
        Ok(QuadraticSurd::new(a, b, c, target.r))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    pub const PHI: QuadraticSurd<i32> = QuadraticSurd::new_raw(1, 1, 2, 5); // 1.618
    pub const PHI_R: QuadraticSurd<i32> = QuadraticSurd::new_raw(-1, 1, 2, 5); // 0.618
    pub const N_PHI: QuadraticSurd<i32> = QuadraticSurd::new_raw(-1, -1, 2, 5); // -1.618
    pub const N_PHI_R: QuadraticSurd<i32> = QuadraticSurd::new_raw(1, -1, 2, 5); // -0.618
    pub const PHI_SQ: QuadraticSurd<i32> = QuadraticSurd::new_raw(3, 1, 2, 5);

    #[test]
    fn new_split_test() {
        assert_eq!(QuadraticSurd::new_split(Ratio::new(1, 2), Ratio::new(1, 2), Ratio::from(5)), PHI);
        assert_eq!(QuadraticSurd::new_split(Ratio::new(-1, 2), Ratio::new(1, 2), Ratio::from(5)), PHI_R);
    }

    #[test]
    fn conversion_test() {
        assert_eq!(PHI.floor().to_i32(), Some(1));
        assert_eq!(PHI_R.floor().to_i32(), Some(0));
        assert_eq!(PHI.to_integer(), Approximation::Approximated(1));
        assert_eq!(PHI_R.to_integer(), Approximation::Approximated(0));
        assert_eq!(N_PHI.floor().to_i32(), Some(-2));
        assert_eq!(N_PHI_R.floor().to_i32(), Some(-1));
        assert_eq!(N_PHI.to_integer(), Approximation::Approximated(-1));
        assert_eq!(N_PHI_R.to_integer(), Approximation::Approximated(0));

        assert!(matches!(PHI.to_f64(),     Some(v) if (v - 1.61803398874989f64).abs() < 1e-10));
        assert!(matches!(PHI_R.to_f64(),   Some(v) if (v - 0.61803398874989f64).abs() < 1e-10));
        assert!(matches!(N_PHI.to_f64(),   Some(v) if (v + 1.61803398874989f64).abs() < 1e-10));
        assert!(matches!(N_PHI_R.to_f64(), Some(v) if (v + 0.61803398874989f64).abs() < 1e-10));
    }

    #[test]
    fn from_xxx_test() {
        assert_eq!(QuadraticSurd::from_sqrt(5).unwrap(), QuadraticSurd::new_raw(0, 1, 1, 5));
        assert!(matches!(QuadraticSurd::from_sqrt(PHI), Err(_)));
        assert!(matches!(QuadraticSurd::from_sqrt(PHI*PHI), Ok(PHI)));

        assert!(matches!(QuadraticSurd::from_root(
            Ratio::from(1), Ratio::from(-1), Ratio::from(-1)), Some(PHI)));
    }

    #[test]
    fn arithmic_test() {
        let sq5 = QuadraticSurd::from_sqrt(5).unwrap();
        assert_eq!(-sq5, QuadraticSurd::new_raw(0, -1, 1, 5));

        // add
        assert_eq!(PHI + PHI_R, sq5);
        assert_eq!(N_PHI + N_PHI_R, -sq5);
        assert!((PHI + N_PHI).is_zero());
        assert!((PHI + N_PHI_R).is_one());

        // sub
        assert!((PHI - PHI).is_zero());
        assert!((PHI - PHI_R).is_one());
        assert!((N_PHI_R - N_PHI).is_one());
        assert_eq!(PHI - N_PHI_R, sq5);
        assert_eq!(N_PHI - PHI_R, -sq5);

        // recip
        assert_eq!(PHI.recip(), PHI_R);
        assert_eq!(N_PHI.recip(), N_PHI_R);

        // conjugate
        assert_eq!(PHI.conjugate(), N_PHI_R);
        assert_eq!(PHI_R.conjugate(), N_PHI);

        // mul
        assert!((PHI * PHI_R).is_one());
        assert!((N_PHI * N_PHI_R).is_one());
        assert_eq!(PHI * PHI, PHI_SQ);
        assert_eq!(N_PHI * N_PHI, PHI_SQ);

        // div
        assert!((PHI / PHI).is_one());
        assert!((N_PHI / N_PHI).is_one());
        assert_eq!(PHI / PHI_R, PHI_SQ);
        assert_eq!(N_PHI / N_PHI_R, PHI_SQ);
    }
}
