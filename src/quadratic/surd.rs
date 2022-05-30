//! Implementation of qudratic irrational numbers

use super::{QuadraticBase, QuadraticOps};
use crate::cont_frac::ContinuedFraction;
use crate::traits::{
    Approximation, Computable, FromSqrt, FromSqrtError, SqrtErrorKind, WithSigned, WithUnsigned,
};
#[cfg(feature = "complex")]
use crate::GaussianInt;
use crate::QuadraticInt;
use core::convert::TryFrom;
use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub};
use num_integer::{sqrt, Integer};
use num_traits::{
    CheckedAdd, CheckedMul, FromPrimitive, NumRef, One, RefNum, Signed, ToPrimitive, Zero,
};
use std::fmt;

use num_rational::Ratio;

/// Underlying representation of a quadratic surd `(a + b*√r) / c` as (a,b,c).
///
/// Note that the representation won't be reduced (normalized) after operations.
/// Therefore, the equality test will requires the coefficients to be completely same,
/// rather than be equal numerically.
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct QuadraticSurdCoeffs<T>(pub T, pub T, pub T);

impl<T> From<(T, T, T)> for QuadraticSurdCoeffs<T> {
    fn from(v: (T, T, T)) -> Self {
        Self(v.0, v.1, v.2)
    }
}

impl<T: QuadraticBase> Add for QuadraticSurdCoeffs<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        if self.2 == rhs.2 {
            Self(self.0 + rhs.0, self.1 + rhs.1, rhs.2)
        } else {
            Self(
                self.0 * &rhs.2 + rhs.0 * &self.2,
                self.1 * &rhs.2 + rhs.1 * &self.2,
                rhs.2 * self.2,
            )
        }
    }
}

impl<T: QuadraticBase> Sub for QuadraticSurdCoeffs<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        if self.2 == rhs.2 {
            Self(self.0 - rhs.0, self.1 - rhs.1, rhs.2)
        } else {
            Self(
                self.0 * &rhs.2 - rhs.0 * &self.2,
                self.1 * &rhs.2 - rhs.1 * &self.2,
                rhs.2 * self.2,
            )
        }
    }
}

impl<T: QuadraticBase> QuadraticOps<Self, &T, Self> for QuadraticSurdCoeffs<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Scalar = Ratio<T>;
    fn mul(self, rhs: Self, discr: &T) -> Self {
        let Self(la, lb, lc) = self;
        let Self(ra, rb, rc) = rhs;
        Self(&la * &ra + &lb * &rb * discr, la * rb + lb * ra, lc * rc)
    }
    fn div(self, rhs: Self, discr: &T) -> Self {
        let Self(la, lb, lc) = self;
        let Self(ra, rb, rc) = rhs;
        let c = lc * (&ra * &ra - &rb * &rb * discr);
        Self(
            &rc * (&la * &ra - &lb * &rb * discr),
            rc * (lb * ra - la * rb),
            c,
        )
    }
    #[inline]
    fn conj(self, _: &T) -> Self {
        Self(self.0, -self.1, self.2)
    }
    #[inline]
    fn norm(self, discr: &T) -> Self::Scalar {
        Ratio::new(
            &self.0 * &self.0 - &self.1 * &self.1 * discr,
            &self.2 * &self.2,
        )
    }
}

impl<T: QuadraticBase> Add<Ratio<T>> for QuadraticSurdCoeffs<T> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Ratio<T>) -> Self {
        let Self(la, lb, lc) = self;
        let (ra, rc) = rhs.into();
        if lc == rc {
            Self(la + ra, lb, lc)
        } else {
            Self(la * &rc + ra * &lc, lb * &rc, lc * rc)
        }
    }
}

impl<T: QuadraticBase> Sub<Ratio<T>> for QuadraticSurdCoeffs<T> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Ratio<T>) -> Self {
        let Self(la, lb, lc) = self;
        let (ra, rc) = rhs.into();
        if lc == rc {
            Self(la - ra, lb, lc)
        } else {
            Self(la * &rc - ra * &lc, lb * &rc, lc * rc)
        }
    }
}

impl<T: QuadraticBase> Mul<Ratio<T>> for QuadraticSurdCoeffs<T> {
    type Output = Self;
    fn mul(self, rhs: Ratio<T>) -> Self {
        let (ra, rc) = rhs.into();
        Self(self.0 * &ra, self.1 * ra, self.2 * rc)
    }
}

impl<T: QuadraticBase> Div<Ratio<T>> for QuadraticSurdCoeffs<T> {
    type Output = Self;
    fn div(self, rhs: Ratio<T>) -> Self {
        let (ra, rc) = rhs.into();
        Self(self.0 * &rc, self.1 * rc, self.2 * ra)
    }
}

impl<T: QuadraticBase> Add<T> for QuadraticSurdCoeffs<T> {
    type Output = Self;
    fn add(self, rhs: T) -> Self::Output {
        Self(self.0 + rhs * &self.2, self.1, self.2)
    }
}

impl<T: QuadraticBase> Sub<T> for QuadraticSurdCoeffs<T> {
    type Output = Self;
    fn sub(self, rhs: T) -> Self::Output {
        Self(self.0 - rhs * &self.2, self.1, self.2)
    }
}

impl<T: QuadraticBase> Mul<T> for QuadraticSurdCoeffs<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        Self(self.0 * &rhs, self.1 * rhs, self.2)
    }
}

impl<T: QuadraticBase> Div<T> for QuadraticSurdCoeffs<T> {
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        Self(self.0, self.1, self.2 * rhs)
    }
}

/// A quadratic number represented as `(a + b*√r) / c`.
/// If the support for complex number is enabled, then this struct can represent
/// any quadratic integers.
///
/// # About the `complex` feature
///
/// Whether enabling the `complex` feature of the crate will lead to some behavior changes
/// of this struct (in a compatible way), some important ones are listed below:
/// - The [new()][QuadraticSurd::new] constructor will panic if `complex` is disabled and the given root is negative
/// - Rounding functions ([fract()][QuadraticSurd::fract], [trunc()][QuadraticSurd::trunc], ..) will round to Gaussian
///   integers instead of rational integers if `complex` is enabled.
/// - Rational approximation ([approximated()][QuadraticSurd::approximated]) will panic if `complex` is enabled and
///   the quadratic surd is complex
/// - [QuadraticSurd::from_sqrt()] won't fail as [SqrtErrorKind::Complex] for negative input if `complex` is enabled
///
#[derive(Hash, Clone, Debug, Copy)]
pub struct QuadraticSurd<T> {
    coeffs: QuadraticSurdCoeffs<T>, // when reduced, coeffs.1 is zero if the surd is rational, coeffs.2 is always positive
    discr: T,                       // when reduced, discriminant is zero if the surd is rational
}

impl<T> QuadraticSurd<T> {
    #[inline]
    pub(crate) const fn new_raw(a: T, b: T, c: T, r: T) -> Self {
        Self {
            coeffs: QuadraticSurdCoeffs(a, b, c),
            discr: r,
        }
    }

    /// Get return-only references to the components `(a, b, c, r)`
    pub const fn parts(&self) -> (&T, &T, &T, &T) {
        (&self.coeffs.0, &self.coeffs.1, &self.coeffs.2, &self.discr)
    }
}

impl<T: Integer> QuadraticSurd<T> {
    /// Determine if the surd is an (rational) integer
    #[inline]
    pub fn is_integer(&self) -> bool {
        self.coeffs.2.is_one() && self.is_rational()
    }

    /// Determine if the surd is a Gaussian integer
    #[inline]
    #[cfg(feature = "complex")]
    pub fn is_gaussint(&self) -> bool {
        self.coeffs.2.is_one()
    }

    /// Determine if the surd is a quadratic integer
    #[inline]
    pub fn is_quadint(&self) -> bool {
        self.coeffs.2.is_one()
    }

    /// Determine if the surd is a rational number
    #[inline]
    pub fn is_rational(&self) -> bool {
        self.coeffs.1.is_zero() || self.discr.is_zero()
    }

    /// Determine if the surd is a complex number
    #[inline]
    #[cfg(feature = "complex")]
    pub fn is_complex(&self) -> bool {
        self.discr < T::zero()
    }

    /// Determine if the quadratic number has no rational part (i.e. a = 0, b != 0)
    #[inline]
    pub fn is_pure(&self) -> bool {
        self.coeffs.0.is_zero() && !self.coeffs.1.is_zero()
    }

    /// Panic when the number is complex and the feature "complex" is disabled
    #[inline(always)]
    fn panic_if_complex(&self) {
        // only need to check when we allow construction of complex surd number
        #[cfg(feature = "complex")]
        if self.discr < T::zero() {
            panic!("operation not supported on complex numbers");
        }
    }
}

impl<T: Integer> From<T> for QuadraticSurd<T> {
    /// Create a `QuadraticSurd` representation of an integer.
    /// The square root base will be zero.
    #[inline]
    fn from(t: T) -> Self {
        Self {
            coeffs: QuadraticSurdCoeffs(t, T::zero(), T::one()),
            discr: T::zero(),
        }
    }
}

impl<T: Integer> From<Ratio<T>> for QuadraticSurd<T> {
    /// Create a `QuadraticSurd` representation of a rational number.
    /// The square root base will be zero.
    #[inline]
    fn from(t: Ratio<T>) -> Self {
        let (a, c) = t.into();
        Self {
            coeffs: QuadraticSurdCoeffs(a, T::zero(), c),
            discr: T::zero(),
        }
    }
}

impl<T: QuadraticBase> QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    // Simplify the surd into normalized form
    fn reduce(&mut self) {
        if self.coeffs.2.is_zero() {
            panic!("denominator is zero");
        }

        // test and reduce if the surd is rational
        let root = sqrt(self.discr.abs());
        if &root * &root == self.discr {
            if self.discr.is_negative() {
                self.coeffs.1 = &self.coeffs.1 * root;
                self.discr = -T::one();
            } else {
                self.coeffs.0 = &self.coeffs.0 + &self.coeffs.1 * root;
                self.coeffs.1 = T::zero();
                self.discr = T::zero();
            }
        }

        if self.coeffs.1.is_zero() || self.discr.is_zero() {
            // shortcut if the surd is rational
            self.discr = T::zero();
            self.coeffs.1 = T::zero();

            let g_ac = self.coeffs.0.gcd(&self.coeffs.2);
            self.coeffs.0 = &self.coeffs.0 / &g_ac;
            self.coeffs.2 = &self.coeffs.2 / g_ac;
        } else {
            // reduce common divisor
            let mut g_ac = self.coeffs.0.gcd(&self.coeffs.2);
            let g_acr = (&g_ac * &g_ac).gcd(&self.discr); // test if squared factor of c can divide r
            let groot = sqrt(g_acr.clone());
            if &groot * &groot == g_acr {
                self.discr = &self.discr / &g_acr;
                self.coeffs.0 = &self.coeffs.0 / &groot;
                self.coeffs.2 = &self.coeffs.2 / &groot;
                g_ac = &g_ac / groot;
            }

            let g_abc = g_ac.gcd(&self.coeffs.1);
            self.coeffs.0 = &self.coeffs.0 / &g_abc;
            self.coeffs.1 = &self.coeffs.1 / &g_abc;
            self.coeffs.2 = &self.coeffs.2 / g_abc;
        }

        // keep denom positive
        if self.coeffs.2 < T::zero() {
            self.coeffs.0 = T::zero() - &self.coeffs.0;
            self.coeffs.1 = T::zero() - &self.coeffs.1;
            self.coeffs.2 = T::zero() - &self.coeffs.2;
        }
    }

    // Try to reduce the root base with a possible factor.
    fn reduce_root_hinted(self, hint: T) -> Self {
        let hint = hint.abs();
        if hint.is_zero() || hint.is_one() || hint.is_negative() {
            return self;
        }

        let (quo, rem) = self.discr.div_rem(&hint);
        if rem.is_zero() {
            // if hint is actually a factor
            let root = sqrt(hint.clone());
            if &root * &root == hint {
                // if hint is a square number
                let g = root.gcd(&self.coeffs.0).gcd(&self.coeffs.2);
                return QuadraticSurd::new_raw(
                    self.coeffs.0 / &g,
                    self.coeffs.1 * root / &g,
                    self.coeffs.2 / g,
                    quo,
                );
            }
        }

        self
    }

    /// Create a surd represented as `(a + b√r)) / c` where `a`, `b`, `c`, `r` are integers.
    ///
    /// # Panics
    /// If `c` is zero or `r` is negative when the `complex` feature is not enabled.
    #[inline]
    pub fn new(a: T, b: T, c: T, r: T) -> Self {
        #[cfg(not(feature = "complex"))]
        if r.is_negative() {
            panic!("Negative root is not supported without the `complex` feature");
        }

        let mut ret = QuadraticSurd::new_raw(a, b, c, r);
        ret.reduce();
        ret
    }

    #[inline]
    pub fn from_coeffs(coeffs: QuadraticSurdCoeffs<T>, r: T) -> Self {
        #[cfg(not(feature = "complex"))]
        if r.is_negative() {
            panic!("Negative root is not supported without the `complex` feature");
        }

        let mut ret = Self { coeffs, discr: r };
        ret.reduce();
        ret
    }

    /// Create a surd represented as `a + b√r` where a, b, r are rationals.
    ///
    /// # Panics
    /// If `r` is negative when the `complex` feature is not enabled.
    #[inline]
    pub fn from_rationals(a: Ratio<T>, b: Ratio<T>, r: Ratio<T>) -> Self {
        #[cfg(not(feature = "complex"))]
        if r.is_negative() {
            panic!("Negative root is not supported without the `complex` feature");
        }

        let surd_r = r.numer() * r.denom();
        let new_b_denom = b.denom() * r.denom();
        let surd_a = a.numer() * &new_b_denom;
        let surd_b = b.numer() * a.denom();
        let surd_c = a.denom() * new_b_denom;
        let mut ret = QuadraticSurd::new_raw(surd_a, surd_b, surd_c, surd_r);
        ret.reduce();
        ret
    }

    /// Get the root of a quadratic equation `ax^2 + bx + c` with rational coefficients.
    /// This method only returns one of the root `(b + √(b^2 - 4ac)) / 2a`, use `conjugate()` for the other root
    /// If there are only complex solutions, then `None` will be returned if the `complex` feature is not enabled.
    #[inline]
    pub fn from_equation(a: Ratio<T>, b: Ratio<T>, c: Ratio<T>) -> Option<Self> {
        // degraded cases
        if a.is_zero() {
            if b.is_zero() {
                return None;
            }
            return Some(Self::from(-c / b));
        }

        let two = T::one() + T::one();
        let four = &two * &two;
        let delta = &b * &b - Ratio::from(four) * &a * c;

        #[cfg(not(feature = "complex"))]
        if delta.is_negative() {
            return None;
        }

        let aa = Ratio::from(two) * a;
        Some(Self::from_rationals(-b * aa.recip(), aa.recip(), delta))
    }

    /// Returns the reciprocal of the surd, TODO: implement as num_traits::Inv trait
    #[inline]
    pub fn recip(self) -> Self {
        let QuadraticSurdCoeffs(a, b, c) = self.coeffs;
        let aa = &a * &a;
        let bb = &b * &b;
        QuadraticSurd::new(-(&c * a), c * b, bb * &self.discr - aa, self.discr)
    }

    /// `.recip()` with reference
    #[inline]
    pub fn recip_ref(&self) -> Self {
        let QuadraticSurdCoeffs(a, b, c) = &self.coeffs;
        let aa = a * a;
        let bb = b * b;
        QuadraticSurd::new(-(c * a), c * b, bb * &self.discr - aa, self.discr.clone())
    }

    /// Return the conjugate of the surd, i.e. `(a - b√r) / c`
    #[inline]
    pub fn conj(self) -> Self {
        Self {
            coeffs: self.coeffs.conj(&self.discr),
            discr: self.discr,
        }
    }

    /// `.conj()` with reference
    #[inline]
    pub fn conj_ref(&self) -> Self {
        self.clone().conj()
    }

    /// Round the surd toward zero. The result will be an integer if the surd is real,
    /// or a Gaussian integer if the surd is complex
    ///
    /// # Panics
    /// if the number is complex, when the `complex` feature is not enabled
    pub fn trunc(self) -> Self {
        let QuadraticSurdCoeffs(a, b, c) = self.coeffs;
        let bsign = b.signum();

        if self.discr.is_negative() {
            let br = bsign * sqrt(&b * &b * -self.discr);
            Self::from_coeffs(QuadraticSurdCoeffs(a / &c, br / c, T::one()), -T::one())
        } else {
            let br = bsign * sqrt(&b * &b * self.discr);
            return Self::from((a + br) / c);
        }
    }

    /// [QuadraticSurd::trunc()] with reference
    #[inline]
    pub fn trunc_ref(&self) -> Self {
        self.clone().trunc()
    }

    /// Get the fractional part of the surd, ensuring `self.trunc() + self.fract() == self`
    ///
    /// # Panics
    /// If the square root base is negative and not -1 (the result won't be representable by this struct)
    #[inline]
    pub fn fract(self) -> Self {
        let trunc = self.trunc_ref();
        self - trunc
    }

    /// [QuadraticSurd::fract()] with reference
    #[inline]
    pub fn fract_ref(&self) -> Self {
        self.clone() - self.trunc_ref()
    }

    /// Get the numerator of the surd
    #[inline]
    pub fn numer(&self) -> Self {
        Self::new_raw(
            self.coeffs.0.clone(),
            self.coeffs.1.clone(),
            T::one(),
            self.discr.clone(),
        )
    }

    /// Get the denumerator of the surd
    #[inline]
    pub fn denom(&self) -> Self {
        Self::from(self.coeffs.2.clone())
    }

    /// Converts to an integer, rounding towards zero
    ///
    /// # Panics
    /// if the number is complex when feature `complex` is not enabled
    #[inline]
    pub fn to_integer(&self) -> Approximation<T> {
        self.panic_if_complex();

        if self.is_integer() {
            Approximation::Exact(self.coeffs.0.clone())
        } else {
            Approximation::Approximated(self.trunc_ref().coeffs.0)
        }
    }

    /// Converts to an Gaussian integer, rounding towards zero
    #[inline]
    #[cfg(feature = "complex")]
    pub fn to_gaussint(&self) -> Approximation<GaussianInt<T>> {
        if self.is_gaussint() {
            Approximation::Exact(GaussianInt::new(
                self.coeffs.0.clone(),
                self.coeffs.1.clone(),
            ))
        } else {
            let trunc = self.trunc_ref();
            Approximation::Approximated(GaussianInt::new(trunc.coeffs.0, trunc.coeffs.1))
        }
    }

    /// Converts to an quadratic integer, rounding towards zero
    #[inline]
    pub fn to_quadint(&self) -> Approximation<QuadraticInt<T>> {
        let QuadraticSurdCoeffs(a, b, c) = &self.coeffs;
        if self.is_quadint() {
            Approximation::Exact(QuadraticInt::new(a.clone(), b.clone(), self.discr.clone()))
        } else {
            Approximation::Approximated(QuadraticInt::new(a / c, b / c, self.discr.clone()))
        }
    }

    /// Converts to a rational, rounding square root towards zero
    ///
    /// # Panics
    /// if the number is complex
    #[inline]
    pub fn to_rational(&self) -> Approximation<Ratio<T>> {
        self.panic_if_complex();

        if self.discr.is_zero() {
            Approximation::Exact(Ratio::new_raw(self.coeffs.0.clone(), self.coeffs.2.clone()))
        } else {
            Approximation::Approximated(Ratio::new_raw(
                &self.coeffs.0 + sqrt(&self.coeffs.1 * &self.coeffs.1 * &self.discr),
                self.coeffs.2.clone(),
            ))
        }
    }

    /// Rounds towards minus infinity
    ///
    /// # Panics
    /// if the number is complex
    pub fn floor(self) -> Self {
        self.panic_if_complex();

        let br = sqrt(&self.coeffs.1 * &self.coeffs.1 * &self.discr);
        let num = if self.coeffs.1 >= T::zero() {
            &self.coeffs.0 + br
        } else {
            &self.coeffs.0 - br - T::one()
        };
        let num = if num >= T::zero() {
            num
        } else {
            num - &self.coeffs.2 + T::one()
        };
        return Self::from(num / &self.coeffs.2);
    }

    /// `.floor()` with reference
    ///
    /// # Panics
    /// if the number is complex
    #[inline]
    pub fn floor_ref(&self) -> Self {
        self.clone().floor()
    }

    fn ceil(self) -> Self { unimplemented!() }

    /// Returns the nearest integer to a number. Round half-way cases away from zero?.
    fn round(self) -> Self { unimplemented!() }

    /// Test if the surd number is positive
    ///
    /// # Panics
    /// if the number is complex
    pub fn is_positive(&self) -> bool {
        self.panic_if_complex();

        if self.coeffs.1.is_zero() {
            self.coeffs.0.is_positive()
        } else if self.coeffs.1.is_positive() {
            if self.coeffs.0.is_positive() {
                true
            } else {
                &self.coeffs.0 * &self.coeffs.0 < &self.coeffs.1 * &self.coeffs.1 * &self.discr
            }
        } else {
            // self.coeffs.1.is_negative()
            if !self.coeffs.0.is_positive() {
                false
            } else {
                &self.coeffs.0 * &self.coeffs.0 > &self.coeffs.1 * &self.coeffs.1 * &self.discr
            }
        }
    }

    /// Test if the surd number is negative
    ///
    /// # Panics
    /// if the number is complex
    pub fn is_negative(&self) -> bool {
        self.panic_if_complex();

        if self.coeffs.1.is_zero() {
            self.coeffs.0.is_negative()
        } else if self.coeffs.1.is_negative() {
            if self.coeffs.0.is_negative() {
                true
            } else {
                &self.coeffs.0 * &self.coeffs.0 < &self.coeffs.1 * &self.coeffs.1 * &self.discr
            }
        } else {
            // self.coeffs.1.is_positive()
            if !self.coeffs.0.is_negative() {
                false
            } else {
                &self.coeffs.0 * &self.coeffs.0 > &self.coeffs.1 * &self.coeffs.1 * &self.discr
            }
        }
    }
}

impl<T: QuadraticBase> PartialEq for QuadraticSurd<T>
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

        let QuadraticSurdCoeffs(la, lb, lc) = &self.coeffs;
        let QuadraticSurdCoeffs(ra, rb, rc) = &other.coeffs;

        // first compare the rational part
        // FIXME: handle overflow like num-rational
        if la * rc != ra * lc {
            return false;
        }
        // then compare the quadratic part
        lb * lb * &self.discr * rc * rc == rb * rb * &other.discr * lc * lc
    }
}

impl<T: QuadraticBase> Eq for QuadraticSurd<T> where for<'r> &'r T: RefNum<T> {}

// TODO: implement PartialOrd

impl<T> Into<(T, T, T, T)> for QuadraticSurd<T> {
    /// Deconstruct the quadratic surd `(a + b√r) / c` into tuple `(a,b,c,r)`
    fn into(self) -> (T, T, T, T) {
        (self.coeffs.0, self.coeffs.1, self.coeffs.2, self.discr)
    }
}

impl<T: QuadraticBase + FromPrimitive + CheckedMul> QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    /// Try to eliminate factors of root that is a squared number.
    ///
    /// It will only try trivial division for several small primes.
    /// For a quadratic surd with large root, consider deconstructing
    /// the surd by `.into()` and then reduce yourself.
    fn reduce_root(&mut self) {
        const SMALL_PRIMES: [u8; 54] = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
            89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
            181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
        ];
        for p in SMALL_PRIMES {
            let p = T::from_u8(p).unwrap();
            if let Some(p2) = p.checked_mul(&p) {
                loop {
                    let (quo, rem) = self.discr.div_rem(&p2);
                    if rem.is_zero() {
                        self.discr = quo;
                        self.coeffs.1 = &self.coeffs.1 * &p;
                        // common factors between a, b, c will be handled by reduce()
                    } else {
                        break;
                    }
                }
            }
        }
    }

    /// Returns a reduced version of self.
    ///
    /// This method will try to eliminate common divisors between a, b, c and
    /// also try to eliminate square factors of r.
    ///
    #[inline]
    pub fn reduced(self) -> Self {
        let mut result = self;
        result.reduce_root();
        result.reduce();
        result
    }
}

impl<
        T: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
        U: QuadraticBase,
    > From<ContinuedFraction<T>> for QuadraticSurd<U>
where
    for<'r> &'r T: RefNum<T>,
    for<'r> &'r U: RefNum<U>,
{
    fn from(s: ContinuedFraction<T>) -> Self {
        // use convergent if it's a rational
        if s.is_rational() {
            return Self::from(s.convergents().last().unwrap());
        }

        // for periodic fraction, assume periodic part x = (ax + b) / (cx + d)
        let mut piter = s.periodic_coeffs().iter().rev();
        let (mut a, mut b) = (T::zero(), T::one());
        let (mut c, mut d) = (T::one(), piter.next().unwrap().clone());

        while let Some(v) = piter.next() {
            let new_c = v * &c + a;
            let new_d = v * &d + b;
            a = c;
            b = d;
            c = new_c;
            d = new_d;
        }

        let psurd = Self::from_equation(
            Ratio::from(c.to_signed()),
            Ratio::from(d.to_signed() - a.to_signed()),
            Ratio::from(-b.to_signed()),
        )
        .unwrap();

        // apply the aperiodic part, assume again result = (ax + b) / (cx + d)
        let surd = if s.aperiodic_coeffs().len() == 1 {
            psurd
        } else {
            let mut aiter = s.aperiodic_coeffs().iter().skip(1).rev();
            let (mut a, mut b) = (T::zero(), T::one());
            let (mut c, mut d) = (T::one(), aiter.next().unwrap().clone());

            while let Some(v) = aiter.next() {
                let new_c = v * &c + a;
                let new_d = v * &d + b;
                a = c;
                b = d;
                c = new_c;
                d = new_d;
            }

            (psurd.clone() * a.to_signed() + b.to_signed())
                / (psurd * c.to_signed() + d.to_signed())
        };
        let surd = surd + s.aperiodic_coeffs().first().unwrap().clone().to_signed();

        // apply sign
        if s.is_negative() {
            -surd
        } else {
            surd
        }
    }
}

// Reference: http://www.numbertheory.org/courses/MP313/lectures/lecture17/page5.html
//            http://www.numbertheory.org/gnubc/surd
// assumes positive surd, parameter `neg` determines the sign of the fraction
fn quadsurd_to_contfrac<T: QuadraticBase + WithUnsigned<Unsigned = U> + AddAssign, U>(
    a: T,
    b: T,
    c: T,
    r: T,
    neg: bool,
) -> ContinuedFraction<U>
where
    for<'r> &'r T: RefNum<T>,
{
    debug_assert!(!c.is_zero() && !r.is_negative());
    debug_assert!(r.sqrt() * r.sqrt() != r);

    // convert to form (p+√d)/q where q|d-p^2
    let mut d = r * &b * &b;
    let (mut p, mut q) = if b.is_negative() { (-a, -c) } else { (a, c) };
    if (&d - &p * &p) % &q != T::zero() {
        d = d * &q * &q;
        p = p * &q;
        q = &q * &q;
    }

    /// Reduction operator for binary quadratic forms
    #[inline]
    fn rho<T: QuadraticBase> (d: &T, rd: &T, p: &mut T, q: &mut T) -> T where
        for<'r> &'r T: RefNum<T>,{
        let a = (rd + &*p).div_floor(&q);
        *p = &a * &*q - &*p;
        *q = (d - &*p * &*p) / &*q; // this step should be exact division
        a
    }

    // find the reduced form and aperiodic coefficients
    let mut a_coeffs: Vec<T> = Vec::new();
    let rd = d.sqrt();
    while a_coeffs.len() == 0 || // ensure that we have a first coefficient
          !(p <= rd && rd < (&p+&q) && (&q-&p) <= rd)
    {
        a_coeffs.push(rho(&d, &rd, &mut p, &mut q));
    }

    // find the periodic coefficients
    let mut p_coeffs: Vec<T> = Vec::new();
    let (init_p, init_q) = (p.clone(), q.clone());
    loop {
        p_coeffs.push(rho(&d, &rd, &mut p, &mut q));
        if p == init_p && q == init_q {
            break;
        }
    }

    ContinuedFraction::new(a_coeffs, p_coeffs, neg)
}

// Although conversion to continued fraction won't fail if complex number is disabled,
// we have to disable From<> to make sure that adding the feature 'complex' won't result
// in break change
impl<T: QuadraticBase + WithUnsigned<Unsigned = U> + AddAssign, U> TryFrom<QuadraticSurd<T>>
    for ContinuedFraction<U>
where
    for<'r> &'r T: RefNum<T>,
{
    type Error = ();

    fn try_from(s: QuadraticSurd<T>) -> Result<Self, ()> {
        if s.discr.is_negative() {
            Err(())
        } else {
            if s.is_negative() {
                Ok(quadsurd_to_contfrac(
                    -s.coeffs.0,
                    -s.coeffs.1,
                    s.coeffs.2,
                    s.discr,
                    true,
                ))
            } else {
                Ok(quadsurd_to_contfrac(
                    s.coeffs.0, s.coeffs.1, s.coeffs.2, s.discr, false,
                ))
            }
        }
    }
}

// TODO: support convert to hurwitz continued fraction based on gaussian integers if `complex` is enabled
//       this should be implemented as From<QuadraticSurd<T>> for HurwitzContinuedFraction (instead of TryFrom)

impl<
        T: QuadraticBase + CheckedAdd + AddAssign + WithUnsigned<Unsigned = U>,
        U: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = T>,
    > Computable<T> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    fn approximated(&self, limit: &T) -> Approximation<Ratio<T>> {
        ContinuedFraction::<U>::try_from(self.clone())
            .expect("only real numbers can be approximated by a rational number")
            .approximated(limit)
    }
}

impl<T: Integer + Signed + fmt::Display + Clone> fmt::Display for QuadraticSurd<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let QuadraticSurdCoeffs(a, b, c) = &self.coeffs;
        if f.alternate() && self.discr.is_negative() {
            // print √-1 as i, √-5 as √5i if alternate flag is set
            // XXX: refactor this function
            let r = -self.discr.clone();
            match (
                a.is_zero(),
                b.is_zero(),
                b.is_one(),
                b == &-T::one(),
                c.is_one(),
                r.is_one(),
            ) {
                (true, true, _, _, _, _) => write!(f, "0"),
                (true, false, true, _, true, true) => write!(f, "i"),
                (true, false, true, _, true, false) => write!(f, "√{}i", r),
                (true, false, true, _, false, true) => write!(f, "i/{}", c),
                (true, false, true, _, false, false) => write!(f, "√{}i/{}", r, c),
                (true, false, false, true, true, true) => write!(f, "-i"),
                (true, false, false, true, true, false) => write!(f, "-√{}i", r),
                (true, false, false, true, false, true) => write!(f, "-i/{}", c),
                (true, false, false, true, false, false) => write!(f, "-√{}i/{}", r, c),
                (true, false, false, false, true, true) => write!(f, "{}i", b),
                (true, false, false, false, true, false) => write!(f, "{}√{}i", b, r),
                (true, false, false, false, false, true) => write!(f, "{}i/{}", b, c),
                (true, false, false, false, false, false) => write!(f, "{}√{}i/{}", b, r, c),
                (false, true, _, _, true, _) => write!(f, "{}", a),
                (false, true, _, _, false, _) => write!(f, "{}/{}", a, c),
                (false, false, true, _, true, true) => write!(f, "{}+i", a),
                (false, false, true, _, true, false) => write!(f, "{}+√{}i", a, r),
                (false, false, false, true, true, true) => write!(f, "{}-i", a),
                (false, false, false, true, true, false) => write!(f, "{}-√{}i", a, r),
                (false, false, false, false, true, true) => {
                    if b.is_negative() {
                        write!(f, "{}{}i", a, b)
                    } else {
                        write!(f, "{}+{}i", a, b)
                    }
                },
                (false, false, false, false, true, false) => {
                    if b.is_negative() {
                        write!(f, "{}{}√{}i", a, b, r)
                    } else {
                        write!(f, "{}+{}√{}i", a, b, r)
                    }
                }
                (false, false, true, _, false, true) => write!(f, "({}+i)/{}", a, c),
                (false, false, true, _, false, false) => write!(f, "({}+√{}i)/{}", a, r, c),
                (false, false, false, true, false, true) => write!(f, "({}-i)/{}", a, c),
                (false, false, false, true, false, false) => write!(f, "({}-√{}i)/{}", a, r, c),
                (false, false, false, false, false, true) => {
                    if b.is_negative() {
                        write!(f, "({}{}i)/{}", a, b, c)
                    } else {
                        write!(f, "({}+{}i)/{}", a, b, c)
                    }
                },
                (false, false, false, false, false, false) => {
                    if b.is_negative() {
                        write!(f, "({}{}√{}i)/{}", a, b, r, c)
                    } else {
                        write!(f, "({}+{}√{}i)/{}", a, b, r, c)
                    }
                }
            }
        } else {
            match (
                a.is_zero(),
                b.is_zero(),
                b.is_one(),
                b == &-T::one(),
                c.is_one(),
            ) {
                (true, true, _, _, _) => write!(f, "0"),
                (true, false, true, _, true) => write!(f, "√{}", self.discr),
                (true, false, true, _, false) => write!(f, "√{}/{}", self.discr, c),
                (true, false, false, true, true) => write!(f, "-√{}", self.discr),
                (true, false, false, true, false) => write!(f, "-√{}/{}", self.discr, c),
                (true, false, false, false, true) => write!(f, "{}√{}", b, self.discr),
                (true, false, false, false, false) => write!(f, "{}√{}/{}", b, self.discr, c),
                (false, true, _, _, true) => write!(f, "{}", a),
                (false, true, _, _, false) => write!(f, "{}/{}", a, c),
                (false, false, true, _, true) => write!(f, "{}+√{}", a, self.discr),
                (false, false, false, true, true) => write!(f, "{}-√{}", a, self.discr),
                (false, false, false, false, true) => {
                    if b.is_negative() {
                        write!(f, "{}{}√{}", a, b, self.discr)
                    } else {
                        write!(f, "{}+{}√{}", a, b, self.discr)
                    }
                }
                (false, false, true, _, false) => write!(f, "({}+√{})/{}", a, self.discr, c),
                (false, false, false, true, false) => {
                    write!(f, "({}-√{})/{}", a, self.discr, c)
                }
                (false, false, false, false, false) => {
                    if b.is_negative() {
                        write!(f, "({}{}√{})/{}", a, b, self.discr, c)
                    } else {
                        write!(f, "({}+{}√{})/{}", a, b, self.discr, c)
                    }
                }
            }
        }
    }
}

// Reduce root base for binary operands. Return None if the two bases
// cannot be matched. This function assumes the root bases on both sides
// are not zero
#[inline]
#[cfg(not(feature = "complex"))]
fn reduce_bin_op<T: QuadraticBase>(
    lhs: QuadraticSurd<T>,
    rhs: QuadraticSurd<T>,
) -> Option<(QuadraticSurd<T>, QuadraticSurd<T>)>
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
    lhs: QuadraticSurd<T>,
    rhs: QuadraticSurd<T>,
) -> Option<(QuadraticSurd<T>, QuadraticSurd<T>)>
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
    lhs: QuadraticSurd<T>,
    rhs: QuadraticSurd<T>,
) -> (QuadraticSurd<T>, QuadraticSurd<T>)
where
    for<'r> &'r T: RefNum<T>,
{
    reduce_bin_op(lhs, rhs).expect("two root bases are not compatible!")
}

macro_rules! forward_binop {
    (impl $imp:ident, $rhs:ty, $method:ident) => {
        impl<T: QuadraticBase> $imp<$rhs> for QuadraticSurd<T>
        where
            for<'r> &'r T: RefNum<T>,
        {
            type Output = Self;
            #[inline]
            fn $method(self, rhs: $rhs) -> Self {
                Self::from_coeffs($imp::$method(self.coeffs, rhs), self.discr)
            }
        }
    };
}

forward_binop!(impl Add, T, add);
forward_binop!(impl Add, Ratio<T>, add);
forward_binop!(impl Sub, T, sub);
forward_binop!(impl Sub, Ratio<T>, sub);
forward_binop!(impl Mul, T, mul);
forward_binop!(impl Mul, Ratio<T>, mul);
forward_binop!(impl Div, T, div);
forward_binop!(impl Div, Ratio<T>, div);

// TODO: implement checked_add, checked_sub, checked_mul, checked_div, checking both overflow and discr dismatch
impl<T: QuadraticBase> Add for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        // shortcuts for trivial cases
        if rhs.is_rational() {
            return self + Ratio::<T>::new(rhs.coeffs.0, rhs.coeffs.2);
        }
        if self.is_rational() {
            return rhs + Ratio::<T>::new(self.coeffs.0, self.coeffs.2);
        }

        // ensure that two operands have compatible bases
        let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
        Self::from_coeffs(lhs.coeffs + rhs.coeffs, lhs.discr)
    }
}

impl<T: QuadraticBase> Sub for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        // shortcuts for trivial cases
        if rhs.is_rational() {
            return self - Ratio::<T>::new(rhs.coeffs.0, rhs.coeffs.2);
        }
        if self.is_rational() {
            return -(rhs - Ratio::<T>::new(self.coeffs.0, self.coeffs.2));
        }

        // ensure that two operands have compatible bases
        let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
        Self::from_coeffs(lhs.coeffs - rhs.coeffs, lhs.discr)
    }
}

impl<T: QuadraticBase> Mul for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        // shortcuts for trivial cases
        if rhs.is_rational() {
            return self * Ratio::<T>::new(rhs.coeffs.0, rhs.coeffs.2);
        }
        if self.is_rational() {
            return rhs * Ratio::<T>::new(self.coeffs.0, self.coeffs.2);
        }
        if self.is_pure() && rhs.is_pure() {
            let gcd = self.discr.gcd(&rhs.discr);
            let discr = (self.discr / &gcd) * (rhs.discr / &gcd);
            return Self::new(
                T::zero(),
                self.coeffs.1 * rhs.coeffs.1 * gcd,
                self.coeffs.2 * rhs.coeffs.2,
                discr,
            );
        }

        // ensure that two operands have compatible bases
        let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
        Self::from_coeffs(
            QuadraticOps::mul(lhs.coeffs, rhs.coeffs, &lhs.discr),
            lhs.discr,
        )
    }
}

impl<T: QuadraticBase> Div for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        // shortcuts for trivial cases
        if rhs.is_rational() {
            return self / Ratio::<T>::new(rhs.coeffs.0, rhs.coeffs.2);
        }
        if self.is_rational() {
            return (rhs / Ratio::<T>::new(self.coeffs.0, self.coeffs.2)).recip();
        }
        if self.is_pure() && rhs.is_pure() {
            let gcd = self.discr.gcd(&rhs.discr);
            let (ld, rd) = (self.discr / &gcd, rhs.discr / gcd);
            return Self::new(
                T::zero(),
                self.coeffs.1 * rhs.coeffs.2,
                self.coeffs.2 * rhs.coeffs.1 * &rd,
                ld * rd,
            );
        }

        // ensure that two operands have compatible bases
        let (lhs, rhs) = reduce_bin_op_unwrap(self, rhs);
        Self::from_coeffs(
            QuadraticOps::div(lhs.coeffs, rhs.coeffs, &lhs.discr),
            lhs.discr,
        )
    }
}

impl<T: QuadraticBase> Neg for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = QuadraticSurd<T>;
    #[inline]
    fn neg(self) -> QuadraticSurd<T> {
        QuadraticSurd::new_raw(-self.coeffs.0, -self.coeffs.1, self.coeffs.2, self.discr)
    }
}

impl<T: QuadraticBase> Zero for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    #[inline]
    fn zero() -> Self {
        Self {
            coeffs: QuadraticSurdCoeffs(T::zero(), T::zero(), T::one()),
            discr: T::zero(),
        }
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.coeffs.0.is_zero() && self.coeffs.1.is_zero()
    }
}

impl<T: QuadraticBase> One for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    #[inline]
    fn one() -> Self {
        Self {
            coeffs: QuadraticSurdCoeffs(T::one(), T::zero(), T::one()),
            discr: T::zero(),
        }
    }
    #[inline]
    fn is_one(&self) -> bool {
        self.coeffs.0.is_one() && self.coeffs.1.is_zero() && self.coeffs.2.is_one()
    }
}

impl<T: Integer + FromPrimitive + Clone> FromPrimitive for QuadraticSurd<T> {
    #[inline]
    fn from_i64(n: i64) -> Option<Self> {
        T::from_i64(n).map(Self::from)
    }

    #[inline]
    fn from_u64(n: u64) -> Option<Self> {
        T::from_u64(n).map(Self::from)
    }

    /// This method depends on [Ratio::from_f64] to convert from a float
    #[inline]
    fn from_f64(f: f64) -> Option<Self> {
        // XXX: this method should be improved when Ratio has a more reasonable API for approximating float
        let frac = Ratio::<i64>::from_f64(f)?;
        let (n, d) = frac.into();
        let frac = Ratio::new(T::from_i64(n)?, T::from_i64(d)?);
        Some(Self::from(frac))
    }
}

impl<T: QuadraticBase + ToPrimitive> ToPrimitive for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    #[inline]
    fn to_i64(&self) -> Option<i64> {
        if self.discr.is_negative() {
            return None;
        }
        match self.to_integer() {
            Approximation::Exact(v) => v.to_i64(),
            Approximation::Approximated(_) => None,
        }
    }

    #[inline]
    fn to_u64(&self) -> Option<u64> {
        if self.discr.is_negative() {
            return None;
        }
        match self.to_integer() {
            Approximation::Exact(v) => v.to_u64(),
            Approximation::Approximated(_) => None,
        }
    }

    #[inline]
    fn to_f64(&self) -> Option<f64> {
        if self.discr < T::zero() {
            return None;
        }
        Some(
            (self.coeffs.0.to_f64()? + self.coeffs.1.to_f64()? * self.discr.to_f64()?.sqrt())
                / self.coeffs.2.to_f64()?,
        )
    }
}

#[cfg(feature = "num-complex")]
mod complex {
    use super::*;
    use num_complex::{Complex32, Complex64};

    impl<T: QuadraticBase + ToPrimitive> QuadraticSurd<T>
    where
        for<'r> &'r T: RefNum<T>,
    {
        pub fn to_complex64(&self) -> Option<Complex64> {
            if self.discr < T::zero() {
                let c = self.coeffs.2.to_f64()?;
                let re = self.coeffs.0.to_f64()? / c;
                let im = self.coeffs.1.to_f64()? * self.discr.to_f64()?.abs().sqrt() / c;
                Some(Complex64::new(re as f64, im as f64))
            } else {
                let re = self.to_f64()?;
                Some(Complex64::new(re as f64, 0f64))
            }
        }

        pub fn to_complex32(&self) -> Option<Complex32> {
            let complex = self.to_complex64()?;
            Some(Complex32::new(complex.re.to_f32()?, complex.im.to_f32()?))
        }
    }
}

impl<T: QuadraticBase> FromSqrt<T> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    #[inline]
    fn from_sqrt(target: T) -> Result<Self, FromSqrtError<T>> {
        #[cfg(not(feature = "complex"))]
        if target.is_negative() {
            return Err(FromSqrtError {
                data: target,
                kind: SqrtErrorKind::Complex,
            });
        }
        Ok(QuadraticSurd::new(T::zero(), T::one(), T::one(), target))
    }
}

impl<T: QuadraticBase + CheckedMul> FromSqrt<Ratio<T>> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    #[inline]
    fn from_sqrt(target: Ratio<T>) -> Result<Self, FromSqrtError<Ratio<T>>> {
        #[cfg(not(feature = "complex"))]
        if target.is_negative() {
            return Err(FromSqrtError {
                data: target,
                kind: SqrtErrorKind::Complex,
            });
        }
        if target.is_integer() {
            let (num, _) = target.into();
            return Ok(QuadraticSurd::new(T::zero(), T::one(), T::one(), num));
        }

        match target.numer().checked_mul(target.denom()) {
            Some(new_r) => Ok(QuadraticSurd::new(
                T::zero(),
                T::one(),
                target.denom().clone(),
                new_r,
            )),
            None => Err(FromSqrtError {
                data: target,
                kind: SqrtErrorKind::Overflow,
            }),
        }
    }
}

impl<T: QuadraticBase + CheckedMul> FromSqrt<QuadraticSurd<T>> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    /// The real part of the output is ensured to be non-negative (if successful).
    #[inline]
    fn from_sqrt(target: QuadraticSurd<T>) -> Result<Self, FromSqrtError<QuadraticSurd<T>>> {
        #[cfg(not(feature = "complex"))]
        if target.is_negative() {
            return Err(FromSqrtError {
                data: target,
                kind: SqrtErrorKind::Complex,
            });
        }

        if target.is_rational() {
            match target.to_rational() {
                Approximation::Exact(v) => {
                    return Self::from_sqrt(v).map_err(|e| FromSqrtError {
                        data: Self::from(e.data),
                        kind: e.kind,
                    })
                }
                _ => unreachable!(),
            };
        }

        // denote a_c = a/c, b_c = b/c
        // suppose (a_c + b_c * sqrt(r))^2 = target = x + y * sqrt(r)
        // let g = a_c/b_c, then g^2 - 2g(x/y) + r = 0
        // result is available only if g is rational

        let QuadraticSurdCoeffs(a, b, c) = target.coeffs;
        let x = Ratio::new(a, c.clone());
        let y = Ratio::new(b, c);
        let x_y = x / &y;
        let delta2 = &x_y * &x_y - &target.discr;

        // reconstruct the original target when error occurred
        #[inline]
        fn reconstruct<T: QuadraticBase>(x_y: Ratio<T>, y: Ratio<T>, r: T, kind: SqrtErrorKind) -> FromSqrtError<QuadraticSurd<T>> where
        for<'r> &'r T: RefNum<T> { 
            FromSqrtError {
                data: QuadraticSurd::from_rationals(x_y * &y, y, Ratio::from(r)),
                kind: kind,
            }
        }

        if delta2.is_negative() {
            // this branch happens only when discr is positive
            return Err(reconstruct(x_y, y, target.discr, SqrtErrorKind::Unrepresentable));
        }
        let delta = match Self::from_sqrt(delta2) {
            Ok(v) => v,
            Err(e) => return Err(reconstruct(x_y, y, target.discr, e.kind))
        };
        let delta = if delta.is_rational() {
            delta.to_rational().value()
        } else {
            return Err(reconstruct(x_y, y, target.discr, SqrtErrorKind::Unrepresentable))
        };

        // from the equation above, y = 2*a_c*b_c => c = sqrt(2ab/y)
        // this function find a tuple of integers (a, b, c) that satisfies the equation, and
        // the output a will be ensured to be non-negative
        fn get_abc<T: QuadraticBase> (g: Ratio<T>, y: &Ratio<T>) -> Option<(T, T, T)> where
        for<'r> &'r T: RefNum<T> {
            let two_ab = (T::one() + T::one()) * g.numer() * g.denom();
            let d = two_ab.gcd(y.numer());
            let scale = y.numer() / d;
            let c2 = two_ab * &scale * &scale * y.denom() / y.numer();
            let (a, b) = (g.numer() * &scale, g.denom() * scale);
            if c2.is_negative() {
                return None;
            }
            let c = c2.sqrt();
            if &c * &c != c2 {
                None
            } else {
                Some((a, b, c))
            }
        }

        let ret = if let Some((a, b, c)) = get_abc(&x_y - &delta, &y) {
            QuadraticSurd::new(a, b, c, target.discr)
        } else if let Some((a, b, c)) = get_abc(&x_y + delta, &y) {
            QuadraticSurd::new(a, b, c, target.discr)
        } else {
            return Err(reconstruct(x_y, y, target.discr, SqrtErrorKind::Unrepresentable))
        };

        debug_assert!(!ret.coeffs.0.is_negative());
        Ok(ret)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::convert::TryFrom;

    pub const PHI: QuadraticSurd<i32> = QuadraticSurd::new_raw(1, 1, 2, 5); // 1.618
    pub const PHI_R: QuadraticSurd<i32> = QuadraticSurd::new_raw(-1, 1, 2, 5); // 0.618
    pub const N_PHI: QuadraticSurd<i32> = QuadraticSurd::new_raw(-1, -1, 2, 5); // -1.618
    pub const N_PHI_R: QuadraticSurd<i32> = QuadraticSurd::new_raw(1, -1, 2, 5); // -0.618
    pub const PHI_SQ: QuadraticSurd<i32> = QuadraticSurd::new_raw(3, 1, 2, 5);

    pub const PHI45: QuadraticSurd<i32> = QuadraticSurd::new_raw(3, 1, 6, 45); // non-reduced version of PHI
    pub const PHI45_R: QuadraticSurd<i32> = QuadraticSurd::new_raw(-3, 1, 6, 45); // non-reduced version of PHI_R

    pub const SQ5: QuadraticSurd<i32> = QuadraticSurd::new_raw(0, 1, 1, 5);
    pub const N_SQ5: QuadraticSurd<i32> = QuadraticSurd::new_raw(0, -1, 1, 5);

    #[test]
    fn creation_test() {
        let coeffs: (i32, i32, i32, i32) = QuadraticSurd::new(2, 6, 2, 5).into();
        assert_eq!(coeffs, (1, 3, 1, 5)); // reduce common divisors
        
        let coeffs: (i32, i32, i32, i32) = QuadraticSurd::new(2, 1, 2, 18).into();
        assert_eq!(coeffs, (2, 1, 2, 18)); // 18 is not trivially reducible

        let coeffs: (i32, i32, i32, i32) = QuadraticSurd::new(3, 1, 3, 18).into();
        assert_eq!(coeffs, (1, 1, 1, 2)); // 18 is reducible with the help of gcd hint

        let coeffs: (i32, i32, i32, i32) = QuadraticSurd::new(3, 1, 3, 9).into();
        assert_eq!(coeffs, (2, 0, 1, 0)); // 9 is a square number
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
        // from_rationals
        assert_eq!(
            QuadraticSurd::from_rationals(Ratio::new(1, 2), Ratio::new(1, 2), Ratio::from(5)),
            PHI
        );
        assert_eq!(
            QuadraticSurd::from_rationals(Ratio::new(-1, 2), Ratio::new(1, 2), Ratio::from(5)),
            PHI_R
        );

        // from_sqrt
        assert_eq!(
            QuadraticSurd::from_sqrt(5).unwrap(),
            QuadraticSurd::new_raw(0, 1, 1, 5)
        );
        assert!(matches!(QuadraticSurd::from_sqrt(PHI), Err(_)));
        assert!(matches!(QuadraticSurd::from_sqrt(PHI * PHI), Ok(v) if v == PHI));

        assert!(matches!(
            QuadraticSurd::from_equation(Ratio::from(1), Ratio::from(-1), Ratio::from(-1)),
            Some(v) if v == PHI
        ));
    }

    #[test]
    fn property_test() {
        // recip
        assert_eq!(PHI.recip(), PHI_R);
        assert_eq!(N_PHI.recip(), N_PHI_R);

        // conjugate
        assert_eq!(PHI.conj(), N_PHI_R);
        assert_eq!(PHI_R.conj(), N_PHI);

        // is_pure
        assert_eq!(PHI.is_pure(), false);
        assert_eq!(SQ5.is_pure(), true);
    }

    #[test]
    fn arithmic_test() {
        // add
        assert_eq!(PHI_R + 1, PHI);
        assert_eq!(PHI45_R + 1, PHI45);
        assert_eq!(N_PHI + 1, N_PHI_R);
        assert_eq!(PHI_R + Ratio::one(), PHI);
        assert_eq!(PHI45_R + Ratio::one(), PHI45);
        assert_eq!(N_PHI + Ratio::one(), N_PHI_R);
        assert_eq!(PHI + PHI_R, SQ5);
        assert_eq!(N_PHI + N_PHI_R, N_SQ5);
        assert!((PHI + N_PHI).is_zero());
        assert!((PHI + N_PHI_R).is_one());

        // sub
        assert_eq!(PHI - 1, PHI_R);
        assert_eq!(PHI45 - 1, PHI45_R);
        assert_eq!(N_PHI_R - 1, N_PHI);
        assert_eq!(PHI - Ratio::one(), PHI_R);
        assert_eq!(PHI45 - Ratio::one(), PHI45_R);
        assert_eq!(N_PHI_R - Ratio::one(), N_PHI);
        assert!((PHI - PHI).is_zero());
        assert!((PHI - PHI_R).is_one());
        assert!((N_PHI_R - N_PHI).is_one());
        assert_eq!(PHI - N_PHI_R, SQ5);
        assert_eq!(N_PHI - PHI_R, -SQ5);

        // mul
        assert_eq!(PHI * 1, PHI);
        assert_eq!(PHI * -1, N_PHI);
        assert_eq!(PHI * Ratio::one(), PHI);
        assert_eq!(PHI * -Ratio::one(), N_PHI);
        assert!((PHI * PHI_R).is_one());
        assert!((N_PHI * N_PHI_R).is_one());
        assert_eq!(PHI * PHI, PHI_SQ);
        assert_eq!(N_PHI * N_PHI, PHI_SQ);

        // div
        assert_eq!(PHI / 1, PHI);
        assert_eq!(PHI / -1, N_PHI);
        assert_eq!(PHI / Ratio::one(), PHI);
        assert_eq!(PHI / -Ratio::one(), N_PHI);
        assert!((PHI / PHI).is_one());
        assert!((N_PHI / N_PHI).is_one());
        assert_eq!(PHI / PHI_R, PHI_SQ);
        assert_eq!(N_PHI / N_PHI_R, PHI_SQ);

        // associativity test
        let three_half = Ratio::new(3, 2);
        assert_eq!(PHI + 5 - PHI, QuadraticSurd::from(5));
        assert_eq!(PHI + three_half - PHI, QuadraticSurd::from(three_half));
        assert_eq!(PHI + SQ5 - PHI, SQ5);
        assert_eq!(PHI + 5 - PHI45, QuadraticSurd::from(5));
        assert_eq!(PHI + three_half - PHI45, QuadraticSurd::from(three_half));
        assert_eq!(PHI * 5 / PHI45, QuadraticSurd::from(5));
        assert_eq!(PHI * three_half / PHI45, QuadraticSurd::from(three_half));

        // mixed
        assert_eq!(PHI * 2 - 1, SQ5);
        assert_eq!(PHI_R * 2 + 1, SQ5);
        assert_eq!(PHI / Ratio::new(1, 2) - 1, SQ5);
        assert_eq!((PHI - Ratio::new(1, 2)) * 2, SQ5);
    }

    #[test]
    fn arithmic_test_diff_base() {
        assert_eq!(PHI45 + PHI_R, SQ5);
        assert_eq!(PHI + PHI45_R, SQ5);
        assert!((PHI - PHI45).is_zero());
        assert!((PHI - PHI45_R).is_one());
        assert!((PHI45 * PHI_R).is_one());
        assert!((PHI * PHI45_R).is_one());
        assert!((PHI45 / PHI).is_one());
        assert!((PHI / PHI45).is_one());
        assert_eq!(PHI + SQ5 - PHI45, SQ5);
        assert_eq!(PHI * SQ5 / PHI45, SQ5);

        let a = QuadraticSurd::new(0, 2, 3, 2);
        let b = QuadraticSurd::new(0, 3, 4, 3);
        let c = QuadraticSurd::new(0, 1, 2, 6);
        assert_eq!(a * b, c);
        assert_eq!(c / a, b);
        assert_eq!(c / b, a);
    }

    #[test]
    fn conversion_between_contfrac() {
        let cf_phi = ContinuedFraction::<u32>::new(vec![1i32], vec![1], false);
        assert_eq!(QuadraticSurd::from(cf_phi.clone()), PHI);
        assert_eq!(ContinuedFraction::try_from(PHI).unwrap(), cf_phi);

        let cf_sq2 = ContinuedFraction::<u32>::new(vec![1i32], vec![2], false);
        let surd_sq2 = QuadraticSurd::new(0, 1, 1, 2);
        assert_eq!(QuadraticSurd::from(cf_sq2.clone()), surd_sq2);
        assert_eq!(ContinuedFraction::try_from(surd_sq2).unwrap(), cf_sq2);

        let cf_n_sq2 = ContinuedFraction::<u32>::new(vec![1i32], vec![2], true);
        let surd_n_sq2 = QuadraticSurd::new(0, -1, 1, 2);
        assert_eq!(QuadraticSurd::from(cf_n_sq2.clone()), surd_n_sq2);
        assert_eq!(ContinuedFraction::try_from(surd_n_sq2).unwrap(), cf_n_sq2);

        let cf_sq2_7 = ContinuedFraction::<u32>::new(vec![0i32, 4], vec![1, 18, 1, 8], false);
        let surd_sq2_7 = QuadraticSurd::new(0, 1, 7, 2);
        assert_eq!(QuadraticSurd::from(cf_sq2_7.clone()), surd_sq2_7);
        assert_eq!(ContinuedFraction::try_from(surd_sq2_7).unwrap(), cf_sq2_7);

        let cf_10_sq2_7 = ContinuedFraction::<u32>::new(vec![1i32, 4], vec![2], false);
        let surd_10_sq2_7 = QuadraticSurd::new(10, -1, 7, 2);
        assert_eq!(QuadraticSurd::from(cf_10_sq2_7.clone()), surd_10_sq2_7);
        assert_eq!(
            ContinuedFraction::try_from(surd_10_sq2_7).unwrap(),
            cf_10_sq2_7
        );
    }

    #[test]
    fn formatting_test() {
        assert_eq!(format!("{}", QuadraticSurd::<i8>::zero()), "0");
        assert_eq!(format!("{}", QuadraticSurd::<i8>::one()), "1");
        assert_eq!(format!("{}", QuadraticSurd::from_sqrt(5).unwrap()), "√5");

        assert_eq!(format!("{}", QuadraticSurd::new(1, 1, 1, 5)), "1+√5");
        assert_eq!(format!("{}", QuadraticSurd::new(1, -1, 1, 5)), "1-√5");
        assert_eq!(format!("{}", QuadraticSurd::new(-1, 1, 1, 5)), "-1+√5");
        assert_eq!(format!("{}", QuadraticSurd::new(-1, -1, 1, 5)), "-1-√5");
        assert_eq!(format!("{}", QuadraticSurd::new(1, 2, 1, 5)), "1+2√5");
        assert_eq!(format!("{}", QuadraticSurd::new(1, -2, 1, 5)), "1-2√5");
        assert_eq!(format!("{}", QuadraticSurd::new(-1, 2, 1, 5)), "-1+2√5");
        assert_eq!(format!("{}", QuadraticSurd::new(-1, -2, 1, 5)), "-1-2√5");
        assert_eq!(format!("{}", QuadraticSurd::new(1, 0, 2, 5)), "1/2");
        assert_eq!(format!("{}", QuadraticSurd::new(-1, 0, 2, 5)), "-1/2");
        assert_eq!(format!("{}", QuadraticSurd::new(1, 1, 2, 5)), "(1+√5)/2");
        assert_eq!(format!("{}", QuadraticSurd::new(1, -1, 2, 5)), "(1-√5)/2");
        assert_eq!(format!("{}", QuadraticSurd::new(-1, 1, 2, 5)), "(-1+√5)/2");
        assert_eq!(format!("{}", QuadraticSurd::new(-1, -1, 2, 5)), "(-1-√5)/2");
        assert_eq!(format!("{}", QuadraticSurd::new(1, 2, 2, 5)), "(1+2√5)/2");
        assert_eq!(format!("{}", QuadraticSurd::new(1, -2, 2, 5)), "(1-2√5)/2");
        assert_eq!(format!("{}", QuadraticSurd::new(-1, 2, 2, 5)), "(-1+2√5)/2");
        assert_eq!(
            format!("{}", QuadraticSurd::new(-1, -2, 2, 5)),
            "(-1-2√5)/2"
        );
    }

    #[test]
    fn trunc_frac_test() {
        assert_eq!(PHI.trunc_ref(), QuadraticSurd::from(1));
        assert_eq!(PHI.fract_ref(), PHI_R);
        assert_eq!(PHI_R.trunc_ref(), QuadraticSurd::zero());
        assert_eq!(PHI_R.fract_ref(), PHI_R);
        assert_eq!(N_PHI.trunc_ref(), QuadraticSurd::from(-1));
        assert_eq!(N_PHI.fract_ref(), N_PHI_R);
        assert_eq!(N_PHI_R.trunc_ref(), QuadraticSurd::from(0));
        assert_eq!(N_PHI_R.fract_ref(), N_PHI_R);
        assert_eq!(PHI_SQ.trunc_ref(), QuadraticSurd::from(2));
        assert_eq!(PHI_SQ.fract_ref(), PHI_R);
        assert_eq!(PHI45.trunc_ref(), QuadraticSurd::from(1));
        assert_eq!(PHI45.fract_ref(), PHI_R);
        assert_eq!(PHI45_R.trunc_ref(), QuadraticSurd::from(0));
        assert_eq!(PHI45_R.fract_ref(), PHI_R);
        assert_eq!(SQ5.trunc_ref(), QuadraticSurd::from(2));
        assert_eq!(SQ5.fract_ref(), QuadraticSurd::new(-2, 1, 1, 5));
        assert_eq!(N_SQ5.trunc_ref(), QuadraticSurd::from(-2));
        assert_eq!(N_SQ5.fract_ref(), QuadraticSurd::new(2, -1, 1, 5));
    }

    #[test]
    fn from_sqrt_test() {
        assert_eq!(
            QuadraticSurd::from_sqrt(5i32).unwrap(),
            QuadraticSurd::new_raw(0, 1, 1, 5)
        );
        assert_eq!(
            QuadraticSurd::from_sqrt(QuadraticSurd::new(3i32, 2, 1, 2)).unwrap(),
            QuadraticSurd::new_raw(1, 1, 1, 2)
        );

        #[cfg(not(feature = "complex"))]
        {
            let err = QuadraticSurd::from_sqrt(-2i32).unwrap_err();
            assert_eq!(
                err.kind,
                SqrtErrorKind::Complex
            );
            assert_eq!(err.data, -2);
            let err = QuadraticSurd::from_sqrt(Ratio::new(-1i32, 2)).unwrap_err();
            assert_eq!(
                err.kind,
                SqrtErrorKind::Complex
            );
            assert_eq!(err.data, Ratio::new(-1, 2));
            let surd = QuadraticSurd::new(-2i32, 1, 1, 2);
            let err = QuadraticSurd::from_sqrt(surd).unwrap_err();
            assert_eq!(
                err.kind,
                SqrtErrorKind::Complex
            );
            assert_eq!(err.data, surd);
        }
    }
}

#[cfg(feature = "complex")]
#[cfg(test)]
mod complex_tests {
    use super::*;

    pub const OMEGA: QuadraticSurd<i32> = QuadraticSurd::new_raw(-1, 1, 2, -3); // -0.5 + 0.866i
    pub const OMEGA3: QuadraticSurd<i32> = QuadraticSurd::new_raw(-3, 3, 2, -3); // -1.5 + 2.598i

    #[test]
    fn trunc_frac_test() {
        assert_eq!(OMEGA.trunc_ref(), QuadraticSurd::zero());
        assert_eq!(OMEGA3.trunc_ref(), QuadraticSurd::new(-1, 2, 1, -1));
    }

    #[test]
    fn from_sqrt_test() {
        assert!(QuadraticSurd::from_sqrt(-2i32).is_ok());
        assert!(QuadraticSurd::from_sqrt(Ratio::new(-1i32, 2)).is_ok());
        assert!(QuadraticSurd::from_sqrt(QuadraticSurd::from(-1i32)).is_ok());

        let sq2i = QuadraticSurd::from_sqrt(-2i32).unwrap();
        let sqhalfi = QuadraticSurd::from_sqrt(Ratio::new(-1i32, 2)).unwrap();
        let i = QuadraticSurd::from_sqrt(QuadraticSurd::from(-1i32)).unwrap();
        assert_eq!(sqhalfi * 2, sq2i);
        assert_eq!(sq2i * i, QuadraticSurd::from_sqrt(2i32).unwrap());
        
        assert_eq!(
            QuadraticSurd::from_sqrt(QuadraticSurd::new(-3i32, 4, 1, -1)).unwrap(),
            QuadraticSurd::new_raw(1, 2, 1, -1)
        );
    }

    #[test]
    fn formatting_test() {
        // test trivial complex cases
        assert_eq!(format!("{:#}", QuadraticSurd::from_sqrt(-1).unwrap()), "i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 1, 1, -1)), "1+i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, -1, 1, -1)), "1-i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, 1, 1, -1)), "-1+i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, -1, 1, -1)), "-1-i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 2, 1, -1)), "1+2i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, -2, 1, -1)), "1-2i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, 2, 1, -1)), "-1+2i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, -2, 1, -1)), "-1-2i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 0, 2, -1)), "1/2");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, 0, 2, -1)), "-1/2");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 1, 2, -1)), "(1+i)/2");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, -1, 2, -1)), "(1-i)/2");
        assert_eq!(
            format!("{:#}", QuadraticSurd::new(-1, 1, 2, -1)),
            "(-1+i)/2"
        );
        assert_eq!(
            format!("{:#}", QuadraticSurd::new(-1, -1, 2, -1)),
            "(-1-i)/2"
        );
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 2, 2, -1)), "(1+2i)/2");
        assert_eq!(
            format!("{:#}", QuadraticSurd::new(1, -2, 2, -1)),
            "(1-2i)/2"
        );
        assert_eq!(
            format!("{:#}", QuadraticSurd::new(-1, 2, 2, -1)),
            "(-1+2i)/2"
        );

        // test non-trivial complex cases
        assert_eq!(format!("{:#}", QuadraticSurd::from_sqrt(-3).unwrap()), "√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 1, 1, -3)), "1+√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, -1, 1, -3)), "1-√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, 1, 1, -3)), "-1+√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, -1, 1, -3)), "-1-√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 2, 1, -3)), "1+2√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, -2, 1, -3)), "1-2√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, 2, 1, -3)), "-1+2√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, -2, 1, -3)), "-1-2√3i");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 0, 2, -3)), "1/2");
        assert_eq!(format!("{:#}", QuadraticSurd::new(-1, 0, 2, -3)), "-1/2");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 1, 2, -3)), "(1+√3i)/2");
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, -1, 2, -3)), "(1-√3i)/2");
        assert_eq!(
            format!("{:#}", QuadraticSurd::new(-1, 1, 2, -3)),
            "(-1+√3i)/2"
        );
        assert_eq!(
            format!("{:#}", QuadraticSurd::new(-1, -1, 2, -3)),
            "(-1-√3i)/2"
        );
        assert_eq!(format!("{:#}", QuadraticSurd::new(1, 2, 2, -3)), "(1+2√3i)/2");
        assert_eq!(
            format!("{:#}", QuadraticSurd::new(1, -2, 2, -3)),
            "(1-2√3i)/2"
        );
        assert_eq!(
            format!("{:#}", QuadraticSurd::new(-1, 2, 2, -3)),
            "(-1+2√3i)/2"
        );
    }
}
