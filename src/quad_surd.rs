use crate::cont_frac::ContinuedFraction;
use crate::traits::{Approximation, Computable, FromSqrt, WithSigned, WithUnsigned};
use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub};
use num_integer::{sqrt, Integer, Roots};
use num_traits::{
    CheckedAdd, CheckedMul, FromPrimitive, NumRef, One, RefNum, Signed, ToPrimitive, Zero,
};
use std::fmt;

use num_rational::Ratio;

/// A helper trait to define valid type that can be used for QuadraticSurd
pub trait QuadraticSurdBase: Integer + NumRef + Clone + Roots + Signed {}
impl<T: Integer + NumRef + Clone + Roots + Signed> QuadraticSurdBase for T {}

/// A type representation quadratic surd number `(a + b*sqrt(r)) / c`.
/// If the support for complex number is enabled, then this struct can represent
/// any quadratic integers.
#[derive(PartialEq, Eq, Hash, Clone, Debug, Copy)]
pub struct QuadraticSurd<T> {
    a: T,
    b: T, // zero when reduced if the surd is a rational number
    c: T, // positive when reduced
    r: T, // zero when reduced if the surd is a rational number
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
            a: t,
            b: T::zero(),
            c: T::one(),
            r: T::zero(),
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
            a,
            b: T::zero(),
            c,
            r: T::zero(),
        }
    }
}

impl<T: QuadraticSurdBase> QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    // TODO: it's possible to support r < 0, but special handling is needed.
    //       And then we can support from_complex, to_complex, is_complex

    fn reduce(&mut self) {
        if self.c.is_zero() {
            panic!("denominator is zero");
        }

        // test if the surd is rational
        let root = sqrt(self.r.clone());
        if &root * &root == self.r {
            self.a = &self.a + &self.b * root;
            self.b = T::zero();
            self.r = T::zero();
        }

        // reduce common divisor
        let mut g_ac = self.a.gcd(&self.c);
        let g_acr = (&g_ac * &g_ac).gcd(&self.r); // test if squared factor of c can divide r
        let groot = sqrt(g_acr.clone());
        if &groot * &groot == g_acr {
            self.r = &self.r / &g_acr;
            self.a = &self.a / &groot;
            self.c = &self.c / &groot;
            g_ac = &g_ac / groot;
        }

        let g_abc = g_ac.gcd(&self.b);
        self.a = &self.a / &g_abc;
        self.b = &self.b / &g_abc;
        self.c = &self.c / g_abc;

        // keep denom positive
        if self.c < T::zero() {
            self.a = T::zero() - &self.a;
            self.b = T::zero() - &self.b;
            self.c = T::zero() - &self.c;
        }
    }

    fn reduce_root_hinted(self, hint: T) -> Self {
        let hint = hint.abs();
        if hint.is_zero() || hint.is_one() {
            return self;
        }

        let (quo, rem) = self.r.div_rem(&hint);
        if rem.is_zero() {
            // if hint is actually a factor
            let root = sqrt(hint.clone());
            if &root * &root == hint {
                // if hint is a square number
                let g = root.gcd(&self.a).gcd(&self.c);
                return QuadraticSurd::new_raw(self.a / &g, self.b * root / &g, self.c / g, quo);
            }
        }

        self
    }

    /// Create a surd represented as `(a + b*sqrt(r)) / c` where `a`, `b`, `c`, `r` are integers.
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

    /// Create a surd represented as `a + b sqrt(r)` where a, b, r are rationals.
    /// 
    /// # Panics
    /// If `r` is negative when the `complex` feature is not enabled.
    #[inline]
    pub fn from_rational(a: Ratio<T>, b: Ratio<T>, r: Ratio<T>) -> Self {
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
    /// This method only returns one of the root `(b + sqrt(b^2 - 4ac)) / 2a`, use `conjugate()` for the other root
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
        Some(Self::from_rational(-b * aa.recip(), aa.recip(), delta))
    }

    /// Returns the reciprocal of the surd
    #[inline]
    pub fn recip(self) -> Self {
        let aa = &self.a * &self.a;
        let bb = &self.b * &self.b;
        QuadraticSurd::new(
            -(&self.c * self.a),
            self.c * self.b,
            bb * &self.r - aa,
            self.r,
        )
    }

    /// `.recip()` with reference
    #[inline]
    pub fn recip_ref(&self) -> Self {
        let aa = &self.a * &self.a;
        let bb = &self.b * &self.b;
        QuadraticSurd::new(
            -(&self.c * &self.a),
            &self.c * &self.b,
            bb * &self.r - aa,
            self.r.clone(),
        )
    }

    /// Return the conjugate of the surd, i.e. (a - b*sqrt(r)) / c
    #[inline]
    pub fn conj(self) -> Self {
        QuadraticSurd {
            a: self.a,
            b: -self.b,
            c: self.c,
            r: self.r,
        }
    }

    /// `.conj()` with reference
    #[inline]
    pub fn conj_ref(&self) -> Self {
        self.clone().conj()
    }

    pub fn trunc(self) -> Self {
        self.panic_if_complex(); // TODO: support complex

        let br = sqrt(&self.b * &self.b * &self.r);
        let num = if self.b >= T::zero() {
            &self.a + br
        } else {
            &self.a - br
        };
        return Self::from(num / &self.c);
    }

    pub fn trunc_ref(&self) -> Self {
        self.clone().trunc()
    }

    pub fn fract(self) -> Self {
        let trunc = self.trunc_ref();
        self - trunc
    }

    pub fn fract_ref(&self) -> Self {
        self.clone() - self.trunc_ref()
    }

    pub fn numer(&self) -> Self {
        Self::new_raw(self.a.clone(), self.b.clone(), T::one(), self.r.clone())
    }

    pub fn denom(&self) -> &T {
        &self.c
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
                &self.a + sqrt(&self.b * &self.b * &self.r),
                self.c.clone(),
            ))
        }
    }

    /// Rounds towards minus infinity
    #[inline]
    pub fn floor(&self) -> Self {
        let br = sqrt(&self.b * &self.b * &self.r);
        let num = if self.b >= T::zero() {
            &self.a + br
        } else {
            &self.a - br - T::one()
        };
        let num = if num >= T::zero() {
            num
        } else {
            num - &self.c + T::one()
        };
        return Self::from(num / &self.c);
    }

    pub fn is_positive(&self) -> bool {
        self.panic_if_complex();

        if self.b.is_zero() {
            self.a.is_positive()
        } else if self.b.is_positive() {
            if self.a.is_positive() {
                true
            } else {
                &self.a * &self.a < &self.b * &self.b * &self.r
            }
        } else {
            // self.b.is_negative()
            if !self.a.is_positive() {
                false
            } else {
                &self.a * &self.a > &self.b * &self.b * &self.r
            }
        }
    }

    pub fn is_negative(&self) -> bool {
        self.panic_if_complex();

        if self.b.is_zero() {
            self.a.is_negative()
        } else if self.b.is_negative() {
            if self.a.is_negative() {
                true
            } else {
                &self.a * &self.a < &self.b * &self.b * &self.r
            }
        } else {
            // self.b.is_positive()
            if !self.a.is_negative() {
                false
            } else {
                &self.a * &self.a > &self.b * &self.b * &self.r
            }
        }
    }
}

impl<T> Into<(T, T, T, T)> for QuadraticSurd<T> {
    /// Deconstruct the quadratic surd `(a + b*sqrt(r)) / c` into tuple `(a,b,c,r)`
    fn into(self) -> (T, T, T, T) {
        (self.a, self.b, self.c, self.r)
    }
}

impl<T: QuadraticSurdBase + FromPrimitive + CheckedMul> QuadraticSurd<T>
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
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
            101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
            197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
        ];
        for p in SMALL_PRIMES {
            let p = T::from_u8(p).unwrap();
            if let Some(p2) = p.checked_mul(&p) {
                loop {
                    let (quo, rem) = self.r.div_rem(&p2);
                    if rem.is_zero() {
                        self.r = quo;
                        self.b = &self.b * &p;
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

fn quadsurd_to_f64<T: Integer + ToPrimitive>(a: T, b: T, c: T, r: T) -> Option<f64> {
    Some((a.to_f64()? + b.to_f64()? * r.to_f64()?.sqrt()) / c.to_f64()?)
}

#[cfg(not(feature = "complex"))]
impl<
        T: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
        U: QuadraticSurdBase,
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
fn quadsurd_to_contfrac<T: QuadraticSurdBase + WithUnsigned<Unsigned = U> + AddAssign, U>(
    a: T,
    b: T,
    c: T,
    r: T,
    neg: bool,
) -> ContinuedFraction<U>
where
    for<'r> &'r T: RefNum<T>,
{
    debug_assert!(!c.is_zero() && r >= T::zero());
    debug_assert!(r.sqrt() * r.sqrt() != r);

    // convert to form (p+sqrt(d))/q where q|d-p^2
    let mut d = r * &b * &b;
    let (mut p, mut q) = if b >= T::zero() { (a, c) } else { (-a, -c) };
    if (&d - &p * &p) % &q != T::zero() {
        d = d * &q * &q;
        p = p * &q;
        q = &q * &q;
    }

    // find the reduced form and aperiodic coefficients
    let mut a_coeffs: Vec<T> = Vec::new();
    let rd = d.sqrt();
    while a_coeffs.len() == 0 || // ensure that we have a first coefficient
          !(p <= rd && rd < (&p+&q) && (&q-&p) <= rd)
    {
        let a = (&rd + &p).div_floor(&q);
        p = &a * &q - p;
        q = (&d - &p * &p) / q;
        a_coeffs.push(a);
    }

    // find the periodic coefficients
    let mut p_coeffs: Vec<T> = Vec::new();
    let (init_p, init_q) = (p.clone(), q.clone());
    loop {
        let a = (&rd + &p).div_floor(&q);
        p = &a * &q - p;
        q = (&d - &p * &p) / q;
        p_coeffs.push(a);
        if p == init_p && q == init_q {
            break;
        }
    }

    ContinuedFraction::new(a_coeffs, p_coeffs, neg)
}

#[cfg(not(feature = "complex"))]
impl<T: QuadraticSurdBase + WithUnsigned<Unsigned = U> + AddAssign, U> From<QuadraticSurd<T>>
    for ContinuedFraction<U>
where
    for<'r> &'r T: RefNum<T>,
{
    fn from(s: QuadraticSurd<T>) -> Self {
        if s.is_negative() {
            quadsurd_to_contfrac(-s.a, -s.b, s.c, s.r, true)
        } else {
            quadsurd_to_contfrac(s.a, s.b, s.c, s.r, false)
        }
    }
}

#[cfg(feature = "complex")]
impl<T: Integer> TryFrom<QuadraticSurd<T>> for ContinuedFraction<T> {
    fn try_from(s: QuadraticSurd<T>) -> Self {
        unimplemented!()
    }
}

impl<
        T: QuadraticSurdBase + CheckedAdd + AddAssign + WithUnsigned<Unsigned = U>,
        U: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = T>,
    > Computable<T> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    #[cfg(not(feature = "complex"))]
    fn approximated(&self, limit: &T) -> Approximation<Ratio<T>> {
        ContinuedFraction::<U>::from(self.clone()).approximated(limit)
    }
}

impl<T: Integer + fmt::Display> fmt::Display for QuadraticSurd<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (
            self.a.is_zero(),
            self.b.is_zero(),
            self.b.is_one(),
            self.c.is_one(),
        ) {
            (true, true, _, _) => write!(f, "0"),
            (true, false, true, true) => write!(f, "√{}", self.r),
            (true, false, false, true) => write!(f, "{}√{}", self.b, self.r),
            (true, false, true, false) => write!(f, "√{}/{}", self.b, self.c),
            (true, false, false, false) => write!(f, "{}√{}/{}", self.b, self.r, self.c),
            (false, true, _, true) => write!(f, "{}", self.a),
            (false, true, _, false) => write!(f, "{}/{}", self.a, self.c),
            (false, false, true, true) => write!(f, "{} + √{}", self.a, self.c),
            (false, false, true, false) => write!(f, "({} + √{})/{}", self.a, self.r, self.c),
            (false, false, false, true) => write!(f, "{} + {}√{}", self.b, self.r, self.c),
            (false, false, false, false) => {
                write!(f, "({} + {}√{})/{}", self.a, self.b, self.r, self.c)
            }
        }
    }
}

// reduce root base for binary operation
#[inline]
fn reduce_bin_op<T: QuadraticSurdBase>(
    lhs: QuadraticSurd<T>,
    rhs: QuadraticSurd<T>,
) -> (QuadraticSurd<T>, QuadraticSurd<T>)
where
    for<'r> &'r T: RefNum<T>,
{
    let result = if lhs.r > rhs.r {
        let hint = &lhs.r / &rhs.r;
        (lhs.reduce_root_hinted(hint), rhs)
    } else if rhs.r > lhs.r {
        let hint = &rhs.r / &lhs.r;
        (lhs, rhs.reduce_root_hinted(hint))
    } else {
        (lhs, rhs)
    };

    if result.0.r != result.1.r {
        panic!("the root base should be the same!");
    }

    result
}

macro_rules! arith_impl {
    (impl $imp:ident, $method:ident) => {
        // Abstracts a/b `op` c/d = (a*lcm/b `op` c*lcm/d)/lcm where lcm = lcm(b,d)
        impl<T: QuadraticSurdBase> $imp<QuadraticSurd<T>> for QuadraticSurd<T>
        where
            for<'r> &'r T: RefNum<T>,
        {
            type Output = QuadraticSurd<T>;
            fn $method(self, rhs: QuadraticSurd<T>) -> QuadraticSurd<T> {
                let (lhs, rhs) = reduce_bin_op(self, rhs);

                if lhs.c == rhs.c {
                    return QuadraticSurd::new(
                        lhs.a.$method(rhs.a),
                        lhs.b.$method(rhs.b),
                        rhs.c,
                        rhs.r,
                    );
                }

                let lcm = lhs.c.lcm(&rhs.c);
                let lhs_r = &lcm / lhs.c;
                let rhs_r = &lcm / rhs.c;
                QuadraticSurd::new(
                    (lhs.a / &lhs_r).$method(rhs.a / &rhs_r),
                    (lhs.b / &lhs_r).$method(rhs.b / &rhs_r),
                    lcm,
                    rhs.r,
                )
            }
        }
        // Abstracts the a/b `op` c/1 = (a*1 `op` b*c) / (b*1) = (a `op` b*c) / b pattern
        impl<T: QuadraticSurdBase> $imp<T> for QuadraticSurd<T>
        where
            for<'r> &'r T: RefNum<T>,
        {
            type Output = QuadraticSurd<T>;
            #[inline]
            fn $method(self, rhs: T) -> QuadraticSurd<T> {
                QuadraticSurd::new(self.a.$method(&self.c * rhs), self.b, self.c, self.r)
            }
        }
    };
}

arith_impl!(impl Add, add);
arith_impl!(impl Sub, sub);

impl<T: QuadraticSurdBase> Mul<T> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = QuadraticSurd<T>;
    #[inline]
    fn mul(self, rhs: T) -> QuadraticSurd<T> {
        let gcd = self.c.gcd(&rhs);
        let rem = rhs / &gcd;
        QuadraticSurd::new(self.a * &rem, self.b * rem, self.c / gcd, self.r)
    }
}

impl<T: QuadraticSurdBase> Mul<QuadraticSurd<T>> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = QuadraticSurd<T>;
    #[inline]
    fn mul(self, rhs: QuadraticSurd<T>) -> QuadraticSurd<T> {
        let (lhs, rhs) = reduce_bin_op(self, rhs);

        let gcd_lr = (&lhs.a).gcd(&lhs.b).gcd(&rhs.c); // gcd between lhs numerator and rhs denominator
        let gcd_rl = (&rhs.a).gcd(&rhs.b).gcd(&rhs.c); // gcd between rhs numerator and lhs denominator

        let (la, lb, lc) = (&lhs.a / &gcd_lr, &lhs.b / &gcd_lr, &rhs.c / &gcd_rl);
        let (ra, rb, rc) = (&rhs.a / &gcd_rl, &rhs.b / &gcd_rl, &lhs.c / &gcd_lr);
        QuadraticSurd::new(
            &la * &ra + &lb * &rb * lhs.r,
            &la * &rb + &lb * &ra,
            lc * rc,
            rhs.r,
        )
    }
}

impl<T: QuadraticSurdBase> Div<T> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = QuadraticSurd<T>;
    #[inline]
    fn div(self, rhs: T) -> QuadraticSurd<T> {
        let gcd = self.a.gcd(&self.b).gcd(&rhs);
        let rem = rhs / &gcd;
        QuadraticSurd::new(self.a / &gcd, self.b / &gcd, self.c * rem, self.r)
    }
}

impl<T: QuadraticSurdBase> Div<QuadraticSurd<T>> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = QuadraticSurd<T>;
    #[inline]
    fn div(self, rhs: QuadraticSurd<T>) -> QuadraticSurd<T> {
        self.mul(rhs.recip())
    }
}

impl<T: QuadraticSurdBase> Neg for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = QuadraticSurd<T>;
    #[inline]
    fn neg(self) -> QuadraticSurd<T> {
        QuadraticSurd::new_raw(-self.a, -self.b, self.c, self.r)
    }
}

impl<T: QuadraticSurdBase> Zero for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    #[inline]
    fn zero() -> Self {
        QuadraticSurd {
            a: T::zero(),
            b: T::zero(),
            c: T::one(),
            r: T::zero(),
        }
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }
}

impl<T: QuadraticSurdBase> One for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    #[inline]
    fn one() -> Self {
        QuadraticSurd {
            a: T::one(),
            b: T::zero(),
            c: T::one(),
            r: T::zero(),
        }
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
where
    for<'r> &'r T: RefNum<T>,
{
    #[inline]
    fn to_i64(&self) -> Option<i64> {
        if self.r < T::zero() {
            return None;
        }
        match self.to_integer() {
            Approximation::Exact(v) => v.to_i64(),
            Approximation::Approximated(_) => None,
        }
    }

    #[inline]
    fn to_u64(&self) -> Option<u64> {
        if self.r < T::zero() {
            return None;
        }
        match self.to_integer() {
            Approximation::Exact(v) => v.to_u64(),
            Approximation::Approximated(_) => None,
        }
    }

    #[inline]
    fn to_f64(&self) -> Option<f64> {
        if self.r < T::zero() {
            return None;
        }
        quadsurd_to_f64(
            self.a.to_i64()?,
            self.b.to_i64()?,
            self.c.to_i64()?,
            self.r.to_i64()?,
        )
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FromSqrtError { // TODO: make this directly an enum
    kind: SqrtErrorKind,
}

#[derive(Copy, Clone, Debug, PartialEq)]
enum SqrtErrorKind {
    Overflow,
    Unrepresentable,
}

impl<T: QuadraticSurdBase> FromSqrt<T> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Error = FromSqrtError;

    #[inline]
    fn from_sqrt(target: T) -> Result<Self, Self::Error> {
        // TODO: complex number check
        let root = sqrt(target.clone());
        if &root * &root == target {
            Ok(Self::from(root))
        } else {
            Ok(QuadraticSurd {
                a: T::zero(),
                b: T::one(),
                c: T::one(),
                r: target,
            })
        }
    }
}

impl<T: QuadraticSurdBase + CheckedMul> FromSqrt<Ratio<T>> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Error = FromSqrtError;

    #[inline]
    fn from_sqrt(target: Ratio<T>) -> Result<Self, Self::Error> {
        if target.is_integer() {
            return Self::from_sqrt(target.to_integer());
        }

        let dd = target.denom().checked_mul(target.denom());
        let nd = target.numer().checked_mul(target.denom());
        match (dd, nd) {
            (Some(new_c), Some(new_r)) => Ok(QuadraticSurd::new(T::zero(), T::one(), new_c, new_r)),
            (_, _) => Err(FromSqrtError {
                kind: SqrtErrorKind::Overflow,
            }),
        }
    }
}

impl<T: QuadraticSurdBase + CheckedMul> FromSqrt<QuadraticSurd<T>> for QuadraticSurd<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Error = FromSqrtError;

    #[inline]
    fn from_sqrt(target: QuadraticSurd<T>) -> Result<Self, Self::Error> {
        if target.is_rational() {
            match target.to_rational() {
                Approximation::Exact(v) => return Self::from_sqrt(v),
                _ => unreachable!(),
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
            return Err(FromSqrtError {
                kind: SqrtErrorKind::Unrepresentable,
            });
        } // TODO: if complex number is enabled, don't report error here
        let sqrt_delta = Self::from_sqrt(delta2)?;
        if !sqrt_delta.is_integer() {
            return Err(FromSqrtError {
                kind: SqrtErrorKind::Unrepresentable,
            });
        }

        let delta2 = match sqrt_delta.to_integer() {
            Approximation::Exact(v) => v,
            _ => unreachable!(),
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

    pub const PHI45: QuadraticSurd<i32> = QuadraticSurd::new_raw(3, 1, 6, 45); // non-reduced version of PHI
    pub const PHI45_R: QuadraticSurd<i32> = QuadraticSurd::new_raw(-3, 1, 6, 45); // non-reduced version of PHI_R

    #[test]
    fn from_rational_test() {
        assert_eq!(
            QuadraticSurd::from_rational(Ratio::new(1, 2), Ratio::new(1, 2), Ratio::from(5)),
            PHI
        );
        assert_eq!(
            QuadraticSurd::from_rational(Ratio::new(-1, 2), Ratio::new(1, 2), Ratio::from(5)),
            PHI_R
        );
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
        assert_eq!(
            QuadraticSurd::from_sqrt(5).unwrap(),
            QuadraticSurd::new_raw(0, 1, 1, 5)
        );
        assert!(matches!(QuadraticSurd::from_sqrt(PHI), Err(_)));
        assert!(matches!(QuadraticSurd::from_sqrt(PHI * PHI), Ok(PHI)));

        assert!(matches!(
            QuadraticSurd::from_equation(Ratio::from(1), Ratio::from(-1), Ratio::from(-1)),
            Some(PHI)
        ));
    }

    #[test]
    fn arithmic_test() {
        let sq5 = QuadraticSurd::from_sqrt(5).unwrap();
        assert_eq!(-sq5, QuadraticSurd::new_raw(0, -1, 1, 5));

        // add
        assert_eq!(PHI + PHI_R, sq5);
        assert_eq!(PHI45 + PHI_R, sq5);
        assert_eq!(PHI + PHI45_R, sq5);
        assert_eq!(N_PHI + N_PHI_R, -sq5);
        assert!((PHI + N_PHI).is_zero());
        assert!((PHI + N_PHI_R).is_one());

        // sub
        assert!((PHI - PHI).is_zero());
        assert!((PHI - PHI45).is_zero());
        assert!((PHI - PHI_R).is_one());
        assert!((PHI - PHI45_R).is_one());
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
        assert!((PHI45 * PHI_R).is_one());
        assert!((PHI * PHI45_R).is_one());
        assert!((N_PHI * N_PHI_R).is_one());
        assert_eq!(PHI * PHI, PHI_SQ);
        assert_eq!(N_PHI * N_PHI, PHI_SQ);

        // div
        assert!((PHI / PHI).is_one());
        assert!((PHI45 / PHI).is_one());
        assert!((PHI / PHI45).is_one());
        assert!((N_PHI / N_PHI).is_one());
        assert_eq!(PHI / PHI_R, PHI_SQ);
        assert_eq!(N_PHI / N_PHI_R, PHI_SQ);
    }

    #[test]
    fn conversion_between_contfrac() {
        let cf_phi = ContinuedFraction::<u32>::new(vec![1], vec![1], false);
        assert_eq!(QuadraticSurd::from(cf_phi.clone()), PHI);
        assert_eq!(ContinuedFraction::from(PHI), cf_phi);

        let cf_sq2 = ContinuedFraction::<u32>::new(vec![1], vec![2], false);
        let surd_sq2 = QuadraticSurd::new(0, 1, 1, 2);
        assert_eq!(QuadraticSurd::from(cf_sq2.clone()), surd_sq2);
        assert_eq!(ContinuedFraction::from(surd_sq2), cf_sq2);

        let cf_n_sq2 = ContinuedFraction::<u32>::new(vec![1], vec![2], true);
        let surd_n_sq2 = QuadraticSurd::new(0, -1, 1, 2);
        assert_eq!(QuadraticSurd::from(cf_n_sq2.clone()), surd_n_sq2);
        assert_eq!(ContinuedFraction::from(surd_n_sq2), cf_n_sq2);

        let cf_sq2_7 = ContinuedFraction::<u32>::new(vec![0, 4], vec![1, 18, 1, 8], false);
        let surd_sq2_7 = QuadraticSurd::new(0, 1, 7, 2);
        assert_eq!(QuadraticSurd::from(cf_sq2_7.clone()), surd_sq2_7);
        assert_eq!(ContinuedFraction::from(surd_sq2_7), cf_sq2_7);

        let cf_10_sq2_7 = ContinuedFraction::<u32>::new(vec![1, 4], vec![2], false);
        let surd_10_sq2_7 = QuadraticSurd::new(10, -1, 7, 2);
        assert_eq!(QuadraticSurd::from(cf_10_sq2_7.clone()), surd_10_sq2_7);
        assert_eq!(ContinuedFraction::from(surd_10_sq2_7), cf_10_sq2_7);
    }
}
