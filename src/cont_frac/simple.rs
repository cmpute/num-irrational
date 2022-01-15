use core::str::FromStr;
use std::ops::{Add, AddAssign, Mul};
use std::mem::swap;
use std::fmt;
use num_traits::{float::FloatCore, Num, One, Zero, Signed, NumRef, RefNum, CheckedAdd, CheckedMul};
use num_integer::{Integer};
use num_rational::Ratio;
use crate::traits::{RationalApproximation, Approximation, WithSigned, WithUnsigned};
use crate::quad_surd::{QuadraticSurd, QuadraticSurdBase};

/// This struct represents a simple continued fraction a0 + 1/(a1 + 1/ (a2 + ...))
/// Where a0, a1, a2 are positive integers
/// It's capable of representing rational numbers and quadratic surds
/// 
/// References:
///     https://pi.math.cornell.edu/~gautam/ContinuedFractions.pdf
///     https://crypto.stanford.edu/pbc/notes/contfrac/
///     http://www.numbertheory.org/continued_fractions.html
///     http://www.numbertheory.org/php/cfrac.html
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuedFraction<T> {
    /// Coefficients of aperiodic part
    a_coeffs: Vec<T>,

    /// Coefficients of periodic part
    p_coeffs: Vec<T>,

    /// Sign of the fraction 
    negative: bool
}

impl<T> ContinuedFraction<T> {
    #[inline]
    pub fn aperiodic_coeffs(&self) -> &[T] {
        &self.a_coeffs[..]
    }

    #[inline]
    pub fn periodic_coeffs(&self) -> &[T] {
        &self.p_coeffs[..]
    }

    #[inline]
    pub fn is_negative(&self) -> bool {
        self.negative
    }

    #[inline]
    pub fn is_rational(&self) -> bool {
        self.p_coeffs.len() == 0
    }

    #[inline]
    pub fn is_integer(&self) -> bool {
        self.a_coeffs.len() == 1 && self.p_coeffs.len() == 0
    }
}

impl<T: Num> ContinuedFraction<T> {
    pub fn new(a_coeffs: Vec<T>, p_coeffs: Vec<T>, negative: bool) -> Self {
        let mut dedup_a = Vec::with_capacity(a_coeffs.len());
        let mut last_zero = false;
        for a in a_coeffs {
            if last_zero {
                if a.is_zero() {
                    continue; // skip consecutive 2 zeros
                }

                last_zero = false;
            } else {
                if a.is_zero() {
                    last_zero = true;
                }
            }

            dedup_a.push(a);
        }

        if dedup_a.len() == 0 && p_coeffs.len() == 0 {
            panic!("at least one coefficient is required!")
        }

        // clear negative sign for 0
        let mut negative = negative;
        if dedup_a.len() == 1 && dedup_a[0] == T::zero() && p_coeffs.len() == 0 && negative {
            negative = false;
        }

        ContinuedFraction { a_coeffs: dedup_a, p_coeffs, negative }
    }
}

#[derive(Debug, Clone)]
pub struct Coefficients<'a, T> {
    a_iter: Option<std::slice::Iter<'a, T>>, // None if aperiodic part has been consumed
    p_ref: &'a Vec<T>,
    p_iter: Option<std::slice::Iter<'a, T>> // None before aperiodic part is consumed, or when periodic part is empty
}

impl<'a, T> Iterator for Coefficients<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(it) = self.a_iter.as_mut() { // in aperiodic part
            match it.next() {
                Some(v) => Some(v),
                None => {
                    self.a_iter = None;
                    if self.p_ref.len() > 0 {
                        let mut new_iter = self.p_ref.iter();
                        let result = new_iter.next();
                        self.p_iter = Some(new_iter);
                        result
                    } else { None }
                }
            }
        } else {
            if let Some(it) = self.p_iter.as_mut() { // in periodic part
                match it.next() {
                    Some(v) => Some(v),
                    None => {
                        let mut new_iter = self.p_ref.iter();
                        let result = new_iter.next();
                        self.p_iter = Some(new_iter);
                        result
                    }
                }
            } else {
                None
            }
        }
    }
}

/// This iterator converts a coefficient iterator of simple continued fraction to
/// coefficients of a general continued fraction. Then it can be consumed by
/// `GeneralContinuedFraction`
pub struct GeneralCoefficients<T> {
    coeffs: T, negative: bool // store the sign
}

impl<'r, I: Iterator<Item = &'r T>, T: 'r + WithSigned<Signed = U> + Clone, U: Signed>
Iterator for GeneralCoefficients<I> {
    type Item = (U, U);

    fn next(&mut self) -> Option<Self::Item> {
        if self.negative {
            self.negative = false; // only apply to the first coefficient
            self.coeffs.next().map(|v| (-v.clone().to_signed(), U::one()))
        } else {
            self.coeffs.next().map(|v| (v.clone().to_signed(), U::one()))
        }
    }
}

pub struct Convergents<'a, T> {
    coeffs: Coefficients<'a, T>,
    pm1: T, // p_(k-1)
    pm2: T, // p_(k-2)
    qm1: T, // q_(k-1)
    qm2: T, // q_(k-2)
    neg: bool // store the sign
}

impl<'a, T: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
         U: Integer + Clone + Signed>
Iterator for Convergents<'a, T> {
    type Item = Ratio<U>;

    fn next(&mut self) -> Option<Self::Item> {
        let a = self.coeffs.next()?;
        // p_k = a_k * p_(k-1) + p_(k-2)
        let p = a.checked_mul(&self.pm1).and_then(|v| v.checked_add(&self.pm2))?;
        // q_k = a_k * q_(k-1) + q_(k-2)
        let q = a.checked_mul(&self.qm1).and_then(|v| v.checked_add(&self.qm2))?;

        swap(&mut self.pm2, &mut self.pm1); // self.pm2 = self.pm1
        swap(&mut self.qm2, &mut self.qm1); // self.qm2 = self.qm1
        self.pm1 = p.clone(); self.qm1 = q.clone();

        let r = Ratio::new(p.to_signed(), q.to_signed());
        if self.neg { Some(-r) } else { Some(r) }
    }
}

pub struct SignedCoefficients<T> {
    coeffs: T, negative: bool // store the sign
}

impl<'r, I: Iterator<Item = &'r T>, T: 'r + WithSigned<Signed = U> + Clone, U: Signed>
Iterator for SignedCoefficients<I> {
    type Item = U;

    fn next(&mut self) -> Option<Self::Item> {
        if self.negative {
            self.negative = false; // only apply to the first coefficient
            self.coeffs.next().map(|v| -v.clone().to_signed())
        } else {
            self.coeffs.next().map(|v| v.clone().to_signed())
        }
    }
}

impl<T> ContinuedFraction<T> {
    /// Returns an iterator of the coefficients in the continued fraction
    /// Note that for a negative number, the coefficients of it's absolute value is returned
    pub fn coeffs(&self) -> Coefficients<T> {
        Coefficients {
            a_iter: Some(self.a_coeffs.iter()),
            p_ref: &self.p_coeffs, p_iter: None
        }
    }

    /// Returns an iterator of generalized coefficients, that can be consumed
    /// by GeneralContinuedFraction
    pub fn generalized(&self) -> GeneralCoefficients<Coefficients<T>> {
        GeneralCoefficients {
            coeffs: self.coeffs(), negative: self.negative
        }
    }

    /// Returns an iterator of the coefficients in the continued fraction.
    /// The first coefficient is signed.
    pub fn coeffs_signed(&self) -> SignedCoefficients<Coefficients<T>> {
        SignedCoefficients {
            coeffs: self.coeffs(), negative: self.negative
        }
    }
}

impl<T: WithSigned<Signed = U> + Clone, U: Signed> ContinuedFraction<T> {
    /// Wrap the continued fraction object as `InfiniteContinuedFraction`
    pub fn expanded(&self) -> InfiniteContinuedFraction<SignedCoefficients<Coefficients<T>>> {
        InfiniteContinuedFraction(self.coeffs_signed())
    }
}

impl<T: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
     U: Integer + Clone + Signed> ContinuedFraction<T> {
    /// Returns an iterator of the convergents. The iterator will stop
    /// if all coefficients are consumed, or numeric overflow happened.
    pub fn convergents(&self) -> Convergents<T> {
        Convergents {
            coeffs: self.coeffs(),
            pm1: T::one(), pm2: T::zero(), qm1: T::zero(), qm2: T::one(),
            neg: self.negative
        }
    }

    #[inline]
    /// This method returns the corresponding rational number if it's rational,
    /// returns the expansion until the first repeating occurence
    pub fn to_rational(&self) -> Approximation<Ratio<U>> {
        if self.is_rational() {
            Approximation::Exact(self.convergents().last().unwrap())
        } else {
            Approximation::Approximated(self.convergents().nth(
                self.a_coeffs.len() + self.p_coeffs.len()).unwrap())
        }
    }
}

impl<T: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
     U: Integer + Clone + Signed + CheckedAdd>
RationalApproximation<U> for ContinuedFraction<T>
{
    fn approx_rational(&self, limit: &U) -> Approximation<Ratio<U>> {
        let within_limit = |v: &U| if v >= &U::zero() { v < limit } else { limit.checked_add(v).unwrap() >= U::zero() };
        let ratio_within_limit = |v: &Ratio<U>| within_limit(v.numer()) && within_limit(v.denom());

        let mut convergents = self.convergents();
        let mut last_conv = convergents.next().unwrap();
        if !ratio_within_limit(&last_conv) { 
            let i = self.a_coeffs.first().unwrap().clone();
            return Approximation::Approximated(Ratio::from(i.to_signed()))
        }
        loop {
            last_conv = match convergents.next() {
                Some(v) => if ratio_within_limit(&v) { v }
                           else { return Approximation::Approximated(last_conv); },
                None => return Approximation::Exact(last_conv)
            }
        }
    }
}


impl<T: fmt::Display> fmt::Display for ContinuedFraction<T>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.negative { write!(f, "-")?; }

        write!(f, "[{}", self.a_coeffs.first().unwrap())?;
        if self.a_coeffs.len() == 1 {
            if self.p_coeffs.len() == 0 {
                return write!(f, "]")
            } else {
                write!(f, "; ")?;
            }
        } else {
            let mut aiter = self.a_coeffs.iter().skip(1);
            write!(f, "; {}", aiter.next().unwrap())?;
            while let Some(v) = aiter.next() {
                write!(f, ", {}", v)?;
            }
            if self.p_coeffs.len() > 0 {
                write!(f, ", ")?;
            }
        }

        if self.p_coeffs.len() > 0 {
            let mut piter = self.p_coeffs.iter();
            write!(f, "({}", piter.next().unwrap())?;
            while let Some(v) = piter.next() {
                write!(f, ", {}", v)?;
            }
            write!(f, ")]")
        } else {
            write!(f, "]")
        }
    }
}

impl<T: Zero> Zero for ContinuedFraction<T> {
    fn zero() -> Self {
        ContinuedFraction { a_coeffs: vec![T::zero()], p_coeffs: Vec::new(), negative: false }
    }

    fn is_zero(&self) -> bool {
        self.a_coeffs.len() == 1 && self.a_coeffs[0].is_zero() && self.p_coeffs.len() == 0
    }
}

impl<T: One + PartialEq> One for ContinuedFraction<T> {
    fn one() -> Self {
        ContinuedFraction { a_coeffs: vec![T::one()], p_coeffs: Vec::new(), negative: false }
    }

    fn is_one(&self) -> bool {
        self.a_coeffs.len() == 1 && self.a_coeffs[0].is_one() && self.p_coeffs.len() == 0
    }
}

// TODO: implement from_float, using FloatCore trait?

impl<T: Num + WithUnsigned<Unsigned = U> + PartialOrd, U> From<T> for ContinuedFraction<U> {
    fn from(t: T) -> Self {
        if t < T::zero() {
            ContinuedFraction { a_coeffs: vec![(T::zero() - t).to_unsigned()], p_coeffs: Vec::new(), negative: true }
        } else {
            ContinuedFraction { a_coeffs: vec![t.to_unsigned()], p_coeffs: Vec::new(), negative: false }
        }
    }
}

impl<T: Integer + Clone + WithUnsigned<Unsigned = U>, U: Zero> From<Ratio<T>> for ContinuedFraction<U> {
    fn from(r: Ratio<T>) -> Self {
        if r.is_zero() { return Self::zero() }

        let mut coeffs = Vec::new();
        let (mut n, mut d) = r.into();

        let negative = if n < T::zero() { n = T::zero() - n; true } else { false };

        if n < d {
            std::mem::swap(&mut n, &mut d);
            coeffs.push(U::zero());
        }

        while !d.is_zero() {
            let (quo, rem) = n.div_rem(&d);
            coeffs.push(quo.to_unsigned());
            n = d; d = rem;
        }

        ContinuedFraction { a_coeffs: coeffs, p_coeffs: Vec::new(), negative }
    }
}

pub struct ParseContFracError {
}

impl<T> FromStr for ContinuedFraction<T> {
    type Err = ParseContFracError;

    /// Parse from standard format (like 355/113 = "[3; 7, 16]", (1+sqrt(5))/2 = "[1; (1)]")
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        unimplemented!()
    }
}

impl<T: Num + WithUnsigned<Unsigned = U> + PartialOrd, U: NumRef + PartialOrd> Add<T> for ContinuedFraction<U>
where for <'r> &'r U: RefNum<U>{
    type Output = Self;

    fn add(self, rhs: T) -> Self {
        let mut new_a = self.a_coeffs;
        let i = new_a.first().unwrap();
        let (new_i, negative) = if self.negative {
            if rhs < T::zero() { // neg + neg
                let rhs = (T::zero() - rhs).to_unsigned();
                (rhs + i, true)
            } else { // neg + pos
                let rhs = rhs.to_unsigned();
                if *i < rhs {
                    (rhs - i, false)
                } else {
                    (i - rhs, true)
                }
            }
        } else {
            if rhs < T::zero() { // pos + neg
                let rhs = (T::zero() - rhs).to_unsigned();
                if *i < rhs {
                    (rhs - i, true)
                } else {
                    (i - rhs, false)
                }
            } else { // pos + pos
                (i + rhs.to_unsigned(), false)
            }
        };
        *new_a.first_mut().unwrap() = new_i;

        let mut result = ContinuedFraction { a_coeffs: new_a, p_coeffs: self.p_coeffs, negative };
        if result.is_zero() { result.negative = false; } // clear negative flag for 0
        result
    }
}

impl<T: QuadraticSurdBase + WithUnsigned<Unsigned = U>,
     U: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = T>>
Mul<T> for ContinuedFraction<U>
where for<'r> &'r T: RefNum<T>, for<'r> &'r U: RefNum<U> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        if self.is_integer() {
            let mut new_a = self.a_coeffs;
            let i = new_a.first().unwrap();
            let (new_i, negative) = if rhs < T::zero() {
                let rhs = (T::zero() - rhs).to_unsigned();
                (rhs * i, !self.negative)
            } else {
                (rhs.to_unsigned() * i, self.negative)
            };
            *new_a.first_mut().unwrap() = new_i;

            let mut result = ContinuedFraction { a_coeffs: new_a, p_coeffs: self.p_coeffs, negative };
            if result.is_zero() { result.negative = false; } // clear negative flag for 0
            result
        } else {
            Self::from(QuadraticSurd::<T>::from(self) * rhs)
        }
    }
}

impl<T> Add<Ratio<T>> for ContinuedFraction<T> {
    type Output = Self;

    fn add(self, rhs: Ratio<T>) -> Self {
        unimplemented!() // TODO: implement by converting to ratio or quadratic surd
    }
}

impl<T> Add<ContinuedFraction<T>> for ContinuedFraction<T> {
    type Output = Self;

    fn add(self, rhs: ContinuedFraction<T>) -> Self {
        // shortcut cases: is_integer, is_rational
        unimplemented!()
    }
}

impl<T> Mul<ContinuedFraction<T>> for ContinuedFraction<T> {
    type Output = Self;

    fn mul(self, rhs: ContinuedFraction<T>) -> Self {
        unimplemented!()
    }
}

/// This trait represents a regular continued fraction with infinite
/// coefficients. All operations here will be done iteratively
#[derive(Clone, Copy)]
pub struct InfiniteContinuedFraction<I: Iterator> (pub I);

impl<I: Iterator<Item = T> + Clone, T: Num + Clone> InfiniteContinuedFraction<I> {
    pub fn generalize(self) -> std::iter::Zip<I, std::iter::Repeat<T>> {
        self.0.zip(std::iter::repeat(T::one()))
    }
}

impl<I: Iterator<Item = T>, T: Integer + Num + NumRef + Clone> InfiniteContinuedFraction<I>
where for <'r> &'r T: RefNum<T> {
    /// This method returns a homographic function result on the fraction
    /// A homographic function is `(ax + b)/(cx + d)`
    pub fn homo(self, a: T, b: T, c: T, d: T)
        -> InfiniteContinuedFraction<HomographicResult<I>> {
        InfiniteContinuedFraction(
            HomographicResult { pm1: a, pm2: b, qm1: c, qm2: d, coeffs: self.0 }
        )
    }

    /// This method returns a bihomographic function result on the fraction
    /// A bihomographic function is `(axy + bx + cy + d)/(exy + fx + gy + h)`
    pub fn bihomo<U: Iterator<Item = T>>(
        self, rhs: InfiniteContinuedFraction<U>,
        a: T, b: T, c: T, d: T,
        e: T, f: T, g: T, h: T)
        -> InfiniteContinuedFraction<BihomographicResult<I, U, T>> {
        InfiniteContinuedFraction(BihomographicResult {
            pm11: a, pm12: b, pm21: c, pm22: d,
            qm11: e, qm12: f, qm21: g, qm22: h,
            x_coeffs: self.0, y_coeffs: rhs.0, down_first: false
        })
    }

    /// Take first N coefficients in the sequence and turn it into a
    /// `ContinuedFraction` object.
    pub fn take<U>(self, count: usize) -> ContinuedFraction<U> {
        unimplemented!()
    }

    /// Take first N coefficients in the sequence and turn it into a
    /// `ContinuedFraction` object with periodic detection.
    pub fn take_periodic<U>(self, count: usize) -> ContinuedFraction<U> {
        unimplemented!()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct HomographicResult<I: Iterator> {
    pm1: I::Item, pm2: I::Item, qm1: I::Item, qm2: I::Item,
    coeffs: I
}

impl<I: Iterator<Item = T>, T: Integer + Num + NumRef + Clone> Iterator for HomographicResult<I>
where for <'r> &'r T: RefNum<T> {
    type Item = T;

    // use the magic table method described in https://crypto.stanford.edu/pbc/notes/contfrac/compute.html
    fn next(&mut self) -> Option<T> {
        loop {
            match self.coeffs.next() {
                Some(v) => {
                    let p = &v * &self.pm1 + &self.pm2;
                    let q = v * &self.qm1 + &self.qm2;
                    
                    let i = if q.is_zero() { T::zero() } else { p.div_floor(&q) };
                    if !q.is_zero() && !self.qm1.is_zero() && i == self.pm1.div_floor(&self.qm1) {
                        let new_qm2 = &self.pm1 - &i * &self.qm1;
                        let new_qm1 = p - &i * &q;
                        self.pm1 = q;
                        swap(&mut self.pm2, &mut self.qm1); // self.pm2 = self.qm1
                        self.qm1 = new_qm1; self.qm2 = new_qm2;
                        break Some(i)
                    } else {
                        swap(&mut self.pm2, &mut self.pm1); // self.pm2 = self.pm1
                        swap(&mut self.qm2, &mut self.qm1); // self.qm2 = self.qm1
                        self.pm1 = p.clone(); self.qm1 = q.clone();
                    }
                },
                None => break None
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct BihomographicResult<X: Iterator<Item = T>, Y: Iterator<Item = T>, T> {
    pm11: T, // p with a_(i-1), b_(j-1)
    pm12: T, // p with a_(i-1), b_(j-2)
    pm21: T, // ..
    pm22: T, // ..
    qm11: T, // q with a_(i-1), b_(j-1)
    qm12: T, // q with a_(i-1), b_(j-2)
    qm21: T, // ..
    qm22: T, // ..
    x_coeffs: X, y_coeffs: Y,
    down_first: bool, // used for letting the table move right and move down iteratively
}

impl<X: Iterator<Item = T>, Y: Iterator<Item = T>, T: Integer + Num + NumRef>
Iterator for BihomographicResult<X, Y, T>
where for <'r> &'r T: RefNum<T> {
    type Item = T;

    // use the magic table method described in https://crypto.stanford.edu/pbc/notes/contfrac/bihom.html
    fn next(&mut self) -> Option<T> {
        loop {
            let i11 = if self.qm11.is_zero() { T::zero() } else { self.pm11.div_floor(&self.qm11) };
            let i22 = if self.qm22.is_zero() { T::zero() } else { self.pm22.div_floor(&self.qm22) };
            if !self.down_first && (self.qm11.is_zero() || self.qm12.is_zero() || self.qm22.is_zero()
                || i11 != self.pm12.div_floor(&self.qm12) || i11 != i22) {

                match self.x_coeffs.next() {
                    Some(v) => {
                        // move down first next time
                        self.down_first = true;

                        let p1 = &v * &self.pm11 + &self.pm21;
                        let q1 = &v * &self.qm11 + &self.qm21;

                        swap(&mut self.pm21, &mut self.pm11); // self.pm21 = self.pm11
                        swap(&mut self.qm21, &mut self.qm11); // self.qm21 = self.qm11
                        self.pm11 = p1; self.qm11 = q1;

                        let p2 = &v * &self.pm12 + &self.pm22;
                        let q2 = v * &self.qm12 + &self.qm22;

                        swap(&mut self.pm22, &mut self.pm12); // self.pm22 = self.pm12
                        swap(&mut self.qm22, &mut self.qm12); // self.qm22 = self.qm12
                        self.pm12 = p2; self.qm12 = q2;
                        continue
                    },
                    None => break None
                }
            }

            self.down_first = false; // clear flag immediately after it takes effect

            if self.qm11.is_zero() || self.qm21.is_zero() || self.qm12.is_zero()
                || i11 != self.pm21.div_floor(&self.qm21) || i11 != i22 {

                match self.y_coeffs.next() {
                    Some(v) => {
                        let p1 = &v * &self.pm11 + &self.pm12;
                        let q1 = &v * &self.qm11 + &self.qm12;

                        swap(&mut self.pm12, &mut self.pm11); // self.pm12 = self.pm11
                        swap(&mut self.qm12, &mut self.qm11); // self.qm12 = self.qm11
                        self.pm11 = p1; self.qm11 = q1;

                        let p2 = &v * &self.pm21 + &self.pm22;
                        let q2 = v * &self.qm21 + &self.qm22;

                        swap(&mut self.pm22, &mut self.pm21); // self.pm22 = self.pm21
                        swap(&mut self.qm22, &mut self.qm21); // self.qm22 = self.qm21
                        self.pm21 = p2; self.qm21 = q2;
                        continue
                    },
                    None => break None
                }
            }

            // now all fractions floor to the same value
            let i = i11;
            let new_qm11 = &self.pm11 - &i * &self.qm11;
            swap(&mut self.pm11, &mut self.qm11); self.qm11 = new_qm11;
            let new_qm12 = &self.pm12 - &i * &self.qm12;
            swap(&mut self.pm12, &mut self.qm12); self.qm12 = new_qm12;
            let new_qm21 = &self.pm21 - &i * &self.qm21;
            swap(&mut self.pm21, &mut self.qm21); self.qm21 = new_qm21;
            let new_qm22 = &self.pm22 - &i * &self.qm22;
            swap(&mut self.pm22, &mut self.qm22); self.qm22 = new_qm22;

            break Some(i)
        }
    }
}

impl<I: Iterator<Item = T>, T: Integer + Num + NumRef + Clone> Add<T> for InfiniteContinuedFraction<I>
where for <'r> &'r T: RefNum<T> {
    type Output = InfiniteContinuedFraction<HomographicResult<I>>;

    fn add(self, rhs: T) -> Self::Output {
        self.homo(T::one(), rhs, T::zero(), T::one())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbols::E;

    #[test]
    fn cont_frac_iter_test() {
        let one = ContinuedFraction::<u32>::new(vec![1], vec![], false);
        assert_eq!(one.coeffs().map(|&v| v).collect::<Vec<_>>(), vec![1]);
        assert_eq!(one.convergents().collect::<Vec<_>>(), vec![Ratio::from(1)]);

        let n_one = ContinuedFraction::<u32>::new(vec![1], vec![], true);
        assert_eq!(n_one.coeffs().map(|&v| v).collect::<Vec<_>>(), vec![1]);
        assert_eq!(n_one.convergents().collect::<Vec<_>>(), vec![Ratio::from(-1)]);

        let sq2 = ContinuedFraction::<u32>::new(vec![1], vec![2], false);
        assert_eq!(sq2.coeffs().take(5).map(|&v| v).collect::<Vec<_>>(), vec![1, 2, 2, 2, 2]);
        assert_eq!(sq2.convergents().take(5).collect::<Vec<_>>(),
            vec![Ratio::from(1), Ratio::new(3, 2), Ratio::new(7, 5), Ratio::new(17, 12), Ratio::new(41, 29)]);

        let n_sq2 = ContinuedFraction::<u32>::new(vec![1], vec![2], true);
        assert_eq!(n_sq2.coeffs().take(5).map(|&v| v).collect::<Vec<_>>(), vec![1, 2, 2, 2, 2]);
        assert_eq!(n_sq2.convergents().take(5).collect::<Vec<_>>(),
            vec![Ratio::from(-1), Ratio::new(-3, 2), Ratio::new(-7, 5), Ratio::new(-17, 12), Ratio::new(-41, 29)]);
    }

    #[test]
    fn cont_frac_conversion_test() {
        assert_eq!(ContinuedFraction::from(Ratio::from(3)), ContinuedFraction::new(vec![3u32], vec![], false));
        assert_eq!(ContinuedFraction::from(Ratio::new(22, 7)), ContinuedFraction::new(vec![3u32, 7], vec![], false));
        assert_eq!(ContinuedFraction::from(Ratio::new(-22, 7)), ContinuedFraction::new(vec![3u32, 7], vec![], true));
        assert_eq!(ContinuedFraction::from(Ratio::new(7, 22)), ContinuedFraction::new(vec![0u32, 3, 7], vec![], false));
        assert_eq!(ContinuedFraction::from(Ratio::new(-7, 22)), ContinuedFraction::new(vec![0u32, 3, 7], vec![], true));
        assert_eq!(ContinuedFraction::from(Ratio::new(355, 113)), ContinuedFraction::new(vec![3u32, 7, 16], vec![], false));
    }

    #[test]
    fn fmt_test() {
       assert_eq!(format!("{}", ContinuedFraction::new(vec![1], vec![], false)), "[1]");
       assert_eq!(format!("{}", ContinuedFraction::new(vec![1, 2, 3], vec![], false)), "[1; 2, 3]");
       assert_eq!(format!("{}", ContinuedFraction::new(vec![1], vec![2], false)), "[1; (2)]");
       assert_eq!(format!("{}", ContinuedFraction::new(vec![1, 2, 3], vec![3, 2], false)), "[1; 2, 3, (3, 2)]");
    }

    #[test]
    fn cont_frac_arithmetic_test() {
        let one = ContinuedFraction::one();
        let n_one = ContinuedFraction::<u32>::new(vec![1], vec![], true);
        let sq2 = ContinuedFraction::new(vec![1], vec![2], false);
        let n_sq2 = ContinuedFraction::new(vec![1], vec![2], true);

        // Add
        assert_eq!(one.clone() + 1, ContinuedFraction::from(2));
        assert_eq!(one.clone() + (-1), ContinuedFraction::zero());
        assert_eq!(n_one.clone() + 1, ContinuedFraction::zero());
        assert_eq!(n_one.clone() + (-1), ContinuedFraction::from(-2));

        assert_eq!(sq2.clone() + 1, ContinuedFraction::new(vec![2u32], vec![2], false));
        assert_eq!(sq2.clone() + (-1), ContinuedFraction::new(vec![0u32], vec![2], false));
        assert_eq!(n_sq2.clone() + 1, ContinuedFraction::new(vec![0u32], vec![2], true));
        assert_eq!(n_sq2.clone() + (-1), ContinuedFraction::new(vec![2u32], vec![2], true));

        // Mul
        assert_eq!(one.clone() * 2, ContinuedFraction::from(2));
        assert_eq!(one.clone() * -2, ContinuedFraction::from(-2));
        assert_eq!(n_one.clone() * 2, ContinuedFraction::from(-2));
        assert_eq!(n_one.clone() * -2, ContinuedFraction::from(2));

        assert_eq!(sq2.clone() * 2, ContinuedFraction::new(vec![2u32], vec![1, 4], false));
        assert_eq!(sq2.clone() * -2, ContinuedFraction::new(vec![2u32], vec![1, 4], true));
        assert_eq!(n_sq2.clone() * 2, ContinuedFraction::new(vec![2u32], vec![1, 4], true));
        assert_eq!(n_sq2.clone() * -2, ContinuedFraction::new(vec![2u32], vec![1, 4], false));
    }

    #[test]
    fn inf_cont_frac_arithmetic_test() {
        let e = E {};
        let e_cf = InfiniteContinuedFraction(e.cfrac::<i64>());
        let ep1_cf = e_cf + 1;
        assert_eq!(ep1_cf.0.take(5).collect::<Vec<_>>(), vec![3,1,2,1,1]);

        let sq2 = ContinuedFraction::<u32>::new(vec![1], vec![2], false);
        let sq2p1 = sq2.clone() + 1;
        assert_eq!((sq2.expanded() + 1).0.take(5).collect::<Vec<_>>(), sq2p1.expanded().0.take(5).collect::<Vec<_>>());
    }
}
