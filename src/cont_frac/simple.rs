//! Implementation of simple continued fractions

use super::block::{Block, DualBlock};
use crate::quadratic::{QuadraticBase, surd::QuadraticSurd};
use crate::traits::{Approximation, Computable, WithSigned, WithUnsigned};
use core::convert::TryFrom;
use core::str::FromStr;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{CheckedAdd, CheckedMul, Num, NumRef, One, RefNum, Signed, Zero};
use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

/// This struct represents a simple continued fraction `a0 + 1/(a1 + 1/ (a2 + ...))`
/// Where a0 is an signed integer, a1, a2, .. are positive integers
/// It's capable of representing rational numbers and quadratic surds
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuedFraction<T> {
    /// Coefficients of aperiodic part
    a_coeffs: Vec<T>,

    /// Coefficients of periodic part
    p_coeffs: Vec<T>,

    /// Sign of the fraction
    negative: bool,
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

impl<U> ContinuedFraction<U> {
    /// Create a continued fraction from
    ///
    /// This function will make sure that if two numbers are equal,
    /// their internal representation in continued fraction will be the same
    pub fn new<T: Num + PartialOrd + AddAssign + WithUnsigned<Unsigned = U>>(
        a_coeffs: Vec<T>,
        p_coeffs: Vec<T>,
        negate: bool,
    ) -> Self {
        if a_coeffs.len() == 0 && p_coeffs.len() == 0 {
            panic!("at least one coefficient is required!");
        }

        let mut negate = negate; // carry the flag when we process through the coefficients
        let mut negative = None; // final sign for the number
        let mut dedup_a: Vec<T> = Vec::with_capacity(a_coeffs.len());
        let mut it = a_coeffs.into_iter();

        // iteratively consume aperiodic coefficients
        loop {
            match it.next() {
                Some(a) => {
                    // apply negate flag
                    let a = if negate { T::zero() - a } else { a };

                    // special cases
                    if dedup_a.len() == 0 {
                        // if first element
                        // make sure the first element is non-negative
                        if a < T::zero() {
                            dedup_a.push(T::zero() - a);
                            negative = Some(true);
                            negate = !negate;
                        } else {
                            if !a.is_zero() {
                                negative = Some(false);
                            }
                            dedup_a.push(a);
                        }
                        continue;
                    } else if a.is_zero() {
                        // if a is zero
                        if dedup_a.last().unwrap().is_zero() {
                            dedup_a.pop(); // remove consecutive 2 zeros
                        } else {
                            dedup_a.push(a); // dealt by the next value
                        }
                        continue;
                    }

                    // make sure all coefficients are non-negative
                    if a < T::zero() {
                        // In this case we use the following formula
                        // [.., a_{n-1}, -a_n, a_{n+1}, ..] = [.., a_{n-1}-1, 1, a_n-1, -[a_{n+1}, ..]]
                        let mut target = a;
                        loop {
                            let am1 = dedup_a.pop().unwrap(); // a_{n-1}
                            let has_am2 = dedup_a.len() != 0; // has a_{n-2}
                            debug_assert!(am1 >= T::zero());

                            if am1.is_one() {
                                if has_am2 {
                                    // [.., a_{n-2}, 1-1, 1, a_n-1, -[a_{n+1}, ..]]
                                    // = [.., a_{n-2}+1, a_n-1, -[a_{n+1}, ..]]
                                    *dedup_a.last_mut().unwrap() += T::one();
                                    dedup_a.push(T::zero() - target - T::one());
                                } else {
                                    // a_{n-1}-1 = 0 will be the first term
                                    dedup_a.push(T::zero());
                                    dedup_a.push(T::one());
                                    dedup_a.push(T::zero() - target - T::one());
                                    negative = Some(false);
                                }
                                negate = !negate;
                                break;
                            } else if am1.is_zero() {
                                if has_am2 {
                                    // [.., a_{n-2}, 0-1, 1, a_n-1, -[a_{n+1}, ..]]
                                    // = [.., a_{n-2}-1, 1, 0, -[1, a_n-1, -[a_{n+1}, ..]]]
                                    // = [.., a_{n-2}-a_n, a_{n+1}, ..]
                                    let am2 = dedup_a.pop().unwrap();
                                    let new_target = am2 + target;
                                    if new_target < T::zero() {
                                        // propagate to previous terms
                                        target = new_target;
                                    } else {
                                        dedup_a.push(new_target);
                                        break;
                                    }
                                } else {
                                    // a_n is the first non-zero coefficient, flip the sign
                                    dedup_a.push(T::zero());
                                    dedup_a.push(T::zero() - target);
                                    negative = Some(true);
                                    negate = !negate;
                                    break;
                                }
                            } else {
                                // normal case
                                debug_assert!(am1 > T::one());
                                dedup_a.push(am1 - T::one());
                                dedup_a.push(T::one());
                                dedup_a.push(T::zero() - target - T::one());
                                negate = !negate;
                                break;
                            }
                        }
                    } else {
                        if dedup_a.last().unwrap().is_zero() && dedup_a.len() > 1 {
                            // combine current value with the value before zero
                            dedup_a.pop();
                            *dedup_a.last_mut().unwrap() += a;
                        } else {
                            debug_assert!(!a.is_zero());
                            if matches!(negative, None) {
                                if a < T::zero() {
                                    println!("Check1!");
                                }
                                negative = Some(a < T::zero());
                            }
                            dedup_a.push(a);
                        }
                    }
                }
                None => {
                    break;
                }
            }
        }

        if dedup_a.len() == 0 && p_coeffs.len() == 0 {
            panic!("no effective coefficient!")
        }

        // when the last term is one or still zero
        if p_coeffs.len() == 0 {
            while dedup_a.len() >= 2 {
                let last_a = dedup_a.last().unwrap();
                if last_a.is_zero() {
                    dedup_a.pop();
                    dedup_a.pop();
                } else if last_a.is_one() {
                    let one = dedup_a.pop().unwrap();
                    *dedup_a.last_mut().unwrap() += one;
                } else {
                    break;
                }
            }

            if dedup_a.len() == 1 && dedup_a[0].is_zero() {
                negative = Some(false);
            }
        }

        if matches!(negative, None) {
            if matches!(p_coeffs.first(), Some(p) if p <= &T::zero()) {
                // TODO: how to handle negative coefficients in the periodic parts
                unimplemented!()
            } else {
                negative = Some(negate);
            }
        }

        // collect the results
        let a_coeffs: Vec<U> = dedup_a.into_iter().map(|v| v.to_unsigned()).collect();
        let p_coeffs: Vec<U> = p_coeffs.into_iter().map(|v| v.to_unsigned()).collect();
        ContinuedFraction {
            a_coeffs,
            p_coeffs,
            negative: negative.unwrap(),
        }
    }
}

/// Iterator of coeffcients in a [ContinuedFraction]
#[derive(Debug, Clone)]
pub struct Coefficients<'a, T> {
    a_iter: Option<std::slice::Iter<'a, T>>, // None if aperiodic part has been consumed
    p_ref: &'a Vec<T>,
    p_iter: Option<std::slice::Iter<'a, T>>, // None before aperiodic part is consumed, or when periodic part is empty
}

impl<'a, T> Iterator for Coefficients<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(it) = self.a_iter.as_mut() {
            // in aperiodic part
            match it.next() {
                Some(v) => Some(v),
                None => {
                    self.a_iter = None;
                    if self.p_ref.len() > 0 {
                        let mut new_iter = self.p_ref.iter();
                        let result = new_iter.next();
                        self.p_iter = Some(new_iter);
                        result
                    } else {
                        None
                    }
                }
            }
        } else {
            if let Some(it) = self.p_iter.as_mut() {
                // in periodic part
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

/// Iterator that converts coefficients of simple continued fraction to
/// coefficients of a general continued fraction. Then it can be consumed by
/// [GeneralContinuedFraction][crate::GeneralContinuedFraction]
pub struct GeneralCoefficients<T> {
    coeffs: T,
    negative: bool, // store the sign
}

impl<'r, I: Iterator<Item = &'r T>, T: 'r + WithSigned<Signed = U> + Clone, U: Signed> Iterator
    for GeneralCoefficients<I>
{
    type Item = (U, U);

    fn next(&mut self) -> Option<Self::Item> {
        if self.negative {
            self.coeffs
                .next()
                .map(|v| (U::one(), -v.clone().to_signed()))
        } else {
            self.coeffs
                .next()
                .map(|v| (U::one(), v.clone().to_signed()))
        }
    }
}

/// Iterator of convergents of a [ContinuedFraction]
pub struct Convergents<'a, T> {
    coeffs: Coefficients<'a, T>,
    block: Block<T>,
    neg: bool, // store the sign
}

impl<
        'a,
        T: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
        U: Integer + Clone + Signed,
    > Iterator for Convergents<'a, T>
{
    type Item = Ratio<U>;

    fn next(&mut self) -> Option<Self::Item> {
        let a = self.coeffs.next()?;
        let (p, q) = self.block.checked_rmove(a.clone())?;
        self.block.update(p.clone(), q.clone());

        let r = Ratio::new(p.to_signed(), q.to_signed());
        if self.neg {
            Some(-r)
        } else {
            Some(r)
        }
    }
}

/// Iterator of coefficients in a [ContinuedFraction] with sign applied. This iterator
/// can be used to construct a [InfiniteContinuedFraction]
pub struct SignedCoefficients<T> {
    coeffs: T,
    negative: bool, // store the sign
}

impl<'r, I: Iterator<Item = &'r T>, T: 'r + WithSigned<Signed = U> + Clone, U: Signed> Iterator
    for SignedCoefficients<I>
{
    type Item = U;

    fn next(&mut self) -> Option<Self::Item> {
        if self.negative {
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
            p_ref: &self.p_coeffs,
            p_iter: None,
        }
    }

    /// Returns an iterator of generalized coefficients, that can be consumed
    /// by GeneralContinuedFraction
    pub fn generalized(&self) -> GeneralCoefficients<Coefficients<T>> {
        GeneralCoefficients {
            coeffs: self.coeffs(),
            negative: self.negative,
        }
    }

    /// Returns an iterator of the coefficients in the continued fraction.
    /// The coefficients will be negative if the number is negative
    pub fn coeffs_signed(&self) -> SignedCoefficients<Coefficients<T>> {
        SignedCoefficients {
            coeffs: self.coeffs(),
            negative: self.negative,
        }
    }
}

impl<T: WithSigned<Signed = U> + Clone, U: Signed> ContinuedFraction<T> {
    /// Wrap the continued fraction object as `InfiniteContinuedFraction`
    pub fn expanded(&self) -> InfiniteContinuedFraction<SignedCoefficients<Coefficients<T>>> {
        InfiniteContinuedFraction(self.coeffs_signed())
    }
}

impl<
        T: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
        U: Integer + Clone + Signed,
    > ContinuedFraction<T>
{
    /// Returns an iterator of the convergents. The iterator will stop
    /// if all coefficients are consumed, or numeric overflow happened.
    pub fn convergents(&self) -> Convergents<T> {
        Convergents {
            coeffs: self.coeffs(),
            block: Block::identity(),
            neg: self.negative,
        }
    }

    /// Converts to an integer by using the first coefficient
    #[inline]
    pub fn to_integer(&self) -> Approximation<U> {
        let a = self.a_coeffs.first().unwrap().clone().to_signed();
        let a = if self.negative { -a } else { a };
        if self.is_integer() {
            Approximation::Exact(a)
        } else {
            Approximation::Approximated(a)
        }
    }

    #[inline]
    /// This method returns the corresponding rational number if it's rational,
    /// returns the expansion until the first repeating occurence
    pub fn to_rational(&self) -> Approximation<Ratio<U>> {
        if self.is_rational() {
            Approximation::Exact(self.convergents().last().unwrap())
        } else {
            Approximation::Approximated(
                self.convergents()
                    .nth(self.a_coeffs.len() + self.p_coeffs.len())
                    .unwrap(),
            )
        }
    }
}

impl<
        T: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
        U: Integer + Clone + Signed + CheckedAdd,
    > Computable<U> for ContinuedFraction<T>
{
    fn approximated(&self, limit: &U) -> Approximation<Ratio<U>> {
        let mut convergents = self.convergents();
        let mut last_conv = convergents.next().unwrap();
        if last_conv.denom() > limit {
            let i = self.a_coeffs.first().unwrap().clone();
            return Approximation::Approximated(Ratio::from(i.to_signed()));
        }
        loop {
            last_conv = match convergents.next() {
                Some(v) => {
                    if v.denom() < limit {
                        v
                    } else {
                        return Approximation::Approximated(last_conv);
                    }
                }
                None => return Approximation::Exact(last_conv),
            }
        }
    }
}

impl<T: fmt::Display> fmt::Display for ContinuedFraction<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.negative {
            write!(f, "-")?;
        }

        write!(f, "[{}", self.a_coeffs.first().unwrap())?;
        if self.a_coeffs.len() == 1 {
            if self.p_coeffs.len() == 0 {
                return write!(f, "]");
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

// trait bound free implementations for several simple functions
impl<T: Zero> ContinuedFraction<T> {
    #[inline]
    fn _zero() -> Self {
        ContinuedFraction {
            a_coeffs: vec![T::zero()],
            p_coeffs: Vec::new(),
            negative: false,
        }
    }

    #[inline]
    fn _is_zero(&self) -> bool {
        self.a_coeffs.len() == 1 && self.a_coeffs[0].is_zero() && self.p_coeffs.len() == 0
    }
}

impl<
        T: QuadraticBase + AddAssign + WithUnsigned<Unsigned = U>,
        U: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = T>,
    > Zero for ContinuedFraction<U>
where
    for<'r> &'r T: RefNum<T>,
    for<'r> &'r U: RefNum<U>,
{
    #[inline]
    fn zero() -> Self {
        Self::_zero()
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self._is_zero()
    }
}

impl<
        T: QuadraticBase + AddAssign + WithUnsigned<Unsigned = U>,
        U: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = T>,
    > One for ContinuedFraction<U>
where
    for<'r> &'r T: RefNum<T>,
    for<'r> &'r U: RefNum<U>,
{
    fn one() -> Self {
        ContinuedFraction {
            a_coeffs: vec![U::one()],
            p_coeffs: Vec::new(),
            negative: false,
        }
    }

    fn is_one(&self) -> bool {
        self.a_coeffs.len() == 1 && self.a_coeffs[0].is_one() && self.p_coeffs.len() == 0
    }
}

// TODO: implement from_float, using FloatCore trait?

impl<T: Num + WithUnsigned<Unsigned = U> + PartialOrd, U> From<T> for ContinuedFraction<U> {
    fn from(t: T) -> Self {
        if t < T::zero() {
            ContinuedFraction {
                a_coeffs: vec![(T::zero() - t).to_unsigned()],
                p_coeffs: Vec::new(),
                negative: true,
            }
        } else {
            ContinuedFraction {
                a_coeffs: vec![t.to_unsigned()],
                p_coeffs: Vec::new(),
                negative: false,
            }
        }
    }
}

impl<T: Integer + Clone + WithUnsigned<Unsigned = U>, U: Zero> From<Ratio<T>>
    for ContinuedFraction<U>
{
    fn from(r: Ratio<T>) -> Self {
        if r.is_zero() {
            return Self::_zero();
        }

        let mut coeffs = Vec::new();
        let (mut n, mut d) = r.into();

        let negative = if n < T::zero() {
            n = T::zero() - n;
            true
        } else {
            false
        };

        if n < d {
            std::mem::swap(&mut n, &mut d);
            coeffs.push(U::zero());
        }

        while !d.is_zero() {
            let (quo, rem) = n.div_rem(&d);
            coeffs.push(quo.to_unsigned());
            n = d;
            d = rem;
        }

        ContinuedFraction {
            a_coeffs: coeffs,
            p_coeffs: Vec::new(),
            negative,
        }
    }
}

// TODO: expose this after implemented
mod parse {
    use super::{ContinuedFraction, FromStr};

    pub enum ParseContFracError {}

    impl<T> FromStr for ContinuedFraction<T> {
        type Err = ParseContFracError;

        /// Parse from standard format (like 355/113 = "[3; 7, 16]", (1+sqrt(5))/2 = "[1; (1)]")
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            unimplemented!()
        }
    }
}

impl<
        T: Num + PartialOrd + AddAssign + Signed + WithUnsigned<Unsigned = U>,
        U: NumRef + PartialOrd + WithSigned<Signed = T>,
    > Add<T> for ContinuedFraction<U>
where
    for<'r> &'r U: RefNum<U>,
{
    type Output = Self;

    fn add(self, rhs: T) -> Self {
        let mut new_a = self.a_coeffs;
        let i = new_a.first().unwrap();

        // find if the signed is flipped after adding
        let (new_i, flipped) = match (self.negative, rhs.is_negative()) {
            (true, true) => {
                // neg + neg
                let rhs = (-rhs).to_unsigned();
                (rhs + i, false)
            }
            (true, false) => {
                // neg + pos
                let rhs = rhs.to_unsigned();
                if i < &rhs {
                    (rhs - i, true)
                } else {
                    (i - rhs, false)
                }
            }
            (false, true) => {
                // pos + neg
                let rhs = (-rhs).to_unsigned();
                if i < &rhs {
                    (rhs - i, true)
                } else {
                    (i - rhs, false)
                }
            }
            (false, false) => {
                // pos + pos
                (i + rhs.to_unsigned(), false)
            }
        };

        if flipped {
            // delegate to the constructor to handle the negative sign
            let mut new_a: Vec<T> = new_a.into_iter().map(|u| u.to_signed()).collect();
            let new_p: Vec<T> = self.p_coeffs.into_iter().map(|u| u.to_signed()).collect();
            *new_a.first_mut().unwrap() = -new_i.to_signed();
            Self::new(new_a, new_p, self.negative)
        } else {
            *new_a.first_mut().unwrap() = new_i;

            let mut result = ContinuedFraction {
                a_coeffs: new_a,
                p_coeffs: self.p_coeffs,
                negative: self.negative,
            };
            if result._is_zero() {
                // clear negative flag for 0
                result.negative = false;
            }
            result
        }
    }
}

impl<
        T: Num + PartialOrd + AddAssign + Signed + WithUnsigned<Unsigned = U>,
        U: NumRef + PartialOrd + WithSigned<Signed = T>,
    > Sub<T> for ContinuedFraction<U>
where
    for<'r> &'r U: RefNum<U>,
{
    type Output = Self;

    fn sub(self, rhs: T) -> Self::Output {
        self.add(-rhs)
    }
}

impl<T: Zero> Neg for ContinuedFraction<T> {
    type Output = ContinuedFraction<T>;

    fn neg(self) -> Self::Output {
        let mut result = self;
        if !result._is_zero() {
            // don't negate when the number is zero
            result.negative = !result.negative;
        }
        result
    }
}

impl<
        T: QuadraticBase + AddAssign + WithUnsigned<Unsigned = U>,
        U: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = T>,
    > Mul<T> for ContinuedFraction<U>
where
    for<'r> &'r T: RefNum<T>,
    for<'r> &'r U: RefNum<U>,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        if self.is_rational() {
            Self::from(self.to_rational().value() * rhs)
        } else {
            Self::try_from(QuadraticSurd::<T>::from(self) * rhs).unwrap()
        }
    }
}

impl<
        T: QuadraticBase + AddAssign + WithUnsigned<Unsigned = U>,
        U: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = T>,
    > Div<T> for ContinuedFraction<U>
where
    for<'r> &'r T: RefNum<T>,
    for<'r> &'r U: RefNum<U>,
{
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        if self.is_rational() {
            Self::from(self.to_rational().value() / rhs)
        } else {
            Self::try_from(QuadraticSurd::<T>::from(self) / rhs).unwrap()
        }
    }
}

macro_rules! impl_binop_for_ratio_surd {
    (impl $imp:ident, $method:ident) => {
        impl<
                T: QuadraticBase + AddAssign + WithUnsigned<Unsigned = U>,
                U: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = T>,
            > $imp<Ratio<T>> for ContinuedFraction<U>
        where
            for<'r> &'r T: RefNum<T>,
            for<'r> &'r U: RefNum<U>,
        {
            type Output = Self;

            fn $method(self, rhs: Ratio<T>) -> Self {
                if self.is_rational() {
                    Self::from(self.to_rational().value().$method(rhs))
                } else {
                    Self::try_from(QuadraticSurd::<T>::from(self).$method(rhs)).unwrap()
                }
            }
        }

        impl<
                T: QuadraticBase + AddAssign + WithUnsigned<Unsigned = U>,
                U: Integer + Clone + NumRef + CheckedAdd + CheckedMul + WithSigned<Signed = T>,
            > $imp<ContinuedFraction<U>> for ContinuedFraction<U>
        where
            for<'r> &'r T: RefNum<T>,
            for<'r> &'r U: RefNum<U>,
        {
            type Output = Self;

            fn $method(self, rhs: ContinuedFraction<U>) -> Self {
                if rhs.is_rational() {
                    self.$method(rhs.to_rational().value())
                } else {
                    Self::try_from(
                        QuadraticSurd::<T>::from(self).$method(QuadraticSurd::<T>::from(rhs)),
                    )
                    .unwrap()
                }
            }
        }
    };
}

impl_binop_for_ratio_surd!(impl Add, add);
impl_binop_for_ratio_surd!(impl Sub, sub);
impl_binop_for_ratio_surd!(impl Mul, mul);
impl_binop_for_ratio_surd!(impl Div, div);

/// Represents a simple continued fraction with infinite
/// coefficients. It's a wrapper of an iterator that returns the continued fraction coefficients.
/// Most operations of this struct will also return an iterator wrapped by this struct.
#[derive(Clone, Copy)]
pub struct InfiniteContinuedFraction<I: Iterator>(pub I);

impl<I: Iterator<Item = T> + Clone, T: Num + Clone> InfiniteContinuedFraction<I> {
    pub fn generalize(self) -> std::iter::Zip<I, std::iter::Repeat<T>> {
        self.0.zip(std::iter::repeat(T::one()))
    }
}

impl<I: Iterator<Item = T>, T: Integer + NumRef + Clone> InfiniteContinuedFraction<I>
where
    for<'r> &'r T: RefNum<T>,
{
    /// This method returns a homographic function result on the fraction
    /// A homographic function is `(ax + b)/(cx + d)`
    pub fn homo(
        self,
        a: T,
        b: T,
        c: T,
        d: T,
    ) -> InfiniteContinuedFraction<HomographicResult<I, T>> {
        InfiniteContinuedFraction(HomographicResult {
            block: Block::new(a, b, c, d),
            coeffs: self.0,
        })
    }

    /// This method returns a bihomographic function result on the fraction
    /// A bihomographic function is `(axy + bx + cy + d)/(exy + fx + gy + h)`
    pub fn bihomo<U: Iterator<Item = T>>(
        self,
        rhs: InfiniteContinuedFraction<U>,
        a: T,
        b: T,
        c: T,
        d: T,
        e: T,
        f: T,
        g: T,
        h: T,
    ) -> InfiniteContinuedFraction<BihomographicResult<I, U, T>> {
        InfiniteContinuedFraction(BihomographicResult {
            block: DualBlock::new(a, b, c, d, e, f, g, h),
            x_coeffs: self.0,
            y_coeffs: rhs.0,
        })
    }

    /// Take first N coefficients in the sequence and turn it into a `ContinuedFraction` object.
    pub fn take<U>(self, count: usize) -> ContinuedFraction<U>
    where
        T: WithUnsigned<Unsigned = U> + AddAssign,
    {
        ContinuedFraction::new(self.0.take(count).collect(), Vec::new(), false)
    }

    /// Take first N coefficients in the sequence and turn it into a
    /// `ContinuedFraction` object with periodic detection.
    fn take_periodic<U>(self, count: usize) -> ContinuedFraction<U> {
        // TODO: detect re-occurence of convergents
        unimplemented!()
    }
}

/// Iterator of [InfiniteContinuedFraction::homo()] result
#[derive(Debug, Clone, Copy)]
pub struct HomographicResult<I: Iterator<Item = T>, T> {
    block: Block<T>,
    coeffs: I,
}

impl<I: Iterator<Item = T>, T: Integer + NumRef> Iterator for HomographicResult<I, T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            match self.block.reduce_recip() {
                Some(i) => break Some(i),
                None => match self.coeffs.next() {
                    Some(v) => self.block.rmove(v),
                    None => break None,
                },
            }
        }
    }
}

/// Iterator of [InfiniteContinuedFraction::bihomo()] result
#[derive(Debug, Clone, Copy)]
pub struct BihomographicResult<X: Iterator<Item = T>, Y: Iterator<Item = T>, T> {
    block: DualBlock<T>,
    x_coeffs: X,
    y_coeffs: Y,
}

impl<X: Iterator<Item = T>, Y: Iterator<Item = T>, T: Integer + NumRef> Iterator
    for BihomographicResult<X, Y, T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            match self.block.reduce_recip() {
                Ok(i) => break Some(i),
                Err((right, down)) => {
                    if right {
                        match self.x_coeffs.next() {
                            Some(v) => self.block.rmove_right(v),
                            None => break None,
                        }
                    }
                    if down {
                        match self.y_coeffs.next() {
                            Some(v) => self.block.rmove_down(v),
                            None => break None,
                        }
                    }
                }
            }
        }
    }
}

impl<I: Iterator<Item = T>, T: Integer + NumRef + Clone> Add<T> for InfiniteContinuedFraction<I>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = InfiniteContinuedFraction<HomographicResult<I, T>>;

    fn add(self, rhs: T) -> Self::Output {
        self.homo(T::one(), rhs, T::zero(), T::one())
    }
}

// TODO: implement sqrt for InfiniteContinuedFraction
// REF: https://crypto.stanford.edu/pbc/notes/contfrac/algebraic.html
//      https://github.com/blynn/frac/blob/master/newton.c#L74

/// This trait provide conversion from iterator of ASCII chars to
/// continued fraction. This can be used for accepting high-precision
/// decimal, or infinite continued fraction representation
trait ParseContinuedFraction {
    // TODO: implement followings and make the trait public
    // fn parse_as_decimals() -> InfiniteContinuedFraction;
    // fn parse_as_cfrac() -> InfiniteContinuedFraction;
}

impl<I: Iterator<Item = u8>> ParseContinuedFraction for I {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbols::E;

    #[test]
    fn cont_frac_creation_test() {
        let cf = ContinuedFraction::new(vec![0, 1, 2], vec![3, 4], false);
        assert_eq!(
            cf,
            ContinuedFraction::new(vec![0, 1, 0, 0, 2], vec![3, 4], false)
        );
        assert_eq!(
            cf,
            ContinuedFraction::new(vec![0, 1, 1, 0, 1], vec![3, 4], false)
        );

        let cf = ContinuedFraction::new(vec![2, 3, 4], vec![], true);
        assert_eq!(
            cf,
            ContinuedFraction::new(vec![2, 1, 0, 2, 4], vec![], true)
        );
        assert_eq!(
            cf,
            ContinuedFraction::new(vec![3, -1, -2, -4], vec![], true)
        );
        assert_eq!(cf, ContinuedFraction::new(vec![-3, 1, 2, 4], vec![], false));
    }

    #[test]
    fn cont_frac_iter_test() {
        let one = ContinuedFraction::<u32>::new(vec![1], vec![], false);
        assert_eq!(one.coeffs().cloned().collect::<Vec<_>>(), vec![1]);
        assert_eq!(one.convergents().collect::<Vec<_>>(), vec![Ratio::from(1)]);

        let n_one = ContinuedFraction::<u32>::new(vec![1], vec![], true);
        assert_eq!(n_one.coeffs().cloned().collect::<Vec<_>>(), vec![1]);
        assert_eq!(
            n_one.convergents().collect::<Vec<_>>(),
            vec![Ratio::from(-1)]
        );

        let sq2 = ContinuedFraction::<u32>::new(vec![1], vec![2], false);
        assert_eq!(
            sq2.coeffs().take(5).cloned().collect::<Vec<_>>(),
            vec![1, 2, 2, 2, 2]
        );
        assert_eq!(
            sq2.convergents().take(5).collect::<Vec<_>>(),
            vec![
                Ratio::from(1),
                Ratio::new(3, 2),
                Ratio::new(7, 5),
                Ratio::new(17, 12),
                Ratio::new(41, 29)
            ]
        );

        let n_sq2 = ContinuedFraction::<u32>::new(vec![1], vec![2], true);
        assert_eq!(
            n_sq2.coeffs().take(5).cloned().collect::<Vec<_>>(),
            vec![1, 2, 2, 2, 2]
        );
        assert_eq!(
            n_sq2.convergents().take(5).collect::<Vec<_>>(),
            vec![
                Ratio::from(-1),
                Ratio::new(-3, 2),
                Ratio::new(-7, 5),
                Ratio::new(-17, 12),
                Ratio::new(-41, 29)
            ]
        );
    }

    #[test]
    fn cont_frac_conversion_test() {
        assert_eq!(
            ContinuedFraction::from(Ratio::from(3)),
            ContinuedFraction::new(vec![3u32], vec![], false)
        );
        assert_eq!(
            ContinuedFraction::from(Ratio::new(22, 7)),
            ContinuedFraction::new(vec![3u32, 7], vec![], false)
        );
        assert_eq!(
            ContinuedFraction::from(Ratio::new(-22, 7)),
            ContinuedFraction::new(vec![3i32, 7], vec![], true)
        );
        assert_eq!(
            ContinuedFraction::from(Ratio::new(7, 22)),
            ContinuedFraction::new(vec![0u32, 3, 7], vec![], false)
        );
        assert_eq!(
            ContinuedFraction::from(Ratio::new(-7, 22)),
            ContinuedFraction::new(vec![0i32, 3, 7], vec![], true)
        );
        assert_eq!(
            ContinuedFraction::from(Ratio::new(355, 113)),
            ContinuedFraction::new(vec![3u32, 7, 16], vec![], false)
        );
    }

    #[test]
    fn fmt_test() {
        assert_eq!(
            format!("{}", ContinuedFraction::new(vec![1], vec![], false)),
            "[1]"
        );
        assert_eq!(
            format!("{}", ContinuedFraction::new(vec![1, 2, 3], vec![], false)),
            "[1; 2, 3]"
        );
        assert_eq!(
            format!("{}", ContinuedFraction::new(vec![1], vec![2], false)),
            "[1; (2)]"
        );
        assert_eq!(
            format!(
                "{}",
                ContinuedFraction::new(vec![1, 2, 3], vec![3, 2], false)
            ),
            "[1; 2, 3, (3, 2)]"
        );
    }

    #[test]
    fn cont_frac_arithmetic_test() {
        // TODO: add more tests to cover all cases

        let one = ContinuedFraction::one();
        let n_one = ContinuedFraction::<u32>::new(vec![1], vec![], true);
        let sq2 = ContinuedFraction::new(vec![1], vec![2], false);
        let n_sq2 = ContinuedFraction::new(vec![1], vec![2], true);

        // Add
        assert_eq!(one.clone() + 1, ContinuedFraction::from(2));
        assert_eq!(one.clone() + (-1), ContinuedFraction::zero());
        assert_eq!(n_one.clone() + 1, ContinuedFraction::zero());
        assert_eq!(n_one.clone() + (-1), ContinuedFraction::from(-2));

        assert_eq!(
            sq2.clone() + 1,
            ContinuedFraction::new(vec![2u32], vec![2], false)
        );
        assert_eq!(
            sq2.clone() + (-1),
            ContinuedFraction::new(vec![0u32], vec![2], false)
        );
        assert_eq!(
            n_sq2.clone() + 1,
            ContinuedFraction::new(vec![0i32], vec![2], true)
        );
        assert_eq!(
            n_sq2.clone() + (-1),
            ContinuedFraction::new(vec![2i32], vec![2], true)
        );

        // Mul
        assert_eq!(one.clone() * 2, ContinuedFraction::from(2));
        assert_eq!(one.clone() * -2, ContinuedFraction::from(-2));
        assert_eq!(n_one.clone() * 2, ContinuedFraction::from(-2));
        assert_eq!(n_one.clone() * -2, ContinuedFraction::from(2));

        assert_eq!(
            sq2.clone() * 2,
            ContinuedFraction::new(vec![2u32], vec![1, 4], false)
        );
        assert_eq!(
            sq2.clone() * -2,
            ContinuedFraction::new(vec![2i32], vec![1, 4], true)
        );
        assert_eq!(
            n_sq2.clone() * 2,
            ContinuedFraction::new(vec![2i32], vec![1, 4], true)
        );
        assert_eq!(
            n_sq2.clone() * -2,
            ContinuedFraction::new(vec![2u32], vec![1, 4], false)
        );
    }

    #[test]
    fn inf_cont_frac_utility_test() {
        let e = E {};

        // take()
        assert_eq!(
            (e.cfrac::<i32>() + (-2)).take(5),
            ContinuedFraction::new(vec![0i32, 1, 2, 2], vec![], false)
        );
        assert_eq!(
            (e.cfrac::<i32>() + (-3)).take(5),
            ContinuedFraction::new(vec![0i32, 3, 1, 1, 4], vec![], true)
        );
    }

    #[test]
    fn inf_cont_frac_arithmetic_test() {
        let e = E {};
        let ep1_cf = e.cfrac::<i32>() + 1;
        assert_eq!(ep1_cf.0.take(5).collect::<Vec<_>>(), vec![3, 1, 2, 1, 1]);

        let sq2 = ContinuedFraction::<u32>::new(vec![1i32], vec![2], false);
        let sq2p1 = sq2.clone() + 1;
        assert_eq!(
            (sq2.expanded() + 1).0.take(5).collect::<Vec<_>>(),
            sq2p1.expanded().0.take(5).collect::<Vec<_>>()
        );
    }
}
