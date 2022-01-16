use core::str::FromStr;
use std::ops::{Add, Mul};
use std::fmt;
use num_traits::{Num, One, Zero, Signed, NumRef, RefNum, CheckedAdd, CheckedMul};
use num_integer::{Integer};
use num_rational::Ratio;
use super::block::{Block, DualBlock};
use crate::traits::{Computable, Approximation, WithSigned, WithUnsigned};
use crate::quad_surd::{QuadraticSurd, QuadraticSurdBase};

/// This struct represents a simple continued fraction `a0 + 1/(a1 + 1/ (a2 + ...))`
/// Where a0, a1, a2 are positive integers
/// It's capable of representing rational numbers and quadratic surds
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
            self.coeffs.next().map(|v| (U::one(), -v.clone().to_signed()))
        } else {
            self.coeffs.next().map(|v| (U::one(), v.clone().to_signed()))
        }
    }
}

pub struct Convergents<'a, T> {
    coeffs: Coefficients<'a, T>, block: Block<T>,
    neg: bool // store the sign
}

impl<'a, T: Integer + Clone + CheckedAdd + CheckedMul + WithSigned<Signed = U>,
         U: Integer + Clone + Signed>
Iterator for Convergents<'a, T> {
    type Item = Ratio<U>;

    fn next(&mut self) -> Option<Self::Item> {
        let a = self.coeffs.next()?;
        let (p, q) = self.block.checked_rmove(a.clone())?;
        self.block.update(p.clone(), q.clone());

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
            coeffs: self.coeffs(), block: Block::identity(),
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
Computable<U> for ContinuedFraction<T>
{
    fn approximated(&self, limit: &U) -> Approximation<Ratio<U>> {
        let mut convergents = self.convergents();
        let mut last_conv = convergents.next().unwrap();
        if last_conv.denom() > limit { 
            let i = self.a_coeffs.first().unwrap().clone();
            return Approximation::Approximated(Ratio::from(i.to_signed()))
        }
        loop {
            last_conv = match convergents.next() {
                Some(v) => if v.denom() < limit { v }
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

        // find if the signed is flipped after adding
        let (new_i, flipped) = if self.negative {
            if rhs < T::zero() { // neg + neg
                let rhs = (T::zero() - rhs).to_unsigned();
                (rhs + i, false)
            } else { // neg + pos
                let rhs = rhs.to_unsigned();
                if *i < rhs {
                    (rhs - i, true)
                } else {
                    (i - rhs, false)
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

        if flipped {
            unimplemented!() // TODO: use InfiniteContinuedFraction
        } else {
            *new_a.first_mut().unwrap() = new_i;

            let mut result = ContinuedFraction { a_coeffs: new_a, p_coeffs: self.p_coeffs, negative: self.negative };
            if result.is_zero() { result.negative = false; } // clear negative flag for 0
            result
        }
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

impl<I: Iterator<Item = T>, T: Integer + NumRef + Clone> InfiniteContinuedFraction<I>
where for <'r> &'r T: RefNum<T> {
    /// This method returns a homographic function result on the fraction
    /// A homographic function is `(ax + b)/(cx + d)`
    pub fn homo(self, a: T, b: T, c: T, d: T)
        -> InfiniteContinuedFraction<HomographicResult<I, T>> {
        InfiniteContinuedFraction(
            HomographicResult { block: Block::new(a, b, c, d), coeffs: self.0 }
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
            block: DualBlock::new(a, b, c, d, e, f, g, h),
            x_coeffs: self.0, y_coeffs: rhs.0,
        })
    }

    /// Take first N coefficients in the sequence and turn it into a
    /// `ContinuedFraction` object. Note that this method only accept
    /// the coefficients with same sign.
    /// TODO: provide simplity method for InfiniteContinuedFraction to allow conversion to this form
    pub fn take<U>(self, count: usize) -> ContinuedFraction<U>
    where T: WithUnsigned<Unsigned = U> {
        let mut it = self.0;
        let mut coeffs = Vec::with_capacity(count);

        // utility function to get abs value and signum
        // signum is None when v is zero
        let abs_signum = |v: T| {
            if v < T::zero() {
                ((T::zero() - v).to_unsigned(), Some(true))
            } else {
                let s = if v.is_zero() { None } else { Some(false) };
                (v.to_unsigned(), s)
            } 
        };

        // parse the first coeff
        let (a1, mut signum) = if let Some(v) = it.next() {
            abs_signum(v)
        } else {
            panic!("there should be at least one coefficient!")
        };
        coeffs.push(a1);

        // collect the remaining coeffs
        let mut count = count - 1;
        while count > 0 {
            match it.next() {
                Some(v) => {
                    let (a, new_signum) = abs_signum(v);
                    match (signum, new_signum) {
                        (Some(s), Some(ns)) => if s != ns { 
                            panic!("all the coefficient should have the same sign");
                        },
                        (None, Some(ns)) => { signum = Some(ns); } // update sign
                        (_, None) => {}
                    };
                    coeffs.push(a);
                }, None => break
            }
            count -= 1;
        }

        ContinuedFraction { a_coeffs: coeffs, p_coeffs: Vec::new(), negative: signum.unwrap() }
    }

    /// Take first N coefficients in the sequence and turn it into a
    /// `ContinuedFraction` object with periodic detection.
    pub fn take_periodic<U>(self, count: usize) -> ContinuedFraction<U> {
        unimplemented!()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct HomographicResult<I: Iterator<Item = T>, T> {
    block: Block<T>, coeffs: I
}

impl<I: Iterator<Item = T>, T: Integer + NumRef> Iterator for HomographicResult<I, T>
where for <'r> &'r T: RefNum<T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            match self.block.reduce_recip() {
                Some(i) => break Some(i),
                None => match self.coeffs.next() {
                    Some(v) => self.block.rmove(v),
                    None => break None
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct BihomographicResult<X: Iterator<Item = T>, Y: Iterator<Item = T>, T> {
    block: DualBlock<T>,
    x_coeffs: X, y_coeffs: Y,
}

impl<X: Iterator<Item = T>, Y: Iterator<Item = T>, T: Integer + NumRef>
Iterator for BihomographicResult<X, Y, T>
where for <'r> &'r T: RefNum<T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            match self.block.reduce_recip() {
                Ok(i) => break Some(i),
                Err((right, down)) => {
                    if right {
                        match self.x_coeffs.next() {
                            Some(v) => self.block.rmove_right(v),
                            None => break None
                        }
                    }
                    if down {
                        match self.y_coeffs.next() {
                            Some(v) => self.block.rmove_down(v),
                            None => break None
                        }
                    }
                }
            }
        }
    }
}

impl<I: Iterator<Item = T>, T: Integer + NumRef + Clone> Add<T> for InfiniteContinuedFraction<I>
where for <'r> &'r T: RefNum<T> {
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
pub trait ParseContinuedFraction {
    // TODO: implement followings
    // fn parse_as_decimals() -> InfiniteContinuedFraction;
    // fn parse_as_cfrac() -> InfiniteContinuedFraction;
}

impl<I: Iterator<Item = u8>> ParseContinuedFraction for I {

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
    fn inf_cont_frac_utility_test() {
        let e = E {};

        // take()
        assert_eq!((e.cfrac::<i32>() + (-2)).take(5), ContinuedFraction::new(vec![0,1,2,1,1], vec![], false));
        assert_eq!((e.cfrac::<i32>() + (-3)).take(5), ContinuedFraction::new(vec![0,3,1,1,4], vec![], true));
    }

    #[test]
    fn inf_cont_frac_arithmetic_test() {
        let e = E {};
        let ep1_cf = e.cfrac::<i32>() + 1;
        assert_eq!(ep1_cf.0.take(5).collect::<Vec<_>>(), vec![3,1,2,1,1]);

        let sq2 = ContinuedFraction::<u32>::new(vec![1], vec![2], false);
        let sq2p1 = sq2.clone() + 1;
        assert_eq!((sq2.expanded() + 1).0.take(5).collect::<Vec<_>>(), sq2p1.expanded().0.take(5).collect::<Vec<_>>());
    }
}
