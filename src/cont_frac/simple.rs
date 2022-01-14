use core::str::FromStr;
use std::ops::{Add, AddAssign};
use std::mem::swap;
use std::fmt;
use num_traits::{float::FloatCore, Num, Signed, NumRef, RefNum, CheckedAdd, CheckedMul};
use num_integer::{Integer};
use num_rational::Ratio;
use crate::traits::{RationalApproximation, Approximation, WithSigned};

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

impl<I: Iterator<Item = T>, T: WithSigned<Signed = U>, U: Signed>
Iterator for GeneralCoefficients<I> {
    type Item = (U, U);

    fn next(&mut self) -> Option<Self::Item> {
        if self.negative {
            self.negative = false; // only apply to the first coefficient
            self.coeffs.next().map(|v| (-v.to_signed(), U::one()))
        } else {
            self.coeffs.next().map(|v| (v.to_signed(), U::one()))
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

        ContinuedFraction { a_coeffs: dedup_a, p_coeffs, negative }
    }

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
    pub fn generalize(&self) -> GeneralCoefficients<Coefficients<T>> {
        GeneralCoefficients {
            coeffs: self.coeffs(), negative: self.negative
        }
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

impl<T> ContinuedFraction<T> {
    /// TODO: implement from_float, from_rational as From traits
    pub fn from_float<U: FloatCore>(f: U) -> Option<Self> {
        unimplemented!()
    }

    /// Return None if bit size of T is not enough
    pub fn from_rational(f: Ratio<T>) -> Option<Self> {
        unimplemented!()
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

impl<T: AddAssign> Add<T> for ContinuedFraction<T> {
    type Output = Self;

    fn add(self, rhs: T) -> Self {
        let mut new_a = self.a_coeffs;
        let i = new_a.first_mut().unwrap();
        *i += rhs;
        ContinuedFraction { a_coeffs: new_a, p_coeffs: self.p_coeffs, negative: self.negative }
    }
}

impl<T> Add<Ratio<T>> for ContinuedFraction<T> {
    type Output = Self;

    fn add(self, rhs: Ratio<T>) -> Self { 
        unimplemented!() // TODO: implement by converting to ratio or quadratic surd
    }
}

// impl<T> Add<ContinuedFraction<T>> for ContinuedFraction<T> {
//     type Output = ContinuedFraction<T>;

//     fn add(self, rhs: ContinuedFraction<T>) -> Self::Output { unimplemented!() }
// }

/// This trait represents a regular continued fraction with infinite
/// coefficients. All operations here will be done iteratively
#[derive(Clone)]
pub struct InfiniteContinuedFraction<T: Iterator> (T);

impl<T: Iterator<Item = U> + Clone, U: Num + Clone> InfiniteContinuedFraction<T> {
    pub fn generalize(&self) -> std::iter::Zip<T, std::iter::Repeat<U>> {
        self.0.zip(std::iter::repeat(U::one()))
    }

    /// This method returns a homographic function result on the fraction
    /// A homographic function is `(ax + b)/(cx + d)`
    pub fn homo(&self, a: T::Item, b: T::Item, c: T::Item, d: T::Item)
        -> InfiniteContinuedFraction<HomographicResult<T>> {
        InfiniteContinuedFraction(
            HomographicResult { pm1: a, pm2: b, qm1: c, qm2: d, coeffs: self.coeffs.clone() }
        )
    }

    /// This method returns a bihomographic function result on the fraction
    /// A bihomographic function is `(axy + bx + cy + d)/(exy + fx + gy + h)`
    pub fn bihomo<U: Iterator>(
        &self, rhs: InfiniteContinuedFraction<U>,
        a: T::Item, b: T::Item, c: T::Item, d: T::Item,
        e: T::Item, f: T::Item, g: T::Item, h: T::Item)
        -> InfiniteContinuedFraction<BihomographicResult<T, U>> {
        InfiniteContinuedFraction {
            coeffs: BihomographicResult {
                a, b, c, d, e, f, g, h,
                x_coeffs: self.coeffs.clone(), y_coeffs: rhs.coeffs
            }, negative: self.negative
        }

    }
}

#[derive(Debug)]
pub struct HomographicResult<T: Iterator> {
    pm1: T::Item, pm2: T::Item, qm1: T::Item, qm2: T::Item,
    coeffs: T
}

impl<T: Iterator<Item = U>, U: Integer + Num + NumRef + Clone> Iterator for HomographicResult<T>
where for <'r> &'r U: RefNum<U> {
    type Item = U;

    fn next(&mut self) -> Option<U> {
        loop {
            match self.coeffs.next() {
                Some(v) => {
                        let p = v.clone() * &self.pm1 + &self.pm2;
                        let q = v * &self.qm1 + &self.qm2;
                        
                        let i = if q.is_zero() { U::zero() } else { p.div_floor(&q) };
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

pub struct BihomographicResult<T: Iterator, U: Iterator> {
    a: T::Item, b: T::Item, c: T::Item, d: T::Item,
    e: T::Item, f: T::Item, g: T::Item, h: T::Item,
    x_coeffs: T, y_coeffs: U
}

// For arithmetics on ContinuedFraction:
// REF: http://inwap.com/pdp10/hbaker/hakmem/cf.html
//      https://www.plover.com/~mjd/cftalk/
//      http://www.idosi.org/aejsr/10(5)15/1.pdf
// TODO: implement arithmetics for InfiniteContinuedFraction
// TODO: implement InfiniteContinuedFraction to ContinuedFraction (with given iteration limit)
// TODO: implement GeneralContinuedFraction to InfiniteContinuedFraction

#[cfg(test)]
mod tests {
    use super::*;

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
    fn fmt_test() {
       assert_eq!(format!("{}", ContinuedFraction::new(vec![1], vec![], false)), "[1]");
       assert_eq!(format!("{}", ContinuedFraction::new(vec![1, 2, 3], vec![], false)), "[1; 2, 3]");
       assert_eq!(format!("{}", ContinuedFraction::new(vec![1], vec![2], false)), "[1; (2)]");
       assert_eq!(format!("{}", ContinuedFraction::new(vec![1, 2, 3], vec![3, 2], false)), "[1; 2, 3, (3, 2)]");
    }
}
