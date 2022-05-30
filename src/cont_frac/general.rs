//! Implementation of general continued fractions

use super::block::Block;
use crate::traits::{Approximation, Computable};
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{CheckedAdd, CheckedMul, FromPrimitive, Num, NumRef, RefNum, Signed, ToPrimitive};

/// This struct defines utility functions for generalized continued fraction number
/// `b_1 + a_2 / (b_2 + a_3 / (b_3 + a_4 / .. ))`.
/// 
/// It can be converted from any iterator that returns a pair of number (a_k, b_k).
/// You need to make sure that a_1 = 1.
/// 
pub trait GeneralContinuedFraction<T: Integer + NumRef>: Iterator<Item = (T, T)>
where
    for<'r> &'r T: RefNum<T>,
{
    /// Compute the convergents of the generalized continued fraction
    fn convergents(self) -> Convergents<Self, T> where Self: Sized {
        Convergents {
            block: Block::identity(),
            g_coeffs: self,
        }
    }

    /// Simplify the generalized continued fraction to an `InfiniteContinuedFraction`
    /// Note that usually this function will generate a simple continued fraction with most
    /// coefficients being positive. If you want to ensure the positiveness, you can either
    /// call collect() on `InfiniteContinuedFraction`, or call simplify() again.
    fn simplify(self) -> Simplified<Self, T> where Self: Sized {
        Simplified {
            block: Block::identity(),
            g_coeffs: self,
        }.into()
    }

    /// Retrieve the decimal representation of the number, as an iterator of digits.
    /// The iterator will stop if the capacity of T is reached
    fn decimals(self) -> DecimalDigits<Self, T> where Self: Sized {
        DecimalDigits {
            block: Block::identity(),
            g_coeffs: self,
        }
    }

    // XXX: we can also implement homo and bihomo function on general continued fraction
    //      however the result will still be an InfiniteContinuedFraction
}
impl<I, T> GeneralContinuedFraction<T> for I where I: Iterator<Item = (T, T)>, T: Integer + NumRef, for<'r> &'r T: RefNum<T>, {}

/// Iterator of convergents of [GeneralContinuedFraction]
#[derive(Debug, Clone)]
pub struct Convergents<I: Iterator<Item = (T, T)> + ?Sized, T> {
    block: Block<T>,
    g_coeffs: I,
}

/// Iterator of [GeneralContinuedFraction::simplify()] result
#[derive(Debug, Clone)]
pub struct Simplified<I: Iterator<Item = (T, T)> + ?Sized, T> {
    block: Block<T>,
    g_coeffs: I,
}

/// Iterator of [GeneralContinuedFraction::decimals()] result
#[derive(Debug, Clone)]
pub struct DecimalDigits<I: Iterator<Item = (T, T)> + ?Sized, T> {
    block: Block<T>,
    g_coeffs: I,
}

impl<I: Iterator<Item = (T, T)>, T: Integer + NumRef + CheckedAdd + CheckedMul + Clone> Iterator
    for Convergents<I, T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Item = Ratio<T>;

    fn next(&mut self) -> Option<Ratio<T>> {
        let (a, b) = self.g_coeffs.next()?;
        let (p, q) = self.block.checked_gmove(a, b)?;
        self.block.update(p.clone(), q.clone());

        Some(Ratio::new(p, q))
    }
}

impl<I: Iterator<Item = (T, T)>, T: Integer + NumRef> Iterator for Simplified<I, T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            match self.block.reduce_recip() {
                Some(i) => break Some(i),
                None => match self.g_coeffs.next() {
                    Some((a, b)) => self.block.gmove(a, b),
                    None => break None,
                },
            }
        }
    }
}

impl<
        I: Iterator<Item = (T, T)>,
        T: Integer + NumRef + CheckedAdd + CheckedMul + FromPrimitive + ToPrimitive,
    > Iterator for DecimalDigits<I, T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        // TODO: print point and handle signed & integer > 10 case
        loop {
            match self.block.reduce_mul(T::from_u8(10).unwrap()) {
                Some(i) => break Some(i.to_u8().unwrap() + b'0'),
                None => match self.g_coeffs.next() {
                    Some((a, b)) => {
                        let (p, q) = self.block.checked_gmove(a, b)?;
                        self.block.update(p, q);
                    }
                    None => break None,
                },
            }
        }
    }
}

impl<I: Iterator<Item = (T, T)> + Clone, T: Integer + NumRef + Clone + CheckedAdd + CheckedMul>
    Computable<T> for I
where
    for<'r> &'r T: RefNum<T>,
{
    fn approximated(&self, limit: &T) -> Approximation<Ratio<T>> {
        let mut convergents = self.clone().convergents();
        let mut last_conv = convergents.next().unwrap();
        if last_conv.denom() > limit {
            return Approximation::Approximated(last_conv);
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

// TODO: implement continued fraction for various functions
// REF: https://crypto.stanford.edu/pbc/notes/contfrac/cheat.html
// List: coth(1/n), tan(1/n), arctan(n), tanh(n), tan(n), log(1+x), exp(2m/n), exp(1/n)

/// Iterator of simple continued fraction coefficients generated by [exp()]
#[derive(Debug, Clone, Copy)]
pub struct ExpCoefficients<T> {
    exponent: T,
    i: T,
}

impl<T: Num + NumRef + Clone + Signed> Iterator for ExpCoefficients<T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Item = (T, T);

    fn next(&mut self) -> Option<Self::Item> {
        let result = if self.i.is_zero() {
            Some((T::one(), T::one()))
        } else if self.i.is_one() {
            Some((self.exponent.clone(), T::one()))
        } else {
            Some((
                -((&self.i - T::one()) * &self.exponent),
                &self.i + &self.exponent,
            ))
        };

        self.i = &self.i + T::one();
        result
    }
}

pub fn exp<T: Num + NumRef + Clone + Signed>(target: T) -> ExpCoefficients<T>
where
    for<'r> &'r T: RefNum<T>,
{
    ExpCoefficients {
        exponent: target,
        i: T::zero(),
    }
}

// TODO: implement operators to caculate HAKMEM Constant (and add as a example)
// https://crypto.stanford.edu/pbc/notes/contfrac/hakmem.html

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbols::{Pi, E};

    #[test]
    fn general_conf_frac_test() {
        let e = E {};
        assert_eq!(e.cfrac::<i32>().take(10).collect::<Vec<_>>(), exp(1i32).simplify().take(10).collect::<Vec<_>>());

        let pi = Pi {};
        assert_eq!(
            pi.gfrac()
                .convergents()
                .skip(1)
                .take(5)
                .map(|c| c.into())
                .collect::<Vec<(i32, i32)>>(),
            vec![(4, 1), (3, 1), (19, 6), (160, 51), (1744, 555)]
        );
        assert_eq!(
            pi.gfrac().simplify().into_iter().take(6).collect::<Vec<u64>>(),
            vec![3, 7, 15, 1, 292, 1]
        );
        assert_eq!(
            pi.gfrac::<u32>().decimals().take(20).collect::<Vec<_>>(),
            b"31415926"
        );
    }

    #[test]
    fn functions_test() {
        assert_eq!(
            exp(1i32).into_iter().take(5).collect::<Vec<_>>(),
            vec![(1i32, 1), (1, 1), (-1, 3), (-2, 4), (-3, 5)]
        )
    }
}
