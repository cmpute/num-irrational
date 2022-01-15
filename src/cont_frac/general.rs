use num_traits::{Num, NumRef, RefNum, Signed};
use num_integer::Integer;
use num_rational::Ratio;
use super::block::Block;
use super::simple::InfiniteContinuedFraction;
use crate::traits::{RationalApproximation, Approximation};

/// This trait defines utility functions for generalized continued fraction number
/// `b_1 + a_2 / (b_2 + a_3 / (b_3 + a_4 / .. ))`. They are available for any
/// iterator that returns a pair of number. The first value will be regarded
/// as a_k while the second value as b_k. You need to make sure that a_1 = 1.
pub trait GeneralContinuedFraction<T: Integer + NumRef> : Iterator<Item = (T, T)>
where for<'r> &'r T: RefNum<T> {
    /// Compute the convergents of the generalized continued fraction
    fn convergents(self) -> GeneralConvergents<Self, T>;

    /// Simplify the generalized continued fraction to an `InfiniteContinuedFraction`
    fn simplify(self) -> InfiniteContinuedFraction<Simplified<Self, T>> where Self: Sized;

    /// Retrieve the decimal representation of the number, as an iterator of digits.
    /// The iterator will stop if the capacity of T is reached
    fn decimals(self) -> DecimalDigits<Self, T>;

    // TODO: we can also implement homo and bihomo function on general continued fraction
    //       however the result will still be an InfiniteContinuedFraction
}

pub struct GeneralConvergents<I: Iterator<Item = (T, T)> + ?Sized, T> {
    block: Block<T>, g_coeffs: I
}

#[derive(Debug, Clone)]
pub struct Simplified<I: Iterator<Item = (T, T)> + ?Sized, T> {
    block: Block<T>, g_coeffs: I
}

pub struct DecimalDigits<I: Iterator<Item = (T, T)> + ?Sized, T> {
    block: Block<T>, g_coeffs: I
}

impl<I: Iterator<Item = (T, T)>, T: Integer + NumRef> Iterator for Simplified<I, T>
where for<'r> &'r T: RefNum<T> {
    type Item = T;

    // use the magic table method described in https://crypto.stanford.edu/pbc/notes/contfrac/nonsimple.html
    fn next(&mut self) -> Option<T> {
        loop {
            match self.g_coeffs.next() {
                Some((a, b)) => {
                    let (p, q) = self.block.gmove(a, b);
                    if let Some(i) = self.block.reduce_integer(p, q) {
                        break Some(i)
                    }
                }, None => break None
            }
        }
    }
}

impl<I: Iterator<Item = (T, T)>, T: Integer + NumRef> GeneralContinuedFraction<T> for I
where for<'r> &'r T: RefNum<T>{
    fn convergents(self) -> GeneralConvergents<I, T> {
        unimplemented!()
    }

    fn simplify(self) -> InfiniteContinuedFraction<Simplified<Self, T>> {
        InfiniteContinuedFraction(Simplified {
            block: Block::identity(), g_coeffs: self
        })
    }

    fn decimals(self) -> DecimalDigits<Self, T> {
        unimplemented!()
    }
}

impl<I: GeneralContinuedFraction<T> + Iterator<Item = (T, T)>, T: Integer + NumRef>
RationalApproximation<T> for I
where for<'r> &'r T: RefNum<T> {
    fn approx_rational(&self, limit: &T) -> Approximation<Ratio<T>> {
        unimplemented!()
    }
}

// TODO: implement continued fraction for various functions
// REF: https://crypto.stanford.edu/pbc/notes/contfrac/cheat.html

#[derive(Debug, Clone, Copy)]
pub struct ExpCoefficients<T> { exponent: T, i: T }

impl<T: Num + NumRef + Clone + Signed> Iterator for ExpCoefficients<T>
where for <'r> &'r T: RefNum<T> {
    type Item = (T, T);

    fn next(&mut self) -> Option<Self::Item> {
        let result = if self.i.is_zero() {
            Some((T::one(), T::one()))
        } else if self.i.is_one() {
            Some((self.exponent.clone(), T::one()))
        } else {
            Some((-((&self.i - T::one()) * &self.exponent), &self.i + &self.exponent))
        };

        self.i = &self.i + T::one();
        result
    }
}

pub fn exp<T: Num + Signed>(target: T) -> ExpCoefficients<T> {
    ExpCoefficients { exponent: target, i: T::zero() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbols::E;

    #[test]
    fn functions_test() {
        assert_eq!(exp(1).take(5).collect::<Vec<_>>(), vec![(1, 1), (1, 1), (-1, 3), (-2, 4), (-3, 5)])
    }

    #[test]
    fn simplify_test() {
        let e = E {};
        assert_eq!(e.cfrac().take(10).collect::<Vec<i32>>(), exp(1).simplify().0.take(10).collect::<Vec<i32>>());
    }
}