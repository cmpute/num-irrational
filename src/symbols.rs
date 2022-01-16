use num_traits::{Num, NumRef};
use std::ops::AddAssign;
use super::cont_frac::InfiniteContinuedFraction;

/// This module contains several predefined irrational math constants

pub struct E { }
pub struct Pi { }

impl E {
    pub fn cfrac<T: Num + NumRef + Clone>(&self) -> InfiniteContinuedFraction<ECoefficients<T>> {
        InfiniteContinuedFraction(ECoefficients { i: T::zero(), m: 0 })
    }
}

#[derive(Debug, Clone, Copy)]
/// The sequence used here is the famous pattern, `e`=[2;1,2,1,1,4,1...]
pub struct ECoefficients<T> { i: T, m: u8 }

impl<T: Num + NumRef + Clone> Iterator for ECoefficients<T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        if self.i.is_zero() {
            self.i = T::one() + T::one();
            Some(T::one() + T::one()) // return 2
        } else {
            let result = match self.m {
                1 => Some(self.i.clone()),
                _ => Some(T::one())
            };

            if self.m == 2 {
                self.m = 0;
                self.i = T::one() + T::one() + &self.i;
            } else {
                self.m += 1;
            }

            result
        }
    }
}

impl Pi {
    /// pi has only generalized continued fraction representation
    pub fn gfrac<T: Num>(&self) -> PiCoefficients<T> {
        PiCoefficients { a: T::zero(), b: T::zero() }
    }
}

#[derive(Debug, Clone, Copy)]
/// The sequence used here is pi = 4/(1+1^2/(3+2^2/(5+..))). The first convergent is 0, but
/// it will converge faster later on.
/// Reference: <https://en.wikipedia.org/wiki/Generalized_continued_fraction#%CF%80>
pub struct PiCoefficients<T> { a: T, b: T }

impl<T: Num + NumRef + Clone + AddAssign> Iterator for PiCoefficients<T> {
    type Item = (T, T);

    fn next(&mut self) -> Option<Self::Item> {
        let two = T::one() + T::one();

        if self.a.is_zero() { // first item
            self.a = T::one();
            Some((T::one(), T::zero())) // return (1, 0)
        } else if self.b.is_zero() { // second item
            self.b = T::one() + two.clone();
            Some((two.clone() + two, T::one())) // return (4, 1)
        } else {
            let result = Some((self.a.clone() * self.a.clone(), self.b.clone()));
            self.a += T::one();
            self.b += two;

            result
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cfrac_test() {
        let e = E {};
        assert_eq!(e.cfrac().0.take(10).collect::<Vec<u32>>(), vec![2u32,1,2,1,1,4,1,1,6,1]);

        let pi = Pi {};
        assert_eq!(pi.gfrac().take(5).collect::<Vec<(u32, u32)>>(),
                   vec![(1,0), (4,1), (1,3), (4,5), (9,7)]);
    }
}
