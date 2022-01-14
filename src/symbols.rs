use num_traits::{Num, NumRef};
use super::cont_frac::{InfiniteContinuedFraction, GeneralContinuedFraction};

/// This module contains several predefined irrational math constants

pub struct E { }
pub struct Pi { }

impl E {
    pub fn cfrac<T: Num>(&self) -> ECoefficients<T> {
        ECoefficients { i: T::zero(), m: 0 }
    }
}

#[derive(Debug, Clone)]
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cfrac_test() {
        let e = E {};
        assert_eq!(e.cfrac().take(10).collect::<Vec<u32>>(), vec![2u32,1,2,1,1,4,1,1,6,1]);

        // let ee = InfiniteContinuedFraction { coeffs: e.cfrac(), negative: false };
        // println!("{:?}", ee.homo(2, 0, 0, 1).coeffs.take(10).collect::<Vec<_>>());

        let em3 = InfiniteContinuedFraction { coeffs: e.cfrac(), negative: false };
        println!("{:?}", em3.homo(1, -3, 0, 1).coeffs.take(30).collect::<Vec<_>>());
    }
}
