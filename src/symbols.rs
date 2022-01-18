//! This module contains several predefined irrational math constants

use num_traits::{One, Num, NumRef, RefNum, Signed, FromPrimitive};
use num_integer::Integer;
use num_rational::Ratio;
use std::ops::{AddAssign};
use super::cont_frac::InfiniteContinuedFraction;

/// √2
pub struct Sqrt2 {}
/// ln2
pub struct Ln2 {}
/// ϕ, the golden ratio
pub struct Phi {}
/// e, Euler's number, the natural logarithmic base
pub struct E {}
/// π, Archimedes' constant
pub struct Pi {}
/// γ, Euler–Mascheroni constant
pub struct Gamma {}
/// G, Catalan's constant
// TODO: https://math.stackexchange.com/questions/3620230/conjectured-continued-fraction-formula-for-catalans-constant
// https://functions.wolfram.com/Constants/Catalan/10/
pub struct G {} 

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

/// The sequence used here is pi = 4/(1+1^2/(3+2^2/(5+..))). The first convergent is 0, but
/// it will converge faster later on.
/// Reference: <https://en.wikipedia.org/wiki/Generalized_continued_fraction#%CF%80>
#[derive(Debug, Clone, Copy)]
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

impl Gamma {
    /// There is no elegant way yet to represent euler constant as a continued
    /// fraction. The formula used here will explode very fast and thus there's
    /// no practical use for it.
    pub fn gfrac<T: Integer + Signed + Clone>(&self) -> GammaCoefficients<T> {
        let two = T::one() + T::one();

        // the initial values for d are starting from n = 2
        GammaCoefficients {
            n: 0, qm1: two, qm2: T::one(), sm1: -Ratio::one(), rm1: Ratio::one()
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GammaCoefficients<T> { n: usize, qm1: T, qm2: T, sm1: Ratio<T>, rm1: Ratio<T> }

impl<T: Integer + Signed + NumRef + Clone + AddAssign + FromPrimitive>
Iterator for GammaCoefficients<T> where for <'r> &'r T: RefNum<T> {
    type Item = (T, T);

    /// Reference: Pilehrood, K. H., & Pilehrood, T. H. (2013). On a continued fraction expansion
    /// for Eulerʼs constant. Journal of Number Theory, 133(2), 769-786.
    fn next(&mut self) -> Option<Self::Item> {
        match self.n {
            0 => {
                self.n = 1;
                return Some((T::one(), T::zero())); // return (1, 0)
            },
            1 => {
                self.n = 2;
                return Some((T::one(), T::from_u8(2).unwrap()));
            },
            _ => {}
        };

        let n = T::from_usize(self.n)?;
        let nm1 = T::from_usize(self.n - 1).unwrap();
        let nm2 = T::from_usize(self.n - 2).unwrap();
        let two = T::from_u8(2).unwrap();

        // Eq(21): q_n = sum_{k=0..n} (binom(n,k)^2 * k!)
        // it satisfies Eq(22): q_{n} = 2n*q_{n-1} - (n-1)^2*q_{n-2} with q0=1, q1=2
        let q = &two * &n * self.qm1.clone() - &nm1 * &nm1 * &self.qm2;

        // Eq(35): s_n = (n-1)^2*s_{n-1} + (n-2)/n*q_{n-1} is a rational
        // s_0 = 1, s_1 = -1
        let s = Ratio::from(&nm1 * &nm1) * &self.sm1 + Ratio::new(nm2.clone(), n.clone()) * &self.qm1;

        // Eq(38): a_n = -s_n/s_{n-1}, b_n=2n+(n-2)/n*q_{n-2}/s_{n-1}
        let a = -&s / &self.sm1;
        let b = Ratio::new(nm2 * &self.qm2, n.clone()) / &self.sm1 + &two * &n;

        // rho_n = n*(n-1)/2*s_{n-1} (n >= 5)
        let rho = match self.n {
            1 | 2 => Ratio::one(),
            3 => Ratio::from(T::from_u8(3).unwrap()),
            4 => Ratio::from(T::from_u8(10).unwrap()),
            _ => {
                Ratio::from(n * nm1 / two) * &self.sm1
            }
        };

        // astar_n = rho_n*rho_{n-1}*a_n, bstar_n = rho_n*b_n
        let astar = (&rho * &self.rm1 * a).to_integer();
        let bstar = (&rho * b).to_integer();

        // store values
        self.n += 1;
        std::mem::swap(&mut self.qm2, &mut self.qm1);
        self.qm1 = q;
        self.sm1 = s;
        self.rm1 = rho;

        Some((astar, bstar))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cont_frac::GeneralContinuedFraction;

    #[test]
    fn cfrac_test() {
        let e = E {};
        assert_eq!(e.cfrac().0.take(10).collect::<Vec<u32>>(), vec![2u32,1,2,1,1,4,1,1,6,1]);

        let pi = Pi {};
        assert_eq!(pi.gfrac().take(5).collect::<Vec<(u32, u32)>>(),
                   vec![(1,0), (4,1), (1,3), (4,5), (9,7)]);
        assert_eq!(pi.gfrac().simplify().0.take(10).collect::<Vec<u64>>(), vec![3,7,15,1,292,1,1,1,2,1]);

        let gamma = Gamma {};
        assert_eq!(gamma.gfrac().take(6).collect::<Vec<(i32, i32)>>(),
                   vec![(1,0), (1,2), (-1,4), (-5,16), (36,59), (-15740, 404)]);
        assert_eq!(gamma.gfrac().simplify().0.take(8).collect::<Vec<i64>>(),
                    vec![0,1,1,2,1,2,1,3]); // [0;1,1,2,1,2,1,3,1,-4,-1,0,..]
    }
}
