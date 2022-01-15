use num_traits::{Num, NumRef, RefNum, Signed};
use num_integer::Integer;
use std::mem::swap;
use super::simple::InfiniteContinuedFraction;

// REF: https://github.com/MostAwesomeDude/continued/blob/master/continued.py
// REF: http://www.numbertheory.org/gnubc/bc_programs.html
/// This trait defines utility functions for general continued fraction number
/// `b_1 + a_2 / (b_2 + a_3 / (b_3 + a_4 / .. ))`. They are available for any
/// iterator that returns a pair of number. The first value will be regarded
/// as a_k while the second value as b_k. You need to make sure that a_1 = 1.
pub trait GeneralContinuedFraction<T: Integer + NumRef> : Iterator<Item = (T, T)>
where for<'r> &'r T: RefNum<T> {
    fn convergents(self) -> GeneralConvergents<Self, T>;
    fn simplify(self) -> InfiniteContinuedFraction<Simplified<Self, T>> where Self: Sized;
    // fn approx_rational(self, limit) -> Ratio<T>, need to implement RationalApproximation
}

pub struct GeneralConvergents<I: Iterator<Item = (T, T)> + ?Sized, T> {
    pm1: T, pm2: T, qm1: T, qm2: T,
    g_coeffs: I
}

#[derive(Debug, Clone)]
pub struct Simplified<I: Iterator<Item = (T, T)> + ?Sized, T> {
    pm1: T, pm2: T, qm1: T, qm2: T,
    g_coeffs: I
}

impl<I: Iterator<Item = (T, T)>, T: Integer + NumRef> Iterator for Simplified<I, T>
where for<'r> &'r T: RefNum<T> {
    type Item = T;

    // use the magic table method described in https://crypto.stanford.edu/pbc/notes/contfrac/nonsimple.html
    fn next(&mut self) -> Option<T> {
        loop {
            match self.g_coeffs.next() {
                Some((a, b)) => {
                    let mut p = &self.pm1 * &b + &self.pm2 * &a;
                    let mut q = &self.qm1 * b + &self.qm2 * a;
                    let g = p.gcd(&q).gcd(&self.pm1).gcd(&self.qm1);
    
                    if g > T::one() {
                        p = p / &g; q = q / &g;
                        self.pm1 = &self.pm1 / &g;
                        self.qm1 = &self.qm1 / g;
                    }

                    if !self.qm1.is_zero() && !q.is_zero() {
                        let (im1, rm1) = self.pm1.div_mod_floor(&self.qm1);
                        let (i, r) = p.div_mod_floor(&q);
                        if im1 == i {
                            self.pm1 = q;
                            swap(&mut self.pm2, &mut self.qm1); // self.pm2 = self.qm1
                            self.qm1 = r; self.qm2 = rm1;
                            break Some(i)
                        }
                    };

                    swap(&mut self.pm2, &mut self.pm1); // self.pm2 = self.pm1
                    swap(&mut self.qm2, &mut self.qm1); // self.qm2 = self.qm1
                    self.pm1 = p; self.qm1 = q;
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
            pm1: T::one(), pm2: T::zero(), qm1: T::zero(), qm2: T::one(), g_coeffs: self
        })
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
        assert_eq!(exp(1).take(5).collect::<Vec<_>>(), vec![(1, 1), (1, 1), (3, -1), (4, -2), (5, -3)])
    }

    #[test]
    fn simplify_test() {
        let e = E {};
        assert_eq!(e.cfrac().take(10).collect::<Vec<i32>>(), exp(1).simplify().0.take(10).collect::<Vec<i32>>());
    }
}