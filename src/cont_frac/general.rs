use num_traits::{Num, NumRef, RefNum, Signed};
use super::simple::InfiniteContinuedFraction;

// REF: https://github.com/MostAwesomeDude/continued/blob/master/continued.py
// REF: http://www.numbertheory.org/gnubc/bc_programs.html
/// This trait defines utility functions for general continued fraction number
/// `a_1 + b_1 / (a_2 + b_2 / (a_3 + b_3 / .. ))`. They are available for any
/// iterator that returns a pair of number. The first value will be regarded
/// as a_k while the second value as b_k.
pub trait GeneralContinuedFraction : Iterator {
    fn convergents(self) -> GeneralConvergents<Self>;
    fn simplify(self) -> Simplified<Self>;
    // fn collect(self) -> ContinuedFraction, need to implement FromIterator
    // fn approx_rational(self, limit) -> Ratio<T>, need to implement RationalApproximation
}

pub struct GeneralConvergents<T: Iterator + ?Sized> {
    pm1: T::Item, pm2: T::Item, qm1: T::Item, qm2: T::Item,
    general_coeffs: T
}

pub struct Simplified<T: ?Sized> { // REF: https://github.com/MostAwesomeDude/continued/blob/master/continued.py#L56
    general_coeffs: T
}

impl<T, U> GeneralContinuedFraction for T
where T: Iterator<Item = (U, U)> {
    fn convergents(self) -> GeneralConvergents<T> {
        unimplemented!()
    }

    fn simplify(self) -> Simplified<T> {
        unimplemented!()
    }
}

// TODO: implement continued fraction for various functions
// REF: https://crypto.stanford.edu/pbc/notes/contfrac/cheat.html 

pub struct ExpCoefficients<T> { exponent: T, i: T }

impl<T: Num + NumRef + Clone + Signed> Iterator for ExpCoefficients<T>
where for <'r> &'r T: RefNum<T> {
    type Item = (T, T);

    fn next(&mut self) -> Option<Self::Item> {
        let result = if self.i.is_zero() {
            Some((T::one(), self.exponent.clone()))
        } else if self.i.is_one() {
            Some((T::one(), -self.exponent.clone() ))
        } else {
            Some((&self.i + &self.exponent, -(&self.i * &self.exponent)))
        };

        self.i = &self.i + T::one();
        result
    }
}

pub fn exp<T: Num + Signed>(target: T) -> ExpCoefficients<T> {
    ExpCoefficients { exponent: target, i: T::zero() }
}
