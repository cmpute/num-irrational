use core::str::FromStr;
use num_traits::float::FloatCore;
#[cfg(feature = "num-rational")]
use num_rational::Ratio;

/// This struct represents a simple continued fraction a0 + 1/(a1 + 1/ (a2 + ...))
/// It's capable of representing rational numbers and quadratic surds
/// REF: https://pi.math.cornell.edu/~gautam/ContinuedFractions.pdf
pub struct ContinuedFraction<T> {
    /// Coefficients of aperiodic part
    a_coeffs: Vec<T>,

    /// Coefficients of periodic part
    p_coeffs: Vec<T>
}

// For arithmetics on ContinuedFraction:
// REF: http://inwap.com/pdp10/hbaker/hakmem/cf.html
//      https://www.plover.com/~mjd/cftalk/
//      http://www.idosi.org/aejsr/10(5)15/1.pdf
// We could implement arithmetics on GeneralContinuedFraction

impl<T> ContinuedFraction<T> {
    /// Return None if bit size of T is not enough
    /// TODO: implement from_float, from_rational, from_quad_surd as TryFrom traits
    pub fn from_float<U: FloatCore>(f: U) -> Option<Self> {
        unimplemented!()
    }

    /// Return None if bit size of T is not enough
    #[cfg(feature = "num-rational")]
    pub fn from_rational<U>(f: U) -> Option<Self> {
        unimplemented!()
    }

    /// Convert the continued fraction to GeneralContinuedFraction
    /// TODO: change to Into/From traits
    pub fn generalize(self) { // -> GeneralContinuedFraction<T, _, _> {
        unimplemented!()
    }

    pub fn is_rational(self) -> bool {
        unimplemented!()
    }
}

pub struct ParseContFracError {
}

impl<T> FromStr for ContinuedFraction<T> {
    type Err = ParseContFracError;

    /// Parse from standard format (like 355/113 = "[3; 7, 16]")
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        unimplemented!()
    }
}

// TODO: quadratic surd (include golden ratio) can be represented by ContinuedFraction

pub struct GeneralContinuedFraction<T, FnA: Fn(u32) -> T, FnB: Fn(u32) -> T> {
    /// Pattern function to generate a series
    a_pattern: FnA,

    /// Pattern function to generate b series
    /// If this function returns 0, then the fraction sequence will terminate
    b_pattern: FnB
}

// TODO: e and pi can be represented by GeneralContinuedFraction

// TODO: implement basic arithmetics for ContinuedFraction (without periodic part)
// TODO: implement arithmetics for GeneralContinuedFraction
