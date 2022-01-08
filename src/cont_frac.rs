
pub struct ContinuedFraction<T> {
    /// Coefficients of aperiodic part
    a_coeffs: Vec<T>,

    /// Coefficients of periodic part
    p_coeffs: Vec<T>
}

// TODO: quadratic surd (include golden ratio) can be represented by ContinuedFraction

pub struct PatternContinuedFraction<T> {
    /// Coefficients of aperiodic part
    a_coeffs: Vec<T>,

    /// Pattern function to generate remaining coefficients
    pattern: Fn(u32) -> T
}

// TODO: e and pi can be represented by ContinuedFraction
