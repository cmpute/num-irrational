# num-irrational

This crate provides representations of irrational numbers within following categories:
- Math constants (`pi`, `e`, etc.)
- [Quadratic Numbers](https://en.wikipedia.org/wiki/Algebraic_number#Examples)
    - Quadratic irrational aka. Quadratic surd
    - Quadratic integer
    - Gaussian integer / Eisenstein integer
- [Continued Fraction](https://en.wikipedia.org/wiki/Continued_fraction)
    - Simple (aka. Regular) continued fraction
    - General continued fraction
    - Hurwitz continued fraction

As you can see, the support for irrational number is not limited in the real field, it also support several
numeric types in the complex field (by enabling the `complex` feature). It's based on the `num` creates.
