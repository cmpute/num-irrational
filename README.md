# num-irrational

This crate provides representations of part of irrational numbers with following categories:
- Math constants (`pi`, `e`, etc.)
- [Quadratic Surd](https://en.wikipedia.org/wiki/Quadratic_irrational_number)
- [Continued Fraction](https://en.wikipedia.org/wiki/Continued_fraction)

It's based on the `num` creates.

# Examples

```rust
use num_irrational::{QuadraticSurd, ContinuedFraction};

let sq2 = QuadraticSurd::from_sqrt(2);
println!("Square root of 2: {}", sq2); // √2

let sq2_approx = sq2.approx_rational(100);
println!("Approximation under 100: {}", sq2_approx); // 99/70

let sq2_fraction = ContinuedFraction::from(sq2);
println!("Continued Fraction: {}", sq2_fraction); // [1; (2)]
```
