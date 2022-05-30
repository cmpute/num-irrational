//! Implementation of infinite continued fractions

use super::block::{Block, DualBlock};
use num_integer::Integer;
use num_traits::{NumRef, RefNum, One};
use core::iter as iter;

/// Represents a simple continued fraction with infinite
/// coefficients. It's a wrapper of an iterator that returns the continued fraction coefficients.
/// Most operations of this struct will also return an iterator wrapped by this struct.
pub trait InfiniteContinuedFraction: Iterator {
    fn generalize(self) -> iter::Zip<Self, core::iter::Repeat<Self::Item>> where Self: Sized, Self::Item: One + Clone {
        self.zip(iter::repeat(Self::Item::one()))
    }

    /// This method returns a homographic function result on the fraction
    /// A homographic function is `(ax + b)/(cx + d)`
    fn homo(
        self,
        a: Self::Item,
        b: Self::Item,
        c: Self::Item,
        d: Self::Item,
    ) -> HomographicResult<Self, Self::Item>
    where Self: Sized, Self::Item: Integer + NumRef + Clone,
    for<'r> &'r Self::Item: RefNum<Self::Item>{
        HomographicResult {
            block: Block::new(a, b, c, d),
            coeffs: self,
        }
    }

    /// This method returns a bihomographic function result on the fraction
    /// A bihomographic function is `(axy + bx + cy + d)/(exy + fx + gy + h)`
    fn bihomo<U: IntoIterator<Item = Self::Item>>(
        self,
        rhs: U,
        a: Self::Item,
        b: Self::Item,
        c: Self::Item,
        d: Self::Item,
        e: Self::Item,
        f: Self::Item,
        g: Self::Item,
        h: Self::Item,
    ) -> BihomographicResult<Self, U::IntoIter, Self::Item>
    where Self: Sized, Self::Item: Integer + NumRef + Clone,
    for<'r> &'r Self::Item: RefNum<Self::Item> {
        BihomographicResult {
            block: DualBlock::new(a, b, c, d, e, f, g, h),
            x_coeffs: self,
            y_coeffs: rhs.into_iter(),
        }
    }
}

impl<T: ?Sized> InfiniteContinuedFraction for T where T: Iterator { }

/// Iterator of [InfiniteContinuedFraction::homo()] result
#[derive(Debug, Clone, Copy)]
pub struct HomographicResult<I: Iterator<Item = T>, T> {
    block: Block<T>,
    coeffs: I,
}

impl<I: Iterator<Item = T>, T: Integer + NumRef> Iterator for HomographicResult<I, T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            match self.block.reduce_recip() {
                Some(i) => break Some(i),
                None => match self.coeffs.next() {
                    Some(v) => self.block.rmove(v),
                    None => break None,
                },
            }
        }
    }
}

/// Iterator of [InfiniteContinuedFraction::bihomo()] result
#[derive(Debug, Clone, Copy)]
pub struct BihomographicResult<X: Iterator<Item = T>, Y: Iterator<Item = T>, T> {
    block: DualBlock<T>,
    x_coeffs: X,
    y_coeffs: Y,
}

impl<X: Iterator<Item = T>, Y: Iterator<Item = T>, T: Integer + NumRef> Iterator
    for BihomographicResult<X, Y, T>
where
    for<'r> &'r T: RefNum<T>,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            match self.block.reduce_recip() {
                Ok(i) => break Some(i),
                Err((right, down)) => {
                    if right {
                        match self.x_coeffs.next() {
                            Some(v) => self.block.rmove_right(v),
                            None => break None,
                        }
                    }
                    if down {
                        match self.y_coeffs.next() {
                            Some(v) => self.block.rmove_down(v),
                            None => break None,
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cont_frac::ContinuedFraction;
    use crate::symbols::E;

    #[test]
    fn arithmetic_test() {
        let e = E {};

        // e - 2
        assert_eq!(
            e.cfrac::<i32>().homo(1, -2, 0, 1).take(5).collect::<ContinuedFraction<_>>(),
            ContinuedFraction::new(vec![0i32, 1, 2, 2], vec![], false)
        );
        // e - 3
        assert_eq!(
            e.cfrac::<i32>().homo(1, -3, 0, 1).take(5).collect::<ContinuedFraction<_>>(),
            ContinuedFraction::new(vec![0i32, 3, 1, 1, 4], vec![], true)
        );
        // e + 1
        assert_eq!(
            e.cfrac::<i32>().homo(1, 1, 0, 1).take(5).collect::<Vec<_>>(),
            vec![3, 1, 2, 1, 1]
        );

        let sq2 = ContinuedFraction::<u32>::new(vec![1i32], vec![2], false);
        let sq2p1 = sq2.clone() + 1;
        assert_eq!(
            sq2.coeffs_signed().homo(1, 1, 0, 1).take(5).collect::<Vec<_>>(),
            sq2p1.coeffs_signed().take(5).collect::<Vec<_>>()
        );
    }
}
