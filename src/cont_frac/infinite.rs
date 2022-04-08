//! Implementation of infinite continued fractions

use super::block::{Block, DualBlock};
use super::simple::ContinuedFraction;
use crate::traits::WithUnsigned;
use num_integer::Integer;
use num_traits::{Num, NumRef, RefNum};
use std::ops::{Add, AddAssign};

// TODO(vnext): move to a separate file
/// Represents a simple continued fraction with infinite
/// coefficients. It's a wrapper of an iterator that returns the continued fraction coefficients.
/// Most operations of this struct will also return an iterator wrapped by this struct.
#[derive(Clone, Copy)]
pub struct InfiniteContinuedFraction<I: Iterator>(pub I);

impl<I: Iterator<Item = T> + Clone, T: Num + Clone> InfiniteContinuedFraction<I> {
    pub fn generalize(self) -> std::iter::Zip<I, std::iter::Repeat<T>> {
        self.0.zip(std::iter::repeat(T::one()))
    }
}

impl<I: Iterator<Item = T>, T: Integer + NumRef + Clone> InfiniteContinuedFraction<I>
where
    for<'r> &'r T: RefNum<T>,
{
    /// This method returns a homographic function result on the fraction
    /// A homographic function is `(ax + b)/(cx + d)`
    pub fn homo(
        self,
        a: T,
        b: T,
        c: T,
        d: T,
    ) -> InfiniteContinuedFraction<HomographicResult<I, T>> {
        InfiniteContinuedFraction(HomographicResult {
            block: Block::new(a, b, c, d),
            coeffs: self.0,
        })
    }

    /// This method returns a bihomographic function result on the fraction
    /// A bihomographic function is `(axy + bx + cy + d)/(exy + fx + gy + h)`
    pub fn bihomo<U: Iterator<Item = T>>(
        self,
        rhs: InfiniteContinuedFraction<U>,
        a: T,
        b: T,
        c: T,
        d: T,
        e: T,
        f: T,
        g: T,
        h: T,
    ) -> InfiniteContinuedFraction<BihomographicResult<I, U, T>> {
        InfiniteContinuedFraction(BihomographicResult {
            block: DualBlock::new(a, b, c, d, e, f, g, h),
            x_coeffs: self.0,
            y_coeffs: rhs.0,
        })
    }

    /// Take first N coefficients in the sequence and turn it into a `ContinuedFraction` object.
    pub fn take<U>(self, count: usize) -> ContinuedFraction<U>
    where
        T: WithUnsigned<Unsigned = U> + AddAssign,
    {
        ContinuedFraction::new(self.0.take(count).collect(), Vec::new(), false)
    }

    /// Take first N coefficients in the sequence and turn it into a
    /// `ContinuedFraction` object with periodic detection.
    fn take_periodic<U>(self, count: usize) -> ContinuedFraction<U> {
        // TODO: detect re-occurence of convergents
        unimplemented!()
    }
}

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

impl<I: Iterator<Item = T>, T: Integer + NumRef + Clone> Add<T> for InfiniteContinuedFraction<I>
where
    for<'r> &'r T: RefNum<T>,
{
    type Output = InfiniteContinuedFraction<HomographicResult<I, T>>;

    fn add(self, rhs: T) -> Self::Output {
        self.homo(T::one(), rhs, T::zero(), T::one())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbols::E;

    #[test]
    fn utility_test() {
        let e = E {};

        // take()
        assert_eq!(
            (e.cfrac::<i32>() + (-2)).take(5),
            ContinuedFraction::new(vec![0i32, 1, 2, 2], vec![], false)
        );
        assert_eq!(
            (e.cfrac::<i32>() + (-3)).take(5),
            ContinuedFraction::new(vec![0i32, 3, 1, 1, 4], vec![], true)
        );
    }

    #[test]
    fn arithmetic_test() {
        let e = E {};
        let ep1_cf = e.cfrac::<i32>() + 1;
        assert_eq!(ep1_cf.0.take(5).collect::<Vec<_>>(), vec![3, 1, 2, 1, 1]);

        let sq2 = ContinuedFraction::<u32>::new(vec![1i32], vec![2], false);
        let sq2p1 = sq2.clone() + 1;
        assert_eq!(
            (sq2.expanded() + 1).0.take(5).collect::<Vec<_>>(),
            sq2p1.expanded().0.take(5).collect::<Vec<_>>()
        );
    }
}
