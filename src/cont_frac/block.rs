use num_integer::Integer;
use num_traits::{CheckedAdd, CheckedMul, NumRef, One, RefNum, Zero};
use std::mem::replace;

/// A block on the magic table for homographic operation computation of continued fractions
/// The method is described in <https://crypto.stanford.edu/pbc/notes/contfrac/compute.html>
#[derive(Debug, Clone, Copy)]
pub struct Block<T> {
    pm1: T, // p_(k-1)
    pm2: T, // p_(k-2)
    qm1: T, // q_(k-1)
    qm2: T, // q_(k-2)
}

impl<T> Block<T> {
    /// create a block that represents (ax + b) / (cx + d)
    pub fn new(a: T, b: T, c: T, d: T) -> Self {
        Block {
            pm1: a,
            pm2: b,
            qm1: c,
            qm2: d,
        }
    }

    /// push the latest convergent to the block
    pub fn update(&mut self, p: T, q: T) {
        self.pm2 = replace(&mut self.pm1, p);
        self.qm2 = replace(&mut self.qm1, q);
    }
}

impl<T: Zero + One> Block<T> {
    /// create a block that represents a identity operation
    pub fn identity() -> Self {
        Block {
            pm1: T::one(),
            pm2: T::zero(),
            qm1: T::zero(),
            qm2: T::one(),
        }
    }
}

impl<T: Integer + NumRef> Block<T>
where
    for<'r> &'r T: RefNum<T>,
{
    /// move with an coefficient from regular continued fraction
    pub fn rmove(&mut self, a: T) {
        let p = &a * &self.pm1 + &self.pm2;
        let q = a * &self.qm1 + &self.qm2;
        self.update(p, q);
    }

    /// move with two coefficients from generalized continued fraction
    pub fn gmove(&mut self, a: T, b: T) {
        let p = &self.pm1 * &b + &self.pm2 * &a;
        let q = &self.qm1 * b + &self.qm2 * a;
        let g = p.gcd(&q).gcd(&self.pm1).gcd(&self.qm1);

        if g > T::one() {
            self.pm1 = &self.pm1 / &g;
            self.qm1 = &self.qm1 / &g;
            self.update(p / &g, q / g);
        } else {
            self.update(p, q);
        };
    }

    /// Check whether we should reduce integer from the convergents
    /// If the convergets has the same integer part i, then return Ok with the integer and remainders
    #[inline]
    pub fn check_integer(&self) -> Result<(T, T, T), ()> {
        if self.qm1.is_zero() || self.qm2.is_zero() {
            return Err(());
        }

        let (im1, rm1) = self.pm1.div_rem(&self.qm1);
        let (im2, rm2) = self.pm2.div_rem(&self.qm2);

        if im1 == im2 {
            Ok((im1, rm1, rm2))
        } else {
            Err(())
        }
    }

    /// extract the integer part if latests two convergents agrees,
    /// and then flip the remaining convergent, this function also
    /// update the convergents if failed
    pub fn reduce_recip(&mut self) -> Option<T> {
        match self.check_integer() {
            Ok((i, rm1, rm2)) => {
                self.pm1 = replace(&mut self.qm1, rm1);
                self.pm2 = replace(&mut self.qm2, rm2);
                Some(i)
            }
            Err(_) => None,
        }
    }

    /// extract the integer part if latest two convergents agrees,
    /// and the mulitply the numerator by some constant. This is used
    /// for decimal digits extraction. This operation is done inplace.
    pub fn reduce_mul(&mut self, base: T) -> Option<T> {
        match self.check_integer() {
            Ok((i, rm1, rm2)) => {
                self.pm1 = rm1 * &base;
                self.pm2 = rm2 * base;
                Some(i)
            }
            Err(_) => None,
        }
    }
}

impl<T: Integer + CheckedAdd + CheckedMul> Block<T> {
    /// Note that update() should be called after checked move
    pub fn checked_rmove(&mut self, a: T) -> Option<(T, T)> {
        let p = a
            .checked_mul(&self.pm1)
            .and_then(|v| v.checked_add(&self.pm2))?;
        let q = a
            .checked_mul(&self.qm1)
            .and_then(|v| v.checked_add(&self.qm2))?;
        Some((p, q))
    }
}

impl<T: Integer + CheckedAdd + CheckedMul + NumRef> Block<T>
where
    for<'r> &'r T: RefNum<T>,
{
    /// Note that update() should be called after checked move
    pub fn checked_gmove(&mut self, a: T, b: T) -> Option<(T, T)> {
        let bpm1 = b.checked_mul(&self.pm1)?;
        let apm2 = a.checked_mul(&self.pm2)?;
        let p = bpm1.checked_add(&apm2)?;
        let bqm1 = b.checked_mul(&self.qm1)?;
        let aqm2 = a.checked_mul(&self.qm2)?;
        let q = bqm1.checked_add(&aqm2)?;

        let g = p.gcd(&q).gcd(&self.pm1).gcd(&self.qm1);
        if g > T::one() {
            self.pm1 = &self.pm1 / &g;
            self.qm1 = &self.qm1 / &g;
            Some((p / &g, q / g))
        } else {
            Some((p, q))
        }
    }
}

/// A block on the magic table for bihomographic operation computation of continued fractions
/// The method is described in <https://crypto.stanford.edu/pbc/notes/contfrac/bihom.html>
#[derive(Debug, Clone, Copy)]
pub struct DualBlock<T> {
    pm11: T, // p with a_(i-1), b_(j-1)
    pm12: T, // p with a_(i-1), b_(j-2)
    pm21: T, // ..
    pm22: T, // ..
    qm11: T, // q with a_(i-1), b_(j-1)
    qm12: T, // q with a_(i-1), b_(j-2)
    qm21: T, // ..
    qm22: T, // ..
}

impl<T> DualBlock<T> {
    /// create a block that represents (axy + bx + cy + d)/(exy + fx + gy + h)
    pub fn new(a: T, b: T, c: T, d: T, e: T, f: T, g: T, h: T) -> Self {
        DualBlock {
            pm11: a,
            pm12: b,
            pm21: c,
            pm22: d,
            qm11: e,
            qm12: f,
            qm21: g,
            qm22: h,
        }
    }

    /// push the latest convergent using x from right to the block
    pub fn update_right(&mut self, p1: T, q1: T, p2: T, q2: T) {
        self.pm21 = replace(&mut self.pm11, p1);
        self.qm21 = replace(&mut self.qm11, q1);

        self.pm22 = replace(&mut self.pm12, p2);
        self.qm22 = replace(&mut self.qm22, q2);
    }

    /// push the latest convergent using y from bottom to the block
    pub fn update_down(&mut self, p1: T, q1: T, p2: T, q2: T) {
        self.pm12 = replace(&mut self.pm11, p1);
        self.qm12 = replace(&mut self.qm11, q1);

        self.pm22 = replace(&mut self.pm21, p2);
        self.qm22 = replace(&mut self.qm21, q2);
    }
}

impl<T: Integer + NumRef> DualBlock<T>
where
    for<'r> &'r T: RefNum<T>,
{
    /// move with an coefficient from the first regular continued fraction (x)
    pub fn rmove_right(&mut self, a: T) {
        let p1 = &a * &self.pm11 + &self.pm21;
        let q1 = &a * &self.qm11 + &self.qm21;
        let p2 = &a * &self.pm12 + &self.pm22;
        let q2 = a * &self.qm12 + &self.qm22;
        self.update_right(p1, q1, p2, q2);
    }

    /// move with an coefficient from the second regular continued fraction (y)
    pub fn rmove_down(&mut self, a: T) {
        let p1 = &a * &self.pm11 + &self.pm12;
        let q1 = &a * &self.qm11 + &self.qm12;
        let p2 = &a * &self.pm21 + &self.pm22;
        let q2 = a * &self.qm21 + &self.qm22;
        self.update_down(p1, q1, p2, q2);
    }

    /// Check whether we should reduce integer from the convergents
    /// Unlike the checking for single block, this is performed inplace (without new convergents)
    /// If the convergets has the same integer part i, then return Ok with the integer and remainders
    /// If not, then return Err with flags about whether to move right or move down
    #[inline]
    pub fn check_integer(&self) -> Result<(T, T, T, T, T), (bool, bool)> {
        if self.qm22.is_zero() {
            return Err((self.qm21.is_zero(), self.qm12.is_zero()));
        }
        if self.qm11.is_zero() {
            return Err((self.qm12.is_zero(), self.qm21.is_zero()));
        }

        let (i11, r11) = self.pm11.div_rem(&self.qm11);
        let (i12, r12) = self.pm12.div_rem(&self.qm12);
        let (i21, r21) = self.pm21.div_rem(&self.qm21);
        let (i22, r22) = self.pm22.div_rem(&self.qm22);

        if i11 == i12 && i11 == i21 && i11 == i22 {
            Ok((i11, r11, r12, r21, r22))
        } else {
            Err((i11 != i12 || i21 != i22, i11 != i21 || i12 != i22))
        }
    }

    /// extract the integer part if all convergents agrees, and the flip
    /// the convergent. this function also move the block if failed
    pub fn reduce_recip(&mut self) -> Result<T, (bool, bool)> {
        match self.check_integer() {
            Ok((i, r11, r12, r21, r22)) => {
                self.pm11 = replace(&mut self.qm11, r11);
                self.pm12 = replace(&mut self.qm12, r12);
                self.pm21 = replace(&mut self.qm21, r21);
                self.pm22 = replace(&mut self.qm22, r22);
                Ok(i)
            }
            Err(f) => Err(f),
        }
    }
}
