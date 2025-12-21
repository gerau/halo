use std::{
    fmt::Display,
    marker::PhantomData,
    ops::{Add, Div, Mul, Sub},
    str::FromStr,
};

use num_bigint::{BigUint, ParseBigIntError, RandBigInt};
use num_traits::ConstZero;
use rand::Rng;

use crate::modulus::{Modulus, OrderP, OrderQ};

#[derive(Debug, PartialEq, PartialOrd, Eq, Ord)]
pub struct Field<M: Modulus> {
    pub number: BigUint,

    modulo: std::marker::PhantomData<M>,
}

impl<M: Modulus> Field<M> {
    pub fn new(number: BigUint) -> Self {
        Self {
            number: number % M::get(),
            modulo: PhantomData,
        }
    }

    pub fn pow_field(&self, power: &Self) -> Self {
        Self {
            number: self.number.modpow(&power.number, M::get()),
            modulo: PhantomData,
        }
    }

    pub fn pow(&self, power: u64) -> Self {
        Self::new(self.number.modpow(&BigUint::from(power), M::get()))
    }

    pub fn get_random<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::new(rng.gen_biguint_below(M::get()))
    }

    pub fn inv(&self) -> Option<Self> {
        self.number.modinv(M::get()).map(Self::new)
    }

    pub fn one() -> Self {
        Self::from(1)
    }

    pub const ZERO: Self = Self {
        number: BigUint::ZERO,
        modulo: PhantomData,
    };
}

impl<M: Modulus> Add for Field<M> {
    type Output = Field<M>;

    fn add(self, rhs: Self) -> Self::Output {
        Field::new(self.number + rhs.number)
    }
}

impl<M: Modulus> Add for &Field<M> {
    type Output = Field<M>;

    fn add(self, rhs: Self) -> Self::Output {
        Field::new(&self.number + &rhs.number)
    }
}

impl<M: Modulus> Add<&Field<M>> for Field<M> {
    type Output = Field<M>;

    fn add(self, rhs: &Self) -> Self::Output {
        Field::new(self.number + &rhs.number)
    }
}

impl<M: Modulus> Add<Field<M>> for &Field<M> {
    type Output = Field<M>;

    fn add(self, rhs: Field<M>) -> Self::Output {
        Field::new(&self.number + &rhs.number)
    }
}

impl<M: Modulus> Mul for Field<M> {
    type Output = Field<M>;

    fn mul(self, rhs: Self) -> Self::Output {
        Field::new(self.number * rhs.number)
    }
}

impl<M: Modulus> Mul for &Field<M> {
    type Output = Field<M>;

    fn mul(self, rhs: Self) -> Self::Output {
        Field::new(&self.number * &rhs.number)
    }
}

impl<M: Modulus> Mul<&Field<M>> for Field<M> {
    type Output = Field<M>;

    fn mul(self, rhs: &Field<M>) -> Self::Output {
        Field::new(self.number * &rhs.number)
    }
}

impl<M: Modulus> Mul<Field<M>> for &Field<M> {
    type Output = Field<M>;

    fn mul(self, rhs: Field<M>) -> Self::Output {
        Field::new(&self.number * &rhs.number)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<M: Modulus> Div for Field<M> {
    type Output = Field<M>;

    fn div(self, rhs: Self) -> Self::Output {
        Field::new(self.number * rhs.number.modinv(M::get()).unwrap())
    }
}

impl<M: Modulus> Sub for Field<M> {
    type Output = Field<M>;

    fn sub(self, rhs: Self) -> Self::Output {
        Field::new(self.number + M::get() - rhs.number)
    }
}

impl<M: Modulus> Sub for &Field<M> {
    type Output = Field<M>;

    fn sub(self, rhs: Self) -> Self::Output {
        Field::new(&self.number + M::get() - &rhs.number)
    }
}

impl<M: Modulus> Sub<&Field<M>> for Field<M> {
    type Output = Field<M>;

    fn sub(self, rhs: &Field<M>) -> Self::Output {
        Field::new(self.number + M::get() - &rhs.number)
    }
}

impl<M: Modulus> Sub<Field<M>> for &Field<M> {
    type Output = Field<M>;

    fn sub(self, rhs: Field<M>) -> Self::Output {
        Field::new(&self.number + M::get() - rhs.number)
    }
}

impl<M: Modulus> Clone for Field<M> {
    fn clone(&self) -> Self {
        Self {
            number: self.number.clone(),
            modulo: PhantomData,
        }
    }
}

impl<M: Modulus> Display for Field<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(format!("({},{})", self.number, M::get()).as_str())?;
        Ok(())
    }
}

impl<M: Modulus> From<u64> for Field<M> {
    fn from(value: u64) -> Self {
        Field::new(BigUint::from(value))
    }
}

impl<M: Modulus> From<i32> for Field<M> {
    fn from(value: i32) -> Self {
        let unsigned: u32 = value.try_into().unwrap();
        Field::new(BigUint::from(unsigned))
    }
}

impl<M: Modulus> FromStr for Field<M> {
    type Err = ParseBigIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Field::new(BigUint::from_str(s)?))
    }
}
#[allow(dead_code)]
pub type FieldP = Field<OrderP>;

#[allow(dead_code)]
pub type FieldQ = Field<OrderQ>;
