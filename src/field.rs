use std::{
    fmt::Display,
    marker::PhantomData,
    ops::{Add, Mul, Sub},
};

use num_bigint::BigUint;

use crate::modulus::{Modulus, OrderP, OrderQ};

#[derive(Debug, PartialEq, PartialOrd, Eq)]
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

    pub fn pow(self, power: Self) -> Self {
        assert!(self.modulo == power.modulo);

        Field {
            number: self.number.modpow(&power.number, M::get()),
            modulo: self.modulo,
        }
    }
}

impl<M: Modulus> Add for Field<M> {
    type Output = Field<M>;

    fn add(self, rhs: Self) -> Self::Output {
        Field::new(self.number + rhs.number)
    }
}

impl<M: Modulus> Mul for Field<M> {
    type Output = Field<M>;

    fn mul(self, rhs: Self) -> Self::Output {
        Field::new(self.number * rhs.number)
    }
}

impl<M: Modulus> Sub for Field<M> {
    type Output = Field<M>;

    fn sub(self, rhs: Self) -> Self::Output {
        Field::new(self.number + M::get() - rhs.number)
    }
}

impl<M: Modulus> Display for Field<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(format!("({},{})", self.number, M::get()).as_str())?;
        Ok(())
    }
}

pub type FieldP = Field<OrderP>;
pub type FieldQ = Field<OrderQ>;
