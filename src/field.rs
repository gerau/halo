use std::ops::{Add, Mul, Sub};

use num_bigint::BigUint;

#[derive(Debug, PartialEq, PartialOrd, Eq)]
struct Field {
    number: BigUint,

    modulo: BigUint,
}

impl Add for Field {
    type Output = Field;

    fn add(self, rhs: Self) -> Self::Output {
        assert!(self.modulo == rhs.modulo);

        let modulo = &self.modulo;

        Field {
            number: (self.number + rhs.number + modulo) % modulo,
            modulo: rhs.modulo,
        }
    }
}

impl Mul for Field {
    type Output = Field;

    fn mul(self, rhs: Self) -> Self::Output {
        assert!(self.modulo == rhs.modulo);

        Field {
            number: (self.number * rhs.number) % self.modulo,
            modulo: rhs.modulo,
        }
    }
}

impl Sub for Field {
    type Output = Field;

    fn sub(self, rhs: Self) -> Self::Output {
        assert!(self.modulo == rhs.modulo);

        Field {
            number: (self.number - rhs.number) % self.modulo,
            modulo: rhs.modulo,
        }
    }
}

impl Field {
    fn one_from(number: &Self) -> Self {
        Field {
            number: BigUint::from(1u32),
            modulo: number.modulo.clone(),
        }
    }

    fn pow(self, power: Self) -> Self {
        assert!(self.modulo == power.modulo);

        Field {
            number: self.number.modpow(&power.number, &power.modulo),
            modulo: self.modulo,
        }
    }
}
