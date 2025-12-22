use crate::{
    field::{self, Field},
    modulus::Modulus,
};

use std::{
    fmt,
    ops::{Mul, Sub},
};
use std::{fmt::Display, ops::Add};

#[derive(Debug, Clone)]
pub struct Polynomial<M: Modulus> {
    pub coefficients: Vec<Field<M>>,

    // This indicates offset from what coefficient is starting.
    pub offset: i64,
}

impl<M: Modulus> Polynomial<M> {
    pub fn evaluate_at(&self, x: Field<M>) -> Field<M> {
        if self.coefficients.is_empty() {
            return Field::<M>::ZERO;
        }

        if x.is_zero() {
            // we don't have constant term, so we return zero
            if self.offset > 0 {
                return Field::<M>::ZERO;
            }

            return self
                .coefficients
                .get((-self.offset).max(0) as usize)
                .unwrap_or(&Field::<M>::ZERO)
                .clone();
        }

        let x_inv = if self.offset < 0 {
            x.inv().expect("X must be invertible")
        } else {
            Field::<M>::ZERO
        };

        self.coefficients
            .iter()
            .enumerate()
            .map(|(i, coeff)| {
                let exponent = i as i64 + self.offset;

                if exponent == 0 {
                    // clone is fine, because it only do this one time
                    coeff.clone()
                } else if exponent > 0 {
                    coeff * &x.pow(exponent as u64)
                } else {
                    coeff * &x_inv.pow(exponent.unsigned_abs())
                }
            })
            .fold(Field::<M>::ZERO, |acc, val| acc + val)
    }

    pub fn shift(&mut self, shift: i64) {
        self.offset += shift;
    }

    pub fn shift_poly(&self, shift: i64) -> Self {
        Self {
            coefficients: self.coefficients.clone(),
            offset: self.offset + shift,
        }
    }

    pub fn empty() -> Self {
        Polynomial {
            coefficients: vec![],
            offset: 0,
        }
    }
}

impl<M: Modulus> Add for &Polynomial<M> {
    type Output = Polynomial<M>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.coefficients.is_empty() {
            return rhs.clone();
        }
        if rhs.coefficients.is_empty() {
            return self.clone();
        }

        let min_offset = self.offset.min(rhs.offset);
        let max_exp_self = self.offset + self.coefficients.len() as i64;
        let max_exp_rhs = rhs.offset + rhs.coefficients.len() as i64;
        let max_exp = max_exp_self.max(max_exp_rhs);

        let new_len = (max_exp - min_offset) as usize;
        let mut new_coeffs = vec![Field::<M>::ZERO; new_len];

        for (i, coeff) in self.coefficients.iter().enumerate() {
            let target_idx = (i as i64 + self.offset - min_offset) as usize;
            new_coeffs[target_idx] = &new_coeffs[target_idx] + coeff;
        }

        for (i, coeff) in rhs.coefficients.iter().enumerate() {
            let target_idx = (i as i64 + rhs.offset - min_offset) as usize;
            new_coeffs[target_idx] = &new_coeffs[target_idx] + coeff;
        }

        Polynomial {
            coefficients: new_coeffs,
            offset: min_offset,
        }
    }
}

impl<M: Modulus> Sub for &Polynomial<M> {
    type Output = Polynomial<M>;

    fn sub(self, rhs: Self) -> Self::Output {
        if self.coefficients.is_empty() {
            return rhs.clone();
        }
        if rhs.coefficients.is_empty() {
            return self.clone();
        }

        let min_offset = self.offset.min(rhs.offset);
        let max_exp_self = self.offset + self.coefficients.len() as i64;
        let max_exp_rhs = rhs.offset + rhs.coefficients.len() as i64;
        let max_exp = max_exp_self.max(max_exp_rhs);

        let new_len = (max_exp - min_offset) as usize;
        let mut new_coeffs = vec![Field::<M>::ZERO; new_len];

        for (i, coeff) in self.coefficients.iter().enumerate() {
            let target_idx = (i as i64 + self.offset - min_offset) as usize;
            new_coeffs[target_idx] = &new_coeffs[target_idx] + coeff;
        }

        for (i, coeff) in rhs.coefficients.iter().enumerate() {
            let target_idx = (i as i64 + rhs.offset - min_offset) as usize;
            new_coeffs[target_idx] = &new_coeffs[target_idx] - coeff;
        }

        Polynomial {
            coefficients: new_coeffs,
            offset: min_offset,
        }
    }
}

impl<M: Modulus> Mul<&Field<M>> for &Polynomial<M> {
    type Output = Polynomial<M>;

    fn mul(self, rhs: &Field<M>) -> Self::Output {
        let new_coeffs = self.coefficients.iter().map(|coeff| coeff * rhs).collect();
        Polynomial {
            coefficients: new_coeffs,
            offset: self.offset,
        }
    }
}

impl<M: Modulus> Mul for &Polynomial<M> {
    type Output = Polynomial<M>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.coefficients.is_empty() || rhs.coefficients.is_empty() {
            return Polynomial {
                coefficients: vec![],
                offset: 0,
            };
        }

        let new_offset = self.offset + rhs.offset;
        let new_len = self.coefficients.len() + rhs.coefficients.len() - 1;
        let mut new_coeffs = vec![Field::<M>::ZERO; new_len];

        for (i, a_i) in self.coefficients.iter().enumerate() {
            for (j, b_j) in rhs.coefficients.iter().enumerate() {
                new_coeffs[i + j] = &new_coeffs[i + j] + (a_i * b_j);
            }
        }

        Polynomial {
            coefficients: new_coeffs,
            offset: new_offset,
        }
    }
}

impl<M: Modulus> Display for Polynomial<M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.coefficients.iter().all(|c| c.is_zero()) {
            return write!(f, "0");
        }

        for (i, coeff) in self.coefficients.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }

            let exponent = i as i64 + self.offset;

            if i != 0 {
                write!(f, " + ")?;
            }

            match exponent {
                0 => write!(f, "{}", coeff)?,
                1 => write!(f, "{} * x", coeff)?,
                _ => write!(f, "{} * x^({})", coeff, exponent)?,
            }
        }

        Ok(())
    }
}

impl<M: Modulus> Polynomial<M> {
    pub fn get_coeff(&self, power: i64) -> Field<M> {
        let index = power - self.offset;
        if index >= 0 && (index as usize) < self.coefficients.len() {
            self.coefficients[index as usize].clone()
        } else {
            Field::<M>::ZERO
        }
    }

    pub fn extract_t_lo(&self, d: usize) -> Polynomial<M> {
        let mut lo_coeffs = Vec::with_capacity(d);

        let d = d as i64;
        for i in 0..d {
            let deg = i - d;
            lo_coeffs.push(self.get_coeff(deg));
        }

        Polynomial {
            coefficients: lo_coeffs,
            offset: -d,
        }
    }

    pub fn extract_t_hi(&self, d: usize) -> Polynomial<M> {
        let mut hi_coeffs = Vec::with_capacity(d);

        for i in 0..d {
            let deg = (i as i64) + 1;
            hi_coeffs.push(self.get_coeff(deg));
        }

        Polynomial {
            coefficients: hi_coeffs,
            offset: 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::field::FieldP;

    use super::*;

    #[test]
    fn test_polynomial_addition() {
        let poly1 = Polynomial {
            coefficients: Vec::from([FieldP::from(1), FieldP::from(2)]),
            offset: -1,
        };

        let poly2 = Polynomial {
            coefficients: Vec::from([FieldP::from(3), FieldP::from(4)]),
            offset: 0,
        };

        let result = &poly1 + &poly2;

        assert_eq!(result.offset, -1);
        assert_eq!(result.coefficients.len(), 3);
        assert_eq!(result.coefficients[0], FieldP::from(1));
        assert_eq!(result.coefficients[1], FieldP::from(5));
        assert_eq!(result.coefficients[2], FieldP::from(4));
    }

    #[test]
    fn test_polynomial_multiplication() {
        let poly1 = Polynomial {
            coefficients: Vec::from([FieldP::from(1), FieldP::from(2)]),
            offset: 1,
        };

        let poly2 = Polynomial {
            coefficients: Vec::from([FieldP::from(3), FieldP::from(4)]),
            offset: -1,
        };

        let result = &poly1 * &poly2;

        assert_eq!(result.offset, 0);
        assert_eq!(result.coefficients.len(), 3);
        assert_eq!(result.coefficients[0], FieldP::from(3));
        assert_eq!(result.coefficients[1], FieldP::from(10));
        assert_eq!(result.coefficients[2], FieldP::from(8));
    }
    #[test]
    fn test_math_and_evaluation() {
        let x = FieldP::from(2);

        let p1 = Polynomial {
            coefficients: Vec::from([FieldP::from(1), FieldP::from(1)]),
            offset: 0,
        };

        let p2 = Polynomial {
            coefficients: Vec::from([FieldP::from(1)]),
            offset: 1,
        };

        let p3 = &p1 * &p2;
        assert_eq!(p3.evaluate_at(x), FieldP::from(6));
    }
}
