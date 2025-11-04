use num_bigint::BigUint;

use std::{
    fmt::Display,
    marker::PhantomData,
    ops::{Add, Mul, Sub},
};

use num_traits::Zero;

use crate::{field::Field, modulus::Modulus};

#[derive(Debug, PartialEq, PartialOrd, Eq, Clone)]
pub struct CurvePoint<M: Modulus> {
    pub x: Field<M>,
    pub y: Field<M>,

    pub is_infinity: bool,

    order: std::marker::PhantomData<M>,
}

impl<M: Modulus> CurvePoint<M> {
    pub fn new(x: Field<M>, y: Field<M>) -> Result<Self, &'static str> {
        let res = Self {
            x,
            y,
            is_infinity: false,
            order: PhantomData,
        };

        if !res.is_on_curve() {
            Err("Point with given coordinates is not on curve")
        } else {
            Ok(res)
        }
    }

    pub fn new_infinity() -> Self {
        Self {
            x: 0.into(),
            y: 0.into(),
            is_infinity: true,
            order: PhantomData,
        }
    }

    fn is_on_curve(&self) -> bool {
        if self.is_infinity {
            return true;
        }

        let x = &self.x;
        let y = &self.y;

        let lhs = x.pow(&3.into()) + 5.into();
        println!("{}", lhs);
        let rhs = y * y;
        println!("{}", rhs);

        lhs == rhs
    }
}

impl<M: Modulus> Add for CurvePoint<M> {
    type Output = CurvePoint<M>;

    fn add(self, rhs: Self) -> Self::Output {
        if rhs.is_infinity {
            return self;
        } else if self.is_infinity {
            return rhs;
        } else if self.x == rhs.x && self.y != rhs.y {
            return Self::new_infinity();
        }

        let two = Field::from(2);
        let three = Field::from(3);

        let s: Field<M> = if self == rhs {
            if self.y.number.is_zero() {
                return Self::new_infinity();
            }
            (&three * &self.x.pow(&two)) / (&two * &self.y)
        } else {
            (&rhs.y - &self.y) / (&rhs.x - &self.x)
        };
        let new_x = &s.pow(&two) - &self.x - &rhs.x;
        let new_y = &s * &(&self.x - &new_x) - &self.y;

        Self::new(new_x, new_y).expect("Addition of points need to be correct")
    }
}

impl<M: Modulus> Display for CurvePoint<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(
            format!(
                "{}\n  x: {}\n  y: {}",
                M::get(),
                self.x.number,
                self.y.number,
            )
            .as_str(),
        )?;
        Ok(())
    }
}

#[cfg(test)]
mod elliptic_curve_tests {
    use std::str::FromStr;

    use crate::modulus::{OrderP, OrderQ};

    use super::*;

    #[test]
    fn addition_test_fp() {
        let a: CurvePoint<OrderP> = CurvePoint::new(
            Field::from_str(
                "2004294055917815597860384095144215206559862545661360726227017460961152747606",
            )
            .unwrap(),
            Field::from_str(
                "28764628678555263360045525070657468757032969501382002203382765555603343529783",
            )
            .unwrap(),
        )
        .unwrap();

        let b: CurvePoint<OrderP> = CurvePoint::new(
            Field::from_str(
                "5471070483644318146131152889433134040510736149701019904800859956377477447406",
            )
            .unwrap(),
            Field::from_str(
                "21706088111463393086001983292662469052283083229778616436834991783518812362858",
            )
            .unwrap(),
        )
        .unwrap();

        let result: CurvePoint<OrderP> = CurvePoint::new(
            Field::from_str(
                "21224905437769693582834397320858778426817179592960913233648874753823773295560",
            )
            .unwrap(),
            Field::from_str(
                "18687433395868625752074020959324078709705853201683474555229546834597093185446",
            )
            .unwrap(),
        )
        .unwrap();

        assert!((a + b) == result)
    }

    #[test]
    fn addition_test_fq() {
        let a: CurvePoint<OrderQ> = CurvePoint::new(
            Field::from_str(
                "15429517692791640599507297261312258517549351511519415466581975689467530435322",
            )
            .unwrap(),
            Field::from_str(
                "8006861636318967171214779219383912529732064645950614224774477786131620663416",
            )
            .unwrap(),
        )
        .unwrap();

        let b: CurvePoint<OrderQ> = CurvePoint::new(
            Field::from_str(
                "4326745279610707763609991096116764230103164341870094283813887583802699972215",
            )
            .unwrap(),
            Field::from_str(
                "6306532159821668765421901642344615381945880463949287658094002319085602426794",
            )
            .unwrap(),
        )
        .unwrap();

        let result: CurvePoint<OrderQ> = CurvePoint::new(
            Field::from_str(
                "26842240685334946436775289762343089862353789842844810776109479930121492577721",
            )
            .unwrap(),
            Field::from_str(
                "25281454545305572737888299705482991177326271983610970161027704840573815321962",
            )
            .unwrap(),
        )
        .unwrap();

        assert!((a + b) == result)
    }

    #[test]
    fn addition_test_with_infinity() {
        let a: CurvePoint<OrderP> = CurvePoint::new(
            Field::from_str(
                "2004294055917815597860384095144215206559862545661360726227017460961152747606",
            )
            .unwrap(),
            Field::from_str(
                "28764628678555263360045525070657468757032969501382002203382765555603343529783",
            )
            .unwrap(),
        )
        .unwrap();

        let infinity: CurvePoint<OrderP> = CurvePoint::new_infinity();

        assert_eq!(a.clone(), a + infinity);
    }
}
