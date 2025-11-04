use std::{
    fmt::Display,
    marker::PhantomData,
    ops::{Add, AddAssign, Mul},
};

use num_traits::Zero;

use crate::{field::Field, modulus::Modulus};

#[derive(Debug, PartialEq, PartialOrd, Eq)]
pub struct CurvePoint<M: Modulus> {
    pub x: Field<M>,
    pub y: Field<M>,

    pub is_infinity: bool,

    order: std::marker::PhantomData<M>,
}

impl<M: Modulus> Clone for CurvePoint<M> {
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.clone(),
            is_infinity: self.is_infinity,
            order: self.order,
        }
    }

    fn clone_from(&mut self, source: &Self) {
        self.x.clone_from(&source.x);
        self.y.clone_from(&source.y);
        self.is_infinity = source.is_infinity;
    }
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
        let rhs = y * y;

        lhs == rhs
    }

    fn set_infinity(&mut self) {
        self.x = 0.into();
        self.y = 0.into();
        self.is_infinity = true;
    }
}

impl<M: Modulus> Add for CurvePoint<M> {
    type Output = CurvePoint<M>;

    fn add(self, rhs: Self) -> Self::Output {
        let (x, y) = (&self.x, &self.y);
        let (x_r, y_r) = (&rhs.x, &rhs.y);

        if rhs.is_infinity {
            return self;
        } else if self.is_infinity {
            return rhs;
        } else if x == x_r && y != y_r {
            return Self::new_infinity();
        }

        let two = Field::from(2);
        let three = Field::from(3);

        let s: Field<M> = if self == rhs {
            if y.number.is_zero() {
                return Self::new_infinity();
            }
            (three * x.pow(&two)) / (&two * y)
        } else {
            (y_r - y) / (x_r - x)
        };
        let new_x = &s.pow(&two) - x - x_r;
        let new_y = &s * &(x - &new_x) - y;

        Self::new(new_x, new_y).expect("Addition of points need to be correct")
    }
}

impl<M: Modulus> Add for &CurvePoint<M> {
    type Output = CurvePoint<M>;

    fn add(self, rhs: Self) -> Self::Output {
        let (x, y) = (&self.x, &self.y);
        let (x_r, y_r) = (&rhs.x, &rhs.y);

        if rhs.is_infinity {
            return self.clone();
        } else if self.is_infinity {
            return rhs.clone();
        } else if x == x_r && y != y_r {
            return Self::Output::new_infinity();
        }

        let two = Field::from(2);
        let three = Field::from(3);

        let s: Field<M> = if self == rhs {
            if y.number.is_zero() {
                return Self::Output::new_infinity();
            }
            (three * x.pow(&two)) / (&two * y)
        } else {
            (y_r - y) / (x_r - x)
        };
        let new_x = &s.pow(&two) - x - x_r;
        let new_y = &s * &(x - &new_x) - y;

        Self::Output::new(new_x, new_y).expect("Addition of points need to be correct")
    }
}

impl<M: Modulus> AddAssign for CurvePoint<M> {
    fn add_assign(&mut self, rhs: Self) {
        let (x, y) = (&self.x, &self.y);
        let (x_r, y_r) = (&rhs.x, &rhs.y);

        if rhs.is_infinity {
            return;
        } else if self.is_infinity {
            self.clone_from(&rhs);
            return;
        } else if x == x_r && y != y_r {
            self.set_infinity();
            return;
        }

        let two = Field::from(2);
        let three = Field::from(3);

        let s: Field<M> = if self == &rhs {
            if y.number.is_zero() {
                self.set_infinity();
                return;
            }
            (three * x.pow(&two)) / (&two * y)
        } else {
            (y_r - y) / (x_r - x)
        };
        let new_x = &s.pow(&two) - x - x_r;
        let new_y = &s * &(x - &new_x) - y;

        self.x = new_x;
        self.y = new_y;
    }
}

impl<M: Modulus> Mul<Field<M>> for CurvePoint<M> {
    type Output = CurvePoint<M>;

    fn mul(self, rhs: Field<M>) -> Self::Output {
        let mut result = Self::new_infinity();
        let mut add = self.clone();

        let bits = rhs.number.bits();

        for i in 0..bits {
            if rhs.number.bit(i) {
                result = &result + &add;
            }
            add = &add + &add;
        }
        result
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

        assert_eq!(a + b, result)
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

        assert_eq!(a + b, result)
    }

    #[test]
    fn doubling_test_fp() {
        let a: CurvePoint<OrderP> = CurvePoint::new(
            Field::from_str(
                "16909226366434337090002029874530773609945429384634010589202308143299843842109",
            )
            .unwrap(),
            Field::from_str(
                "22289652916292038233802578916041176193424404491330669531128447027442819547009",
            )
            .unwrap(),
        )
        .unwrap();

        let result: CurvePoint<OrderP> = CurvePoint::new(
            Field::from_str(
                "26627901895311929412061073520912504423663864690969978810618273611899328804957",
            )
            .unwrap(),
            Field::from_str(
                "5645391259068071962401054088667169471848861916041555603144023837132996877952",
            )
            .unwrap(),
        )
        .unwrap();

        assert_eq!(a.clone() + a, result)
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

    #[test]
    fn multiply() {
        let a: CurvePoint<OrderP> = CurvePoint::new(
            Field::from_str(
                "19422280595109069951862896517863616574418686534058761073670180784027971855802",
            )
            .unwrap(),
            Field::from_str(
                "11090450506703525262630446804833734397652482233395951028154717035371611573726",
            )
            .unwrap(),
        )
        .unwrap();

        let scalar: Field<OrderP> = Field::from_str(
            "1563750830861819851520880065691358712351730162896247138866720840651323457775",
        )
        .unwrap();

        dbg!(&scalar);

        let result: CurvePoint<OrderP> = CurvePoint::new(
            Field::from_str(
                "24900093740427460719317843302751314869260252172138935506494925835061423875240",
            )
            .unwrap(),
            Field::from_str(
                "10690880907742924436758301380673321005651802536989903562762736696002004031756",
            )
            .unwrap(),
        )
        .unwrap();

        assert_eq!(result, a * scalar);
    }
}
