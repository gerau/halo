use lazy_static::lazy_static;
use num_bigint::BigUint;
use std::str::FromStr;

pub trait Modulus: PartialEq + Clone {
    fn get() -> &'static BigUint;
}

lazy_static! {
    static ref P: BigUint = BigUint::from_str(
        "28948022309329048855892746252171976963363056481941647379679742748393362948097"
    )
    .unwrap();
    static ref Q: BigUint = BigUint::from_str(
        "28948022309329048855892746252171976963363056481941560715954676764349967630337"
    )
    .unwrap();
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct OrderP;

impl Modulus for OrderP {
    fn get() -> &'static BigUint {
        &P
    }
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct OrderQ;

impl Modulus for OrderQ {
    fn get() -> &'static BigUint {
        &Q
    }
}
