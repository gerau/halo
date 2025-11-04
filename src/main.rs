use crate::field::FieldP;

mod elliptic_curve;
mod field;
mod modulus;

fn main() {
    let a = FieldP::new(123123123123u64.into());
    let b = FieldP::new(1231123u64.into());

    println!("{}", a.pow(b));
}
