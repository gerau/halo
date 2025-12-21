use itertools::izip;

use crate::{
    field::{self, Field},
    modulus::Modulus,
};

// TODO: refactor, properly check bound

pub fn evaluate_r<M: Modulus>(
    x: &Field<M>,
    y: &Field<M>,
    n: usize,
    a: &[Field<M>],
    b: &[Field<M>],
    c: &[Field<M>],
) -> Field<M> {
    assert_eq!(n, a.len());
    assert_eq!(n, b.len());
    assert_eq!(n, c.len());

    let xy = x * y;
    let z = xy.inv().expect("XY must be invertible");
    let z_n = z.pow(n as u64);

    let sum_a: Field<M> = std::iter::successors(Some(xy.clone()), |prev| Some(&xy * prev))
        .take(n)
        .zip(a.iter())
        .map(|(power, coeff_a)| coeff_a * power)
        .fold(Field::<M>::ZERO, |acc, val| acc + val);

    let sum_bc: Field<M> = std::iter::successors(Some(z.clone()), |prev| Some(&z * prev))
        .take(n)
        .zip(b.iter().zip(c.iter()))
        .map(|(pwr_z, (coeff_b, coeff_c))| pwr_z * (coeff_b + coeff_c * &z_n))
        .fold(Field::<M>::ZERO, |acc, val| acc + val);

    sum_a + sum_bc
}

pub fn evaluate_s<M: Modulus>(
    x: &Field<M>,
    y: &Field<M>,
    n: usize,
    q: usize,
    u: &[impl AsRef<[Field<M>]>],
    v: &[impl AsRef<[Field<M>]>],
    w: &[impl AsRef<[Field<M>]>],
) -> Field<M> {
    let y_powers: Vec<Field<M>> = std::iter::successors(Some(y.clone()), |prev| Some(y * prev))
        .take(q)
        .collect();

    let (u_i, v_i, w_i): (Vec<_>, Vec<_>, Vec<_>) = izip!(u, v, w)
        .map(|(row_u, row_v, row_w)| {
            let dot_product = |row: &dyn AsRef<[Field<M>]>| {
                row.as_ref()
                    .iter()
                    .zip(&y_powers)
                    .map(|(coeff, pwr)| coeff * pwr)
                    .fold(Field::<M>::ZERO, |acc, val| acc + val)
            };

            (dot_product(row_u), dot_product(row_v), dot_product(row_w))
        })
        .collect();
    let z = x.inv().expect("X must be invertible");

    let x_n = x.pow(n as u64);

    let sum_a: Field<M> = std::iter::successors(Some(z.clone()), |prev| Some(&z * prev))
        .take(n)
        .zip(u_i.iter())
        .map(|(power, coeff_u)| coeff_u * power)
        .fold(Field::<M>::ZERO, |acc, val| acc + val);

    let sum_bc: Field<M> = std::iter::successors(Some(x.clone()), |prev| Some(x * prev))
        .take(n)
        .zip(v_i.iter().zip(w_i.iter()))
        .map(|(pwr_z, (coeff_v, coeff_w))| pwr_z * (coeff_v + coeff_w * &x_n))
        .fold(Field::<M>::ZERO, |acc, val| acc + val);

    sum_a + sum_bc
}

pub fn evaluate_ss<M: Modulus + std::fmt::Debug>(
    x: &Field<M>,
    y: &Field<M>,
    n: usize,
    q: usize,
    u: &[impl AsRef<[Field<M>]>],
    v: &[impl AsRef<[Field<M>]>],
    w: &[impl AsRef<[Field<M>]>],
) -> Field<M> {
    let z = y.inv().expect("Y must be invertible");
    let y_s: Vec<Field<M>> = std::iter::successors(Some(y.clone()), |prev| Some(y * prev))
        .take(n)
        .collect();
    let y_inv: Vec<Field<M>> = std::iter::successors(Some(z.clone()), |prev| Some(&z * prev))
        .take(n)
        .collect();

    let x_n: Vec<Field<M>> =
        std::iter::successors(Some(x.pow(n as u64 + 1)), |prev| Some(x * prev))
            .take(n)
            .collect();

    let y_n = y.pow(n as u64);
    let lhs = y_n * evaluate_s(x, y, n, q, u, v, w);

    let rhs = izip!(y_s, y_inv, x_n)
        .map(|(y, y_inv, x_s)| (y + y_inv) * x_s)
        .fold(Field::<M>::ZERO, |acc, val| acc + val);

    lhs - rhs
}

/// We would refactor this soon, right now it's only need to work properly
#[allow(clippy::too_many_arguments)]
pub fn evaluate_t<M: Modulus + std::fmt::Debug>(
    x: &Field<M>,
    y: &Field<M>,
    n: usize,
    q: usize,
    u: &[impl AsRef<[Field<M>]>],
    v: &[impl AsRef<[Field<M>]>],
    w: &[impl AsRef<[Field<M>]>],
    k: &[Field<M>],
    a: &[Field<M>],
    b: &[Field<M>],
    c: &[Field<M>],
) -> Field<M> {
    let lhs = evaluate_r(x, &Field::<M>::from(1), n, a, b, c);

    let lhs = lhs * (evaluate_r(x, y, n, a, b, c) + evaluate_ss(x, y, n, q, u, v, w));

    let rhs = y.pow(n as u64) * evaluate_poly(q, k, y);

    lhs - rhs
}

#[allow(clippy::too_many_arguments)]
pub fn calculate_k_vector<M: Modulus>(
    q: usize,
    a: &[Field<M>],
    b: &[Field<M>],
    c: &[Field<M>],
    u: &[impl AsRef<[Field<M>]>],
    v: &[impl AsRef<[Field<M>]>],
    w: &[impl AsRef<[Field<M>]>],
) -> Vec<Field<M>> {
    (0..q)
        .map(|q_i| {
            let u_i = u[q_i].as_ref();
            let v_i = v[q_i].as_ref();
            let w_i = w[q_i].as_ref();

            let dot_product = |coeffs: &[Field<M>], row: &[Field<M>]| {
                coeffs
                    .iter()
                    .zip(row.iter())
                    .map(|(coeff_i, val_i)| coeff_i * val_i)
                    .fold(Field::<M>::ZERO, |acc, val| acc + val)
            };

            let term_a = dot_product(a, u_i);
            let term_b = dot_product(b, v_i);
            let term_c = dot_product(c, w_i);

            term_a + term_b + term_c
        })
        .collect()
}

pub fn evaluate_poly<M: Modulus>(n: usize, poly: &[Field<M>], x: &Field<M>) -> Field<M> {
    assert_eq!(n, poly.len());

    std::iter::successors(Some(x.clone()), |prev| Some(x * prev))
        .take(n)
        .zip(poly.iter())
        .map(|(power, coeff_a)| coeff_a * power)
        .fold(Field::<M>::ZERO, |acc, val| acc + val)
}

#[cfg(test)]
mod utils_tests {
    use std::str::FromStr;

    use num_bigint::BigUint;

    use crate::{
        field::Field,
        modulus::{self, OrderP},
        utils::{evaluate_r, evaluate_s, evaluate_ss},
    };

    type FieldP = Field<modulus::OrderP>;

    #[test]
    fn evaluate_r_test() {
        let a = [FieldP::from(2), FieldP::from(3)];
        let b = [FieldP::from(4), FieldP::from(2)];

        let c = [FieldP::from(8), FieldP::from(6)];

        let r = evaluate_r(&FieldP::from(3), &FieldP::from(10), 2, &a, &b, &c);

        assert_eq!(
            r,
            FieldP::from_str(
                "8866028729065393948478499105841146159290165158213639216352728321613484282871"
            )
            .unwrap()
        );
    }

    #[test]
    fn evaluate_s_test() {
        let u = vec![
            [FieldP::from(2), FieldP::from(2)],
            [FieldP::from(2), FieldP::from(2)],
        ];
        let r = evaluate_s(
            &FieldP::from(3),
            &FieldP::from(10),
            2,
            2,
            &u.clone(),
            &u.clone(),
            &u.clone(),
        );

        assert_eq!(
            r,
            FieldP::from_str(
                "6432893846517566412420610278260439325191790329320366084373276166309636237186"
            )
            .unwrap()
        );
    }

    #[test]
    fn evaluate_ss_test() {
        let u = vec![
            [FieldP::from(2), FieldP::from(2)],
            [FieldP::from(2), FieldP::from(2)],
        ];
        let r = evaluate_ss(
            &FieldP::from(3),
            &FieldP::from(10),
            2,
            2,
            &u.clone(),
            &u.clone(),
            &u.clone(),
        );

        assert_eq!(
            r,
            FieldP::from_str(
                "1511730053931628106918843415391203241420070727390286029827719899082767150916"
            )
            .unwrap()
        );
    }

    #[test]
    fn evaluate_t_test() {
        let u = vec![
            [FieldP::from(2), FieldP::from(2)],
            [FieldP::from(2), FieldP::from(2)],
        ];
        let r = evaluate_ss(
            &FieldP::from(3),
            &FieldP::from(10),
            2,
            2,
            &u.clone(),
            &u.clone(),
            &u.clone(),
        );

        assert_eq!(
            r,
            FieldP::from_str(
                "1511730053931628106918843415391203241420070727390286029827719899082767150916"
            )
            .unwrap()
        );
    }
}
