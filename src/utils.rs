use itertools::izip;

use crate::{
    field::{self, Field},
    modulus::Modulus,
    poly::Polynomial,
};

fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>>
where
    T: Clone,
{
    assert!(!v.is_empty());
    (0..v[0].len())
        .map(|i| v.iter().map(|inner| inner[i].clone()).collect::<Vec<T>>())
        .collect()
}

pub struct ProoverHelper<M: Modulus> {
    pub n: usize,
    pub q: usize,

    pub a: Vec<Field<M>>,
    pub b: Vec<Field<M>>,
    pub c: Vec<Field<M>>,

    // NOTE: we assume that this matricies is QxN, and not NxQ, because it simpler to work with
    pub u: Vec<Vec<Field<M>>>,
    pub v: Vec<Vec<Field<M>>>,
    pub w: Vec<Vec<Field<M>>>,
    pub k: Vec<Field<M>>,
}

impl<M: Modulus> ProoverHelper<M> {
    pub fn new(
        n: usize,
        q: usize,
        a: Vec<Field<M>>,
        b: Vec<Field<M>>,
        c: Vec<Field<M>>,
        u: Vec<Vec<Field<M>>>,
        v: Vec<Vec<Field<M>>>,
        w: Vec<Vec<Field<M>>>,
    ) -> Self {
        let k = calculate_k_vector(q, &a, &b, &c, &u, &v, &w);
        Self {
            n,
            q,
            a,
            b,
            c,
            u: transpose(u),
            v: transpose(v),
            w: transpose(w),
            k,
        }
    }

    pub fn r_evaluate_at_x(&self, x: &Field<M>) -> Polynomial<M> {
        let z = x.inv().expect("value must be invertible");
        let z_n = z.pow(self.n as u64 + 1);

        let coeff_a: Vec<Field<M>> = std::iter::successors(Some(x.clone()), |prev| Some(x * prev))
            .take(self.n)
            .zip(self.a.iter())
            .map(|(x, a)| x * a)
            .collect();

        let poly_a = Polynomial {
            coefficients: coeff_a,
            offset: 1,
        };

        let mut coeff_b: Vec<Field<M>> =
            std::iter::successors(Some(z.clone()), |prev| Some(&z * prev))
                .take(self.n)
                .zip(self.b.iter())
                .map(|(x, b)| x * b)
                .collect();

        coeff_b.reverse();
        let poly_b = Polynomial {
            coefficients: coeff_b,
            offset: -(self.n as i64),
        };

        let mut coeff_c: Vec<Field<M>> =
            std::iter::successors(Some(z_n.clone()), |prev| Some(&z * prev))
                .take(self.n)
                .zip(self.c.iter())
                .map(|(x, c)| x * c)
                .collect();

        coeff_c.reverse();

        let poly_c = Polynomial {
            coefficients: coeff_c,
            offset: -(2 * self.n as i64),
        };

        &(&poly_a + &poly_b) + &poly_c
    }

    pub fn evaluate_s_at_x(&self, x: &Field<M>) -> Polynomial<M> {
        let z = x.inv().expect("X must be invertible");
        let x_n = x.pow(self.n as u64 + 1);

        let poly_a = std::iter::successors(Some(z.clone()), |prev| Some(&z * prev))
            .take(self.n)
            .zip(self.u.iter())
            .map(|(x, a)| {
                let poly_a = Polynomial {
                    coefficients: (*a).to_vec(),
                    offset: 1,
                };
                &poly_a * &x
            })
            .fold(Polynomial::<M>::empty(), |acc, val| &acc + &val);

        let poly_b = std::iter::successors(Some(x.clone()), |prev| Some(x * prev))
            .take(self.n)
            .zip(self.v.iter())
            .map(|(x, a)| {
                let poly_b = Polynomial {
                    coefficients: (*a).to_vec(),
                    offset: 1,
                };
                &poly_b * &x
            })
            .fold(Polynomial::<M>::empty(), |acc, val| &acc + &val);

        let poly_c = std::iter::successors(Some(x_n.clone()), |prev| Some(x * prev))
            .take(self.n)
            .zip(self.w.iter())
            .map(|(x, a)| {
                let poly_a = Polynomial {
                    coefficients: (*a).to_vec(),
                    offset: 1,
                };
                &poly_a * &x
            })
            .fold(Polynomial::<M>::empty(), |acc, val| &acc + &val);

        &(&poly_a + &poly_b) + &poly_c
    }

    pub fn evaluate_s_at_y(&self, y: &Field<M>) -> Polynomial<M> {
        let y_powers: Vec<Field<M>> = std::iter::successors(Some(y.clone()), |prev| Some(y * prev))
            .take(self.q)
            .collect();

        let (mut u, v, w): (Vec<Field<M>>, Vec<Field<M>>, Vec<Field<M>>) =
            izip!(&self.u, &self.v, &self.w)
                .map(|(row_u, row_v, row_w)| {
                    let dot_product = |coeffs: &[Field<M>], row: &[Field<M>]| {
                        coeffs
                            .iter()
                            .zip(row.iter())
                            .map(|(coeff_i, val_i)| coeff_i * val_i)
                            .fold(Field::<M>::ZERO, |acc, val| acc + val)
                    };

                    (
                        dot_product(&y_powers, row_u),
                        dot_product(&y_powers, row_v),
                        dot_product(&y_powers, row_w),
                    )
                })
                .collect();
        u.reverse();
        let poly_a = Polynomial {
            coefficients: u,
            offset: -(self.n as i64),
        };

        let poly_b = Polynomial {
            coefficients: v,
            offset: 1,
        };

        let poly_c = Polynomial {
            coefficients: w,
            offset: self.n as i64 + 1,
        };

        &(&poly_a + &poly_b) + &poly_c
    }
    pub fn evaluate_ss_at_x(&self, x: &Field<M>) -> Polynomial<M> {
        let mut eval_s = self.evaluate_s_at_x(x);
        println!("s: {eval_s}");

        eval_s.shift(self.n as i64);

        let mut x_n: Vec<Field<M>> =
            std::iter::successors(Some(x.pow(self.n as u64 + 1)), |prev| Some(x * prev))
                .take(self.n)
                .collect();

        let poly_y = Polynomial {
            coefficients: x_n.clone(),
            offset: 1,
        };

        x_n.reverse();

        let poly_y_inv = Polynomial {
            coefficients: x_n,
            offset: -(self.n as i64),
        };

        &eval_s - &(&poly_y_inv + &poly_y)
    }

    pub fn evaluate_ss_at_y(&self, y: &Field<M>) -> Polynomial<M> {
        let mut eval_s = &self.evaluate_s_at_y(y) * &y.pow(self.n as u64);

        let y_powers: Vec<Field<M>> = std::iter::successors(Some(y.clone()), |prev| Some(y * prev))
            .take(self.n)
            .collect();

        let y_inv = y.inv().expect("Y must be invertible");

        let y_inv_powers: Vec<Field<M>> =
            std::iter::successors(Some(y_inv.clone()), |prev| Some(&y_inv * prev))
                .take(self.n)
                .collect();

        let coeff = izip!(y_powers, y_inv_powers)
            .map(|(y, y_inv)| y + y_inv)
            .collect();
        let poly = Polynomial {
            coefficients: coeff,
            offset: self.n as i64 + 1,
        };

        &eval_s - &poly
    }

    pub fn evaluate_t_at_x(&self, x: &Field<M>) -> Polynomial<M> {
        let r = self.r_evaluate_at_x(x);

        let eval_ss = self.evaluate_ss_at_x(x);

        let r_x1 = r.evaluate_at(Field::<M>::one());

        let k = Polynomial {
            coefficients: self.k.iter().map(|x| x.inv_additive()).collect(),
            offset: 1 + self.n as i64,
        };

        let mult = &(&eval_ss + &r) * &r_x1;
        &(mult) + &k
    }

    pub fn evaluate_t_at_y(&self, y: &Field<M>) -> Polynomial<M> {
        let r = self.r_evaluate_at_x(y);
        let eval_ss = self.evaluate_ss_at_y(y);

        let r_1 = self.r_evaluate_at_x(&Field::one());

        let y = std::iter::successors(Some(y.pow(self.n as u64 + 1).clone()), |prev| {
            Some(y * prev)
        })
        .take(self.q)
        .zip(&self.k)
        .map(|(y, k)| y)
        .fold(Field::<M>::ZERO, |acc, val| acc + val);

        let y_poly = Polynomial {
            coefficients: vec![y],
            offset: 0,
        };

        &(&(&r + &eval_ss) * &r_1) + &y_poly
    }
}

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

pub fn evaluate_s<M: Modulus + std::fmt::Debug>(
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

    let mut u_i = vec![Field::<M>::ZERO; n];
    let mut v_i = vec![Field::<M>::ZERO; n];
    let mut w_i = vec![Field::<M>::ZERO; n];

    izip!(y_powers, u, v, w).for_each(|(pwr_y, row_u, row_v, row_w)| {
        let update = |target: &mut [Field<M>], row: &dyn AsRef<[Field<M>]>| {
            target
                .iter_mut()
                .zip(row.as_ref().iter())
                .for_each(|(acc, coeff)| {
                    *acc += (coeff * &pwr_y);
                });
        };

        update(&mut u_i, row_u);
        update(&mut v_i, row_v);
        update(&mut w_i, row_w);
    });
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
mod poly_test {
    use std::sync::LazyLock;

    use crate::{
        field::Field,
        modulus::{self, OrderP},
        utils::{ProoverHelper, transpose},
    };

    const N: usize = 4;
    const Q: usize = 3;

    type FieldP = Field<modulus::OrderP>;

    pub static TEST_DATA: LazyLock<ProoverHelper<OrderP>> = LazyLock::new(|| {
        ProoverHelper::new(
            N,
            Q,
            Vec::from([854, 114, 212, 987].map(FieldP::from)),
            Vec::from([345, 678, 123, 888].map(FieldP::from)),
            Vec::from([555, 444, 333, 222].map(FieldP::from)),
            Vec::from([
                Vec::from([123, 456, 789, 101].map(FieldP::from)),
                Vec::from([202, 303, 404, 505].map(FieldP::from)),
                Vec::from([606, 707, 808, 909].map(FieldP::from)),
            ]),
            Vec::from([
                Vec::from([111, 222, 333, 444].map(FieldP::from)),
                Vec::from([555, 666, 777, 888].map(FieldP::from)),
                Vec::from([999, 100, 200, 300].map(FieldP::from)),
            ]),
            Vec::from([
                Vec::from([150, 250, 350, 450].map(FieldP::from)),
                Vec::from([550, 650, 750, 850].map(FieldP::from)),
                Vec::from([950, 105, 205, 305].map(FieldP::from)),
            ]),
        )
    });

    #[test]
    fn s_test() {
        println!("t1 {}\n", TEST_DATA.evaluate_t_at_x(&FieldP::from(228)));
        println!("t2 {}\n", TEST_DATA.evaluate_t_at_y(&FieldP::from(1337)));
    }
}

#[cfg(test)]
mod utils_tests {
    use std::str::FromStr;

    use num_bigint::BigUint;

    use crate::{
        field::Field,
        modulus::{self, OrderP},
        utils::{calculate_k_vector, evaluate_r, evaluate_s, evaluate_ss, evaluate_t},
    };

    use std::sync::LazyLock;

    const N: usize = 4;
    const Q: usize = 3;

    pub static TEST_DATA: LazyLock<TestData> = LazyLock::new(|| TestData {
        x: FieldP::from(228),
        y: FieldP::from(1337),

        a: [854, 114, 212, 987].map(FieldP::from),
        b: [345, 678, 123, 888].map(FieldP::from),
        c: [555, 444, 333, 222].map(FieldP::from),
        u: [
            [123, 456, 789, 101].map(FieldP::from),
            [202, 303, 404, 505].map(FieldP::from),
            [606, 707, 808, 909].map(FieldP::from),
        ],
        v: [
            [111, 222, 333, 444].map(FieldP::from),
            [555, 666, 777, 888].map(FieldP::from),
            [999, 100, 200, 300].map(FieldP::from),
        ],
        w: [
            [150, 250, 350, 450].map(FieldP::from),
            [550, 650, 750, 850].map(FieldP::from),
            [950, 105, 205, 305].map(FieldP::from),
        ],
    });

    pub struct TestData {
        pub x: FieldP,
        pub y: FieldP,

        pub a: [FieldP; N],
        pub b: [FieldP; N],
        pub c: [FieldP; N],
        pub u: [[FieldP; N]; Q],
        pub v: [[FieldP; N]; Q],
        pub w: [[FieldP; N]; Q],
    }

    type FieldP = Field<modulus::OrderP>;

    #[test]
    fn evaluate_r_test() {
        let r = evaluate_r(
            &TEST_DATA.x,
            &TEST_DATA.y,
            N,
            &TEST_DATA.a,
            &TEST_DATA.b,
            &TEST_DATA.c,
        );

        assert_eq!(
            r,
            FieldP::from_str(
                "9883121794604549144343819521516913311653873111753789353038581158791770311053"
            )
            .unwrap()
        );
    }

    #[test]
    fn evaluate_s_test() {
        let s = evaluate_s(
            &TEST_DATA.x,
            &TEST_DATA.y,
            N,
            Q,
            &TEST_DATA.u,
            &TEST_DATA.v,
            &TEST_DATA.w,
        );

        assert_eq!(
            s,
            FieldP::from_str(
                "26379642808924966178563811802021192142145621114491300471681632055362572754274"
            )
            .unwrap()
        );
    }

    #[test]
    fn evaluate_ss_test() {
        let t = evaluate_ss(
            &TEST_DATA.x,
            &TEST_DATA.y,
            N,
            Q,
            &TEST_DATA.u,
            &TEST_DATA.v,
            &TEST_DATA.w,
        );

        assert_eq!(
            t,
            FieldP::from_str(
                "5783717380048885944422759096669824392116921961451678152134778150967934883412"
            )
            .unwrap()
        );
    }

    #[test]
    fn evaluate_k_test() {
        let k: Vec<FieldP> = calculate_k_vector(
            Q,
            &TEST_DATA.a,
            &TEST_DATA.b,
            &TEST_DATA.c,
            &TEST_DATA.u,
            &TEST_DATA.v,
            &TEST_DATA.w,
        );

        assert_eq!(
            k,
            [
                FieldP::from(1458723),
                FieldP::from(3350571),
                FieldP::from(3079901)
            ]
        );
    }

    #[test]
    fn evaluate_t_test() {
        let k: Vec<FieldP> = calculate_k_vector(
            Q,
            &TEST_DATA.a,
            &TEST_DATA.b,
            &TEST_DATA.c,
            &TEST_DATA.u,
            &TEST_DATA.v,
            &TEST_DATA.w,
        );
        let t = evaluate_t(
            &TEST_DATA.x,
            &TEST_DATA.y,
            N,
            Q,
            &TEST_DATA.u,
            &TEST_DATA.v,
            &TEST_DATA.w,
            &k,
            &TEST_DATA.a,
            &TEST_DATA.b,
            &TEST_DATA.c,
        );

        assert_eq!(
            t,
            FieldP::from_str(
                "9577296933606094075357504016059020464592551159956904365253909461538337282914"
            )
            .unwrap()
        );
    }
}
