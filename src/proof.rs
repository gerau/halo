use std::default;

use num_bigint::BigUint;
use rand::Rng;
use rand::prelude::*;

use crate::poly::Polynomial;
use crate::utils::OpeningData;
use crate::utils::ProoverHelper;
use crate::utils::compute_g_poly;
use crate::utils::compute_h_at_y;
use crate::{
    elliptic_curve::CurvePoint,
    field::{self, Field},
    modulus::Modulus,
};

use sha3::{
    Shake128,
    digest::{ExtendableOutput, Update, XofReader},
};

pub enum Data<'a, M: Modulus> {
    Field(&'a Field<M>),
    Point(&'a CurvePoint<M>),
}

pub fn sample_challenge<M: Modulus>(inputs: &[Data<M>]) -> Field<M> {
    let mut hasher = Shake128::default();

    for input in inputs {
        match input {
            Data::Field(f) => hasher.update(&f.number.to_bytes_be()),
            Data::Point(p) => {
                let x_bytes = &p.x.number.to_bytes_be();
                let mut y_bytes = p.y.number.to_bytes_be();
                y_bytes.extend_from_slice(&x_bytes);

                hasher.update(&y_bytes);
            }
        }
    }

    let mut reader = hasher.finalize_xof();

    read_field_from_xof(&mut reader)
}

pub fn sample_two_challenges<M: Modulus>(inputs: &[Data<M>]) -> (Field<M>, Field<M>) {
    let mut hasher = Shake128::default();

    for input in inputs {
        match input {
            Data::Field(f) => hasher.update(&f.number.to_bytes_be()),
            Data::Point(p) => {
                let x_bytes = &p.x.number.to_bytes_be();
                let mut y_bytes = p.y.number.to_bytes_be();
                y_bytes.extend_from_slice(&x_bytes);

                hasher.update(&y_bytes);
            }
        }
    }

    let mut reader = hasher.finalize_xof();

    let c1 = read_field_from_xof(&mut reader);
    let c2 = read_field_from_xof(&mut reader); // XOF continues generating new bytes

    (c1, c2)
}

fn read_field_from_xof<M: Modulus>(reader: &mut impl XofReader) -> Field<M> {
    let mut buf = [0u8; 64];
    reader.read(&mut buf);

    Field::<M>::new(BigUint::from_bytes_le(&buf))
}

fn inner_product<M: Modulus>(a: &[Field<M>], b: &[Field<M>]) -> Field<M> {
    a.iter()
        .zip(b)
        .map(|(x, y)| x * y)
        .fold(Field::<M>::ZERO, |acc, val| acc + val)
}

fn split_vec<T: Clone>(vec: &[T]) -> (Vec<T>, Vec<T>) {
    let mid = vec.len() / 2;
    (vec[0..mid].to_vec(), vec[mid..].to_vec())
}

fn evaluation_step<M: Modulus>(
    a: &[Field<M>],
    b: &[Field<M>],
    g: &[CurvePoint<M>],
    h_blind: &CurvePoint<M>,
    u: &CurvePoint<M>,
    rng: &mut impl rand::Rng,
) -> ((CurvePoint<M>, CurvePoint<M>), (Field<M>, Field<M>)) {
    let (a_l, a_r) = split_vec(a);
    let (b_l, b_r) = split_vec(b);
    let (g_l, g_r) = split_vec(g);

    let c_l = inner_product(&a_l, &b_r);
    let c_r = inner_product(&a_r, &b_l);

    let l_blind = Field::<M>::get_random(rng);
    let r_blind = Field::<M>::get_random(rng);

    let commitment_l = CurvePoint::msm(&g_r, &a_l);
    let l = (commitment_l + u * &c_l) + h_blind * &l_blind;

    let commitment_r = CurvePoint::msm(&g_l, &a_r);
    let r = (commitment_r + u * &c_r) + h_blind * &r_blind;

    ((l, r), (l_blind, r_blind))
}

fn update_vectors<M: Modulus>(
    a: Vec<Field<M>>,
    b: Vec<Field<M>>,
    g: Vec<CurvePoint<M>>,
    u: &Field<M>,
) -> (Vec<Field<M>>, Vec<Field<M>>, Vec<CurvePoint<M>>) {
    let u_inv = u.inv().expect("Challenge u must be invertible");

    let (a_l, a_r) = split_vec(&a);
    let (b_l, b_r) = split_vec(&b);
    let (g_l, g_r) = split_vec(&g);

    let a_new: Vec<Field<M>> = a_l
        .iter()
        .zip(&a_r)
        .map(|(al, ar)| al * &u_inv + ar * u)
        .collect();

    let b_new: Vec<Field<M>> = b_l
        .iter()
        .zip(&b_r)
        .map(|(bl, br)| bl * u + br * &u_inv)
        .collect();

    let g_new: Vec<CurvePoint<M>> = g_l
        .iter()
        .zip(&g_r)
        .map(|(gl, gr)| gl * &u_inv + gr * &u)
        .collect::<Vec<_>>();

    (a_new, b_new, g_new)
}

type Commit<M> = CurvePoint<M>;

#[derive(Debug, Default, Clone)]
pub struct ProverParams<M: Modulus> {
    pub g: Vec<CurvePoint<M>>,
    pub h: CurvePoint<M>,
}

impl<M: Modulus> ProverParams<M> {
    pub fn new(k: usize) -> Self {
        let mut rng = rand::thread_rng();
        let d = 1 << k;
        ProverParams {
            g: (0..d).map(|_| CurvePoint::get_random(&mut rng)).collect(),
            h: CurvePoint::get_random(&mut rng),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct ProverState<M: Modulus> {
    pub z1: Field<M>,
    pub z2: Field<M>,
    pub z3: Field<M>,
    pub z4: Field<M>,

    pub delta_r: Field<M>,
    pub delta_tlo: Field<M>,
    pub delta_thi: Field<M>,
    pub delta_h: Field<M>,
    pub r_x1_poly: Polynomial<M>,
    pub s_xy_old_poly: Polynomial<M>,
    pub s_xy_poly: Polynomial<M>,
    pub s_xy_new_poly: Polynomial<M>,

    pub t_lo_poly: Polynomial<M>,
    pub t_hi_poly: Polynomial<M>,
    pub g_x_poly: Polynomial<M>,
    pub h_x_poly: Polynomial<M>,

    pub k_y_poly: Polynomial<M>,
    pub s_at_x_poly_y: Polynomial<M>,

    pub v: [Field<M>; 13],

    pub r: Commit<M>,
    pub t_low: Commit<M>,
    pub t_high: Commit<M>,
    pub s: Commit<M>,
    pub s_new: Commit<M>,
    pub c: Commit<M>,
    pub h: Commit<M>,
}

#[derive(Clone, Debug)]
pub struct ZkProof<M: Modulus> {
    pub l_vec: Vec<CurvePoint<M>>,
    pub r_vec: Vec<CurvePoint<M>>,
    pub r_zk: CurvePoint<M>,
    pub z1_zk: Field<M>,
    pub z2_zk: Field<M>,
    pub final_g: CurvePoint<M>,
    pub open_value: Field<M>,
}

#[derive(Debug, Default)]
pub struct Prover<M: Modulus> {
    pub n: usize,
    pub d: usize,
    pub k: usize,

    pub state: ProverState<M>,
    sigma: ProverParams<M>,
}

fn generate_polynom_for_prove<M: Modulus>(
    z1: &Field<M>,
    z4: &Field<M>,
    groups: &[Vec<&Polynomial<M>>],
) -> Polynomial<M> {
    let mut grand_total = Polynomial::<M>::empty();
    let mut z4_pow = Field::<M>::one();

    for group in groups {
        let mut group_sum = Polynomial::<M>::empty();
        let mut z1_pow = Field::<M>::one();

        for poly in group {
            let term = *poly * &z1_pow;
            group_sum = &group_sum + &term;
            z1_pow = &z1_pow * z1;
        }

        let group_term = &group_sum * &z4_pow;
        grand_total = &grand_total + &group_term;

        z4_pow = &z4_pow * z4;
    }

    grand_total
}

fn compute_r_blind<M: Modulus>(
    l_blinds: &[Field<M>],
    r_blinds: &[Field<M>],
    u_vec: &[Field<M>],
    initial_blind: &Field<M>,
) -> Field<M> {
    let mut r_total = initial_blind.clone();

    for ((l_b, r_b), u) in l_blinds.iter().zip(r_blinds).zip(u_vec) {
        let u_sq = u * u;
        let u_inv = u.inv().unwrap();
        let u_inv_sq = &u_inv * &u_inv;

        r_total = r_total + (l_b * &u_sq) + (r_b * &u_inv_sq);
    }
    r_total
}

impl<M: Modulus + Default + std::fmt::Debug> Prover<M> {
    pub fn setup(sigma: ProverParams<M>, k: usize) -> Self {
        let d = 1 << k;

        let n = d / 4;

        Self {
            n,
            d,
            k,
            sigma,
            ..Default::default()
        }
    }

    pub fn commit(&self, poly: &Polynomial<M>, delta: &Field<M>) -> Commit<M> {
        let padded_poly = poly.pad_to_zero_offset(self.d);

        let lhs = self
            .sigma
            .g
            .iter()
            .zip(padded_poly.coefficients)
            .map(|(point, scalar)| point.clone() * scalar)
            .fold(CurvePoint::<M>::new_infinity(), |acc, val| acc + val);

        lhs + &self.sigma.h.clone() * delta
    }

    pub fn generate_commitments(
        &mut self,
        proover: &ProoverHelper<M>,
        y_old: &Field<M>,
        s_old_poly: &Polynomial<M>,
        u_old: &[Field<M>],
    ) {
        let mut rng = rand::thread_rng();

        self.state.delta_r = Field::<M>::get_random(&mut rng);

        self.state.r_x1_poly = proover
            .evaluate_r_at_x(&Field::<M>::from(1))
            .shift_poly((3 * proover.n - 1) as i64);

        self.state.r = self.commit(&self.state.r_x1_poly, &self.state.delta_r);

        let y: Field<M> = sample_challenge(&[Data::Point(&self.state.r)]);

        self.state.delta_tlo = Field::<M>::get_random(&mut rng);
        self.state.delta_thi = Field::<M>::get_random(&mut rng);

        let t_poly = proover.evaluate_t_at_y(&y);

        self.state.t_lo_poly = t_poly.extract_t_lo(self.d);
        self.state.t_hi_poly = t_poly.extract_t_hi(self.d);

        self.state.t_low = self.commit(&self.state.t_lo_poly, &self.state.delta_tlo);
        self.state.t_high = self.commit(&self.state.t_hi_poly, &self.state.delta_thi);

        self.state.s_xy_poly = proover.evaluate_s_at_y(&y).shift_poly(self.n as i64);

        self.state.s_xy_old_poly = s_old_poly.clone();

        self.state.s = self.commit(&self.state.s_xy_poly, &Field::<M>::get_random(&mut rng));

        let x: Field<M> = sample_challenge(&[
            Data::Point(&self.state.t_low),
            Data::Point(&self.state.t_high),
            Data::Point(&self.state.s),
        ]);

        println!(
            "tlo: {}\n thi: {}\n t: {}",
            self.state.t_lo_poly, self.state.t_hi_poly, t_poly
        );
        println!("t_poly: {}", t_poly.evaluate_at(&x));
        let x_n = x.pow(self.n as u64);
        let s_at_x = proover.evaluate_s_at_x(&x);
        self.state.s_at_x_poly_y = &s_at_x * &x_n;

        self.state.k_y_poly = Polynomial {
            coefficients: proover.k.clone(),
            offset: 1,
        };

        self.state.c = self.commit(&self.state.s_at_x_poly_y, &Field::<M>::get_random(&mut rng));

        let y_new: Field<M> = sample_challenge(&[Data::Point(&self.state.c)]);

        self.state.s_xy_new_poly = proover.evaluate_s_at_y(&y_new).shift_poly(self.n as i64);

        self.state.s_new =
            self.commit(&self.state.s_xy_new_poly, &Field::<M>::get_random(&mut rng));

        self.state.g_x_poly = compute_g_poly(u_old);

        self.state.v[0] = self.state.r_x1_poly.evaluate_at(&x);
        self.state.v[1] = self.state.s_at_x_poly_y.evaluate_at(y_old);
        self.state.v[2] = self.state.s_at_x_poly_y.evaluate_at(&y);
        self.state.v[3] = self.state.s_at_x_poly_y.evaluate_at(&y_new);
        self.state.v[4] = self.state.t_lo_poly.evaluate_at(&x);
        self.state.v[5] = self.state.t_hi_poly.evaluate_at(&x);
        self.state.v[6] = self.state.g_x_poly.evaluate_at(&x);
        self.state.v[7] = self.state.k_y_poly.evaluate_at(&y);

        self.state.v[8] =
            proover.evaluate_r_at_x(&x).evaluate_at(&y) * &(&x * &y).pow((3 * self.n - 1) as u64);

        let (z1, z2) = sample_two_challenges(&[
            Data::Point(&self.state.s_new),
            Data::Field(&self.state.v[0]),
            Data::Field(&self.state.v[1]),
            Data::Field(&self.state.v[2]),
            Data::Field(&self.state.v[3]),
            Data::Field(&self.state.v[4]),
            Data::Field(&self.state.v[5]),
            Data::Field(&self.state.v[6]),
            Data::Field(&self.state.v[7]),
            Data::Field(&self.state.v[8]),
        ]);
        self.state.z1 = z1;
        self.state.z2 = z2;

        let p_poly =
            proover.compute_p_poly(&y, &y_new, &self.state.z1, u_old, &self.state.s_xy_old_poly);

        let q_poly = &self.state.k_y_poly + &(&self.state.s_at_x_poly_y * &self.state.z1);

        let openings = vec![
            OpeningData {
                polynomial: p_poly,
                point: x.clone(),
                values: vec![
                    self.state.v[0].clone(),
                    self.state.s_xy_old_poly.evaluate_at(&x) * &x_n,
                    self.state.v[2].clone(),
                    self.state.v[3].clone(),
                    self.state.v[4].clone(),
                    self.state.v[5].clone(),
                    self.state.v[6].clone(),
                ],
            },
            OpeningData {
                polynomial: q_poly.pad_to_zero_offset(self.d),
                point: y.clone(),
                values: vec![self.state.v[7].clone(), self.state.v[2].clone()],
            },
            OpeningData {
                polynomial: self.state.r_x1_poly.pad_to_zero_offset(self.d),
                point: &x * &y,
                values: vec![self.state.v[8].clone()],
            },
            OpeningData {
                polynomial: self.state.s_at_x_poly_y.pad_to_zero_offset(self.d),
                point: y_old.clone(),
                values: vec![self.state.v[1].clone()],
            },
            OpeningData {
                polynomial: self.state.s_at_x_poly_y.pad_to_zero_offset(self.d),
                point: y_new.clone(),
                values: vec![self.state.v[3].clone()],
            },
        ];

        self.state.h_x_poly = compute_h_at_y(&openings, &self.state.z2, &self.state.z1);
        self.state.delta_h = Field::<M>::get_random(&mut rng);

        self.state.h = self.commit(&self.state.h_x_poly, &self.state.delta_h);

        self.state.z3 = sample_challenge(&[Data::Point(&self.state.h)]);

        self.state.v[9] = proover
            .compute_p_poly(&y, &y_new, &self.state.z1, u_old, &self.state.s_xy_old_poly)
            .evaluate_at(&self.state.z3);

        self.state.v[10] = proover
            .compute_q_poly(&x, &self.state.z1)
            .evaluate_at(&self.state.z3);

        self.state.v[11] = self.state.r_x1_poly.evaluate_at(&self.state.z3);

        self.state.v[12] = self.state.s_at_x_poly_y.evaluate_at(&self.state.z3);

        self.state.z4 = sample_challenge(&[
            Data::Field(&self.state.v[9]),
            Data::Field(&self.state.v[10]),
            Data::Field(&self.state.v[11]),
            Data::Field(&self.state.v[12]),
        ]);

        println!(
            "z: {}, {}, {}, {}",
            self.state.z1, self.state.z2, self.state.z3, self.state.z4
        );
    }

    pub fn compute_open_value_at_z3(&self) -> Field<M> {
        let v = &self.state.v;
        let z4 = &self.state.z4;
        let z3 = &self.state.z3;

        let h_at_z3 = self.state.h_x_poly.evaluate_at(z3);

        let mut val = v[9].clone();
        let mut z4_pow = z4.clone();

        val = val + &v[10] * &z4_pow;
        z4_pow = &z4_pow * z4;

        val = val + &v[11] * &z4_pow;
        z4_pow = &z4_pow * z4;

        val = val + &v[12] * &z4_pow;
        z4_pow = &z4_pow * z4;

        val = val + h_at_z3 * &z4_pow;

        val
    }

    pub fn zero_knowledge_opening(&mut self, u: &CurvePoint<M>) -> ZkProof<M> {
        let mut rng = rand::thread_rng();

        let groups: Vec<Vec<&Polynomial<M>>> = vec![
            vec![
                &self.state.r_x1_poly,
                &self.state.s_xy_old_poly,
                &self.state.s_xy_poly,
                &self.state.s_xy_new_poly,
                &self.state.t_lo_poly,
                &self.state.t_hi_poly,
                &self.state.g_x_poly,
            ],
            vec![&self.state.k_y_poly, &self.state.s_at_x_poly_y],
            vec![&self.state.r_x1_poly],
            vec![&self.state.s_at_x_poly_y],
            vec![&self.state.h_x_poly],
        ];

        let poly_prove = generate_polynom_for_prove(&self.state.z1, &self.state.z4, &groups);

        let z1 = &self.state.z1;
        let z4 = &self.state.z4;

        let blind_for_prove = self.state.delta_r.clone()
            + z1.pow(4) * &self.state.delta_tlo
            + z1.pow(5) * &self.state.delta_thi
            + z4.pow(2) * &self.state.delta_r
            + z4.pow(4) * &self.state.delta_h;

        let mut a_vec = poly_prove.coefficients.clone();

        if a_vec.len() < self.d {
            a_vec.resize(self.d, Field::<M>::ZERO);
        }

        let mut b_vec = Vec::with_capacity(self.d);
        let mut z3_pow = Field::<M>::one();
        for _ in 0..self.d {
            b_vec.push(z3_pow.clone());
            z3_pow = &z3_pow * &self.state.z3;
        }

        let mut g_vec = self.sigma.g.clone();
        if g_vec.len() > self.d {
            g_vec.truncate(self.d);
        }

        let h_blind = &self.sigma.h;

        let mut l_vec = Vec::new();
        let mut r_vec = Vec::new();
        let mut l_blind_vec = Vec::new();
        let mut r_blind_vec = Vec::new();
        let mut u_vec = Vec::new();

        for _ in 0..self.k {
            let ((l_point, r_point), (l_blind, r_blind)) =
                evaluation_step(&a_vec, &b_vec, &g_vec, h_blind, u, &mut rng);

            l_vec.push(l_point.clone());
            r_vec.push(r_point.clone());
            l_blind_vec.push(l_blind);
            r_blind_vec.push(r_blind);

            let u_j: Field<M> = sample_challenge(&[Data::Point(&l_point), Data::Point(&r_point)]);
            u_vec.push(u_j.clone());

            let (next_a, next_b, next_g) = update_vectors(a_vec, b_vec, g_vec, &u_j);
            a_vec = next_a;
            b_vec = next_b;
            g_vec = next_g;
        }

        let r_blind = compute_r_blind(&l_blind_vec, &r_blind_vec, &u_vec, &blind_for_prove);

        let d_blind = Field::<M>::get_random(&mut rng);
        let s_blind = Field::<M>::get_random(&mut rng);

        let final_a = &a_vec[0]; // Not used in commitment, but is the secret value
        let final_b = &b_vec[0];
        let final_g = &g_vec[0];

        let base_point = final_g + &(u * final_b);
        let r_zk_point = (&base_point * &d_blind) + (h_blind * &s_blind);

        let c: Field<M> = sample_challenge(&[Data::Point(&r_zk_point)]);

        let z1_zk = (final_a * &c) + &d_blind;

        let z2_zk = (&r_blind * &c) + &s_blind;

        let open_value = self.compute_open_value_at_z3();

        ZkProof {
            l_vec,
            r_vec,
            r_zk: r_zk_point,
            z1_zk,
            z2_zk,
            final_g: final_g.clone(),
            open_value,
        }
    }
}

#[derive(Debug, Clone)]
pub struct VerifierProof<M: Modulus> {
    pub r_comm: Commit<M>,
    pub s_comm: Commit<M>,
    pub s_new_comm: Commit<M>,
    pub t_lo_comm: Commit<M>,
    pub t_hi_comm: Commit<M>,
    pub c_comm: Commit<M>,
    pub h_comm: Commit<M>,

    pub values: Vec<Field<M>>,

    pub zk_proof: ZkProof<M>,
}

pub struct Verifier<M: Modulus> {
    pub sigma: ProverParams<M>,
    pub n: usize,
    pub d: usize,
    pub k: usize,
}

impl<M: Modulus + std::fmt::Debug> Verifier<M> {
    pub fn new(sigma: ProverParams<M>, d: usize) -> Self {
        Self {
            sigma,
            d,
            n: d / 4,
            k: (d as f64).log2() as usize,
        }
    }

    pub fn check_g(&self, v7: &Field<M>, x: &Field<M>, u_old: &[Field<M>]) -> bool {
        let expected_v7 = compute_g_poly(u_old).evaluate_at(x);

        *v7 == expected_v7
    }

    pub fn check_polynomial_identity(&self, v: &[Field<M>], x: &Field<M>, y: &Field<M>) -> bool {
        let n = self.n as u64;
        let d = self.d as u64;
        let exp_3n_1 = (3 * self.n - 1) as u64;
        let x_inv = x.inv().unwrap();

        let x_n = x.pow(n);
        let y_n = y.pow(n);
        let xy = x * y;
        let xy_inv = xy.inv().unwrap();

        let y_inv = y.inv().unwrap();
        let mut y_i = y.clone();
        let mut y_inv_i = y.inv().unwrap();
        let mut x_pow = &x_n * x;

        let mut sum = Field::<M>::ZERO;
        for _ in 0..self.n {
            let term = (&y_i + &y_inv_i) * &x_pow;
            sum += term;

            y_i = y_i * y;
            y_inv_i = y_inv_i * &y_inv;
            x_pow = x_pow * x;
        }

        let term_v1 = &v[0] * x_inv.pow(exp_3n_1);

        let term_v9 = &v[8] * &xy_inv.pow(exp_3n_1);

        let term_v3 = &v[2] * &y_n * &x_inv.pow(n);

        let bracket = term_v9 + term_v3 - sum;

        let term_v8 = &v[7] * &y_n;

        let rhs = (term_v1 * bracket) - term_v8;

        let x_d_inv = x.inv().unwrap().pow(d);
        let term_v5 = &v[4] * &x_d_inv;

        let term_v6 = &v[5] * x;

        let lhs = term_v5.clone() + term_v6.clone();

        lhs == rhs
    }

    pub fn verify_commitment_check(
        &self,
        proof: &VerifierProof<M>,
        s_old_comm: &CurvePoint<M>,
        k_comm: &CurvePoint<M>,
        u_point: &CurvePoint<M>,
        x: &Field<M>,
        y: &Field<M>,
        y_old: &Field<M>,
        y_new: &Field<M>,
        z1: &Field<M>,
        z2: &Field<M>,
        z3: &Field<M>,
        z4: &Field<M>,
        c: &Field<M>,
        u_j_vec: &[Field<M>],
    ) -> bool {
        let z1_2 = z1 * z1;
        let z1_3 = &z1_2 * z1;
        let z1_4 = &z1_3 * z1;
        let z1_5 = &z1_4 * z1;
        let z1_6 = &z1_5 * z1;

        let z4_2 = z4 * z4;
        let z4_3 = &z4_2 * z4;
        let z4_4 = &z4_3 * z4;

        let mut p_comb = proof.r_comm.clone();
        p_comb += s_old_comm * z1;
        p_comb += &proof.s_comm * &z1_2;
        p_comb += &proof.s_new_comm * &z1_3;
        p_comb += &proof.t_lo_comm * &z1_4;
        p_comb += &proof.t_hi_comm * &z1_5;

        p_comb += k_comm.clone() * z4.clone();
        let z4_z1 = z4 * z1;
        p_comb += proof.c_comm.clone() * z4_z1;

        p_comb += proof.r_comm.clone() * z4_2.clone();
        p_comb += proof.c_comm.clone() * z4_3.clone();
        p_comb += proof.h_comm.clone() * z4_4.clone();

        // v = v10 + z4*v11 + z4^2*v12 + z4^3*v13 + z4^4 * [H_term]

        let v = &proof.values;

        let z2_2 = z2 * z2;
        let z2_3 = &z2_2 * z2;
        let z2_4 = &z2_3 * z2;

        //  (v10 - (v1 + z1*v2 + ... + z1^6*v7)) / (z3 - x)
        let mut p_val = v[0].clone();
        p_val += &v[1] * z1;
        p_val += &v[2] * &z1_2;
        p_val += &v[3] * &z1_3;
        p_val += &v[4] * &z1_4;
        p_val += &v[5] * &z1_5;
        p_val += &v[6] * &z1_6;

        let term_p_num = &v[9] - &p_val; // v10 - P(x)
        let term_p_den = z3 - x;
        let term_p = term_p_num * term_p_den.inv().unwrap();

        // z2 * (v11 - (v8 + z1*v3)) / (z3 - y)
        let q_val = &v[7] + &(&v[2] * z1); // v8 + z1*v3
        let term_q_num = &v[10] - &q_val; // v11 - Q(y)
        let term_q_den = z3 - y;
        let term_q = (term_q_num * term_q_den.inv().unwrap()) * z2;

        // z2^2 * (v12 - v9) / (z3 - xy)
        let term_r_num = &v[11] - &v[8]; // v12 - v9
        let term_r_den = z3 - &(x * y);
        let term_r = (term_r_num * term_r_den.inv().unwrap()) * &z2_2;

        // z2^3 * (v13 - v2) / (z3 - y_old)
        let term_c1_num = &v[12] - &v[1]; // v13 - v2
        let term_c1_den = z3 - y_old;
        let term_c1 = (term_c1_num * term_c1_den.inv().unwrap()) * &z2_3;

        //z2^4 * (v13 - v4) / (z3 - y_new)
        let term_c2_num = &v[12] - &v[3];
        let term_c2_den = z3 - y_new;
        let term_c2 = (term_c2_num * term_c2_den.inv().unwrap()) * &z2_4;

        let h_term = term_p + term_q + term_r + term_c1 + term_c2;

        let mut expected_val = v[9].clone();
        expected_val += &v[10] * z4;
        expected_val += &v[11] * &z4_2;
        expected_val += &v[12] * &z4_3;
        expected_val += h_term * &z4_4;

        let p_stroke = p_comb + u_point.clone() * expected_val;

        let mut q_ipa = p_stroke;
        for ((l, r), u) in proof
            .zk_proof
            .l_vec
            .iter()
            .zip(&proof.zk_proof.r_vec)
            .zip(u_j_vec)
        {
            let u_sq = u * u;
            let u_inv_sq = u_sq.inv().unwrap();
            q_ipa = q_ipa + l.clone() * u_sq + r.clone() * u_inv_sq;
        }

        let s_vec = generate_s_vector(u_j_vec, self.d);
        let g_folded = CurvePoint::msm(&self.sigma.g, &s_vec);

        let mut b_val = Field::<M>::ZERO;
        let mut z3_pow = Field::<M>::one();
        for s in &s_vec {
            b_val += s * &z3_pow;
            z3_pow = z3_pow * z3;
        }

        let lhs = q_ipa * c.clone() + proof.zk_proof.r_zk.clone();

        let term_base = g_folded + u_point.clone() * b_val;
        let rhs = term_base * proof.zk_proof.z1_zk.clone()
            + self.sigma.h.clone() * proof.zk_proof.z2_zk.clone();

        lhs == rhs
    }
    pub fn reconstruct_challenges(
        &self,
        proof: &VerifierProof<M>,
        l_vec: &[CurvePoint<M>],
        r_vec: &[CurvePoint<M>],
        r_zk: &CurvePoint<M>,
    ) -> (
        Field<M>,
        Field<M>,
        Field<M>,
        Field<M>,
        Field<M>,
        Field<M>,
        Field<M>,
        Vec<Field<M>>,
        Field<M>,
    ) {
        let y: Field<M> = sample_challenge(&[Data::Point(&proof.r_comm)]);

        let x: Field<M> = sample_challenge(&[
            Data::Point(&proof.t_lo_comm),
            Data::Point(&proof.t_hi_comm),
            Data::Point(&proof.s_comm),
        ]);

        let y_new: Field<M> = sample_challenge(&[Data::Point(&proof.c_comm)]);

        let (z1, z2) = sample_two_challenges(&[
            Data::Point(&proof.s_new_comm),
            Data::Field(&proof.values[0]),
            Data::Field(&proof.values[1]),
            Data::Field(&proof.values[2]),
            Data::Field(&proof.values[3]),
            Data::Field(&proof.values[4]),
            Data::Field(&proof.values[5]),
            Data::Field(&proof.values[6]),
            Data::Field(&proof.values[7]),
            Data::Field(&proof.values[8]),
        ]);

        let z3: Field<M> = sample_challenge(&[Data::Point(&proof.h_comm)]);

        let z4: Field<M> = sample_challenge(&[
            Data::Field(&proof.values[9]),
            Data::Field(&proof.values[10]),
            Data::Field(&proof.values[11]),
            Data::Field(&proof.values[12]),
        ]);

        let mut u_vec = Vec::new();
        for (l, r) in l_vec.iter().zip(r_vec.iter()) {
            let u = sample_challenge(&[Data::Point(l), Data::Point(r)]);
            u_vec.push(u);
        }

        let c: Field<M> = sample_challenge(&[Data::Point(r_zk)]);

        (y, x, y_new, z1, z2, z3, z4, u_vec, c)
    }
    pub fn verify_proof(
        &self,
        proof: &VerifierProof<M>,
        g_old: &CurvePoint<M>,
        u_old: &[Field<M>],
        s_old_comm: &CurvePoint<M>,
        y_old: &Field<M>,
        k_comm: &CurvePoint<M>,
        u_point: &CurvePoint<M>,
    ) -> bool {
        let l_vec = &proof.zk_proof.l_vec;
        let r_vec = &proof.zk_proof.r_vec;
        let r_zk = &proof.zk_proof.r_zk;

        let (y, x, y_new, z1, z2, z3, z4, u_j_vec, c) =
            self.reconstruct_challenges(proof, l_vec, r_vec, r_zk);

        let g_check = self.check_g(&proof.values[6], &x, u_old);

        let poly_check = self.check_polynomial_identity(&proof.values, &x, &y);

        let commitment_check = self.verify_commitment_check(
            proof, s_old_comm, k_comm, u_point, &x, &y, y_old, &y_new, &z1, &z2, &z3, &z4, &c,
            &u_j_vec,
        );

        poly_check && g_check && commitment_check
    }
}

//something go wrong here ^()
fn generate_s_vector<M: Modulus>(u_vec: &[Field<M>], n: usize) -> Vec<Field<M>> {
    let k = u_vec.len();

    let mut s = vec![Field::<M>::one(); n];

    for (i, u) in u_vec.iter().enumerate() {
        let bit_mask = 1 << (k - 1 - i);
        let u_inv = u.inv().unwrap();

        for j in 0..n {
            if (j & bit_mask) != 0 {
                s[j] = &s[j] * u;
            } else {
                s[j] = &s[j] * &u_inv;
            }
        }
    }
    s
}
