use std::default;

use rand::Rng;
use rand::prelude::*;

use crate::poly::Polynomial;
use crate::utils::ProoverHelper;
use crate::{
    elliptic_curve::CurvePoint,
    field::{self, Field},
    modulus::Modulus,
};

pub fn sample_challenge<M: Modulus, R: Rng + Sized>(rng: &mut R) -> Field<M> {
    Field::get_random(rng)
}

type Commit<M> = CurvePoint<M>;

#[derive(Debug, Default)]
struct ProofParams<M: Modulus> {
    pub g: Vec<CurvePoint<M>>,
    pub h: CurvePoint<M>,
}

#[derive(Debug, Default)]
struct Proof<M: Modulus> {
    d: usize,
    n: usize,

    state: [Field<M>; 13],
    t_low: Commit<M>,
    t_high: Commit<M>,
    s: Commit<M>,
    s_new: Commit<M>,
    c: Commit<M>,

    sigma: ProofParams<M>,
}

impl<M: Modulus + Default> Proof<M> {
    pub fn setup(d: usize) -> Self {
        let mut rng = rand::thread_rng();
        let sigma = ProofParams {
            g: (0..d).map(|_| CurvePoint::get_random(&mut rng)).collect(),
            h: CurvePoint::get_random(&mut rng),
        };

        let n = d / 4;

        Self {
            n,
            d,
            sigma,
            ..Default::default()
        }
    }

    pub fn commit(&self, poly: Polynomial<M>, delta: Field<M>) -> Commit<M> {
        let lhs = self
            .sigma
            .g
            .iter()
            .zip(poly.coefficients)
            .map(|(point, scalar)| point.clone() * scalar)
            .fold(CurvePoint::<M>::new_infinity(), |acc, val| acc + val);

        lhs + self.sigma.h.clone() * delta
    }

    pub fn generate_commitments(&self, proover: ProoverHelper<M>) {
        let mut rng = rand::thread_rng();

        let delta_r = Field::<M>::get_random(&mut rng);

        let r = self.commit(
            proover
                .r_evaluate_at_x(&Field::<M>::from(1))
                .shift_poly((3 * proover.n - 1) as i64),
            delta_r,
        );

        let y: Field<M> = sample_challenge(&mut rng);
    }
}
