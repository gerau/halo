use std::str::FromStr;

use sha3::digest::typenum::Quot;

use crate::{
    elliptic_curve::CurvePoint,
    field::Field,
    modulus::{Modulus, OrderP},
    poly::Polynomial,
    proof::{Prover, ProverParams, Verifier, VerifierProof},
    utils::{ProoverHelper, compute_g_poly},
};

mod elliptic_curve;
mod field;
mod modulus;
mod poly;
mod proof;
mod utils;

type FieldP = Field<OrderP>;

fn main() {
    let a = Vec::from([854, 114, 212, 987].map(FieldP::from));
    let b = Vec::from([345, 678, 123, 888].map(FieldP::from));
    let c: Vec<FieldP> = a.iter().zip(&b).map(|(a, b)| a * b).collect();

    let u_mat = Vec::from([
        Vec::from([123, 456, 789, 101].map(FieldP::from)),
        Vec::from([202, 303, 404, 505].map(FieldP::from)),
        Vec::from([606, 707, 808, 909].map(FieldP::from)),
    ]);
    let v_mat = Vec::from([
        Vec::from([111, 222, 333, 444].map(FieldP::from)),
        Vec::from([555, 666, 777, 888].map(FieldP::from)),
        Vec::from([999, 100, 200, 300].map(FieldP::from)),
    ]);
    let w_mat = Vec::from([
        Vec::from([150, 250, 350, 450].map(FieldP::from)),
        Vec::from([550, 650, 750, 850].map(FieldP::from)),
        Vec::from([950, 105, 205, 305].map(FieldP::from)),
    ]);

    let u_old = vec![
        FieldP::from_str("64754676296675759630").unwrap(),
        FieldP::from_str("73765879184746879147").unwrap(),
        FieldP::from_str("97773392004576823364").unwrap(),
        FieldP::from_str("60165719114422929560").unwrap(),
    ];

    let k_log = 4;
    let d = 1 << k_log;
    let n = d / 4;

    let sigma = ProverParams::<OrderP>::new(k_log);

    let mut prover = Prover::setup(sigma.clone(), k_log);
    let verifier = Verifier::new(sigma.clone(), d);

    let y_old = FieldP::from(123456789);

    let init_helper = ProoverHelper::new(
        n,
        3,
        vec![],
        vec![],
        vec![],
        u_mat.clone(),
        v_mat.clone(),
        w_mat.clone(),
    );

    let s_old_poly = init_helper.evaluate_s_at_y(&y_old).shift_poly(n as i64);
    let s_old_comm = prover.commit(&s_old_poly, &FieldP::ZERO);

    let g_old_poly = compute_g_poly(&u_old);
    let g_old_comm = prover.commit(&g_old_poly, &FieldP::ZERO);

    let helper = ProoverHelper::new(
        n,
        3,
        a.clone(),
        b.clone(),
        c.clone(),
        u_mat.clone(),
        v_mat.clone(),
        w_mat.clone(),
    );

    let k_poly = Polynomial {
        coefficients: helper.k.clone(),
        offset: 1,
    };
    let k_comm = prover.commit(&k_poly, &FieldP::ZERO);

    let u_point = CurvePoint::get_random(&mut rand::thread_rng());

    prover.generate_commitments(&helper, &y_old, &s_old_poly, &u_old);

    let zk_proof = prover.zero_knowledge_opening(&u_point);

    let proof_package = VerifierProof {
        r_comm: prover.state.r.clone(),
        s_comm: prover.state.s.clone(),
        s_new_comm: prover.state.s_new.clone(),
        t_lo_comm: prover.state.t_low.clone(),
        t_hi_comm: prover.state.t_high.clone(),
        c_comm: prover.state.c.clone(),

        h_comm: prover.state.h,

        values: prover.state.v.to_vec(),

        zk_proof: zk_proof.clone(),
    };

    let is_valid = verifier.verify_proof(
        &proof_package,
        &g_old_comm,
        &u_old,
        &s_old_comm,
        &y_old,
        &k_comm,
        &u_point,
    );

    println!("--------------------------------------------------");
    if is_valid {
        println!("VERIFICATION SUCCESSFUL");
    } else {
        println!("VERIFICATION FAILED");
    }
    println!("--------------------------------------------------");
}
