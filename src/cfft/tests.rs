// +--------------------------+
// |                          |
// |    Copyright (c) 2023    |
// |       Keith Cullen       |
// |                          |
// +--------------------------+

use super::*;
use std::f64::consts::PI;

fn compute_snr(s1: &[f64], s2: &[f64]) -> f64
{
    let mut sss: f64 = 0.0;
    let mut ses: f64 = 0.0;
    let mut err: f64;

    for i in 0..s1.len() {
        sss += s1[i] * s1[i];
        err  = s1[i] - s2[i];
        ses += err * err;
    }
    10.0 * (sss / ses).log10()
}

#[test]
fn forward() {
    const THRESH_SNR: f64 = 300.0;
    const SIZE: usize = 8;
    const PHI: f64 = PI / 3.5;
    const W: f64 = PI / 177.0;
    let mut x1: [f64; 2 * SIZE] = [0.0; 2 * SIZE];
    let x2: [f64; 2 * SIZE] = [
        4.586201271849956740,
        0.000000000000000000,
        0.056022113303054955,
        -0.141206857527545360,
        0.057745057812326424,
        -0.058467550303321847,
        0.058040536646275523,
        -0.024216472082145848,
        0.058101727496598254,
        0.000000000000000000,
        0.058040536646275509,
        0.024216472082145848,
        0.057745057812326410,
        0.058467550303321847,
        0.056022113303054941,
        0.141206857527545360];
    let ctx = Ctx::new(SIZE);
    for i in 0..SIZE {
        x1[2 * i] = (W * i as f64 + PHI).cos();
    }
    ctx.fwd(&mut x1[..]);
    let snr = compute_snr(&x1[..], &x2[..]);
    if snr < THRESH_SNR {
        panic!("SNR below threshold");
    }
}

#[test]
fn forward_and_backward() {
    const THRESH_SNR: f64 = 300.0;
    const SIZE: usize = 8;
    const PHI: f64 = PI / 3.5;
    const W: f64 = PI / 177.0;
    let mut x1: [f64; 2 * SIZE] = [0.0; 2 * SIZE];
    let mut x2: [f64; 2 * SIZE] = [0.0; 2 * SIZE];
    let ctx = Ctx::new(SIZE);
    for i in 0..SIZE {
        x1[2 * i] = (W * i as f64 + PHI).cos();
        x2[2 * i] = x1[2 * i];
    }
    ctx.fwd(&mut x1[..]);
    ctx.bwd(&mut x1[..]);
    for i in 0..SIZE*2 {
        x1[i] /= SIZE as f64;
    }
    let snr = compute_snr(&x1[..], &x2[..]);
    if snr < THRESH_SNR {
        panic!("SNR below threshold");
    }
}
