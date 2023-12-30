// +--------------------------+
// |                          |
// |    Copyright (c) 2023    |
// |       Keith Cullen       |
// |                          |
// +--------------------------+

use super::*;
use std::f64::consts::PI;

pub struct Ctx {
    cfft_ctx: cfft::Ctx,
    size: usize,
    size_sin_tbl: usize,
    sin_tbl: Vec<f64>,
}

impl Ctx {
    pub fn new(size: usize) -> Ctx {
        let mut ctx = Ctx {
            cfft_ctx: cfft::Ctx::new(size >> 1),
            size: size,
            size_sin_tbl: 0,
            sin_tbl: Vec::new(),
        };
        ctx.gen_sin_tbl();
        ctx
    }

    fn gen_sin_tbl(&mut self) {
        self.size_sin_tbl = (self.size >> 2) + 1;
        self.sin_tbl.resize(self.size_sin_tbl, 0.0);
        let w = (2.0 * PI) / self.size as f64;
        for i in 0..self.size_sin_tbl {
            let x: f64 = w * i as f64;
            self.sin_tbl[i] = x.sin();
        }
    }

    pub fn fwd(&self, x: &mut [f64]) {
        self.run(x, -1.0);
    }

    pub fn bwd(&self, x: &mut [f64]) {
        self.run(x, 1.0);
    }

    pub fn run(&self, x: &mut [f64], isign: f64) {
        let c2 = isign * 0.5;
        let c1 = 0.5;
        let mut cos_index = self.size_sin_tbl - 1;
        let mut sin_index = 0;
    
        if isign < 0.0 {
            self.cfft_ctx.run(x, -1.0);
        }
        for k in 1..((self.size >> 2) + 1) {
            cos_index -= 1;
            sin_index += 1;
            let wr =         self.sin_tbl[cos_index];
            let wi = isign * self.sin_tbl[sin_index];

            let k0 = k << 1;
            let k1 = k0 + 1;
            let k2 = self.size - k0;
            let k3 = k2 + 1;
    
            let t1 = c1 * (x[k0] + x[k2]);
            let t2 = c2 * (x[k0] - x[k2]);
            let t3 = c2 * (x[k1] + x[k3]);
            let t4 = c1 * (x[k1] - x[k3]);
    
            let wrt2 = wr * t2;
            let wrt3 = wr * t3;
            let wit2 = wi * t2;
            let wit3 = wi * t3;
    
            x[k0] =  t1 - wrt3 - wit2;
            x[k1] =  t4 + wrt2 - wit3;
            x[k2] =  t1 + wrt3 + wit2;
            x[k3] = -t4 + wrt2 - wit3;
        }
        if isign < 0.0
        {
            let xt = x[0];
            x[0] = x[0] + x[1];
            x[1] = xt   - x[1];
        }
        else
        {
            let xt = x[0];
            x[0] = 0.5 * (xt + x[1]);
            x[1] = 0.5 * (xt - x[1]);
            self.cfft_ctx.run(x, 1.0);
        }
    }

    pub fn unpack(&self, y: &mut [f64], x: &[f64])
    {
        for i in 0..self.size {
            y[i] = x[i];
        }
        y[self.size] = y[1];
        y[self.size + 1] = 0.0;
        y[1] = 0.0;
        for i in 1..(self.size >> 1) {
            let j = self.size - i;
            y[j << 1] = x[i << 1];
            y[(j << 1) + 1] = -x[(i << 1) + 1];
        }
    }

    pub fn pack(&self, y: &mut [f64], x: &[f64])
    {
        for i in 0..self.size {
            y[i] = x[i];
        }
        y[1] = x[self.size];
    }
}

#[cfg(test)]
mod tests;
