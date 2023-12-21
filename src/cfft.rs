// +--------------------------+
// |                          |
// |    Copyright (c) 2023    |
// |       Keith Cullen       |
// |                          |
// +--------------------------+

use std::f64::consts::PI;

const CFFT_MIN_SIZE: usize = 4;
const CFFT_MAX_SIZE: usize = 1 << 16;

pub struct Ctx {
    size: usize,
    bit_rev_tbl: Vec<usize>,
    sin_tbl: Vec<f64>,
}

impl Ctx {
    pub fn new(size: usize, ) -> Ctx {
        validate_size(size);
        let mut ctx = Ctx {
            size: size,
            bit_rev_tbl: Vec::new(),
            sin_tbl: Vec::new(),
        };
        ctx.gen_bit_rev_tbl();
        ctx.gen_sin_tbl();
        ctx
    }

    fn gen_bit_rev_tbl(&mut self)
    {
        self.bit_rev_tbl.resize(self.size, 0);
        let mut i: usize = 1;
        let mut j: usize = 1;
        for k in 0..self.size {
            if j > i {
                self.bit_rev_tbl[k] = j - 1;
            } else {
                self.bit_rev_tbl[k] = 0;
            }
            let mut m = self.size;
            while m >= 2 && j > m {
                j -= m;
                m >>= 1;
            }
            j += m;
            i += 2;
        }
    }

    fn gen_sin_tbl(&mut self) {
        self.sin_tbl.resize(3 * (self.size >> 2), 0.0);
        let size_sin_tbl = 3 * (self.size >> 2);
        for i in 0..size_sin_tbl {
            let x: f64 = ((2.0 * PI) / self.size as f64) * i as f64;
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
        let n: usize = self.size << 1;
    
        for i in (0..n).step_by(2) {
            let bri = self.bit_rev_tbl[i >> 1];
            if bri != 0 {
                let xt = x[bri];
                x[bri] = x[i];
                x[i] = xt;
    
                let xt = x[bri + 1];
                x[bri + 1] = x[i + 1];
                x[i + 1] = xt;
            }
        }
        /* stage 1 */
        for i in (0..n).step_by(4) {
            let j = i + 2;
            let tr = x[j];
            let ti = x[j + 1];
            x[j] = x[i] - tr;
            x[j + 1] = x[i + 1] - ti;
            x[i] += tr;
            x[i + 1] += ti;
        }
        /* stage 2 */
        for i in (0..n).step_by(8) {
            let j = i + 4;
            let tr = x[j];
            let ti = x[j + 1];
            x[j] = x[i] - tr;
            x[j + 1] = x[i + 1] - ti;
            x[i] += tr;
            x[i + 1] += ti;
        }
        for i in (2..n).step_by(8) {
            let j = i + 4;
            let tr = -isign * x[j + 1];
            let ti =  isign * x[j];
            x[j] = x[i] - tr;
            x[j + 1] = x[i + 1] - ti;
            x[i] += tr;
            x[i + 1] += ti;
        }
        /* remaining stages */
        let cos_offset = self.size >> 2;   /* PI / 2 */
        let mut sin_step = self.size >> 3; /* start with stage 3 */
        let mut mmax: usize = 8;
        while mmax <= self.size {
            for m in (0..mmax).step_by(2) {
                let sin_index = (m >> 1) * sin_step;
                let wr =         self.sin_tbl[sin_index + cos_offset];
                let wi = isign * self.sin_tbl[sin_index];
                for i in (m..n).step_by(mmax << 1) {
                    let j = i + mmax;
                    let tr = wr * x[j] - wi * x[j + 1];
                    let ti = wr * x[j + 1] + wi * x[j];
                    x[j] = x[i] - tr;
                    x[j + 1] = x[i + 1] - ti;
                    x[i] += tr;
                    x[i + 1] += ti;
                }
            }
            sin_step >>= 1;
            mmax <<= 1;
        }
    }
}

fn validate_size(size: usize) {
    let mut i: usize = CFFT_MIN_SIZE;
    while i <= CFFT_MAX_SIZE {
        if i == size {
            return;
        }
        i *= 2;
    }
    panic!("FFT size must be an integer power of 2")
}
