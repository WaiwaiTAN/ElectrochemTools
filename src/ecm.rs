use anyhow::{Result, bail};
use num_complex::Complex;
use serde::Serialize;
use std::f64::consts::PI;

pub trait EquivalentCircuitModel {
    fn impedance(&self, freq_hz: &[f64]) -> Vec<Complex<f64>>;
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct RqrParams {
    pub rs: f64,
    pub rct: f64,
    pub q: f64,
    pub n: f64,
}

impl RqrParams {
    pub fn validate(&self) -> Result<()> {
        if self.rs <= 0.0 || !self.rs.is_finite() {
            bail!("Rs must be finite and > 0");
        }
        if self.rct <= 0.0 || !self.rct.is_finite() {
            bail!("Rct must be finite and > 0");
        }
        if self.q <= 0.0 || !self.q.is_finite() {
            bail!("Q must be finite and > 0");
        }
        if self.n <= 0.0 || self.n > 1.0 || !self.n.is_finite() {
            bail!("n must be finite and in (0, 1]");
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RqrModel {
    pub params: RqrParams,
}

impl EquivalentCircuitModel for RqrModel {
    fn impedance(&self, freq_hz: &[f64]) -> Vec<Complex<f64>> {
        freq_hz
            .iter()
            .map(|&freq| rqr_impedance_at(freq, self.params))
            .collect()
    }
}

pub fn rqr_impedance_at(freq_hz: f64, params: RqrParams) -> Complex<f64> {
    let omega = 2.0 * PI * freq_hz;
    let mag = omega.powf(params.n);
    let phase = params.n * PI / 2.0;
    let jw_n = Complex::new(mag * phase.cos(), mag * phase.sin());
    let y_cpe = jw_n * params.q;
    let y_rct = Complex::new(1.0 / params.rct, 0.0);
    Complex::new(params.rs, 0.0) + Complex::new(1.0, 0.0) / (y_rct + y_cpe)
}
