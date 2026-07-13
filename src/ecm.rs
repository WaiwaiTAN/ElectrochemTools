use anyhow::{Result, bail};
use num_complex::Complex;
use serde::Serialize;
use std::f64::consts::PI;
use std::str::FromStr;

pub trait EquivalentCircuitModel {
    fn impedance(&self, freq_hz: &[f64]) -> Vec<Complex<f64>>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum BranchKind {
    Rc,
    Rq,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct EcmModelSpec {
    pub branches: Vec<BranchKind>,
    pub warburg: bool,
}

impl EcmModelSpec {
    pub const SUPPORTED: &'static str =
        "R_QR, R_CR, R_QR_CR, R_CR_CR, R_QR_QR, and each with suffix _W";

    pub fn canonical_name(&self) -> String {
        let mut name = String::from("R");
        for branch in &self.branches {
            name.push('_');
            name.push_str(match branch {
                BranchKind::Rc => "CR",
                BranchKind::Rq => "QR",
            });
        }
        if self.warburg {
            name.push_str("_W");
        }
        name
    }

    pub fn parameter_labels(&self) -> Vec<String> {
        let mut labels = vec!["Rs".to_string()];
        for (index, branch) in self.branches.iter().enumerate() {
            let number = index + 1;
            labels.push(format!("R{number}"));
            match branch {
                BranchKind::Rc => labels.push(format!("C{number}")),
                BranchKind::Rq => {
                    labels.push(format!("Q{number}"));
                    labels.push(format!("n{number}"));
                }
            }
        }
        if self.warburg {
            labels.push("sigma_w".to_string());
        }
        labels
    }
}

impl FromStr for EcmModelSpec {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        // Parentheses are accepted so the literature-style R_(QR)_(CR) notation
        // and the compact CLI spelling R_QR_CR select the same topology.
        let mut normalized = value
            .trim()
            .to_ascii_uppercase()
            .replace(['(', ')', ' '], "")
            .replace('-', "_");
        while normalized.contains("__") {
            normalized = normalized.replace("__", "_");
        }
        let warburg = normalized.ends_with("_W");
        if warburg {
            normalized.truncate(normalized.len() - 2);
        }
        let branches = match normalized.as_str() {
            "R_QR" => vec![BranchKind::Rq],
            "R_CR" => vec![BranchKind::Rc],
            "R_QR_CR" => vec![BranchKind::Rq, BranchKind::Rc],
            "R_CR_CR" => vec![BranchKind::Rc, BranchKind::Rc],
            "R_QR_QR" => vec![BranchKind::Rq, BranchKind::Rq],
            _ => bail!(
                "unsupported model '{}'; supported models: {}",
                value,
                Self::SUPPORTED
            ),
        };
        Ok(Self { branches, warburg })
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
#[serde(tag = "element", rename_all = "lowercase")]
pub enum ReactiveElementParams {
    Capacitor { c: f64 },
    Cpe { q: f64, n: f64 },
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct ParallelBranchParams {
    pub r: f64,
    #[serde(flatten)]
    pub reactive: ReactiveElementParams,
}

#[derive(Debug, Clone, Serialize)]
pub struct EcmParams {
    pub rs: f64,
    pub branches: Vec<ParallelBranchParams>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub warburg_sigma: Option<f64>,
}

impl EcmParams {
    pub fn validate_for(&self, model: &EcmModelSpec) -> Result<()> {
        validate_positive("Rs", self.rs)?;
        if self.branches.len() != model.branches.len() {
            bail!(
                "model {} requires {} relaxation branch(es), but {} were provided",
                model.canonical_name(),
                model.branches.len(),
                self.branches.len()
            );
        }
        for (index, (actual, expected)) in self.branches.iter().zip(&model.branches).enumerate() {
            let number = index + 1;
            validate_positive(&format!("R{number}"), actual.r)?;
            match (actual.reactive, expected) {
                (ReactiveElementParams::Capacitor { c }, BranchKind::Rc) => {
                    validate_positive(&format!("C{number}"), c)?;
                }
                (ReactiveElementParams::Cpe { q, n }, BranchKind::Rq) => {
                    validate_positive(&format!("Q{number}"), q)?;
                    if !n.is_finite() || n <= 0.0 || n > 1.0 {
                        bail!("n{number} must be finite and in (0, 1]");
                    }
                }
                _ => bail!("branch {number} parameters do not match the selected model"),
            }
        }
        match (model.warburg, self.warburg_sigma) {
            (true, Some(sigma)) => validate_positive("sigma_w", sigma)?,
            (true, None) => bail!("model {} requires sigma_w", model.canonical_name()),
            (false, Some(_)) => bail!("sigma_w was supplied for a model without Warburg diffusion"),
            (false, None) => {}
        }
        Ok(())
    }

    pub fn values(&self) -> Vec<f64> {
        let mut values = vec![self.rs];
        for branch in &self.branches {
            values.push(branch.r);
            match branch.reactive {
                ReactiveElementParams::Capacitor { c } => values.push(c),
                ReactiveElementParams::Cpe { q, n } => {
                    values.push(q);
                    values.push(n);
                }
            }
        }
        if let Some(sigma) = self.warburg_sigma {
            values.push(sigma);
        }
        values
    }
}

#[derive(Debug, Clone)]
pub struct EcmModel {
    pub spec: EcmModelSpec,
    pub params: EcmParams,
}

impl EquivalentCircuitModel for EcmModel {
    fn impedance(&self, freq_hz: &[f64]) -> Vec<Complex<f64>> {
        freq_hz
            .iter()
            .map(|&freq| ecm_impedance_at(freq, &self.spec, &self.params))
            .collect()
    }
}

pub fn ecm_impedance_at(freq_hz: f64, model: &EcmModelSpec, params: &EcmParams) -> Complex<f64> {
    let omega = 2.0 * PI * freq_hz;
    let mut impedance = Complex::new(params.rs, 0.0);
    for (kind, branch) in model.branches.iter().zip(&params.branches) {
        let y_reactive = match (kind, branch.reactive) {
            (BranchKind::Rc, ReactiveElementParams::Capacitor { c }) => {
                Complex::new(0.0, omega * c)
            }
            (BranchKind::Rq, ReactiveElementParams::Cpe { q, n }) => {
                let magnitude = omega.powf(n);
                let phase = n * PI / 2.0;
                Complex::new(magnitude * phase.cos(), magnitude * phase.sin()) * q
            }
            _ => unreachable!("validated ECM parameters must match their branch kinds"),
        };
        let y_resistor = Complex::new(1.0 / branch.r, 0.0);
        impedance += Complex::new(1.0, 0.0) / (y_resistor + y_reactive);
    }
    if model.warburg {
        let sigma = params
            .warburg_sigma
            .expect("validated Warburg model must have sigma_w");
        let scale = sigma / omega.sqrt();
        impedance += Complex::new(scale, -scale);
    }
    impedance
}

fn validate_positive(label: &str, value: f64) -> Result<()> {
    if !value.is_finite() || value <= 0.0 {
        bail!("{label} must be finite and > 0");
    }
    Ok(())
}

// Backward-compatible R(QR) API used by existing library clients.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct RqrParams {
    pub rs: f64,
    pub rct: f64,
    pub q: f64,
    pub n: f64,
}

impl RqrParams {
    pub fn validate(&self) -> Result<()> {
        self.to_ecm().validate_for(&EcmModelSpec {
            branches: vec![BranchKind::Rq],
            warburg: false,
        })
    }

    fn to_ecm(self) -> EcmParams {
        EcmParams {
            rs: self.rs,
            branches: vec![ParallelBranchParams {
                r: self.rct,
                reactive: ReactiveElementParams::Cpe {
                    q: self.q,
                    n: self.n,
                },
            }],
            warburg_sigma: None,
        }
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
    ecm_impedance_at(
        freq_hz,
        &EcmModelSpec {
            branches: vec![BranchKind::Rq],
            warburg: false,
        },
        &params.to_ecm(),
    )
}
