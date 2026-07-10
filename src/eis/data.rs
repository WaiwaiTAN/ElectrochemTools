use crate::eis::format::EisFormat;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ImaginaryConvention {
    NegativeCapacitive,
    PositiveCapacitive,
    Mixed,
    ZeroOrUnknown,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EisPoint {
    pub frequency_hz: f64,
    pub z_real_ohm: f64,
    pub z_imag_ohm: f64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PreprocessRecord {
    pub operation: String,
    pub detail: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EisMetadata {
    pub source_path: Option<PathBuf>,
    pub source_format: EisFormat,
    pub imaginary_convention: ImaginaryConvention,
    pub original_point_count: usize,
    pub preprocessing: Vec<PreprocessRecord>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EisSpectrum {
    pub points: Vec<EisPoint>,
    pub metadata: EisMetadata,
}
