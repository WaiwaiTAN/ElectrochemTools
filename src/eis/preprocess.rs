use crate::eis::data::{EisSpectrum, ImaginaryConvention, PreprocessRecord};
use crate::eis::validation::detect_imaginary_convention;
use clap::ValueEnum;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, ValueEnum)]
#[serde(rename_all = "kebab-case")]
pub enum ImagSignPolicy {
    Preserve,
    Flip,
    NegativeCapacitive,
    PositiveCapacitive,
}

pub fn apply_imag_sign(spectrum: &mut EisSpectrum, policy: ImagSignPolicy) {
    let before = detect_imaginary_convention(&spectrum.points);
    let flip = match policy {
        ImagSignPolicy::Preserve => false,
        ImagSignPolicy::Flip => true,
        ImagSignPolicy::NegativeCapacitive => before == ImaginaryConvention::PositiveCapacitive,
        ImagSignPolicy::PositiveCapacitive => before == ImaginaryConvention::NegativeCapacitive,
    };
    if flip {
        for point in &mut spectrum.points {
            point.z_imag_ohm = -point.z_imag_ohm;
        }
    }
    spectrum.metadata.imaginary_convention = detect_imaginary_convention(&spectrum.points);
    spectrum.metadata.preprocessing.push(PreprocessRecord {
        operation: "imaginary-sign".to_string(),
        detail: format!("policy={policy:?}; flipped={flip}"),
    });
}
