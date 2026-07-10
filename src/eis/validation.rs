use crate::eis::data::{EisPoint, ImaginaryConvention};
use anyhow::{Result, bail};

pub fn validate_point(point: &EisPoint) -> Result<()> {
    if !point.frequency_hz.is_finite() {
        bail!("frequency is not finite");
    }
    if point.frequency_hz <= 0.0 {
        bail!("frequency must be strictly positive");
    }
    if !point.z_real_ohm.is_finite() {
        bail!("real impedance is not finite");
    }
    if !point.z_imag_ohm.is_finite() {
        bail!("imaginary impedance is not finite");
    }
    Ok(())
}

pub fn detect_imaginary_convention(points: &[EisPoint]) -> ImaginaryConvention {
    let positive = points.iter().filter(|point| point.z_imag_ohm > 0.0).count();
    let negative = points.iter().filter(|point| point.z_imag_ohm < 0.0).count();
    match (positive, negative) {
        (0, 0) => ImaginaryConvention::ZeroOrUnknown,
        (p, 0) if p > 0 => ImaginaryConvention::PositiveCapacitive,
        (0, n) if n > 0 => ImaginaryConvention::NegativeCapacitive,
        (p, n) if p > n => ImaginaryConvention::PositiveCapacitive,
        (p, n) if n > p => ImaginaryConvention::NegativeCapacitive,
        _ => ImaginaryConvention::Mixed,
    }
}
