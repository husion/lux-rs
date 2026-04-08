use crate::error::{LuxError, LuxResult};
use crate::spectrum::{getwld, Spectrum};

struct SpectralMismatchContext {
    spacing: Vec<f64>,
    target: Vec<f64>,
    relative_detectors: Vec<Vec<f64>>,
}

pub fn spectral_mismatch_f1prime(
    detector: &Spectrum,
    calibration_illuminant: &Spectrum,
    target_responsivity: &Spectrum,
) -> LuxResult<f64> {
    Ok(spectral_mismatch_f1primes(
        &Spectrum::new(
            detector.wavelengths().to_vec(),
            vec![detector.values().to_vec()],
        )?,
        calibration_illuminant,
        target_responsivity,
    )?[0])
}

pub fn spectral_mismatch_f1primes(
    detectors: &Spectrum,
    calibration_illuminant: &Spectrum,
    target_responsivity: &Spectrum,
) -> LuxResult<Vec<f64>> {
    let context = build_context(
        detectors,
        calibration_illuminant,
        target_responsivity,
        detectors.wavelengths(),
    )?;
    let target_integral = weighted_sum(&context.target, &context.spacing);
    if target_integral == 0.0 {
        return Err(LuxError::InvalidInput(
            "target responsivity integral must be non-zero",
        ));
    }

    Ok(context
        .relative_detectors
        .iter()
        .map(|relative_detector| {
            weighted_sum_abs_difference(relative_detector, &context.target, &context.spacing)
                / target_integral
        })
        .collect())
}

pub fn spectral_mismatch_correction_factor(
    measured_source: &Spectrum,
    detector: &Spectrum,
    calibration_illuminant: &Spectrum,
    target_responsivity: &Spectrum,
) -> LuxResult<f64> {
    Ok(spectral_mismatch_correction_factors(
        &Spectrum::new(
            measured_source.wavelengths().to_vec(),
            vec![measured_source.values().to_vec()],
        )?,
        &Spectrum::new(
            detector.wavelengths().to_vec(),
            vec![detector.values().to_vec()],
        )?,
        calibration_illuminant,
        target_responsivity,
    )?[0][0])
}

pub fn spectral_mismatch_correction_factors(
    measured_sources: &Spectrum,
    detectors: &Spectrum,
    calibration_illuminant: &Spectrum,
    target_responsivity: &Spectrum,
) -> LuxResult<Vec<Vec<f64>>> {
    let context = build_context(
        detectors,
        calibration_illuminant,
        target_responsivity,
        measured_sources.wavelengths(),
    )?;

    measured_sources
        .spectra()
        .iter()
        .map(|source| {
            let numerator = weighted_dot(source, &context.target, &context.spacing);
            context
                .relative_detectors
                .iter()
                .map(|relative_detector| {
                    let denominator = weighted_dot(source, relative_detector, &context.spacing);
                    if denominator == 0.0 {
                        Err(LuxError::InvalidInput(
                            "spectral mismatch denominator must be non-zero",
                        ))
                    } else {
                        Ok(numerator / denominator)
                    }
                })
                .collect()
        })
        .collect()
}

fn build_context(
    detectors: &Spectrum,
    calibration_illuminant: &Spectrum,
    target_responsivity: &Spectrum,
    wavelengths: &[f64],
) -> LuxResult<SpectralMismatchContext> {
    let spacing = getwld(wavelengths)?;
    let target = target_responsivity.interpolate_linear(wavelengths)?;
    let calibration = calibration_illuminant.interpolate_linear(wavelengths)?;
    let detectors = detectors.cie_interp_linear(wavelengths, true)?;

    let calibration_weights: Vec<f64> = calibration
        .values()
        .iter()
        .zip(spacing.iter())
        .map(|(value, dl)| value * dl)
        .collect();
    let target_weighted_sum = dot(target.values(), &calibration_weights);
    if target_weighted_sum == 0.0 {
        return Err(LuxError::InvalidInput(
            "target responsivity under calibration illuminant must be non-zero",
        ));
    }

    let relative_detectors = detectors
        .spectra()
        .iter()
        .map(|detector| {
            let detector_weighted_sum = dot(detector, &calibration_weights);
            if detector_weighted_sum == 0.0 {
                Err(LuxError::InvalidInput(
                    "detector responsivity under calibration illuminant must be non-zero",
                ))
            } else {
                let scale = target_weighted_sum / detector_weighted_sum;
                Ok(detector.iter().map(|value| value * scale).collect())
            }
        })
        .collect::<LuxResult<Vec<Vec<f64>>>>()?;

    Ok(SpectralMismatchContext {
        spacing,
        target: target.values().to_vec(),
        relative_detectors,
    })
}

fn weighted_dot(lhs: &[f64], rhs: &[f64], weights: &[f64]) -> f64 {
    lhs.iter()
        .zip(rhs.iter())
        .zip(weights.iter())
        .map(|((lhs, rhs), weight)| lhs * rhs * weight)
        .sum()
}

fn weighted_sum(values: &[f64], weights: &[f64]) -> f64 {
    values
        .iter()
        .zip(weights.iter())
        .map(|(value, weight)| value * weight)
        .sum()
}

fn weighted_sum_abs_difference(lhs: &[f64], rhs: &[f64], weights: &[f64]) -> f64 {
    lhs.iter()
        .zip(rhs.iter())
        .zip(weights.iter())
        .map(|((lhs, rhs), weight)| (lhs - rhs).abs() * weight)
        .sum()
}

fn dot(lhs: &[f64], rhs: &[f64]) -> f64 {
    lhs.iter().zip(rhs.iter()).map(|(lhs, rhs)| lhs * rhs).sum()
}
