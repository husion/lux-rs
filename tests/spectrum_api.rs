mod common;

use common::{observer_1931, spectrum_400_420, spectrum_555_556};
use lux_rs::{getwld, getwlr, SpectralMatrix, Spectrum, SpectrumNormalization, WavelengthGrid};

#[test]
fn grid_matches_luxpy_style_range() {
    let wl = getwlr(WavelengthGrid::new(360.0, 365.0, 1.0).unwrap()).unwrap();
    assert_eq!(wl, vec![360.0, 361.0, 362.0, 363.0, 364.0, 365.0]);
}

#[test]
fn unequal_spacing_matches_luxpy_formula() {
    let dl = getwld(&[400.0, 410.0, 430.0]).unwrap();
    assert_eq!(dl, vec![10.0, 15.0, 20.0]);
}

#[test]
fn validates_spectrum_lengths() {
    let result = Spectrum::new(vec![400.0], vec![1.0, 2.0]);
    assert!(result.is_err());
}

#[test]
fn linearly_interpolates_and_extrapolates() {
    let spectrum = spectrum_400_420();
    let resampled = spectrum
        .interpolate_linear(&[395.0, 405.0, 420.0, 425.0])
        .unwrap();
    assert_eq!(resampled.values(), &[0.5, 1.5, 3.0, 3.5]);
}

#[test]
fn constructs_spectral_matrix() {
    let matrix = SpectralMatrix::new(
        vec![400.0, 410.0, 420.0],
        vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]],
    )
    .unwrap();
    assert_eq!(matrix.spectrum_count(), 2);
    assert_eq!(matrix.wavelength_count(), 3);
}

#[test]
fn rejects_spectral_matrix_with_bad_row_length() {
    let result = SpectralMatrix::new(
        vec![400.0, 410.0, 420.0],
        vec![vec![1.0, 2.0], vec![4.0, 5.0, 6.0]],
    );
    assert!(result.is_err());
}

#[test]
fn interpolates_spectral_matrix_linearly() {
    let matrix = SpectralMatrix::new(
        vec![400.0, 410.0, 420.0],
        vec![vec![1.0, 2.0, 3.0], vec![2.0, 3.0, 4.0]],
    )
    .unwrap();
    let resampled = matrix
        .cie_interp_linear(&[395.0, 405.0, 420.0, 425.0], true)
        .unwrap();
    assert_eq!(resampled.spectra()[0], vec![0.5, 1.5, 3.0, 3.5]);
    assert_eq!(resampled.spectra()[1], vec![1.5, 2.5, 4.0, 4.5]);
}

#[test]
fn clips_negative_values_when_requested() {
    let spectrum = Spectrum::new(vec![400.0, 410.0], vec![1.0, 0.0]).unwrap();
    let resampled = spectrum.cie_interp_linear(&[420.0], false).unwrap();
    assert_eq!(resampled.values(), &[0.0]);
}

#[test]
fn normalizes_spectrum_to_max() {
    let spectrum = spectrum_400_420();
    let normalized = spectrum
        .normalize(SpectrumNormalization::Max(2.0), None)
        .unwrap();
    assert_eq!(normalized.values(), &[2.0 / 3.0, 4.0 / 3.0, 2.0]);
}

#[test]
fn normalizes_spectrum_to_area() {
    let spectrum = spectrum_400_420();
    let normalized = spectrum
        .normalize(SpectrumNormalization::Area(1.0), None)
        .unwrap();
    assert_eq!(normalized.values(), &[1.0 / 60.0, 2.0 / 60.0, 3.0 / 60.0]);
}

#[test]
fn normalizes_spectrum_to_lambda() {
    let spectrum = spectrum_400_420();
    let normalized = spectrum
        .normalize(SpectrumNormalization::Lambda(410.0), None)
        .unwrap();
    assert_eq!(normalized.values(), &[0.5, 1.0, 1.5]);
}

#[test]
fn normalizes_spectrum_to_photometric_power() {
    let observer = observer_1931();
    let spectrum = spectrum_555_556();
    let normalized = spectrum
        .normalize(SpectrumNormalization::Photometric(1000.0), Some(&observer))
        .unwrap();
    assert!((normalized.values()[0] - 0.732_114_734_022_806_9).abs() < 1e-12);
    assert!((normalized.values()[1] - 0.732_114_734_022_806_9).abs() < 1e-12);
}

#[test]
fn one_row_matrix_normalization_preserves_numeric_baselines() {
    let observer = observer_1931();
    let matrix = SpectralMatrix::new(vec![555.0, 556.0], vec![vec![1.0, 1.0]]).unwrap();

    let normalized = matrix
        .normalize_each(&[SpectrumNormalization::Photometric(1000.0)], Some(&observer))
        .unwrap();

    assert_eq!(normalized.spectrum_count(), 1);
    assert!((normalized.spectra()[0][0] - 0.732_114_734_022_806_9).abs() < 1e-12);
    assert!((normalized.spectra()[0][1] - 0.732_114_734_022_806_9).abs() < 1e-12);
}

#[test]
fn normalizes_each_spectrum_in_matrix() {
    let matrix = SpectralMatrix::new(
        vec![400.0, 410.0, 420.0],
        vec![vec![1.0, 2.0, 3.0], vec![2.0, 4.0, 6.0]],
    )
    .unwrap();
    let normalized = matrix
        .normalize_each(
            &[
                SpectrumNormalization::Max(2.0),
                SpectrumNormalization::Area(1.0),
            ],
            None,
        )
        .unwrap();
    assert_eq!(normalized.spectra()[0], vec![2.0 / 3.0, 4.0 / 3.0, 2.0]);
    assert_eq!(
        normalized.spectra()[1],
        vec![2.0 / 120.0, 4.0 / 120.0, 6.0 / 120.0]
    );
}
