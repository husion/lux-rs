use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

use lux_rs::{
    individual_observer_cmf, individual_observer_default_std_devs,
    individual_observer_lms_to_xyz_matrix, IndividualObserverParameters,
};

fn parse_vec(value: &str) -> Vec<f64> {
    value
        .split(',')
        .map(|item| item.parse::<f64>().unwrap())
        .collect()
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    let diff = (actual - expected).abs();
    assert!(
        diff <= tolerance,
        "expected {expected}, got {actual}, diff {diff}, tolerance {tolerance}"
    );
}

fn assert_vec_close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len(), "length mismatch");
    for (&actual, &expected) in actual.iter().zip(expected.iter()) {
        assert_close(actual, expected, tolerance);
    }
}

fn flatten_sample_columns(matrix: &lux_rs::Spectrum, idxs: &[usize]) -> Vec<f64> {
    let mut sampled = Vec::new();
    for &idx in idxs {
        sampled.push(matrix.wavelengths()[idx]);
        for spectrum in matrix.spectra() {
            sampled.push(spectrum[idx]);
        }
    }
    sampled
}

fn flatten_sample_values(spectrum: &lux_rs::Spectrum, idxs: &[usize]) -> Vec<f64> {
    let mut sampled = Vec::new();
    for &idx in idxs {
        sampled.push(spectrum.wavelengths()[idx]);
        sampled.push(spectrum.values()[idx]);
    }
    sampled
}

fn load_python_baselines() -> HashMap<String, String> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let python = root.join("luxpy/.venv/bin/python");
    let script = root.join("tests/python_ref/baseline_indvcmf.py");

    let output = Command::new(python)
        .env("MPLCONFIGDIR", "/tmp/mpl")
        .arg(script)
        .output()
        .expect("failed to run Python baseline script");

    assert!(
        output.status.success(),
        "python baseline script failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| {
            let (key, value) = line
                .split_once('=')
                .unwrap_or_else(|| panic!("invalid baseline line: {line}"));
            (key.to_string(), value.to_string())
        })
        .collect()
}

#[test]
fn exposes_default_std_devs() {
    let std_devs = individual_observer_default_std_devs();
    assert_eq!(std_devs.lens_density, 19.1);
    assert_eq!(std_devs.macular_density, 37.2);
    assert_eq!(std_devs.cone_density, [17.9, 17.9, 14.7]);
    assert_eq!(std_devs.cone_peak_shift, [4.0, 3.0, 2.5]);
}

#[test]
fn clamps_lms_to_xyz_matrix_to_supported_field_sizes() {
    assert_eq!(
        individual_observer_lms_to_xyz_matrix(1.0),
        individual_observer_lms_to_xyz_matrix(2.0)
    );
    assert_eq!(
        individual_observer_lms_to_xyz_matrix(12.0),
        individual_observer_lms_to_xyz_matrix(10.0)
    );
}

#[test]
fn rejects_non_zero_peak_shift_in_current_slice() {
    let parameters = IndividualObserverParameters {
        cone_peak_shift: [1.0, 0.0, 0.0],
        ..Default::default()
    };
    let error = individual_observer_cmf(parameters).unwrap_err();
    assert_eq!(
        error.to_string(),
        "invalid input: non-zero cone peak shifts are not yet supported in the Rust indvcmf slice"
    );
}

#[test]
fn matches_python_indvcmf_baselines() {
    let baselines = load_python_baselines();
    let sample_idxs = [0usize, 1, 33, 46, 47, 78];
    let trans_sample_idxs = [0usize, 33, 78];

    assert_vec_close(
        &individual_observer_lms_to_xyz_matrix(2.0)
            .into_iter()
            .flatten()
            .collect::<Vec<_>>(),
        &parse_vec(&baselines["indvcmf_matrix_2"]),
        1e-12,
    );
    assert_vec_close(
        &individual_observer_lms_to_xyz_matrix(5.0)
            .into_iter()
            .flatten()
            .collect::<Vec<_>>(),
        &parse_vec(&baselines["indvcmf_matrix_5"]),
        1e-12,
    );
    assert_vec_close(
        &individual_observer_lms_to_xyz_matrix(10.0)
            .into_iter()
            .flatten()
            .collect::<Vec<_>>(),
        &parse_vec(&baselines["indvcmf_matrix_10"]),
        1e-12,
    );

    let default_observer =
        individual_observer_cmf(IndividualObserverParameters::default()).unwrap();
    assert_eq!(default_observer.lms.wavelength_count(), 79);
    assert_eq!(default_observer.xyz.wavelength_count(), 79);
    assert_vec_close(
        &flatten_sample_columns(&default_observer.lms, &sample_idxs),
        &parse_vec(&baselines["indvcmf_default_lms_samples"]),
        1e-9,
    );
    assert_vec_close(
        &flatten_sample_columns(&default_observer.xyz, &sample_idxs),
        &parse_vec(&baselines["indvcmf_default_xyz_samples"]),
        1e-9,
    );

    let varied = individual_observer_cmf(IndividualObserverParameters {
        age: 60.0,
        field_size: 2.0,
        lens_density_variation: 15.0,
        macular_density_variation: -10.0,
        cone_density_variation: [5.0, -7.0, 3.0],
        cone_peak_shift: [0.0, 0.0, 0.0],
        allow_negative_xyz_values: false,
    })
    .unwrap();
    assert_vec_close(
        &flatten_sample_columns(&varied.lms, &sample_idxs),
        &parse_vec(&baselines["indvcmf_varied_lms_samples"]),
        1e-9,
    );
    assert_vec_close(
        &flatten_sample_columns(&varied.xyz, &sample_idxs),
        &parse_vec(&baselines["indvcmf_varied_xyz_samples"]),
        1e-9,
    );
    assert_vec_close(
        &flatten_sample_values(&varied.lens_transmission, &trans_sample_idxs),
        &parse_vec(&baselines["indvcmf_varied_lens_samples"]),
        1e-9,
    );
    assert_vec_close(
        &flatten_sample_values(&varied.macular_transmission, &trans_sample_idxs),
        &parse_vec(&baselines["indvcmf_varied_macular_samples"]),
        1e-9,
    );
}
