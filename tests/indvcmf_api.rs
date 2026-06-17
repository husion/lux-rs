mod common;

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

use lux_rs::{
    individual_observer_categorical_observers, individual_observer_cmf,
    individual_observer_cmf_aicom_plus, individual_observer_cmf_stockman2023,
    individual_observer_cmf_with_source, individual_observer_default_std_devs,
    individual_observer_generate, individual_observer_generate_population,
    individual_observer_lms_to_xyz_matrix, individual_observer_lms_to_xyz_matrix_stockman2023,
    individual_observer_monte_carlo, individual_observer_monte_carlo_parameters,
    individual_observer_us_census_age_distribution, IndividualObserverCategoricalOptions,
    IndividualObserverDataSource, IndividualObserverMonteCarloOptions,
    IndividualObserverParameters, IndividualObserverPopulationRequest,
    IndividualObserverPopulationStrategy, IndividualObserverRequest,
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
    let python = common::resolve_python_executable(&root);
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

    let shifted = individual_observer_cmf(IndividualObserverParameters {
        age: 45.0,
        field_size: 5.0,
        cone_peak_shift: [1.5, -0.75, 0.5],
        ..Default::default()
    })
    .unwrap();
    assert_vec_close(
        &flatten_sample_columns(&shifted.lms, &sample_idxs),
        &parse_vec(&baselines["indvcmf_shifted_lms_samples"]),
        1e-9,
    );
    assert_vec_close(
        &flatten_sample_columns(&shifted.xyz, &sample_idxs),
        &parse_vec(&baselines["indvcmf_shifted_xyz_samples"]),
        1e-9,
    );
}

#[test]
fn supports_cietc197_data_source_for_single_observer() {
    let observer = individual_observer_cmf_with_source(
        IndividualObserverParameters::default(),
        IndividualObserverDataSource::CieTc197,
    )
    .unwrap();
    assert!(observer.lms.wavelength_count() > 1000);
    assert_eq!(observer.lms.wavelengths().first().copied(), Some(390.0));
    assert_eq!(observer.lms.wavelengths().last().copied(), Some(830.0));
    assert_ne!(
        observer.lms_to_xyz_matrix,
        individual_observer_lms_to_xyz_matrix(IndividualObserverParameters::default().field_size)
    );
}

#[test]
fn cietc197_matrix_fit_is_finite_and_changes_xyz_output() {
    let params = IndividualObserverParameters::default();
    let cietc = individual_observer_cmf_with_source(params, IndividualObserverDataSource::CieTc197)
        .unwrap();
    let asano_xyz = lux_rs::individual_observer_lms_to_xyz(
        &cietc.lms,
        params.field_size,
        params.allow_negative_xyz_values,
    )
    .unwrap();

    for row in cietc.lms_to_xyz_matrix {
        for value in row {
            assert!(value.is_finite());
        }
    }
    let diff = cietc
        .xyz
        .spectra()
        .iter()
        .zip(asano_xyz.spectra().iter())
        .flat_map(|(lhs, rhs)| lhs.iter().zip(rhs.iter()))
        .map(|(lhs, rhs)| (lhs - rhs).abs())
        .sum::<f64>();
    assert!(diff > 0.0);
}

#[test]
fn monte_carlo_parameter_sampling_is_seeded_and_bounded() {
    let options = IndividualObserverMonteCarloOptions {
        n_observers: 8,
        field_size: 10.0,
        age_pool: vec![22.0, 30.0, 45.0],
        seed: 7,
        ..Default::default()
    };
    let first = individual_observer_monte_carlo_parameters(&options).unwrap();
    let second = individual_observer_monte_carlo_parameters(&options).unwrap();
    assert_eq!(first, second);
    for params in &first {
        assert!(options.age_pool.contains(&params.age));
        assert!(params.lens_density_variation >= -100.0);
        assert!(params.macular_density_variation >= -100.0);
        assert!(params
            .cone_density_variation
            .iter()
            .all(|value| *value >= -100.0));
    }
}

#[test]
fn monte_carlo_generation_returns_population() {
    let population = individual_observer_monte_carlo(IndividualObserverMonteCarloOptions {
        n_observers: 3,
        field_size: 5.0,
        age_pool: vec![32.0, 40.0],
        seed: 1234,
        data_source: IndividualObserverDataSource::Asano,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(population.parameters.len(), 3);
    assert_eq!(population.cmfs.len(), 3);
    assert!(population
        .cmfs
        .iter()
        .all(|cmf| cmf.lms.wavelength_count() == 79));
}

#[test]
fn categorical_observer_generation_matches_requested_count() {
    let population = individual_observer_categorical_observers(
        4,
        2.0,
        IndividualObserverDataSource::Asano,
        false,
    )
    .unwrap();
    assert_eq!(population.parameters.len(), 4);
    assert_eq!(population.cmfs.len(), 4);
}

#[test]
fn us_census_age_distribution_is_available() {
    let ages = individual_observer_us_census_age_distribution().unwrap();
    assert!(!ages.is_empty());
    assert!(ages.iter().all(|age| *age >= 10.0 && *age <= 70.0));
}

#[test]
fn stockman2023_scheme_generates_observer() {
    let observer = individual_observer_cmf_stockman2023(IndividualObserverParameters {
        age: 32.0,
        field_size: 2.0,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(observer.lms.wavelength_count(), 79);
    assert_eq!(observer.xyz.wavelength_count(), 79);
    assert_eq!(
        observer.lms_to_xyz_matrix,
        individual_observer_lms_to_xyz_matrix_stockman2023(2.0)
    );
}

#[test]
fn stockman2023_log_shift_changes_lms_shape() {
    let base =
        individual_observer_cmf_stockman2023(IndividualObserverParameters::default()).unwrap();
    let shifted = individual_observer_cmf_stockman2023(IndividualObserverParameters {
        cone_peak_shift: [2.0, 0.0, 0.0],
        ..Default::default()
    })
    .unwrap();
    let base_l = &base.lms.spectra()[0];
    let shifted_l = &shifted.lms.spectra()[0];
    let diff = base_l
        .iter()
        .zip(shifted_l.iter())
        .map(|(lhs, rhs)| (lhs - rhs).abs())
        .sum::<f64>();
    assert!(diff > 0.0);
}

#[test]
fn stockman2023_lms_to_xyz_matrix_matches_paper_endpoints() {
    let m2 = individual_observer_lms_to_xyz_matrix_stockman2023(2.0);
    let m10 = individual_observer_lms_to_xyz_matrix_stockman2023(10.0);
    let expected_2 = [
        [1.947_354_69, -1.414_451_23, 0.364_763_27],
        [0.689_902_72, 0.348_321_89, 0.0],
        [0.0, 0.0, 1.934_853_43],
    ];
    let expected_10 = [
        [1.939_864_43, -1.346_643_59, 0.430_449_35],
        [0.692_839_32, 0.349_675_67, 0.0],
        [0.0, 0.0, 2.146_879_45],
    ];
    for row in 0..3 {
        for col in 0..3 {
            assert_close(m2[row][col], expected_2[row][col], 1e-12);
            assert_close(m10[row][col], expected_10[row][col], 1e-12);
        }
    }
}

#[test]
fn unified_single_observer_request_matches_legacy_entrypoint() {
    let parameters = IndividualObserverParameters {
        age: 40.0,
        field_size: 5.0,
        cone_peak_shift: [1.0, -0.5, 0.2],
        ..Default::default()
    };
    let legacy =
        individual_observer_cmf_with_source(parameters, IndividualObserverDataSource::Stockman2023)
            .unwrap();
    let unified = individual_observer_generate(IndividualObserverRequest {
        model: IndividualObserverDataSource::Stockman2023,
        parameters,
    })
    .unwrap();
    assert_eq!(legacy, unified);
}

#[test]
fn unified_population_request_supports_model_and_strategy() {
    let population = individual_observer_generate_population(IndividualObserverPopulationRequest {
        model: IndividualObserverDataSource::CieTc197,
        strategy: IndividualObserverPopulationStrategy::Categorical(
            IndividualObserverCategoricalOptions {
                n_categories: 3,
                field_size: 2.0,
                allow_negative_xyz_values: false,
            },
        ),
    })
    .unwrap();
    assert_eq!(population.cmfs.len(), 3);
    assert!(population
        .cmfs
        .iter()
        .all(|cmf| cmf.lms.wavelength_count() > 1000));
}

#[test]
fn aicom_plus_model_generates_valid_observer() {
    let observer = individual_observer_cmf_aicom_plus(IndividualObserverParameters {
        age: 42.0,
        field_size: 6.0,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(observer.lms.wavelength_count(), 79);
    assert_eq!(observer.xyz.wavelength_count(), 79);
}

#[test]
fn aicom_plus_differs_from_asano_under_same_parameters() {
    let parameters = IndividualObserverParameters {
        age: 42.0,
        field_size: 6.0,
        ..Default::default()
    };
    let asano =
        individual_observer_cmf_with_source(parameters, IndividualObserverDataSource::Asano)
            .unwrap();
    let aicom =
        individual_observer_cmf_with_source(parameters, IndividualObserverDataSource::AicomPlus)
            .unwrap();
    let diff = aicom
        .xyz
        .spectra()
        .iter()
        .zip(asano.xyz.spectra().iter())
        .flat_map(|(lhs, rhs)| lhs.iter().zip(rhs.iter()))
        .map(|(lhs, rhs)| (lhs - rhs).abs())
        .sum::<f64>();
    assert!(diff > 0.0);
}

#[test]
fn test_individual_observer_cmf_from_measured() {
    use lux_rs::{
        individual_observer_cmf_from_measured, IndividualObserverMeasuredParameters,
        LMS_TO_XYZ_2DEG_FIXED, LMS_TO_XYZ_10DEG_FIXED,
    };

    let params = IndividualObserverMeasuredParameters {
        lshift: 0.0,
        mshift: 0.0,
        sshift: 0.0,
        lod: 0.38,
        mod_: 0.38,
        sod: 0.2,
        mac: 0.35,
        lens: 1.7649,
        field_size: 2.0,
    };

    let wavelengths = vec![400.0, 500.0, 600.0, 700.0];
    let cmf = individual_observer_cmf_from_measured(&wavelengths, params, None).unwrap();

    assert_eq!(cmf.lms.wavelengths(), &wavelengths);
    assert_eq!(cmf.xyz.wavelengths(), &wavelengths);
    assert_eq!(cmf.lms_to_xyz_matrix, LMS_TO_XYZ_2DEG_FIXED);

    // Verify 10-degree field size uses LMS_TO_XYZ_10DEG_FIXED
    let params_10 = IndividualObserverMeasuredParameters {
        field_size: 10.0,
        ..params
    };
    let cmf_10 = individual_observer_cmf_from_measured(&wavelengths, params_10, None).unwrap();
    assert_eq!(cmf_10.lms_to_xyz_matrix, LMS_TO_XYZ_10DEG_FIXED);

    // Verify error on empty wavelengths
    let err = individual_observer_cmf_from_measured(&[], params, None);
    assert!(err.is_err());
}

