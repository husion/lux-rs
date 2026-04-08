mod common;

use common::{matrix_555_556, observer_1931, spectrum_555_556};
use lux_rs::{spd_to_ler, spd_to_power, spd_to_xyz, PowerType, SpectralMatrix, Spectrum};

#[test]
fn computes_radiometric_power() {
    let spectrum = Spectrum::new(vec![400.0, 410.0, 420.0], vec![1.0, 2.0, 3.0]).unwrap();
    let power = spd_to_power(&spectrum, PowerType::Radiometric, None).unwrap();
    assert_eq!(power, 60.0);
}

#[test]
fn computes_quantal_power() {
    let spectrum = Spectrum::new(vec![500.0, 510.0], vec![1.0, 1.0]).unwrap();
    let power = spd_to_power(&spectrum, PowerType::Quantal, None).unwrap();
    assert!(power.is_finite());
    assert!(power > 0.0);
}

#[test]
fn computes_photometric_power() {
    let observer = observer_1931();
    let spectrum = spectrum_555_556();
    let power = spd_to_power(&spectrum, PowerType::Photometric, Some(&observer)).unwrap();
    assert!(power > 680.0);
}

#[test]
fn photometric_power_requires_observer() {
    let spectrum = Spectrum::new(vec![555.0, 556.0], vec![1.0, 1.0]).unwrap();
    let result = spd_to_power(&spectrum, PowerType::Photometric, None);
    assert!(result.is_err());
}

#[test]
fn computes_relative_xyz() {
    let observer = observer_1931();
    let spectrum = spectrum_555_556();
    let xyz = spd_to_xyz(&spectrum, &observer, true).unwrap();
    assert!((xyz[0] - 52.021_027_306_606_52).abs() < 1e-9);
    assert!((xyz[1] - 100.0).abs() < 1e-12);
    assert!((xyz[2] - 0.552_719_552_355_926_3).abs() < 1e-9);
}

#[test]
fn one_row_batch_photometry_preserves_numeric_baselines() {
    let observer = observer_1931();
    let spectra = SpectralMatrix::new(vec![555.0, 556.0], vec![vec![1.0, 1.0]]).unwrap();

    let xyz = spectra.spd_to_xyz(&observer, true).unwrap();
    assert_eq!(xyz.len(), 1);
    assert!((xyz[0][0] - 52.021_027_306_606_52).abs() < 1e-9);
    assert!((xyz[0][1] - 100.0).abs() < 1e-12);
    assert!((xyz[0][2] - 0.552_719_552_355_926_3).abs() < 1e-9);

    let ler = spectra.spd_to_ler(&observer).unwrap();
    assert_eq!(ler.len(), 1);
    assert!((ler[0] - 682.953_062_906_7).abs() < 1e-9);
}

#[test]
fn computes_relative_xyz_for_multiple_spectra() {
    let observer = observer_1931();
    let spectra = matrix_555_556();
    let xyz = spectra.spd_to_xyz(&observer, true).unwrap();
    assert_eq!(xyz.len(), 2);
    assert!((xyz[0][0] - 52.021_027_306_606_52).abs() < 1e-9);
    assert!((xyz[1][1] - 100.0).abs() < 1e-12);
}

#[test]
fn computes_ler() {
    let observer = observer_1931();
    let spectrum = spectrum_555_556();
    let ler = spd_to_ler(&spectrum, &observer).unwrap();
    assert!((ler - 682.953_062_906_7).abs() < 1e-9);
}

#[test]
fn computes_ler_for_multiple_spectra() {
    let observer = observer_1931();
    let spectra = matrix_555_556();
    let ler = spectra.spd_to_ler(&observer).unwrap();
    assert_eq!(ler.len(), 2);
    assert!((ler[0] - 682.953_062_906_7).abs() < 1e-9);
    assert!((ler[1] - 682.953_062_906_7).abs() < 1e-9);
}
