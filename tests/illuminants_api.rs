mod common;

use common::{grid_360_365, grid_380_385};
use lux_rs::{
    blackbody, cct_to_xyz, cri_ref, daylightlocus, daylightphase, standard_illuminant,
    standard_illuminant_names, xyz_to_cct, Observer, WavelengthGrid,
};

#[test]
fn computes_relative_blackbody_spectrum() {
    let spectrum = blackbody(6500.0, Some(grid_360_365()), None, true).unwrap();

    assert_eq!(
        spectrum.wavelengths(),
        &[360.0, 361.0, 362.0, 363.0, 364.0, 365.0]
    );
    assert!((spectrum.values()[0] - 0.995_951_160_513_301_6).abs() < 1e-12);
    assert!((spectrum.values()[5] - 1.011_458_643_342_760_8).abs() < 1e-12);
}

#[test]
fn computes_absolute_blackbody_spectrum() {
    let spectrum = blackbody(
        6500.0,
        Some(WavelengthGrid::new(560.0, 560.0, 1.0).unwrap()),
        None,
        false,
    )
    .unwrap();

    assert!((spectrum.values()[0] - 42_340_048_320_714.19).abs() < 1e-3);
}

#[test]
fn normalizes_relative_blackbody_to_560_nm() {
    let spectrum = blackbody(
        6500.0,
        Some(WavelengthGrid::new(560.0, 560.0, 1.0).unwrap()),
        None,
        true,
    )
    .unwrap();

    assert!((spectrum.values()[0] - 1.0).abs() < 1e-12);
}

#[test]
fn rejects_non_positive_blackbody_inputs() {
    let result = blackbody(0.0, None, Some(1.0), true);
    assert!(result.is_err());

    let result = blackbody(6500.0, None, Some(0.0), true);
    assert!(result.is_err());
}

#[test]
fn computes_daylight_locus_for_6500k() {
    let [x_d, y_d] = daylightlocus(6500.0, false, false).unwrap();
    assert!((x_d - 0.312_778_876_194_811_15).abs() < 1e-12);
    assert!((y_d - 0.329_183_498_498_040_96).abs() < 1e-12);
}

#[test]
fn computes_daylightphase_spectrum() {
    let spectrum = daylightphase(6500.0, Some(grid_360_365()), false, false, None).unwrap();

    assert!((spectrum.values()[0] - 0.944_787_686_958_161_2).abs() < 1e-12);
    assert!((spectrum.values()[5] - 1.0).abs() < 1e-12);
}

#[test]
fn daylightphase_uses_blackbody_below_4000k() {
    let spectrum = daylightphase(3500.0, Some(grid_360_365()), false, false, None).unwrap();

    assert!((spectrum.values()[0] - 0.154_169_092_129_982_6).abs() < 1e-12);
    assert!((spectrum.values()[5] - 0.168_260_026_564_076_65).abs() < 1e-12);
}

#[test]
fn computes_default_cri_ref_spectra() {
    let spectra = cri_ref(&[3000.0, 6500.0], Some(grid_360_365())).unwrap();

    assert_eq!(spectra.spectrum_count(), 2);
    assert!((spectra.spectra()[0][0] - 0.078_162_762_564_806_27).abs() < 1e-12);
    assert!((spectra.spectra()[0][5] - 0.087_559_911_695_420_67).abs() < 1e-12);
    assert!((spectra.spectra()[1][0] - 0.944_787_686_958_161_2).abs() < 1e-12);
    assert!((spectra.spectra()[1][5] - 1.0).abs() < 1e-12);
}

#[test]
fn converts_cct_to_xyz() {
    let xyz = cct_to_xyz(6500.0, Observer::Cie1931_2).unwrap();
    assert!((xyz[0] - 96.878_415_094_956_67).abs() < 1e-6);
    assert!((xyz[1] - 100.0).abs() < 1e-9);
    assert!((xyz[2] - 112.116_528_133_993_16).abs() < 1e-6);
}

#[test]
fn converts_xyz_to_cct() {
    let (cct, duv) = xyz_to_cct([100.0, 100.0, 100.0], Observer::Cie1931_2).unwrap();
    assert!((cct - 5_455.485_887_350_497).abs() < 1.0);
    assert!((duv - (-0.004_423_324_748_595_847)).abs() < 1e-4);
}

#[test]
fn loads_standard_illuminant_a() {
    let spectrum = standard_illuminant("A", Some(grid_360_365())).unwrap();
    assert!((spectrum.values()[0] - 6.144_62).abs() < 1e-12);
    assert!((spectrum.values()[5] - 6.947_2).abs() < 1e-12);
}

#[test]
fn loads_standard_illuminant_f4() {
    let spectrum = standard_illuminant("F4", Some(grid_380_385())).unwrap();
    assert!((spectrum.values()[0] - 0.57).abs() < 1e-12);
    assert!((spectrum.values()[5] - 0.7).abs() < 1e-12);
}

#[test]
fn loads_standard_illuminant_led_b1() {
    let spectrum = standard_illuminant("LED_B1", Some(grid_380_385())).unwrap();
    assert!((spectrum.values()[0] - 0.0).abs() < 1e-12);
    assert!((spectrum.values()[5] - 0.01).abs() < 1e-12);
}

#[test]
fn loads_nominal_daylight_illuminants() {
    let spectrum = standard_illuminant("D50", Some(grid_360_365())).unwrap();
    assert!((spectrum.values()[0] - 0.940_694_581_416_273_4).abs() < 1e-12);
    assert!((spectrum.values()[5] - 1.0).abs() < 1e-12);
}

#[test]
fn exposes_standard_illuminant_names() {
    assert_eq!(
        standard_illuminant_names(),
        &[
            "A", "D50", "D55", "D65", "D75", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9",
            "F10", "F11", "F12", "LED_B1", "LED_B2", "LED_B3", "LED_B4", "LED_B5", "LED_BH1",
            "LED_RGB1", "LED_V1", "LED_V2",
        ]
    );
}

#[test]
fn resolves_standard_illuminant_aliases() {
    let canonical = standard_illuminant("D65", None).unwrap();

    for alias in ["d65", "CIE D65", "cie-d65", " cie_d65 "] {
        let spectrum = standard_illuminant(alias, None).expect(alias);
        assert_eq!(spectrum, canonical);
    }

    let f4 = standard_illuminant("F4", None).unwrap();
    assert_eq!(standard_illuminant("cie f4", None).unwrap(), f4);

    let led = standard_illuminant("LED_B1", None).unwrap();
    assert_eq!(standard_illuminant("led-b1", None).unwrap(), led);
    assert_eq!(standard_illuminant("CIELED_B1", None).unwrap(), led);
}

#[test]
fn rejects_unknown_illuminant_aliases() {
    assert!(standard_illuminant("D93", None).is_err());
    assert!(standard_illuminant("LED_X1", None).is_err());
    assert!(standard_illuminant("", None).is_err());
}
