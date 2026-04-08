mod common;

use common::{observer_1931, WHITE_D65, WHITE_E, XYZ_BRIGHT, XYZ_SAMPLE, XYZ_SAMPLE_ALT};
use core::str::FromStr;
use lux_rs::CamSurround as ModelCamSurround;
use lux_rs::{
    cam16_viewing_conditions, ciecam02_viewing_conditions, delta_e, delta_e_cie76,
    delta_e_ciede2000, get_cie_mesopic_adaptation, lab_to_xyz, lms_to_xyz, luv_to_xyz, srgb_to_xyz,
    vlbar_cie_mesopic, xyz_to_lab, xyz_to_lms, xyz_to_luv, xyz_to_srgb, xyz_to_yuv, xyz_to_yxy,
    yuv_to_xyz, yxy_to_xyz, CamUcsType, DeltaEFormula, Observer, Tristimulus,
};

#[test]
fn loads_standard_observer() {
    let observer = observer_1931();
    assert_eq!(observer.wavelengths.first().copied(), Some(360.0));
    assert_eq!(observer.wavelengths.last().copied(), Some(830.0));
    assert_eq!(observer.wavelengths.len(), 471);
}

#[test]
fn lists_canonical_observer_names() {
    let names: Vec<_> = Observer::all()
        .iter()
        .map(|observer| observer.name())
        .collect();
    assert_eq!(
        names,
        vec!["1931_2", "1964_10", "2006_2", "2006_10", "2015_2", "2015_10"]
    );
}

#[test]
fn parses_observer_aliases() {
    assert_eq!(Observer::from_name("1931_2").unwrap(), Observer::Cie1931_2);
    assert_eq!(
        Observer::from_name("CIE 1964 10").unwrap(),
        Observer::Cie1964_10
    );
    assert_eq!(
        Observer::from_name("cie-2006_2").unwrap(),
        Observer::Cie2006_2
    );
    assert_eq!(
        Observer::from_name("2006-10").unwrap(),
        Observer::Cie2006_10
    );
    assert_eq!(Observer::from_name("2015 2").unwrap(), Observer::Cie2015_2);
    assert_eq!(
        Observer::from_name("2015_10").unwrap(),
        Observer::Cie2015_10
    );
    assert_eq!(Observer::from_name("1931").unwrap(), Observer::Cie1931_2);
    assert_eq!(Observer::from_name("1964").unwrap(), Observer::Cie1964_10);
}

#[test]
fn formats_and_parses_observers_via_str_traits() {
    let observer = Observer::Cie2006_10;
    assert_eq!(observer.to_string(), "2006_10");
    assert_eq!(Observer::from_str("CIE2006-10").unwrap(), observer);
}

#[test]
fn rejects_unknown_observer_aliases() {
    assert!(Observer::from_name("cie 2012 2").is_err());
    assert!(Observer::from_str("").is_err());
}

#[test]
fn loads_cie2006_standard_observers() {
    let observer_2 = Observer::Cie2006_2.standard().unwrap();
    assert_eq!(observer_2.wavelengths.first().copied(), Some(390.0));
    assert_eq!(observer_2.wavelengths.last().copied(), Some(830.0));
    assert_eq!(observer_2.wavelengths.len(), 441);

    let observer_10 = Observer::Cie2006_10.standard().unwrap();
    assert_eq!(observer_10.wavelengths.first().copied(), Some(390.0));
    assert_eq!(observer_10.wavelengths.last().copied(), Some(830.0));
    assert_eq!(observer_10.wavelengths.len(), 441);
}

#[test]
fn loads_cie2015_standard_observers() {
    let observer_2 = Observer::Cie2015_2.standard().unwrap();
    assert_eq!(observer_2.wavelengths.first().copied(), Some(390.0));
    assert_eq!(observer_2.wavelengths.last().copied(), Some(830.0));
    assert_eq!(observer_2.wavelengths.len(), 441);

    let observer_10 = Observer::Cie2015_10.standard().unwrap();
    assert_eq!(observer_10.wavelengths.first().copied(), Some(390.0));
    assert_eq!(observer_10.wavelengths.last().copied(), Some(830.0));
    assert_eq!(observer_10.wavelengths.len(), 441);
}

#[test]
fn exposes_xyzbar() {
    let xyzbar = Observer::Cie1931_2.xyzbar().unwrap();
    assert_eq!(xyzbar.wavelength_count(), 471);
    assert_eq!(xyzbar.spectrum_count(), 3);
}

#[test]
fn exposes_vlbar_and_k() {
    let (vl, k) = Observer::Cie1931_2.vlbar().unwrap();
    assert_eq!(vl.wavelengths().len(), 471);
    assert_eq!(vl.values()[195], 1.0);
    assert_eq!(k, 683.002);
}

#[test]
fn exposes_vlbar_and_k_for_cie2006_observers() {
    let (vl_2, k_2) = Observer::Cie2006_2.vlbar().unwrap();
    assert_eq!(vl_2.wavelengths().len(), 441);
    assert!(vl_2.values().iter().all(|value| value.is_finite()));
    assert_eq!(k_2, 683.358);

    let (vl_10, k_10) = Observer::Cie2006_10.vlbar().unwrap();
    assert_eq!(vl_10.wavelengths().len(), 441);
    assert!(vl_10.values().iter().all(|value| value.is_finite()));
    assert_eq!(k_10, 683.144);
}

#[test]
fn exposes_vlbar_and_k_for_cie2015_observers() {
    let (vl_2, k_2) = Observer::Cie2015_2.vlbar().unwrap();
    assert_eq!(vl_2.wavelengths().len(), 441);
    assert!(vl_2.values().iter().all(|value| value.is_finite()));
    assert_eq!(k_2, 683.358);

    let (vl_10, k_10) = Observer::Cie2015_10.vlbar().unwrap();
    assert_eq!(vl_10.wavelengths().len(), 441);
    assert!(vl_10.values().iter().all(|value| value.is_finite()));
    assert_eq!(k_10, 683.144);
}

#[test]
fn interpolates_xyzbar_linearly() {
    let xyzbar = Observer::Cie1931_2
        .xyzbar_linear(&[554.5, 555.0, 555.5, 556.0])
        .unwrap();
    assert!((xyzbar.spectra()[0][0] - 0.504_010_7).abs() < 1e-9);
    assert!((xyzbar.spectra()[1][1] - 1.0).abs() < 1e-12);
    assert!((xyzbar.spectra()[2][3] - 0.005_303_6).abs() < 1e-9);
}

#[test]
fn computes_cie_mesopic_adaptation_from_s_p_ratio() {
    let (lmes, m_values) = get_cie_mesopic_adaptation(&[1.0], None, Some(&[1.0])).unwrap();
    assert!((lmes[0] - 1.0).abs() < 1e-12);
    assert!((m_values[0] - 0.767).abs() < 1e-12);
}

#[test]
fn computes_mesopic_luminous_efficiency_curve() {
    let mesopic = vlbar_cie_mesopic(&[0.5, 1.0], None).unwrap();
    assert_eq!(mesopic.curves.spectrum_count(), 2);
    assert_eq!(mesopic.curves.wavelength_count(), 471);
    assert!((mesopic.k_mesopic[0] - 974.322_396_576_319_4).abs() < 1e-9);
    assert!((mesopic.k_mesopic[1] - 683.0).abs() < 1e-12);
    assert!((mesopic.curves.spectra()[0][195] - 0.837_061_500_974_263_2).abs() < 1e-9);
    assert!((mesopic.curves.spectra()[1][195] - 1.0).abs() < 1e-12);
}

#[test]
fn converts_xyz_to_yxy() {
    let yxy = xyz_to_yxy(XYZ_SAMPLE);
    assert!((yxy[0] - 0.5).abs() < 1e-12);
    assert!((yxy[1] - 0.25).abs() < 1e-12);
    assert!((yxy[2] - 0.5).abs() < 1e-12);
}

#[test]
fn converts_yxy_to_xyz() {
    let xyz = yxy_to_xyz([0.5, 0.25, 0.5]);
    assert!((xyz[0] - 0.25).abs() < 1e-12);
    assert!((xyz[1] - 0.5).abs() < 1e-12);
    assert!((xyz[2] - 0.25).abs() < 1e-12);
}

#[test]
fn converts_xyz_to_yuv() {
    let yuv = xyz_to_yuv([0.25, 0.5, 0.25]);
    assert!((yuv[0] - 0.5).abs() < 1e-12);
    assert!((yuv[1] - 0.117_647_058_823_529_41).abs() < 1e-12);
    assert!((yuv[2] - 0.529_411_764_705_882_4).abs() < 1e-12);
}

#[test]
fn converts_yuv_to_xyz() {
    let xyz = yuv_to_xyz([0.5, 0.117_647_058_823_529_41, 0.529_411_764_705_882_4]);
    assert!((xyz[0] - 0.25).abs() < 1e-12);
    assert!((xyz[1] - 0.5).abs() < 1e-12);
    assert!((xyz[2] - 0.25).abs() < 1e-12);
}

#[test]
fn converts_xyz_to_lab() {
    let lab = xyz_to_lab(XYZ_SAMPLE, WHITE_E);
    assert!((lab[0] - 100.0).abs() < 1e-12);
    assert!((lab[1] + 103.149_737_007_950_17).abs() < 1e-9);
    assert!((lab[2] - 41.259_894_803_180_07).abs() < 1e-9);
}

#[test]
fn converts_lab_to_xyz() {
    let xyz = lab_to_xyz(
        [100.0, -103.149_737_007_950_17, 41.259_894_803_180_07],
        WHITE_E,
    );
    assert!((xyz[0] - 0.25).abs() < 1e-9);
    assert!((xyz[1] - 0.5).abs() < 1e-12);
    assert!((xyz[2] - 0.25).abs() < 1e-9);
}

#[test]
fn converts_xyz_to_luv() {
    let luv = xyz_to_luv(XYZ_SAMPLE, WHITE_E);
    assert!((luv[0] - 100.0).abs() < 1e-12);
    assert!((luv[1] + 120.743_034_055_727_58).abs() < 1e-9);
    assert!((luv[2] - 72.445_820_433_436_54).abs() < 1e-9);
}

#[test]
fn converts_luv_to_xyz() {
    let xyz = luv_to_xyz(
        [100.0, -120.743_034_055_727_58, 72.445_820_433_436_54],
        WHITE_E,
    );
    assert!((xyz[0] - 0.25).abs() < 1e-9);
    assert!((xyz[1] - 0.5).abs() < 1e-12);
    assert!((xyz[2] - 0.25).abs() < 1e-9);
}

#[test]
fn computes_delta_e_cie76() {
    let xyz1 = lab_to_xyz([50.0, 2.5, -80.0], WHITE_D65);
    let xyz2 = lab_to_xyz([50.0, 0.0, -82.5], WHITE_D65);
    let delta = delta_e_cie76(xyz1, xyz2, WHITE_D65);
    assert!((delta - 3.535_533_905_932_737_8).abs() < 1e-12);
}

#[test]
fn computes_delta_e_ciede2000() {
    let xyz1 = lab_to_xyz([50.0, 2.6772, -79.7751], WHITE_D65);
    let xyz2 = lab_to_xyz([50.0, 0.0, -82.7485], WHITE_D65);
    let delta = delta_e_ciede2000(xyz1, xyz2, WHITE_D65);
    assert!((delta - 2.042_459_680_156_574).abs() < 1e-12);
}

#[test]
fn computes_delta_e_with_formula_dispatch() {
    let xyz1 = lab_to_xyz([50.0, 2.6772, -79.7751], WHITE_D65);
    let xyz2 = lab_to_xyz([50.0, 0.0, -82.7485], WHITE_D65);
    let delta = delta_e(xyz1, xyz2, WHITE_D65, DeltaEFormula::Ciede2000);
    assert!((delta - 2.042_459_680_156_574).abs() < 1e-12);
}

#[test]
fn computes_delta_e_from_xyz() {
    let lab1 = [50.0, 2.5, -80.0];
    let lab2 = [50.0, 0.0, -82.5];
    let xyz1 = lab_to_xyz(lab1, WHITE_D65);
    let xyz2 = lab_to_xyz(lab2, WHITE_D65);
    let delta = delta_e(xyz1, xyz2, WHITE_D65, DeltaEFormula::Cie76);
    assert!((delta - 3.535_533_905_932_737_8).abs() < 1e-9);
}

#[test]
fn batch_chromaticity_transforms_match_scalar_versions() {
    let xyz = [XYZ_SAMPLE, XYZ_SAMPLE_ALT];
    assert_eq!(
        Tristimulus::new(xyz.to_vec()).xyz_to_yxy().into_vec(),
        vec![xyz_to_yxy(xyz[0]), xyz_to_yxy(xyz[1])]
    );
    let yxy = [
        [0.5, 0.25, 0.5],
        [0.3, 0.222_222_222_222_222_2, 0.333_333_333_333_333_3],
    ];
    assert_eq!(
        Tristimulus::new(yxy.to_vec()).yxy_to_xyz().into_vec(),
        vec![yxy_to_xyz(yxy[0]), yxy_to_xyz(yxy[1])]
    );
    assert_eq!(
        Tristimulus::new(xyz.to_vec()).xyz_to_yuv().into_vec(),
        vec![xyz_to_yuv(xyz[0]), xyz_to_yuv(xyz[1])]
    );
    let yuv = [
        [0.5, 0.117_647_058_823_529_41, 0.529_411_764_705_882_4],
        [0.3, 0.129_032_258_064_516_13, 0.435_483_870_967_741_94],
    ];
    assert_eq!(
        Tristimulus::new(yuv.to_vec()).yuv_to_xyz().into_vec(),
        vec![yuv_to_xyz(yuv[0]), yuv_to_xyz(yuv[1])]
    );
}

#[test]
fn batch_color_space_transforms_match_scalar_versions() {
    let xyz = [XYZ_SAMPLE, XYZ_SAMPLE_ALT];
    let xyz_set = Tristimulus::new(xyz.to_vec());
    let lab = xyz_set.xyz_to_lab(WHITE_E).into_vec();
    assert_eq!(
        lab,
        vec![xyz_to_lab(xyz[0], WHITE_E), xyz_to_lab(xyz[1], WHITE_E)]
    );
    assert_eq!(
        Tristimulus::new(lab.clone())
            .lab_to_xyz(WHITE_E)
            .into_vec(),
        vec![lab_to_xyz(lab[0], WHITE_E), lab_to_xyz(lab[1], WHITE_E)]
    );
    let luv = xyz_set.xyz_to_luv(WHITE_E).into_vec();
    assert_eq!(
        luv,
        vec![xyz_to_luv(xyz[0], WHITE_E), xyz_to_luv(xyz[1], WHITE_E)]
    );
    assert_eq!(
        Tristimulus::new(luv.clone())
            .luv_to_xyz(WHITE_E)
            .into_vec(),
        vec![luv_to_xyz(luv[0], WHITE_E), luv_to_xyz(luv[1], WHITE_E)]
    );
}

#[test]
fn batch_lms_and_srgb_transforms_match_scalar_versions() {
    let xyz = [XYZ_SAMPLE, XYZ_BRIGHT];
    let lms_many = Tristimulus::new(vec![xyz[0], xyz[0]])
        .xyz_to_lms(Observer::Cie1931_2)
        .unwrap()
        .into_vec();
    assert_eq!(
        lms_many,
        vec![xyz_to_lms(xyz[0], Observer::Cie1931_2).unwrap(); 2]
    );
    let lms_input = [lms_many[0], lms_many[1]];
    assert_eq!(
        Tristimulus::new(lms_input.to_vec())
            .lms_to_xyz(Observer::Cie1931_2)
            .unwrap()
            .into_vec(),
        vec![
            lms_to_xyz(lms_input[0], Observer::Cie1931_2).unwrap(),
            lms_to_xyz(lms_input[1], Observer::Cie1931_2).unwrap()
        ]
    );
    assert_eq!(
        Tristimulus::new(vec![xyz[1], xyz[1]])
            .xyz_to_srgb(2.4, -0.055, true)
            .into_vec(),
        vec![xyz_to_srgb(xyz[1], 2.4, -0.055, true); 2]
    );
    let rgb = [[64.0, 128.0, 192.0], [32.0, 64.0, 96.0]];
    assert_eq!(
        Tristimulus::new(rgb.to_vec())
            .srgb_to_xyz(2.4, -0.055, true)
            .into_vec(),
        vec![
            srgb_to_xyz(rgb[0], 2.4, -0.055, true),
            srgb_to_xyz(rgb[1], 2.4, -0.055, true)
        ]
    );
}

#[test]
fn tristimulus_wrapper_matches_scalar_transforms() {
    let xyz = Tristimulus::new(vec![XYZ_SAMPLE]);
    assert_eq!(xyz.xyz_to_yxy().values(), &[xyz_to_yxy(XYZ_SAMPLE)]);
    assert_eq!(
        xyz.xyz_to_lab(WHITE_E).values(),
        &[xyz_to_lab(XYZ_SAMPLE, WHITE_E)]
    );
    assert_eq!(
        xyz.xyz_to_lms(Observer::Cie1931_2).unwrap().values(),
        &[xyz_to_lms(XYZ_SAMPLE, Observer::Cie1931_2).unwrap()]
    );
    let cam16_conditions = cam16_viewing_conditions(
        [95.047, 100.0, 108.883],
        None,
        100.0,
        20.0,
        ModelCamSurround::Average,
        Some(1.0),
        None,
    )
    .unwrap();
    let cam16 = xyz.cam16_forward(cam16_conditions).unwrap();
    let cam16_scalar = lux_rs::cam16_forward(XYZ_SAMPLE, cam16_conditions).unwrap();
    assert_eq!(cam16[0], cam16_scalar);
    let cam16_ucs = xyz
        .cam16_ucs_forward(cam16_conditions, CamUcsType::Ucs)
        .unwrap();
    let cam16_ucs_scalar =
        lux_rs::cam16_ucs_forward(XYZ_SAMPLE, cam16_conditions, CamUcsType::Ucs).unwrap();
    assert_eq!(cam16_ucs[0], cam16_ucs_scalar);
    let cam16_back = Tristimulus::new(vec![[
        cam16_ucs[0].j_prime,
        cam16_ucs[0].a_prime,
        cam16_ucs[0].b_prime,
    ]])
    .cam16_ucs_inverse(cam16_conditions, CamUcsType::Ucs)
    .unwrap();
    assert!((cam16_back.values()[0][0] - 0.25).abs() < 1e-10);
    assert!((cam16_back.values()[0][1] - 0.5).abs() < 1e-10);
    assert!((cam16_back.values()[0][2] - 0.25).abs() < 1e-10);
}

#[test]
fn tristimulus_set_wrapper_matches_batch_transforms() {
    let xyz = Tristimulus::new(vec![XYZ_SAMPLE, XYZ_SAMPLE_ALT]);
    assert_eq!(
        xyz.xyz_to_yxy().values(),
        vec![xyz_to_yxy([0.25, 0.5, 0.25]), xyz_to_yxy([0.2, 0.3, 0.4])]
    );
    assert_eq!(
        xyz.xyz_to_lab([0.5, 0.5, 0.5]).values(),
        vec![
            xyz_to_lab([0.25, 0.5, 0.25], [0.5, 0.5, 0.5]),
            xyz_to_lab([0.2, 0.3, 0.4], [0.5, 0.5, 0.5])
        ]
    );
    assert_eq!(
        xyz.xyz_to_lms(Observer::Cie1931_2).unwrap().values(),
        vec![
            xyz_to_lms([0.25, 0.5, 0.25], Observer::Cie1931_2).unwrap(),
            xyz_to_lms([0.2, 0.3, 0.4], Observer::Cie1931_2).unwrap()
        ]
    );
    let ciecam02_conditions = ciecam02_viewing_conditions(
        [95.047, 100.0, 108.883],
        None,
        100.0,
        20.0,
        ModelCamSurround::Average,
        Some(1.0),
        None,
    )
    .unwrap();
    let cam_many = xyz.ciecam02_forward(ciecam02_conditions).unwrap();
    let cam_scalar = xyz
        .values()
        .iter()
        .copied()
        .map(|value| lux_rs::ciecam02_forward(value, ciecam02_conditions).unwrap())
        .collect::<Vec<_>>();
    assert_eq!(cam_many, cam_scalar);
    let cam_ucs_many = xyz
        .ciecam02_ucs_forward(ciecam02_conditions, CamUcsType::Ucs)
        .unwrap();
    let cam_ucs_scalar = xyz
        .values()
        .iter()
        .copied()
        .map(|value| {
            lux_rs::ciecam02_ucs_forward(value, ciecam02_conditions, CamUcsType::Ucs).unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(cam_ucs_many, cam_ucs_scalar);
    let ucs_triplets = Tristimulus::new(
        cam_ucs_many
            .iter()
            .map(|value| [value.j_prime, value.a_prime, value.b_prime])
            .collect(),
    );
    let xyz_back = ucs_triplets
        .ciecam02_ucs_inverse(ciecam02_conditions, CamUcsType::Ucs)
        .unwrap();
    assert!((xyz_back.values()[0][0] - 0.25).abs() < 1e-10);
    assert!((xyz_back.values()[0][1] - 0.5).abs() < 1e-10);
    assert!((xyz_back.values()[0][2] - 0.25).abs() < 1e-10);
}

#[test]
fn one_row_batch_color_transforms_preserve_numeric_baselines() {
    let xyz = TristimulusSet::new(vec![XYZ_SAMPLE]);

    let yxy = xyz.xyz_to_yxy().into_vec();
    assert_eq!(yxy.len(), 1);
    assert!((yxy[0][0] - 0.5).abs() < 1e-12);
    assert!((yxy[0][1] - 0.25).abs() < 1e-12);
    assert!((yxy[0][2] - 0.5).abs() < 1e-12);

    let lab = xyz.xyz_to_lab(WHITE_E).into_vec();
    assert_eq!(lab.len(), 1);
    assert!((lab[0][0] - 100.0).abs() < 1e-12);
    assert!((lab[0][1] + 103.149_737_007_950_17).abs() < 1e-9);
    assert!((lab[0][2] - 41.259_894_803_180_07).abs() < 1e-9);

    let lms = xyz.xyz_to_lms(Observer::Cie1931_2).unwrap().into_vec();
    assert_eq!(lms.len(), 1);
    assert!((lms[0][0] - 0.422_247_5).abs() < 1e-12);
    assert!((lms[0][1] - 0.545_850_000_000_000_1).abs() < 1e-12);
    assert!((lms[0][2] - 0.25).abs() < 1e-12);
}

#[test]
fn one_row_batch_delta_e_preserves_numeric_baseline() {
    let left = TristimulusSet::new(vec![lab_to_xyz([50.0, 2.6772, -79.7751], WHITE_D65)]);
    let right = TristimulusSet::new(vec![lab_to_xyz([50.0, 0.0, -82.7485], WHITE_D65)]);

    let result = left
        .delta_e(&right, WHITE_D65, DeltaEFormula::Ciede2000)
        .unwrap();

    assert_eq!(result.len(), 1);
    assert!((result[0] - 2.042_459_680_156_574).abs() < 1e-12);
}

#[test]
fn tristimulus_set_delta_e_matches_pairwise_scalar_computation() {
    let left = Tristimulus::new(vec![
        lab_to_xyz([50.0, 2.5, -80.0], WHITE_D65),
        lab_to_xyz([50.0, 2.6772, -79.7751], WHITE_D65),
    ]);
    let right = Tristimulus::new(vec![
        lab_to_xyz([50.0, 0.0, -82.5], WHITE_D65),
        lab_to_xyz([50.0, 0.0, -82.7485], WHITE_D65),
    ]);
    let result = left
        .delta_e(&right, WHITE_D65, DeltaEFormula::Ciede2000)
        .unwrap();
    assert_eq!(result.len(), 2);
    assert!((result[1] - 2.042_459_680_156_574).abs() < 1e-12);
}

#[test]
fn converts_xyz_to_lms_for_1931() {
    let lms = xyz_to_lms([0.25, 0.5, 0.25], Observer::Cie1931_2).unwrap();
    assert!((lms[0] - 0.422_247_5).abs() < 1e-12);
    assert!((lms[1] - 0.545_850_000_000_000_1).abs() < 1e-12);
    assert!((lms[2] - 0.25).abs() < 1e-12);
}

#[test]
fn converts_lms_to_xyz_for_1931() {
    let xyz = lms_to_xyz(
        [0.422_247_5, 0.545_850_000_000_000_1, 0.25],
        Observer::Cie1931_2,
    )
    .unwrap();
    assert!((xyz[0] - 0.25).abs() < 1e-12);
    assert!((xyz[1] - 0.5).abs() < 1e-12);
    assert!((xyz[2] - 0.25).abs() < 1e-12);
}

#[test]
fn converts_xyz_to_lms_for_1964() {
    let lms = xyz_to_lms([0.25, 0.5, 0.25], Observer::Cie1964_10).unwrap();
    assert!((lms[0] - 0.461_241_798_178_5).abs() < 1e-12);
    assert!((lms[1] - 0.516_002_575_170_388).abs() < 1e-12);
    assert!((lms[2] - 0.116_448_084_684_028_24).abs() < 1e-12);
}

#[test]
fn converts_xyz_to_lms_and_back_for_2006_observers() {
    for observer in [Observer::Cie2006_2, Observer::Cie2006_10] {
        let lms = xyz_to_lms(XYZ_SAMPLE, observer).unwrap();
        assert!(lms.iter().all(|value| value.is_finite()));
        let xyz = lms_to_xyz(lms, observer).unwrap();
        assert!((xyz[0] - XYZ_SAMPLE[0]).abs() < 1e-12);
        assert!((xyz[1] - XYZ_SAMPLE[1]).abs() < 1e-12);
        assert!((xyz[2] - XYZ_SAMPLE[2]).abs() < 1e-12);
    }
}

#[test]
fn converts_xyz_to_lms_and_back_for_2015_observers() {
    for observer in [Observer::Cie2015_2, Observer::Cie2015_10] {
        let lms = xyz_to_lms(XYZ_SAMPLE, observer).unwrap();
        assert!(lms.iter().all(|value| value.is_finite()));
        let xyz = lms_to_xyz(lms, observer).unwrap();
        assert!((xyz[0] - XYZ_SAMPLE[0]).abs() < 1e-12);
        assert!((xyz[1] - XYZ_SAMPLE[1]).abs() < 1e-12);
        assert!((xyz[2] - XYZ_SAMPLE[2]).abs() < 1e-12);
    }
}

#[test]
fn exposes_xyz_to_lms_matrix_for_1931() {
    let matrix = Observer::Cie1931_2.xyz_to_lms_matrix().unwrap();
    assert_eq!(matrix[0], [0.38971, 0.68898, -0.07868]);
    assert_eq!(matrix[2], [0.0, 0.0, 1.0]);
}

#[test]
fn converts_xyz_to_srgb() {
    let rgb = xyz_to_srgb([20.0, 21.0, 22.0], 2.4, -0.055, true);
    assert!((rgb[0] - 127.932_633_053_083_4).abs() < 1e-9);
    assert!((rgb[1] - 126.171_697_951_843_17).abs() < 1e-9);
    assert!((rgb[2] - 123.804_791_369_705).abs() < 1e-9);
}

#[test]
fn converts_srgb_to_xyz() {
    let xyz = srgb_to_xyz([64.0, 128.0, 192.0], 2.4, -0.055, true);
    assert!((xyz[0] - 19.344_430_750_022_802).abs() < 1e-9);
    assert!((xyz[1] - 20.332_127_014_120_942).abs() < 1e-9);
    assert!((xyz[2] - 52.763_974_844_108_34).abs() < 1e-9);
}
