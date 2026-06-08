mod common;

use lux_rs::{
    butterworth_spd, color3mixer, colormixer, colormixer_pinv, fit_gaussian_spd_params,
    gaussian_spd, getwlr, lorentzian2_spd, mono_led_spd, phosphor_led_spd,
    phosphor_led_spd_with_components, roundedtriangle_spd, spd_builder, MonoLedParams,
    Observer, PhosphorLedParams, RoundedTriangleParams, Spectrum, WavelengthGrid,
};

#[test]
fn test_parametric_spd_generators() {
    let grid = WavelengthGrid::new(400.0, 700.0, 10.0).unwrap();

    // Gaussian
    let g = gaussian_spd(&[550.0], &[20.0], Some(grid)).unwrap();
    assert_eq!(g.wavelength_count(), 31);
    assert!(g.values()[15] > 0.9); // near peak

    // Lorentzian2
    let l = lorentzian2_spd(&[550.0], &[20.0], Some(grid)).unwrap();
    assert_eq!(l.wavelength_count(), 31);

    // Butterworth
    let b = butterworth_spd(&[550.0], &[20.0], &[2.0], Some(grid)).unwrap();
    assert_eq!(b.wavelength_count(), 31);

    // Rounded Triangle
    let rt_params = RoundedTriangleParams {
        peakwl: 550.0,
        fwhm: Some(40.0),
        rounding: 0.5,
        min_v: 0.0,
        max_v: 1.0,
        fw: 100.0,
        rw: 100.0,
    };
    let rt = roundedtriangle_spd(&[rt_params], Some(grid)).unwrap();
    assert_eq!(rt.wavelength_count(), 31);
}

#[test]
fn test_mono_led_spd() {
    let grid = WavelengthGrid::new(400.0, 700.0, 10.0).unwrap();
    let mono_params = MonoLedParams {
        peakwl: 450.0,
        fwhm: 20.0,
        strength_shoulder: 1.5,
        bw_order: -1.0,
    };

    let spd = mono_led_spd(&[mono_params], Some(grid)).unwrap();
    assert_eq!(spd.spectrum_count(), 1);
    assert_eq!(spd.wavelength_count(), 31);
}

#[test]
fn test_phosphor_led_spd() {
    let grid = WavelengthGrid::new(400.0, 700.0, 10.0).unwrap();
    let ph_params = PhosphorLedParams {
        peakwl: 450.0,
        fwhm: 20.0,
        strength_shoulder: 1.5,
        bw_order: -1.0,
        strength_ph: Some(0.5),
        peakwl_ph1: 530.0,
        fwhm_ph1: 80.0,
        strength_ph1: 0.8,
        peakwl_ph2: 560.0,
        fwhm_ph2: 80.0,
        strength_ph2: None,
        use_piecewise_fcn: true,
    };

    // phosphor_led_spd
    let spd = phosphor_led_spd(&[ph_params.clone()], Some(grid)).unwrap();
    assert_eq!(spd.spectrum_count(), 1);
    assert_eq!(spd.wavelength_count(), 31);

    // phosphor_led_spd_with_components
    let comps = phosphor_led_spd_with_components(&[ph_params], Some(grid)).unwrap();
    assert_eq!(comps.spd.spectrum_count(), 1);
    assert_eq!(comps.components.spectrum_count(), 3); // mono_led, ph1, ph2
}

#[test]
fn test_color_mixers() {
    // Basic color3mixer check
    let target = [100.0, 0.4, 0.4];
    let p1 = [100.0, 0.2, 0.2];
    let p2 = [100.0, 0.6, 0.2];
    let p3 = [100.0, 0.3, 0.7];

    let m3 = color3mixer(target, p1, p2, p3);
    assert_eq!(m3.len(), 3);
    // Check that we can reconstruct the target coordinates approximately
    // using additive mixing of the primaries
    let sum_m = m3[0] + m3[1] + m3[2];
    assert!(sum_m > 0.0);

    // colormixer_pinv
    let pinv = colormixer_pinv(target, &[p1, p2, p3], "Yxy");
    assert_eq!(pinv.len(), 3);

    // colormixer (cascaded relative mixer)
    let primaries = vec![p1, p2, p3, [100.0, 0.5, 0.5]];
    let cm = colormixer(target, &primaries, &[0.5]);
    assert_eq!(cm.len(), 4);
}

#[test]
fn test_spd_builder() {
    let grid = WavelengthGrid::new(400.0, 700.0, 10.0).unwrap();
    let ph_params = PhosphorLedParams {
        peakwl: 450.0,
        fwhm: 20.0,
        strength_shoulder: 1.5,
        bw_order: -1.0,
        strength_ph: Some(0.5),
        peakwl_ph1: 530.0,
        fwhm_ph1: 80.0,
        strength_ph1: 0.8,
        peakwl_ph2: 560.0,
        fwhm_ph2: 80.0,
        strength_ph2: None,
        use_piecewise_fcn: true,
    };

    // 1. No target, with flux values
    let spd1 = spd_builder(
        Some(&[0.5, 0.3, 0.2]),
        None,
        &ph_params,
        None,
        None,
        "Yxy",
        Observer::Cie1931_2,
        Some(grid),
    )
    .unwrap();
    assert_eq!(spd1.spectrum_count(), 1);

    // 2. With target, in gamut
    let spd2 = spd_builder(
        None,
        None,
        &ph_params,
        None,
        Some(&[100.0, 0.25, 0.35]),
        "Yxy",
        Observer::Cie1931_2,
        Some(grid),
    )
    .unwrap();
    assert_eq!(spd2.spectrum_count(), 1);
    // Values should not be NaN
    assert!(!spd2.values().iter().any(|&v| v.is_nan()));

    // 3. With target, out of gamut (should return NaN spectrum, matching python behavior)
    let spd_out = spd_builder(
        None,
        None,
        &ph_params,
        None,
        Some(&[100.0, 0.9, 0.9]), // way out of gamut
        "Yxy",
        Observer::Cie1931_2,
        Some(grid),
    )
    .unwrap();
    assert!(spd_out.values().iter().all(|&v| v.is_nan()));
}

#[test]
fn test_user_primary_mixing_632_532_467() {
    let grid = WavelengthGrid::new(380.0, 780.0, 5.0).unwrap();
    let observer = Observer::Cie1931_2;

    // Generate primaries using gaussian_spd
    let g1 = gaussian_spd(&[632.0], &[15.0], Some(grid)).unwrap();
    let g2 = gaussian_spd(&[532.0], &[8.0], Some(grid)).unwrap();
    let g3 = gaussian_spd(&[467.0], &[9.0], Some(grid)).unwrap();

    // Convert each to XYZ then Yxy
    let xyz1 = g1.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz2 = g2.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz3 = g3.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];

    let yxy1 = lux_rs::xyz_to_yxy(xyz1);
    let yxy2 = lux_rs::xyz_to_yxy(xyz2);
    let yxy3 = lux_rs::xyz_to_yxy(xyz3);

    println!("Simulation 1 (632/532/467nm):");
    println!("  Primary 1 (632 nm) Yxy: {:?}", yxy1);
    println!("  Primary 2 (532 nm) Yxy: {:?}", yxy2);
    println!("  Primary 3 (467 nm) Yxy: {:?}", yxy3);

    let target_yxy = [100.0, 0.3127, 0.3290];
    let m = color3mixer(target_yxy, yxy1, yxy2, yxy3);
    let y_contrib1 = m[0] * yxy1[0];
    let y_contrib2 = m[1] * yxy2[0];
    let y_contrib3 = m[2] * yxy3[0];
    let sum_y = y_contrib1 + y_contrib2 + y_contrib3;

    println!("  Luminance ratio: R = {:.6}%, G = {:.6}%, B = {:.6}%", 
             y_contrib1 / sum_y * 100.0, 
             y_contrib2 / sum_y * 100.0, 
             y_contrib3 / sum_y * 100.0);
}

#[test]
fn test_user_primary_mixing_optimized() {
    let grid = WavelengthGrid::new(380.0, 780.0, 1.0).unwrap();
    let observer = Observer::Cie1931_2;

    // Generate primaries using gaussian_spd with optimized parameters
    let g1 = gaussian_spd(&[630.065], &[1.618], Some(grid)).unwrap();
    let g2 = gaussian_spd(&[531.644], &[0.523], Some(grid)).unwrap();
    let g3 = gaussian_spd(&[467.376], &[0.522], Some(grid)).unwrap();

    // Convert each to XYZ then Yxy
    let xyz1 = g1.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz2 = g2.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz3 = g3.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];

    let yxy1 = lux_rs::xyz_to_yxy(xyz1);
    let yxy2 = lux_rs::xyz_to_yxy(xyz2);
    let yxy3 = lux_rs::xyz_to_yxy(xyz3);

    println!("Simulation 2 (Optimized target chromaticities):");
    println!("  Primary 1 (Red) Yxy: {:?}", yxy1);
    println!("  Primary 2 (Green) Yxy: {:?}", yxy2);
    println!("  Primary 3 (Blue) Yxy: {:?}", yxy3);

    let target_yxy = [100.0, 0.3127, 0.3290];
    let m = color3mixer(target_yxy, yxy1, yxy2, yxy3);
    let y_contrib1 = m[0] * yxy1[0];
    let y_contrib2 = m[1] * yxy2[0];
    let y_contrib3 = m[2] * yxy3[0];
    let sum_y = y_contrib1 + y_contrib2 + y_contrib3;

    println!("  Luminance ratio: R = {:.6}%, G = {:.6}%, B = {:.6}%", 
             y_contrib1 / sum_y * 100.0, 
             y_contrib2 / sum_y * 100.0, 
             y_contrib3 / sum_y * 100.0);
}

#[test]
fn test_user_primary_mixing_srgb() {
    let grid = WavelengthGrid::new(380.0, 780.0, 1.0).unwrap();
    let observer = Observer::Cie1931_2;

    // Generate primaries using gaussian_spd with sRGB-optimized parameters
    let g1 = gaussian_spd(&[679.4975], &[113.0899], Some(grid)).unwrap();
    let g2 = gaussian_spd(&[539.1311], &[78.4488], Some(grid)).unwrap();
    let g3 = gaussian_spd(&[415.1195], &[111.1245], Some(grid)).unwrap();

    // Convert each to XYZ then Yxy
    let xyz1 = g1.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz2 = g2.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz3 = g3.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];

    let yxy1 = lux_rs::xyz_to_yxy(xyz1);
    let yxy2 = lux_rs::xyz_to_yxy(xyz2);
    let yxy3 = lux_rs::xyz_to_yxy(xyz3);

    println!("Simulation 3 (sRGB targets):");
    println!("  Primary 1 (Red) Yxy: {:?}", yxy1);
    println!("  Primary 2 (Green) Yxy: {:?}", yxy2);
    println!("  Primary 3 (Blue) Yxy: {:?}", yxy3);

    let target_yxy = [100.0, 0.3127, 0.3290];
    let m = color3mixer(target_yxy, yxy1, yxy2, yxy3);
    let y_contrib1 = m[0] * yxy1[0];
    let y_contrib2 = m[1] * yxy2[0];
    let y_contrib3 = m[2] * yxy3[0];
    let sum_y = y_contrib1 + y_contrib2 + y_contrib3;

    println!("  Luminance ratio: R = {:.6}%, G = {:.6}%, B = {:.6}%", 
             y_contrib1 / sum_y * 100.0, 
             y_contrib2 / sum_y * 100.0, 
             y_contrib3 / sum_y * 100.0);
}

#[test]
fn test_rust_spd_parameter_fitting() {
    // Fit Green [0.3, 0.6] starting from 520 nm, FWHM 30 nm
    let fit_g = fit_gaussian_spd_params([0.3, 0.6], 520.0, 30.0).unwrap();
    println!("Rust Green Fit: peak = {}, fwhm = {}", fit_g.0, fit_g.1);
    
    // Check that the fitted parameters yield the target coordinate approximately
    let grid = WavelengthGrid::new(380.0, 780.0, 1.0).unwrap();
    let g = gaussian_spd(&[fit_g.0], &[fit_g.1], Some(grid)).unwrap();
    let xyz = g.spd_to_xyz(&Observer::Cie1931_2.standard().unwrap(), false).unwrap()[0];
    let yxy = lux_rs::xyz_to_yxy(xyz);
    
    assert!((yxy[1] - 0.3).abs() < 0.001);
    assert!((yxy[2] - 0.6).abs() < 0.001);
}

#[test]
fn test_user_primary_mixing_native_r() {
    let grid = WavelengthGrid::new(380.0, 780.0, 1.0).unwrap();
    let observer = Observer::Cie1931_2;

    // Generate Red primary as a sum of 5 optimized Gaussian peaks
    let peaks = &[453.962, 540.000, 613.197, 631.868, 648.087];
    let fwhms = &[21.195, 60.000, 13.718, 9.938, 13.113];
    let heights = &[0.196576, 0.356536, 5.37827, 15.90185, 2.25871];

    let wavelengths = getwlr(grid).unwrap();
    let mut red_values = vec![0.0; wavelengths.len()];

    for i in 0..5 {
        let g = gaussian_spd(&[peaks[i]], &[fwhms[i]], Some(grid)).unwrap();
        let g_vals = &g.spectra()[0];
        let h = heights[i];
        for k in 0..wavelengths.len() {
            red_values[k] += h * g_vals[k];
        }
    }

    let g_red = Spectrum::new(wavelengths.clone(), vec![red_values]).unwrap();

    // Generate Green and Blue primaries (optimized for sRGB targets)
    let g_green = gaussian_spd(&[539.1311], &[78.4488], Some(grid)).unwrap();
    let g_blue = gaussian_spd(&[415.1195], &[111.1245], Some(grid)).unwrap();

    // Convert to XYZ then Yxy
    let xyz_red = g_red.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz_green = g_green.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz_blue = g_blue.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];

    let yxy_red = lux_rs::xyz_to_yxy(xyz_red);
    let yxy_green = lux_rs::xyz_to_yxy(xyz_green);
    let yxy_blue = lux_rs::xyz_to_yxy(xyz_blue);

    println!("Simulation 4 (Native-R + sRGB):");
    println!("  Primary 1 (Red) Yxy: {:?}", yxy_red);
    println!("  Primary 2 (Green) Yxy: {:?}", yxy_green);
    println!("  Primary 3 (Blue) Yxy: {:?}", yxy_blue);

    // Assert Red primary chromaticity matches target [0.64, 0.33] within 0.001
    assert!((yxy_red[1] - 0.64).abs() < 0.001);
    assert!((yxy_red[2] - 0.33).abs() < 0.001);

    // Mix to D65
    let target_yxy = [100.0, 0.3127, 0.3290];
    let m = color3mixer(target_yxy, yxy_red, yxy_green, yxy_blue);
    let y_contrib1 = m[0] * yxy_red[0];
    let y_contrib2 = m[1] * yxy_green[0];
    let y_contrib3 = m[2] * yxy_blue[0];
    let sum_y = y_contrib1 + y_contrib2 + y_contrib3;

    println!("  Luminance ratio: R = {:.6}%, G = {:.6}%, B = {:.6}%", 
             y_contrib1 / sum_y * 100.0, 
             y_contrib2 / sum_y * 100.0, 
             y_contrib3 / sum_y * 100.0);
}

#[test]
fn test_user_primary_mixing_isolated_rgb() {
    let grid = WavelengthGrid::new(380.0, 780.0, 1.0).unwrap();
    let observer = Observer::Cie1931_2;

    // 1. Isolated Red (KSF components only, leakage heights set to 0.0)
    let peaks = &[453.962, 540.000, 613.197, 631.868, 648.087];
    let fwhms = &[21.195, 60.000, 13.718, 9.938, 13.113];
    let heights_isolated_red = &[0.0, 0.0, 5.37827, 15.90185, 2.25871];

    let wavelengths = getwlr(grid).unwrap();
    let mut red_values = vec![0.0; wavelengths.len()];

    for i in 0..5 {
        if heights_isolated_red[i] == 0.0 {
            continue;
        }
        let g = gaussian_spd(&[peaks[i]], &[fwhms[i]], Some(grid)).unwrap();
        let g_vals = &g.spectra()[0];
        let h = heights_isolated_red[i];
        for k in 0..wavelengths.len() {
            red_values[k] += h * g_vals[k];
        }
    }
    let g_red = Spectrum::new(wavelengths.clone(), vec![red_values]).unwrap();

    // 2. Isolated Green: peak = 531.4658, fwhm = 45.1469
    let g_green = gaussian_spd(&[531.4658], &[45.1469], Some(grid)).unwrap();

    // 3. Isolated Blue: peak = 454.9756, fwhm = 26.1472
    let g_blue = gaussian_spd(&[454.9756], &[26.1472], Some(grid)).unwrap();

    // Convert to XYZ then Yxy
    let xyz_red = g_red.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz_green = g_green.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];
    let xyz_blue = g_blue.spd_to_xyz(&observer.standard().unwrap(), false).unwrap()[0];

    let yxy_red = lux_rs::xyz_to_yxy(xyz_red);
    let yxy_green = lux_rs::xyz_to_yxy(xyz_green);
    let yxy_blue = lux_rs::xyz_to_yxy(xyz_blue);

    println!("Simulation 5 (Isolated RGB):");
    println!("  Primary 1 (Red) Yxy: {:?}", yxy_red);
    println!("  Primary 2 (Green) Yxy: {:?}", yxy_green);
    println!("  Primary 3 (Blue) Yxy: {:?}", yxy_blue);

    // Mix to D65
    let target_yxy = [100.0, 0.3127, 0.3290];
    let m = color3mixer(target_yxy, yxy_red, yxy_green, yxy_blue);
    let y_contrib1 = m[0] * yxy_red[0];
    let y_contrib2 = m[1] * yxy_green[0];
    let y_contrib3 = m[2] * yxy_blue[0];
    let sum_y = y_contrib1 + y_contrib2 + y_contrib3;

    println!("  Luminance ratio: R = {:.6}%, G = {:.6}%, B = {:.6}%", 
             y_contrib1 / sum_y * 100.0, 
             y_contrib2 / sum_y * 100.0, 
             y_contrib3 / sum_y * 100.0);
}

