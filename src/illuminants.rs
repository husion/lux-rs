use crate::color::Observer;
use crate::error::{LuxError, LuxResult};
use crate::photometry::spd_to_xyz;
use crate::spectrum::{getwlr, SpectralMatrix, Spectrum, WavelengthGrid};

const DEFAULT_WAVELENGTH_GRID: WavelengthGrid = WavelengthGrid {
    start: 360.0,
    end: 830.0,
    step: 1.0,
};
const BLACKBODY_C1: f64 = 3.741_771_85e-16;
const BLACKBODY_C2: f64 = 1.438_8e-2;
const DEFAULT_REFRACTIVE_INDEX: f64 = 1.0;
const DAYLIGHT_WAVELENGTHS_10NM: [f64; 48] = [
    360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0, 430.0, 440.0, 450.0, 460.0, 470.0, 480.0,
    490.0, 500.0, 510.0, 520.0, 530.0, 540.0, 550.0, 560.0, 570.0, 580.0, 590.0, 600.0, 610.0,
    620.0, 630.0, 640.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, 710.0, 720.0, 730.0, 740.0,
    750.0, 760.0, 770.0, 780.0, 790.0, 800.0, 810.0, 820.0, 830.0,
];
const DAYLIGHT_S0: [f64; 48] = [
    61.5, 68.8, 63.4, 65.8, 94.8, 104.8, 105.9, 96.8, 113.9, 125.6, 125.5, 121.3, 121.3, 113.5,
    113.1, 110.8, 106.5, 108.8, 105.3, 104.4, 100.0, 96.0, 95.1, 89.1, 90.5, 90.3, 88.4, 84.0,
    85.1, 81.9, 82.6, 84.9, 81.3, 71.9, 74.3, 76.4, 63.3, 71.7, 77.0, 65.2, 47.7, 68.6, 65.0, 66.0,
    61.0, 53.3, 58.9, 61.9,
];
const DAYLIGHT_S1: [f64; 48] = [
    38.0, 42.4, 38.5, 35.0, 43.4, 46.3, 43.9, 37.1, 36.7, 35.9, 32.6, 27.9, 24.3, 20.1, 16.2, 13.2,
    8.6, 6.1, 4.2, 1.9, 0.0, -1.6, -3.5, -3.5, -5.8, -7.2, -8.6, -9.5, -10.9, -10.7, -12.0, -14.0,
    -13.6, -12.0, -13.3, -12.9, -10.6, -11.6, -12.2, -10.2, -7.8, -11.2, -10.4, -10.6, -9.7, -8.3,
    -9.3, -9.8,
];
const DAYLIGHT_S2: [f64; 48] = [
    5.3, 6.1, 3.0, 1.2, -1.1, -0.5, -0.7, -1.2, -2.6, -2.9, -2.8, -2.6, -2.6, -1.8, -1.5, -1.3,
    -1.2, -1.0, -0.5, -0.3, 0.0, 0.2, 0.5, 2.1, 3.2, 4.1, 4.7, 5.1, 6.7, 7.3, 8.6, 9.8, 10.2, 8.3,
    9.6, 8.5, 7.0, 7.6, 8.0, 6.7, 5.2, 7.4, 6.8, 7.0, 6.4, 5.5, 6.1, 6.5,
];
const DAYLIGHT_I: f64 = 0.0241;
const DAYLIGHT_J: f64 = 0.2562;
const DAYLIGHT_K: f64 = -0.7341;
const DAYLIGHT_I1: f64 = -1.3515;
const DAYLIGHT_J1: f64 = -1.7703;
const DAYLIGHT_K1: f64 = 5.9114;
const DAYLIGHT_I2: f64 = 0.03;
const DAYLIGHT_J2: f64 = -31.4424;
const DAYLIGHT_K2: f64 = 30.0717;
const CIERA_MIX_RANGE: [f64; 2] = [5000.0, 5000.0];
const CCT_SEARCH_MIN: f64 = 1000.0;
const CCT_SEARCH_MAX: f64 = 100_000.0;
const CCT_SEARCH_SAMPLES: usize = 256;
const CIE_A_CSV: &str = include_str!("../data/spds/CIE_A.csv");
const CIE_D65_CSV: &str = include_str!("../data/spds/CIE_D65.csv");
const CIE_F_SERIES_CSV: &str = include_str!("../data/spds/CIE_F_1to12_1nm.csv");
const CIE_LED_SERIES_CSV: &str = include_str!("../data/spds/CIE_LED_B1toB5_BH1_RGB1_V1_V2.csv");
const LED_NAMES: [&str; 9] = [
    "LED_B1", "LED_B2", "LED_B3", "LED_B4", "LED_B5", "LED_BH1", "LED_RGB1", "LED_V1", "LED_V2",
];

pub fn standard_illuminant_names() -> &'static [&'static str] {
    &[
        "A", "D50", "D55", "D65", "D75", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9",
        "F10", "F11", "F12", "LED_B1", "LED_B2", "LED_B3", "LED_B4", "LED_B5", "LED_BH1",
        "LED_RGB1", "LED_V1", "LED_V2",
    ]
}

pub fn blackbody(
    cct: f64,
    wavelength_grid: Option<WavelengthGrid>,
    refractive_index: Option<f64>,
    relative: bool,
) -> LuxResult<Spectrum> {
    if !cct.is_finite() || cct <= 0.0 {
        return Err(LuxError::InvalidInput(
            "cct must be finite and strictly positive",
        ));
    }

    let n = refractive_index.unwrap_or(DEFAULT_REFRACTIVE_INDEX);
    if !n.is_finite() || n <= 0.0 {
        return Err(LuxError::InvalidInput(
            "refractive index must be finite and strictly positive",
        ));
    }

    let wavelengths = getwlr(wavelength_grid.unwrap_or(DEFAULT_WAVELENGTH_GRID))?;
    let normalization = if relative {
        planck_spectral_radiance(560.0, cct, n)
    } else {
        1.0
    };

    let values = wavelengths
        .iter()
        .map(|&wavelength| planck_spectral_radiance(wavelength, cct, n) / normalization)
        .collect();

    Spectrum::new(wavelengths, values)
}

pub fn daylightlocus(
    cct: f64,
    force_daylight_below_4000k: bool,
    cct_is_nominal: bool,
) -> LuxResult<[f64; 2]> {
    if !cct.is_finite() || cct <= 0.0 {
        return Err(LuxError::InvalidInput(
            "cct must be finite and strictly positive",
        ));
    }
    if cct < 4000.0 && !force_daylight_below_4000k {
        return Err(LuxError::InvalidInput(
            "daylight locus is undefined below 4000 K unless forced",
        ));
    }

    let cct = if cct_is_nominal {
        cct * (1.4388 / 1.4380)
    } else {
        cct
    };

    let x_d = if cct < 7000.0 {
        -4.607 * (1e3 / cct).powi(3)
            + 2.9678 * (1e3 / cct).powi(2)
            + 0.09911 * (1e3 / cct)
            + 0.244_063
    } else {
        -2.0064 * (1e3 / cct).powi(3)
            + 1.9018 * (1e3 / cct).powi(2)
            + 0.24748 * (1e3 / cct)
            + 0.23704
    };
    let y_d = -3.0 * x_d.powi(2) + 2.87 * x_d - 0.275;

    Ok([x_d, y_d])
}

pub fn daylightphase(
    cct: f64,
    wavelength_grid: Option<WavelengthGrid>,
    cct_is_nominal: bool,
    force_daylight_below_4000k: bool,
    refractive_index: Option<f64>,
) -> LuxResult<Spectrum> {
    let grid = wavelength_grid.unwrap_or(DEFAULT_WAVELENGTH_GRID);
    if cct < 4000.0 && !force_daylight_below_4000k {
        return blackbody(cct, Some(grid), refractive_index, true);
    }

    let wavelengths = getwlr(grid)?;
    let [x_d, y_d] = daylightlocus(cct, force_daylight_below_4000k, cct_is_nominal)?;
    let denominator = DAYLIGHT_I + DAYLIGHT_J * x_d + DAYLIGHT_K * y_d;
    let mut m1 = (DAYLIGHT_I1 + DAYLIGHT_J1 * x_d + DAYLIGHT_K1 * y_d) / denominator;
    let mut m2 = (DAYLIGHT_I2 + DAYLIGHT_J2 * x_d + DAYLIGHT_K2 * y_d) / denominator;
    if cct_is_nominal {
        m1 = (m1 * 1_000.0).round() / 1_000.0;
        m2 = (m2 * 1_000.0).round() / 1_000.0;
    }

    let s0 = interpolate_daylight_component(&DAYLIGHT_S0, &wavelengths);
    let s1 = interpolate_daylight_component(&DAYLIGHT_S1, &wavelengths);
    let s2 = interpolate_daylight_component(&DAYLIGHT_S2, &wavelengths);

    let mut values: Vec<f64> = s0
        .iter()
        .zip(s1.iter())
        .zip(s2.iter())
        .map(|((s0, s1), s2)| s0 + m1 * s1 + m2 * s2)
        .collect();

    let idx_560 = wavelengths
        .iter()
        .enumerate()
        .min_by(|(_, lhs), (_, rhs)| (*lhs - 560.0).abs().total_cmp(&(*rhs - 560.0).abs()))
        .map(|(index, _)| index)
        .ok_or(LuxError::EmptyInput)?;
    let normalization = values[idx_560];
    for value in &mut values {
        *value /= normalization;
    }

    Spectrum::new(wavelengths, values)
}

pub fn cri_ref(ccts: &[f64], wavelength_grid: Option<WavelengthGrid>) -> LuxResult<SpectralMatrix> {
    if ccts.is_empty() {
        return Err(LuxError::EmptyInput);
    }

    let grid = wavelength_grid.unwrap_or(DEFAULT_WAVELENGTH_GRID);
    let wavelengths = getwlr(grid)?;
    let mut spectra = Vec::with_capacity(ccts.len());

    for &cct in ccts {
        let spectrum = if cct < CIERA_MIX_RANGE[0] {
            blackbody(cct, Some(grid), None, true)?
        } else {
            daylightphase(cct, Some(grid), false, false, None)?
        };
        spectra.push(spectrum.values().to_vec());
    }

    SpectralMatrix::new(wavelengths, spectra)
}

pub fn standard_illuminant(
    name: &str,
    wavelength_grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    let key = name.trim().to_ascii_uppercase();
    match key.as_str() {
        "A" => load_single_illuminant_csv(CIE_A_CSV, wavelength_grid),
        "D50" => daylightphase(5000.0, wavelength_grid, true, false, None),
        "D55" => daylightphase(5500.0, wavelength_grid, true, false, None),
        "D65" => load_single_illuminant_csv(CIE_D65_CSV, wavelength_grid),
        "D75" => daylightphase(7500.0, wavelength_grid, true, false, None),
        _ if key.starts_with('F') => {
            let index = key[1..]
                .parse::<usize>()
                .map_err(|_| LuxError::UnsupportedObserver("unknown standard illuminant"))?;
            if (1..=12).contains(&index) {
                load_series_illuminant_csv(CIE_F_SERIES_CSV, index, wavelength_grid)
            } else {
                Err(LuxError::UnsupportedObserver("unknown standard illuminant"))
            }
        }
        _ => {
            if let Some(index) = LED_NAMES.iter().position(|candidate| *candidate == key) {
                load_series_illuminant_csv(CIE_LED_SERIES_CSV, index + 1, wavelength_grid)
            } else {
                Err(LuxError::UnsupportedObserver("unknown standard illuminant"))
            }
        }
    }
}

pub fn cct_to_xyz(cct: f64, observer: Observer) -> LuxResult<[f64; 3]> {
    let spectrum = blackbody(cct, None, None, true)?;
    let observer = observer.standard()?;
    spd_to_xyz(&spectrum, &observer, true)
}

pub fn xyz_to_cct(xyz: [f64; 3], observer: Observer) -> LuxResult<(f64, f64)> {
    if xyz.iter().any(|value| !value.is_finite() || *value < 0.0) {
        return Err(LuxError::InvalidInput(
            "xyz values must be finite and non-negative",
        ));
    }

    let target_uv = xyz_to_uv1960(xyz);
    let (mut best_cct, mut best_distance) = coarse_cct_search(target_uv, observer)?;

    let coarse_step = best_cct * 0.05;
    let mut low = (best_cct - coarse_step).max(CCT_SEARCH_MIN);
    let mut high = (best_cct + coarse_step).min(CCT_SEARCH_MAX);
    if low == high {
        low = (best_cct * 0.9).max(CCT_SEARCH_MIN);
        high = (best_cct * 1.1).min(CCT_SEARCH_MAX);
    }

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let inv_phi = 1.0 / phi;
    let mut c = high - (high - low) * inv_phi;
    let mut d = low + (high - low) * inv_phi;
    let mut fc = cct_distance_squared(target_uv, c, observer)?;
    let mut fd = cct_distance_squared(target_uv, d, observer)?;

    for _ in 0..80 {
        if (high - low) < 1e-6 {
            break;
        }
        if fc < fd {
            high = d;
            d = c;
            fd = fc;
            c = high - (high - low) * inv_phi;
            fc = cct_distance_squared(target_uv, c, observer)?;
        } else {
            low = c;
            c = d;
            fc = fd;
            d = low + (high - low) * inv_phi;
            fd = cct_distance_squared(target_uv, d, observer)?;
        }
    }

    let candidates = [(c, fc), (d, fd), (best_cct, best_distance)];
    let (refined_cct, refined_distance) = candidates
        .into_iter()
        .min_by(|lhs, rhs| lhs.1.total_cmp(&rhs.1))
        .ok_or(LuxError::InvalidInput("failed to refine cct search"))?;
    best_cct = refined_cct;
    best_distance = refined_distance;

    let best_uv = cct_to_uv1960(best_cct, observer)?;
    let delta_cct = (best_cct * 1e-4).max(0.1);
    let uv_low = cct_to_uv1960((best_cct - delta_cct).max(CCT_SEARCH_MIN), observer)?;
    let uv_high = cct_to_uv1960((best_cct + delta_cct).min(CCT_SEARCH_MAX), observer)?;
    let tangent = [uv_high[0] - uv_low[0], uv_high[1] - uv_low[1]];
    let offset = [target_uv[0] - best_uv[0], target_uv[1] - best_uv[1]];
    let cross = tangent[0] * offset[1] - tangent[1] * offset[0];
    let duv = if cross == 0.0 {
        0.0
    } else {
        -cross.signum() * best_distance.sqrt()
    };

    Ok((best_cct, duv))
}

fn planck_spectral_radiance(wavelength_nm: f64, cct: f64, refractive_index: f64) -> f64 {
    let wavelength_m = wavelength_nm * 1e-9;
    let exponent = BLACKBODY_C2 / (refractive_index * wavelength_m * (cct + f64::EPSILON));
    (1.0 / std::f64::consts::PI)
        * BLACKBODY_C1
        * wavelength_m.powi(-5)
        * refractive_index.powi(-2)
        * (exponent.exp() - 1.0).powi(-1)
}

fn coarse_cct_search(target_uv: [f64; 2], observer: Observer) -> LuxResult<(f64, f64)> {
    let mut best = None;
    let min_ln = CCT_SEARCH_MIN.ln();
    let max_ln = CCT_SEARCH_MAX.ln();

    for index in 0..CCT_SEARCH_SAMPLES {
        let t = index as f64 / (CCT_SEARCH_SAMPLES - 1) as f64;
        let cct = (min_ln + (max_ln - min_ln) * t).exp();
        let distance = cct_distance_squared(target_uv, cct, observer)?;
        match best {
            None => best = Some((cct, distance)),
            Some((_, best_distance)) if distance < best_distance => best = Some((cct, distance)),
            _ => {}
        }
    }

    best.ok_or(LuxError::InvalidInput("failed to search cct"))
}

fn cct_distance_squared(target_uv: [f64; 2], cct: f64, observer: Observer) -> LuxResult<f64> {
    let uv = cct_to_uv1960(cct, observer)?;
    let du = target_uv[0] - uv[0];
    let dv = target_uv[1] - uv[1];
    Ok(du * du + dv * dv)
}

fn cct_to_uv1960(cct: f64, observer: Observer) -> LuxResult<[f64; 2]> {
    Ok(xyz_to_uv1960(cct_to_xyz(cct, observer)?))
}

fn xyz_to_uv1960(xyz: [f64; 3]) -> [f64; 2] {
    let denominator = xyz[0] + 15.0 * xyz[1] + 3.0 * xyz[2];
    if denominator == 0.0 {
        [0.0, 0.0]
    } else {
        [4.0 * xyz[0] / denominator, 6.0 * xyz[1] / denominator]
    }
}

fn load_single_illuminant_csv(
    csv: &str,
    wavelength_grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    let spectrum = parse_single_illuminant_csv(csv)?;
    match wavelength_grid {
        Some(grid) => spectrum.cie_interp_linear(&getwlr(grid)?, true),
        None => Ok(spectrum),
    }
}

fn load_series_illuminant_csv(
    csv: &str,
    column_index: usize,
    wavelength_grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    let spectrum = parse_series_illuminant_csv(csv, column_index)?;
    match wavelength_grid {
        Some(grid) => spectrum.cie_interp_linear(&getwlr(grid)?, true),
        None => Ok(spectrum),
    }
}

fn parse_single_illuminant_csv(csv: &str) -> LuxResult<Spectrum> {
    let mut wavelengths = Vec::new();
    let mut values = Vec::new();

    for line in csv.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let mut parts = trimmed.split(',');
        let wavelength = parts
            .next()
            .ok_or(LuxError::ParseError("missing wavelength"))?
            .trim()
            .parse::<f64>()
            .map_err(|_| LuxError::ParseError("invalid wavelength"))?;
        let value = parts
            .next()
            .ok_or(LuxError::ParseError("missing illuminant value"))?
            .trim()
            .parse::<f64>()
            .map_err(|_| LuxError::ParseError("invalid illuminant value"))?;
        wavelengths.push(wavelength);
        values.push(value);
    }

    Spectrum::new(wavelengths, values)
}

fn parse_series_illuminant_csv(csv: &str, column_index: usize) -> LuxResult<Spectrum> {
    let mut wavelengths = Vec::new();
    let mut values = Vec::new();

    for line in csv.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let parts: Vec<&str> = trimmed.split(',').collect();
        let wavelength = parts
            .first()
            .ok_or(LuxError::ParseError("missing wavelength"))?
            .trim()
            .parse::<f64>()
            .map_err(|_| LuxError::ParseError("invalid wavelength"))?;
        let value = parts
            .get(column_index)
            .ok_or(LuxError::ParseError("missing illuminant series value"))?
            .trim()
            .parse::<f64>()
            .map_err(|_| LuxError::ParseError("invalid illuminant series value"))?;
        wavelengths.push(wavelength);
        values.push(value);
    }

    Spectrum::new(wavelengths, values)
}

fn interpolate_daylight_component(source: &[f64; 48], target_wavelengths: &[f64]) -> Vec<f64> {
    target_wavelengths
        .iter()
        .map(|&target| interpolate_linear(&DAYLIGHT_WAVELENGTHS_10NM, source, target))
        .collect()
}

fn interpolate_linear(wavelengths: &[f64], values: &[f64], target: f64) -> f64 {
    if target <= wavelengths[0] {
        return linear_segment(wavelengths[0], values[0], wavelengths[1], values[1], target);
    }
    if target >= wavelengths[wavelengths.len() - 1] {
        let last = wavelengths.len() - 1;
        return linear_segment(
            wavelengths[last - 1],
            values[last - 1],
            wavelengths[last],
            values[last],
            target,
        );
    }

    let idx = wavelengths.partition_point(|wavelength| *wavelength < target);
    if wavelengths[idx] == target {
        values[idx]
    } else {
        linear_segment(
            wavelengths[idx - 1],
            values[idx - 1],
            wavelengths[idx],
            values[idx],
            target,
        )
    }
}

fn linear_segment(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    y0 + (y1 - y0) * ((x - x0) / (x1 - x0))
}

#[cfg(test)]
mod tests {
    use super::{
        blackbody, cct_to_xyz, cri_ref, daylightlocus, daylightphase, standard_illuminant,
        standard_illuminant_names, xyz_to_cct,
    };
    use crate::color::Observer;
    use crate::spectrum::WavelengthGrid;

    #[test]
    fn computes_relative_blackbody_spectrum() {
        let spectrum = blackbody(
            6500.0,
            Some(WavelengthGrid::new(360.0, 365.0, 1.0).unwrap()),
            None,
            true,
        )
        .unwrap();

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
        let spectrum = daylightphase(
            6500.0,
            Some(WavelengthGrid::new(360.0, 365.0, 1.0).unwrap()),
            false,
            false,
            None,
        )
        .unwrap();

        assert!((spectrum.values()[0] - 0.944_787_686_958_161_2).abs() < 1e-12);
        assert!((spectrum.values()[5] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn daylightphase_uses_blackbody_below_4000k() {
        let spectrum = daylightphase(
            3500.0,
            Some(WavelengthGrid::new(360.0, 365.0, 1.0).unwrap()),
            false,
            false,
            None,
        )
        .unwrap();

        assert!((spectrum.values()[0] - 0.154_169_092_129_982_6).abs() < 1e-12);
        assert!((spectrum.values()[5] - 0.168_260_026_564_076_65).abs() < 1e-12);
    }

    #[test]
    fn computes_default_cri_ref_spectra() {
        let spectra = cri_ref(
            &[3000.0, 6500.0],
            Some(WavelengthGrid::new(360.0, 365.0, 1.0).unwrap()),
        )
        .unwrap();

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
        assert!((cct - 5455.485_887_350_497).abs() < 1.0);
        assert!((duv - (-0.004_423_324_748_595_847)).abs() < 1e-4);
    }

    #[test]
    fn loads_standard_illuminant_a() {
        let spectrum =
            standard_illuminant("A", Some(WavelengthGrid::new(360.0, 365.0, 1.0).unwrap()))
                .unwrap();
        assert!((spectrum.values()[0] - 6.144_62).abs() < 1e-12);
        assert!((spectrum.values()[5] - 6.947_2).abs() < 1e-12);
    }

    #[test]
    fn loads_standard_illuminant_f4() {
        let spectrum =
            standard_illuminant("F4", Some(WavelengthGrid::new(380.0, 385.0, 1.0).unwrap()))
                .unwrap();
        assert!((spectrum.values()[0] - 0.57).abs() < 1e-12);
        assert!((spectrum.values()[5] - 0.7).abs() < 1e-12);
    }

    #[test]
    fn loads_standard_illuminant_led_b1() {
        let spectrum = standard_illuminant(
            "LED_B1",
            Some(WavelengthGrid::new(380.0, 385.0, 1.0).unwrap()),
        )
        .unwrap();
        assert!((spectrum.values()[0] - 0.0).abs() < 1e-12);
        assert!((spectrum.values()[5] - 0.01).abs() < 1e-12);
    }

    #[test]
    fn loads_nominal_daylight_illuminants() {
        let spectrum =
            standard_illuminant("D50", Some(WavelengthGrid::new(360.0, 365.0, 1.0).unwrap()))
                .unwrap();
        assert!((spectrum.values()[0] - 0.940_694_581_416_273_4).abs() < 1e-12);
        assert!((spectrum.values()[5] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn exposes_standard_illuminant_names() {
        assert!(standard_illuminant_names().contains(&"A"));
        assert!(standard_illuminant_names().contains(&"D65"));
        assert!(standard_illuminant_names().contains(&"F4"));
        assert!(standard_illuminant_names().contains(&"LED_B1"));
    }
}
