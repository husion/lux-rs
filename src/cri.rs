use crate::cam::{ciecam02_viewing_conditions, xyz_to_jab_cam02ucs, CamSurround};
use crate::color::{
    cat_apply, xyz_to_yuv, xyz_to_yxy, yxy_to_xyz, CatTransform, Observer, TristimulusObserver,
};
use crate::error::{LuxError, LuxResult};
use crate::illuminants::{blackbody, cri_ref, daylightphase, xyz_to_cct};
use crate::photometry::integrate_xyz;
use crate::spectrum::{getwlr, SpectralMatrix, Spectrum, WavelengthGrid};

const CIE_RA_SAMPLE_CSV: &str = include_str!("../data/rfls/CIE_13_3_1995_R14.dat");
const CIE_RF_SAMPLE_CSV: &str = include_str!("../data/rfls/CIE224_2017_R99_1nm.dat");
const CIE_RA_SAMPLE_COUNT: usize = 8;
const CIE_RF_SAMPLE_COUNT: usize = 99;
const CIE_RA_ROUNDING_DIGITS: f64 = 10_000.0;
const CIE_RA_SCALE_FACTOR: f64 = 4.6;
const CIE_RF_SCALE_FACTOR: f64 = 6.73;
const CIE_RF_GRID: WavelengthGrid = WavelengthGrid {
    start: 380.0,
    end: 780.0,
    step: 1.0,
};
const CIE_RF_HUE_BINS: usize = 8;
const CIE_RF_ADAPTING_LUMINANCE: f64 = 100.0;
const CIE_RF_BACKGROUND_LUMINANCE: f64 = 20.0;

#[derive(Debug, Clone, PartialEq)]
pub struct CieRaResult {
    pub ra: f64,
    pub ri: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CieRfResult {
    pub cct: f64,
    pub rf: f64,
    pub rg: f64,
    pub dei: Vec<f64>,
    pub rfi: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tm30HueBin {
    pub hue_center_rad: f64,
    pub sample_count: usize,
    pub test_jab: [f64; 3],
    pub reference_jab: [f64; 3],
    pub mean_delta_e: f64,
    pub local_fidelity: f64,
    pub test_chroma: f64,
    pub reference_chroma: f64,
    pub chroma_shift: f64,
    pub chroma_shift_ratio: f64,
    pub test_hue_rad: f64,
    pub reference_hue_rad: f64,
    pub hue_shift_rad: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tm30Result {
    pub rf: CieRfResult,
    pub hue_bins: Vec<Tm30HueBin>,
    pub test_gamut_area: f64,
    pub reference_gamut_area: f64,
}

pub fn spd_to_ciera(spectrum: &Spectrum) -> LuxResult<f64> {
    Ok(spd_to_ciera_result(spectrum)?.ra)
}

pub fn spd_to_ciera_special(spectrum: &Spectrum) -> LuxResult<Vec<f64>> {
    Ok(spd_to_ciera_result(spectrum)?.ri)
}

pub fn spd_to_ciera_result(spectrum: &Spectrum) -> LuxResult<CieRaResult> {
    let observer = Observer::Cie1931_2.standard()?;
    let samples = cie_ra_samples()?;
    let test_white = relative_white_xyz(spectrum, &observer)?;
    let (cct, _) = xyz_to_cct(test_white, Observer::Cie1931_2)?;
    let source_grid = wavelength_grid_from_spectrum(spectrum)?;
    let reference = cri_ref(&[cct], Some(source_grid))?;
    let reference_spectrum = Spectrum::new(
        reference.wavelengths().to_vec(),
        reference.spectra()[0].clone(),
    )?;

    let test_xyz = sample_xyz_relative_to_white(spectrum, &samples, &observer)?;
    let reference_xyz = sample_xyz_relative_to_white(&reference_spectrum, &samples, &observer)?;
    let reference_white = relative_white_xyz(&reference_spectrum, &observer)?;

    let rounded_test_white = round_xyz_by_yxy(test_white);
    let rounded_reference_white = round_xyz_by_yxy(reference_white);

    let adapted_test_xyz: LuxResult<Vec<[f64; 3]>> = test_xyz
        .into_iter()
        .map(|xyz| {
            cat_apply(
                round_xyz_by_yxy(xyz),
                rounded_test_white,
                rounded_reference_white,
                CatTransform::Judd1945,
                1.0,
            )
        })
        .collect();
    let adapted_test_xyz = adapted_test_xyz?;
    let adapted_test_white = cat_apply(
        rounded_test_white,
        rounded_test_white,
        rounded_reference_white,
        CatTransform::Judd1945,
        1.0,
    )?;

    let mut ri = Vec::with_capacity(CIE_RA_SAMPLE_COUNT);
    for (test, reference) in adapted_test_xyz.into_iter().zip(reference_xyz.into_iter()) {
        let delta_e = cie_uvw_delta_e(
            test,
            adapted_test_white,
            round_xyz_by_yxy(reference),
            rounded_reference_white,
        );
        ri.push(100.0 - CIE_RA_SCALE_FACTOR * delta_e);
    }

    let rounded_sum: f64 = ri.iter().map(|value| value.round()).sum();
    Ok(CieRaResult {
        ra: rounded_sum / ri.len() as f64,
        ri,
    })
}

pub fn spds_to_ciera(spectra: &SpectralMatrix) -> LuxResult<Vec<f64>> {
    Ok(spds_to_ciera_result(spectra)?
        .into_iter()
        .map(|result| result.ra)
        .collect())
}

pub fn spds_to_ciera_result(spectra: &SpectralMatrix) -> LuxResult<Vec<CieRaResult>> {
    spectra
        .spectra()
        .iter()
        .map(|values| Spectrum::new(spectra.wavelengths().to_vec(), values.clone()))
        .map(|spectrum| spectrum.and_then(|s| spd_to_ciera_result(&s)))
        .collect()
}

pub fn spds_to_ciera_special(spectra: &SpectralMatrix) -> LuxResult<Vec<Vec<f64>>> {
    spectra
        .spectra()
        .iter()
        .map(|values| Spectrum::new(spectra.wavelengths().to_vec(), values.clone()))
        .map(|spectrum| spectrum.and_then(|s| spd_to_ciera_special(&s)))
        .collect()
}

pub fn spd_to_cierf(spectrum: &Spectrum) -> LuxResult<f64> {
    Ok(spd_to_cierf_result(spectrum)?.rf)
}

pub fn spd_to_iesrf(spectrum: &Spectrum) -> LuxResult<f64> {
    spd_to_cierf(spectrum)
}

pub fn spd_to_cierg(spectrum: &Spectrum) -> LuxResult<f64> {
    Ok(spd_to_cierf_result(spectrum)?.rg)
}

pub fn spd_to_iesrg(spectrum: &Spectrum) -> LuxResult<f64> {
    spd_to_cierg(spectrum)
}

pub fn spd_to_cierf_special(spectrum: &Spectrum) -> LuxResult<Vec<f64>> {
    Ok(spd_to_cierf_result(spectrum)?.rfi)
}

pub fn spd_to_iesrf_special(spectrum: &Spectrum) -> LuxResult<Vec<f64>> {
    spd_to_cierf_special(spectrum)
}

pub fn spd_to_cierf_result(spectrum: &Spectrum) -> LuxResult<CieRfResult> {
    Ok(spd_to_tm30_result(spectrum)?.rf)
}

pub fn spd_to_iesrf_result(spectrum: &Spectrum) -> LuxResult<CieRfResult> {
    spd_to_cierf_result(spectrum)
}

pub fn spd_to_tm30_result(spectrum: &Spectrum) -> LuxResult<Tm30Result> {
    let wavelengths = getwlr(CIE_RF_GRID)?;
    let source = spectrum.cie_interp_linear(&wavelengths, false)?;
    let samples = cie_rf_samples()?.cie_interp_linear(&wavelengths, false)?;
    let observer_xyz = Observer::Cie1964_10.standard()?;
    let observer_cct = Observer::Cie1931_2.standard()?;

    let test_white_cct = source.spd_to_xyz(&observer_cct, true)?;
    let (cct, _) = xyz_to_cct(test_white_cct, Observer::Cie1931_2)?;
    let reference = cierf_reference(cct, CIE_RF_GRID)?;

    let test_xyz = sample_xyz_relative_to_white(&source, &samples, &observer_xyz)?;
    let reference_xyz = sample_xyz_relative_to_white(&reference, &samples, &observer_xyz)?;
    let test_white = source.spd_to_xyz(&observer_xyz, true)?;
    let reference_white = reference.spd_to_xyz(&observer_xyz, true)?;

    let test_conditions = ciecam02_viewing_conditions(
        test_white,
        None,
        CIE_RF_ADAPTING_LUMINANCE,
        CIE_RF_BACKGROUND_LUMINANCE,
        CamSurround::Average,
        Some(1.0),
        Some(CatTransform::Cat02),
    )?;
    let reference_conditions = ciecam02_viewing_conditions(
        reference_white,
        None,
        CIE_RF_ADAPTING_LUMINANCE,
        CIE_RF_BACKGROUND_LUMINANCE,
        CamSurround::Average,
        Some(1.0),
        Some(CatTransform::Cat02),
    )?;

    let jab_test: LuxResult<Vec<[f64; 3]>> = test_xyz
        .iter()
        .copied()
        .map(|xyz| xyz_to_jab_cam02ucs(xyz, test_conditions))
        .collect();
    let jab_reference: LuxResult<Vec<[f64; 3]>> = reference_xyz
        .iter()
        .copied()
        .map(|xyz| xyz_to_jab_cam02ucs(xyz, reference_conditions))
        .collect();
    let jab_test = jab_test?;
    let jab_reference = jab_reference?;

    let dei: Vec<f64> = jab_test
        .iter()
        .zip(jab_reference.iter())
        .map(|(test, reference)| euclidean3(*test, *reference))
        .collect();
    let dea = mean(&dei)?;
    let rfi: Vec<f64> = dei
        .iter()
        .copied()
        .map(|de| log_scale(de, CIE_RF_SCALE_FACTOR))
        .collect();
    let rf = log_scale(dea, CIE_RF_SCALE_FACTOR);

    let hue_summary = tm30_hue_bin_summary(&jab_test, &jab_reference, CIE_RF_HUE_BINS)?;
    let rg = 100.0 * hue_summary.test_area / hue_summary.reference_area;

    Ok(Tm30Result {
        rf: CieRfResult {
            cct,
            rf,
            rg,
            dei,
            rfi,
        },
        hue_bins: hue_summary.bins,
        test_gamut_area: hue_summary.test_area,
        reference_gamut_area: hue_summary.reference_area,
    })
}

pub fn spd_to_ies_tm30_result(spectrum: &Spectrum) -> LuxResult<Tm30Result> {
    spd_to_tm30_result(spectrum)
}

pub fn spds_to_cierf(spectra: &SpectralMatrix) -> LuxResult<Vec<f64>> {
    Ok(spds_to_cierf_result(spectra)?
        .into_iter()
        .map(|result| result.rf)
        .collect())
}

pub fn spds_to_cierf_result(spectra: &SpectralMatrix) -> LuxResult<Vec<CieRfResult>> {
    spectra
        .spectra()
        .iter()
        .map(|values| Spectrum::new(spectra.wavelengths().to_vec(), values.clone()))
        .map(|spectrum| spectrum.and_then(|s| spd_to_cierf_result(&s)))
        .collect()
}

pub fn spds_to_iesrf(spectra: &SpectralMatrix) -> LuxResult<Vec<f64>> {
    spds_to_cierf(spectra)
}

pub fn spds_to_iesrf_result(spectra: &SpectralMatrix) -> LuxResult<Vec<CieRfResult>> {
    spds_to_cierf_result(spectra)
}

pub fn spds_to_tm30_result(spectra: &SpectralMatrix) -> LuxResult<Vec<Tm30Result>> {
    spectra
        .spectra()
        .iter()
        .map(|values| Spectrum::new(spectra.wavelengths().to_vec(), values.clone()))
        .map(|spectrum| spectrum.and_then(|s| spd_to_tm30_result(&s)))
        .collect()
}

pub fn spds_to_ies_tm30_result(spectra: &SpectralMatrix) -> LuxResult<Vec<Tm30Result>> {
    spds_to_tm30_result(spectra)
}

pub fn spds_to_cierg(spectra: &SpectralMatrix) -> LuxResult<Vec<f64>> {
    spectra
        .spectra()
        .iter()
        .map(|values| Spectrum::new(spectra.wavelengths().to_vec(), values.clone()))
        .map(|spectrum| spectrum.and_then(|s| spd_to_cierg(&s)))
        .collect()
}

pub fn spds_to_iesrg(spectra: &SpectralMatrix) -> LuxResult<Vec<f64>> {
    spds_to_cierg(spectra)
}

pub fn spds_to_cierf_special(spectra: &SpectralMatrix) -> LuxResult<Vec<Vec<f64>>> {
    spectra
        .spectra()
        .iter()
        .map(|values| Spectrum::new(spectra.wavelengths().to_vec(), values.clone()))
        .map(|spectrum| spectrum.and_then(|s| spd_to_cierf_special(&s)))
        .collect()
}

pub fn spds_to_iesrf_special(spectra: &SpectralMatrix) -> LuxResult<Vec<Vec<f64>>> {
    spds_to_cierf_special(spectra)
}

fn cie_ra_samples() -> LuxResult<SpectralMatrix> {
    sample_matrix_from_csv(CIE_RA_SAMPLE_CSV, CIE_RA_SAMPLE_COUNT)
}

fn cie_rf_samples() -> LuxResult<SpectralMatrix> {
    sample_matrix_from_csv(CIE_RF_SAMPLE_CSV, CIE_RF_SAMPLE_COUNT)
}

fn sample_matrix_from_csv(csv: &str, sample_count: usize) -> LuxResult<SpectralMatrix> {
    let mut wavelengths = Vec::new();
    let mut spectra = vec![Vec::new(); sample_count];

    for line in csv.lines() {
        let row: Vec<f64> = line
            .trim()
            .split(',')
            .map(|value| {
                value
                    .parse::<f64>()
                    .map_err(|_| LuxError::ParseError("invalid CIE Ra sample data"))
            })
            .collect::<LuxResult<Vec<_>>>()?;
        if row.len() < sample_count + 1 {
            return Err(LuxError::ParseError("invalid CRI sample row width"));
        }
        wavelengths.push(row[0]);
        for index in 0..sample_count {
            spectra[index].push(row[index + 1]);
        }
    }

    SpectralMatrix::new(wavelengths, spectra)
}

fn cierf_reference(cct: f64, grid: WavelengthGrid) -> LuxResult<Spectrum> {
    let wavelengths = getwlr(grid)?;
    if cct < 4000.0 {
        return blackbody(cct, Some(grid), None, true);
    }
    if cct >= 5000.0 {
        return daylightphase(cct, Some(grid), false, false, None);
    }

    let blackbody = blackbody(cct, Some(grid), None, true)?;
    let daylight = daylightphase(cct, Some(grid), false, false, None)?;
    let observer = Observer::Cie1931_2.standard()?;
    let blackbody_y = blackbody.spd_to_xyz(&observer, false)?[1];
    let daylight_y = daylight.spd_to_xyz(&observer, false)?[1];
    let tb = 4000.0;
    let te = 5000.0;
    let c_bb = clamp01((te - cct) / (te - tb));
    let c_dl = clamp01((cct - tb) / (te - tb));

    let mixed: Vec<f64> = blackbody
        .values()
        .iter()
        .zip(daylight.values().iter())
        .map(|(bb, dl)| (100.0 * bb / blackbody_y) * c_bb + (100.0 * dl / daylight_y) * c_dl)
        .collect();
    let index_560 = wavelengths
        .iter()
        .enumerate()
        .min_by(|(_, left), (_, right)| (*left - 560.0).abs().total_cmp(&(*right - 560.0).abs()))
        .map(|(index, _)| index)
        .ok_or(LuxError::EmptyInput)?;
    let normalization = mixed[index_560];
    Spectrum::new(
        wavelengths,
        mixed
            .into_iter()
            .map(|value| value / normalization)
            .collect(),
    )
}

fn wavelength_grid_from_spectrum(spectrum: &Spectrum) -> LuxResult<WavelengthGrid> {
    let wavelengths = spectrum.wavelengths();
    if wavelengths.len() < 2 {
        return Err(LuxError::InvalidGridSpec);
    }
    let step = wavelengths[1] - wavelengths[0];
    for pair in wavelengths.windows(2) {
        if (pair[1] - pair[0] - step).abs() > 1e-9 {
            return Err(LuxError::InvalidInput(
                "CIE Ra currently requires an equal-step wavelength grid",
            ));
        }
    }
    WavelengthGrid::new(wavelengths[0], *wavelengths.last().unwrap(), step)
}

fn relative_white_xyz(spectrum: &Spectrum, observer: &TristimulusObserver) -> LuxResult<[f64; 3]> {
    spectrum.spd_to_xyz(observer, true)
}

fn sample_xyz_relative_to_white(
    source: &Spectrum,
    samples: &SpectralMatrix,
    observer: &TristimulusObserver,
) -> LuxResult<Vec<[f64; 3]>> {
    let wavelengths = source.wavelengths();
    let reflectances = samples.cie_interp_linear(wavelengths, false)?;
    let x_bar = observer.x_bar_spectrum()?.interpolate_linear(wavelengths)?;
    let y_bar = observer.vl_spectrum()?.interpolate_linear(wavelengths)?;
    let z_bar = observer.z_bar_spectrum()?.interpolate_linear(wavelengths)?;
    let white_raw = integrate_xyz(
        source,
        x_bar.values(),
        y_bar.values(),
        z_bar.values(),
        1.0,
        false,
    )?;
    let scale = 100.0 / white_raw[1];

    reflectances
        .spectra()
        .iter()
        .map(|reflectance| {
            let reflected = Spectrum::new(
                wavelengths.to_vec(),
                source
                    .values()
                    .iter()
                    .zip(reflectance.iter())
                    .map(|(spd, reflectance)| spd * reflectance)
                    .collect(),
            )?;
            let raw = integrate_xyz(
                &reflected,
                x_bar.values(),
                y_bar.values(),
                z_bar.values(),
                1.0,
                false,
            )?;
            Ok([raw[0] * scale, raw[1] * scale, raw[2] * scale])
        })
        .collect()
}

fn round_xyz_by_yxy(xyz: [f64; 3]) -> [f64; 3] {
    let mut yxy = xyz_to_yxy(xyz);
    yxy[1] = (yxy[1] * CIE_RA_ROUNDING_DIGITS).round() / CIE_RA_ROUNDING_DIGITS;
    yxy[2] = (yxy[2] * CIE_RA_ROUNDING_DIGITS).round() / CIE_RA_ROUNDING_DIGITS;
    yxy_to_xyz(yxy)
}

fn xyz_to_wuv(xyz: [f64; 3], white: [f64; 3]) -> [f64; 3] {
    let yuv = xyz_to_yuv(xyz);
    let white_yuv = xyz_to_yuv(white);
    let w = 25.0 * yuv[0].cbrt() - 17.0;
    [
        w,
        13.0 * w * (yuv[1] - white_yuv[1]),
        13.0 * w * (yuv[2] - white_yuv[2]) * (2.0 / 3.0),
    ]
}

fn cie_uvw_delta_e(
    test_xyz: [f64; 3],
    test_white: [f64; 3],
    ref_xyz: [f64; 3],
    ref_white: [f64; 3],
) -> f64 {
    let test = xyz_to_wuv(test_xyz, test_white);
    let reference = xyz_to_wuv(ref_xyz, ref_white);
    ((test[0] - reference[0]).powi(2)
        + (test[1] - reference[1]).powi(2)
        + (test[2] - reference[2]).powi(2))
    .sqrt()
}

fn log_scale(delta_e: f64, scale_factor: f64) -> f64 {
    10.0 * (((100.0 - scale_factor * delta_e) / 10.0).exp() + 1.0).ln()
}

fn euclidean3(left: [f64; 3], right: [f64; 3]) -> f64 {
    ((left[0] - right[0]).powi(2) + (left[1] - right[1]).powi(2) + (left[2] - right[2]).powi(2))
        .sqrt()
}

fn mean(values: &[f64]) -> LuxResult<f64> {
    if values.is_empty() {
        return Err(LuxError::EmptyInput);
    }
    Ok(values.iter().sum::<f64>() / values.len() as f64)
}

struct Tm30HueSummary {
    bins: Vec<Tm30HueBin>,
    test_area: f64,
    reference_area: f64,
}

fn tm30_hue_bin_summary(
    jab_test: &[[f64; 3]],
    jab_reference: &[[f64; 3]],
    hue_bins: usize,
) -> LuxResult<Tm30HueSummary> {
    if jab_test.len() != jab_reference.len() || jab_test.is_empty() {
        return Err(LuxError::InvalidInput(
            "test and reference CAM sample sets must be non-empty and aligned",
        ));
    }

    let edges: Vec<f64> = (0..=hue_bins)
        .map(|index| (index as f64) * std::f64::consts::TAU / hue_bins as f64)
        .collect();
    let mut bins = Vec::with_capacity(hue_bins);

    for bin_index in 0..hue_bins {
        let low = edges[bin_index];
        let high = edges[bin_index + 1];
        let mut test_sum = [0.0; 3];
        let mut reference_sum = [0.0; 3];
        let mut delta_e_sum = 0.0;
        let mut count = 0usize;

        for (test, reference) in jab_test.iter().zip(jab_reference.iter()) {
            let hue = positive_hue_angle(reference[1], reference[2]);
            let in_bin = if bin_index + 1 == hue_bins {
                hue >= low && hue <= high
            } else {
                hue >= low && hue < high
            };
            if in_bin {
                for axis in 0..3 {
                    test_sum[axis] += test[axis];
                    reference_sum[axis] += reference[axis];
                }
                delta_e_sum += euclidean3(*test, *reference);
                count += 1;
            }
        }

        if count == 0 {
            return Err(LuxError::InvalidInput(
                "CIE Rg encountered an empty hue bin",
            ));
        }

        let test_jab = [
            test_sum[0] / count as f64,
            test_sum[1] / count as f64,
            test_sum[2] / count as f64,
        ];
        let reference_jab = [
            reference_sum[0] / count as f64,
            reference_sum[1] / count as f64,
            reference_sum[2] / count as f64,
        ];
        let mean_delta_e = delta_e_sum / count as f64;
        let test_chroma = ab_radius(test_jab);
        let reference_chroma = ab_radius(reference_jab);
        let chroma_shift = test_chroma - reference_chroma;
        let test_hue_rad = positive_hue_angle(test_jab[1], test_jab[2]);
        let reference_hue_rad = positive_hue_angle(reference_jab[1], reference_jab[2]);
        bins.push(Tm30HueBin {
            hue_center_rad: 0.5 * (low + high),
            sample_count: count,
            test_jab,
            reference_jab,
            mean_delta_e,
            local_fidelity: log_scale(mean_delta_e, CIE_RF_SCALE_FACTOR),
            test_chroma,
            reference_chroma,
            chroma_shift,
            chroma_shift_ratio: test_chroma / reference_chroma,
            test_hue_rad,
            reference_hue_rad,
            hue_shift_rad: signed_hue_delta(test_hue_rad, reference_hue_rad),
        });
    }

    let test_points: Vec<[f64; 3]> = bins.iter().map(|bin| bin.test_jab).collect();
    let reference_points: Vec<[f64; 3]> = bins.iter().map(|bin| bin.reference_jab).collect();
    Ok(Tm30HueSummary {
        bins,
        test_area: polygon_area_ab(&test_points),
        reference_area: polygon_area_ab(&reference_points),
    })
}

fn polygon_area_ab(points: &[[f64; 3]]) -> f64 {
    let mut area = 0.0;
    for index in 0..points.len() {
        let next = (index + 1) % points.len();
        area += points[index][1] * points[next][2] - points[index][2] * points[next][1];
    }
    area.abs() * 0.5
}

fn ab_radius(jab: [f64; 3]) -> f64 {
    (jab[1].powi(2) + jab[2].powi(2)).sqrt()
}

fn positive_hue_angle(a: f64, b: f64) -> f64 {
    let angle = b.atan2(a);
    if angle < 0.0 {
        angle + std::f64::consts::TAU
    } else {
        angle
    }
}

fn signed_hue_delta(test_hue_rad: f64, reference_hue_rad: f64) -> f64 {
    let mut delta = test_hue_rad - reference_hue_rad;
    while delta <= -std::f64::consts::PI {
        delta += std::f64::consts::TAU;
    }
    while delta > std::f64::consts::PI {
        delta -= std::f64::consts::TAU;
    }
    delta
}

fn clamp01(value: f64) -> f64 {
    value.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::{
        spd_to_ciera, spd_to_ciera_result, spd_to_cierf, spd_to_cierf_result, spd_to_cierg,
        spd_to_ies_tm30_result, spd_to_iesrf, spd_to_iesrf_result, spd_to_iesrf_special,
        spd_to_iesrg, spd_to_tm30_result, spds_to_ciera, spds_to_ciera_result, spds_to_cierf,
        spds_to_cierf_result, spds_to_ies_tm30_result, spds_to_iesrf, spds_to_iesrf_result,
        spds_to_iesrf_special, spds_to_iesrg, spds_to_tm30_result,
    };
    use crate::illuminants::standard_illuminant;
    use crate::spectrum::{SpectralMatrix, Spectrum};

    #[test]
    fn ciera_d65_is_near_perfect() {
        let d65 = standard_illuminant("D65", None).unwrap();
        let result = spd_to_ciera_result(&d65).unwrap();
        assert_eq!(result.ri.len(), 8);
        assert!(result.ra > 99.0);
    }

    #[test]
    fn ciera_discriminates_f4_from_d65() {
        let d65 = standard_illuminant("D65", None).unwrap();
        let f4 = standard_illuminant("F4", None).unwrap();
        assert!(spd_to_ciera(&d65).unwrap() > spd_to_ciera(&f4).unwrap());
    }

    #[test]
    fn ciera_batch_api_matches_scalar_path() {
        let grid = Some(crate::spectrum::WavelengthGrid::new(380.0, 780.0, 1.0).unwrap());
        let d65 = standard_illuminant("D65", grid).unwrap();
        let f4 = standard_illuminant("F4", grid).unwrap();
        let spectra = SpectralMatrix::new(
            d65.wavelengths().to_vec(),
            vec![d65.values().to_vec(), f4.values().to_vec()],
        )
        .unwrap();
        let batch = spds_to_ciera(&spectra).unwrap();
        assert_eq!(batch.len(), 2);
        assert!((batch[0] - spd_to_ciera(&d65).unwrap()).abs() < 1e-9);
        assert!((batch[1] - spd_to_ciera(&f4).unwrap()).abs() < 1e-9);
    }

    #[test]
    fn ciera_requires_uniform_wavelength_grid() {
        let spectrum = Spectrum::new(vec![400.0, 410.0, 425.0], vec![1.0, 1.0, 1.0]).unwrap();
        assert!(spd_to_ciera(&spectrum).is_err());
    }

    #[test]
    fn cierf_d65_is_near_perfect_and_rg_is_neutral() {
        let d65 = standard_illuminant("D65", None).unwrap();
        let result = spd_to_cierf_result(&d65).unwrap();
        assert!(result.cct > 6000.0);
        assert_eq!(result.dei.len(), 99);
        assert_eq!(result.rfi.len(), 99);
        assert!(result.rf > 99.0);
        assert!(result.rg > 99.0);
    }

    #[test]
    fn cierf_discriminates_f4_from_d65() {
        let d65 = standard_illuminant("D65", None).unwrap();
        let f4 = standard_illuminant("F4", None).unwrap();
        assert!(spd_to_cierf(&d65).unwrap() > spd_to_cierf(&f4).unwrap());
        assert!(spd_to_cierg(&d65).unwrap() > spd_to_cierg(&f4).unwrap());
    }

    #[test]
    fn cierf_batch_api_matches_scalar_path() {
        let grid = Some(crate::spectrum::WavelengthGrid::new(380.0, 780.0, 1.0).unwrap());
        let d65 = standard_illuminant("D65", grid).unwrap();
        let f4 = standard_illuminant("F4", grid).unwrap();
        let spectra = SpectralMatrix::new(
            d65.wavelengths().to_vec(),
            vec![d65.values().to_vec(), f4.values().to_vec()],
        )
        .unwrap();
        let batch = spds_to_cierf(&spectra).unwrap();
        assert_eq!(batch.len(), 2);
        assert!((batch[0] - spd_to_cierf(&d65).unwrap()).abs() < 1e-9);
        assert!((batch[1] - spd_to_cierf(&f4).unwrap()).abs() < 1e-9);
    }

    #[test]
    fn ciera_batch_result_matches_scalar_path() {
        let grid = Some(crate::spectrum::WavelengthGrid::new(380.0, 780.0, 1.0).unwrap());
        let d65 = standard_illuminant("D65", grid).unwrap();
        let f4 = standard_illuminant("F4", grid).unwrap();
        let spectra = SpectralMatrix::new(
            d65.wavelengths().to_vec(),
            vec![d65.values().to_vec(), f4.values().to_vec()],
        )
        .unwrap();
        let batch = spds_to_ciera_result(&spectra).unwrap();
        assert_eq!(batch.len(), 2);
        assert_eq!(batch[0], spd_to_ciera_result(&d65).unwrap());
        assert_eq!(batch[1], spd_to_ciera_result(&f4).unwrap());
    }

    #[test]
    fn cierf_batch_result_matches_scalar_path() {
        let grid = Some(crate::spectrum::WavelengthGrid::new(380.0, 780.0, 1.0).unwrap());
        let d65 = standard_illuminant("D65", grid).unwrap();
        let f4 = standard_illuminant("F4", grid).unwrap();
        let spectra = SpectralMatrix::new(
            d65.wavelengths().to_vec(),
            vec![d65.values().to_vec(), f4.values().to_vec()],
        )
        .unwrap();
        let batch = spds_to_cierf_result(&spectra).unwrap();
        assert_eq!(batch.len(), 2);
        assert_eq!(batch[0], spd_to_cierf_result(&d65).unwrap());
        assert_eq!(batch[1], spd_to_cierf_result(&f4).unwrap());
    }

    #[test]
    fn ies_aliases_match_cie224_scalar_results() {
        let d65 = standard_illuminant("D65", None).unwrap();
        assert!((spd_to_iesrf(&d65).unwrap() - spd_to_cierf(&d65).unwrap()).abs() < 1e-12);
        assert!((spd_to_iesrg(&d65).unwrap() - spd_to_cierg(&d65).unwrap()).abs() < 1e-12);
        assert_eq!(
            spd_to_iesrf_special(&d65).unwrap(),
            super::spd_to_cierf_special(&d65).unwrap()
        );
        assert_eq!(
            spd_to_iesrf_result(&d65).unwrap(),
            spd_to_cierf_result(&d65).unwrap()
        );
    }

    #[test]
    fn ies_aliases_match_cie224_batch_results() {
        let grid = Some(crate::spectrum::WavelengthGrid::new(380.0, 780.0, 1.0).unwrap());
        let d65 = standard_illuminant("D65", grid).unwrap();
        let f4 = standard_illuminant("F4", grid).unwrap();
        let spectra = SpectralMatrix::new(
            d65.wavelengths().to_vec(),
            vec![d65.values().to_vec(), f4.values().to_vec()],
        )
        .unwrap();

        assert_eq!(
            spds_to_iesrf(&spectra).unwrap(),
            spds_to_cierf(&spectra).unwrap()
        );
        assert_eq!(
            spds_to_iesrg(&spectra).unwrap(),
            super::spds_to_cierg(&spectra).unwrap()
        );
        assert_eq!(
            spds_to_iesrf_special(&spectra).unwrap(),
            super::spds_to_cierf_special(&spectra).unwrap()
        );
        assert_eq!(
            spds_to_iesrf_result(&spectra).unwrap(),
            spds_to_cierf_result(&spectra).unwrap()
        );
    }

    #[test]
    fn tm30_result_exposes_hue_bin_summary() {
        let d65 = standard_illuminant("D65", None).unwrap();
        let result = spd_to_tm30_result(&d65).unwrap();
        assert_eq!(result.hue_bins.len(), 8);
        assert!(result.test_gamut_area > 0.0);
        assert!(result.reference_gamut_area > 0.0);
        assert!(
            (result.rf.rg - 100.0 * result.test_gamut_area / result.reference_gamut_area).abs()
                < 1e-9
        );
    }

    #[test]
    fn tm30_hue_bins_expose_local_metrics() {
        let d65 = standard_illuminant("D65", None).unwrap();
        let result = spd_to_tm30_result(&d65).unwrap();
        for bin in result.hue_bins {
            assert!(bin.sample_count > 0);
            assert!(bin.mean_delta_e >= 0.0);
            assert!(bin.local_fidelity > 0.0);
            assert!(bin.test_chroma >= 0.0);
            assert!(bin.reference_chroma > 0.0);
            assert!((bin.chroma_shift - (bin.test_chroma - bin.reference_chroma)).abs() < 1e-12);
            assert!(
                (bin.chroma_shift_ratio - (bin.test_chroma / bin.reference_chroma)).abs() < 1e-12
            );
            assert!((0.0..=std::f64::consts::TAU).contains(&bin.test_hue_rad));
            assert!((0.0..=std::f64::consts::TAU).contains(&bin.reference_hue_rad));
            assert!((-std::f64::consts::PI..=std::f64::consts::PI).contains(&bin.hue_shift_rad));
        }
    }

    #[test]
    fn tm30_local_fidelity_is_high_for_d65() {
        let d65 = standard_illuminant("D65", None).unwrap();
        let result = spd_to_tm30_result(&d65).unwrap();
        let min_local_fidelity = result
            .hue_bins
            .iter()
            .map(|bin| bin.local_fidelity)
            .fold(f64::INFINITY, f64::min);
        assert!(min_local_fidelity > 95.0);
    }

    #[test]
    fn tm30_hue_shift_is_small_for_d65() {
        let d65 = standard_illuminant("D65", None).unwrap();
        let max_abs_hue_shift = spd_to_tm30_result(&d65)
            .unwrap()
            .hue_bins
            .iter()
            .map(|bin| bin.hue_shift_rad.abs())
            .fold(0.0, f64::max);
        assert!(max_abs_hue_shift < 0.05);
    }

    #[test]
    fn tm30_alias_and_batch_results_match() {
        let grid = Some(crate::spectrum::WavelengthGrid::new(380.0, 780.0, 1.0).unwrap());
        let d65 = standard_illuminant("D65", grid).unwrap();
        let f4 = standard_illuminant("F4", grid).unwrap();
        let spectra = SpectralMatrix::new(
            d65.wavelengths().to_vec(),
            vec![d65.values().to_vec(), f4.values().to_vec()],
        )
        .unwrap();

        assert_eq!(
            spd_to_ies_tm30_result(&d65).unwrap(),
            spd_to_tm30_result(&d65).unwrap()
        );
        assert_eq!(
            spds_to_ies_tm30_result(&spectra).unwrap(),
            spds_to_tm30_result(&spectra).unwrap()
        );
    }
}
