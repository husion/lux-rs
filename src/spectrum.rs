use crate::color::TristimulusObserver;
use crate::cri::{
    spd_to_ciera, spd_to_ciera_result, spd_to_ciera_special, spd_to_cierf, spd_to_cierf_result,
    spd_to_cierf_special, spd_to_cierg, spd_to_ies_tm30_result, spd_to_iesrf, spd_to_iesrf_result,
    spd_to_iesrf_special, spd_to_iesrg, spd_to_tm30_result, spds_to_ciera, spds_to_ciera_result,
    spds_to_ciera_special, spds_to_cierf, spds_to_cierf_result, spds_to_cierf_special,
    spds_to_cierg, spds_to_ies_tm30_result, spds_to_iesrf, spds_to_iesrf_result,
    spds_to_iesrf_special, spds_to_iesrg, spds_to_tm30_result, CieRaResult, CieRfResult,
    Tm30Result,
};
use crate::error::{LuxError, LuxResult};
use crate::photometry::{spd_to_power, PowerType};
use crate::spectral_mismatch::{
    spectral_mismatch_correction_factor, spectral_mismatch_correction_factors,
    spectral_mismatch_f1prime, spectral_mismatch_f1primes,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WavelengthGrid {
    pub start: f64,
    pub end: f64,
    pub step: f64,
}

impl WavelengthGrid {
    pub fn new(start: f64, end: f64, step: f64) -> LuxResult<Self> {
        if !start.is_finite() || !end.is_finite() || !step.is_finite() || step <= 0.0 || end < start
        {
            return Err(LuxError::InvalidGridSpec);
        }
        Ok(Self { start, end, step })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SingleSpectrum {
    wavelengths: Vec<f64>,
    values: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Spectrum {
    wavelengths: Vec<f64>,
    spectra: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectrumNormalization {
    Max(f64),
    Area(f64),
    Lambda(f64),
    Radiometric(f64),
    Photometric(f64),
    Quantal(f64),
}

trait IntoSpectrumRows {
    fn into_rows(self) -> Vec<Vec<f64>>;
}

impl IntoSpectrumRows for Vec<f64> {
    fn into_rows(self) -> Vec<Vec<f64>> {
        vec![self]
    }
}

impl IntoSpectrumRows for Vec<Vec<f64>> {
    fn into_rows(self) -> Vec<Vec<f64>> {
        self
    }
}

impl SingleSpectrum {
    pub fn new(wavelengths: Vec<f64>, values: Vec<f64>) -> LuxResult<Self> {
        if wavelengths.is_empty() || values.is_empty() {
            return Err(LuxError::EmptyInput);
        }
        if wavelengths.len() != values.len() {
            return Err(LuxError::MismatchedLengths {
                wavelengths: wavelengths.len(),
                values: values.len(),
            });
        }
        if wavelengths
            .windows(2)
            .any(|pair| !(pair[1].is_finite() && pair[0].is_finite() && pair[1] > pair[0]))
        {
            return Err(LuxError::NonMonotonicWavelengths);
        }
        Ok(Self {
            wavelengths,
            values,
        })
    }

    pub fn wavelengths(&self) -> &[f64] {
        &self.wavelengths
    }

    pub fn values(&self) -> &[f64] {
        &self.values
    }

    pub fn spacing(&self) -> LuxResult<Vec<f64>> {
        getwld(&self.wavelengths)
    }

    pub fn interpolate_linear(&self, target_wavelengths: &[f64]) -> LuxResult<Self> {
        if target_wavelengths.is_empty() {
            return Err(LuxError::EmptyInput);
        }
        if target_wavelengths
            .windows(2)
            .any(|pair| !(pair[1].is_finite() && pair[0].is_finite() && pair[1] > pair[0]))
        {
            return Err(LuxError::NonMonotonicWavelengths);
        }

        let mut values = Vec::with_capacity(target_wavelengths.len());
        for &target in target_wavelengths {
            values.push(self.interpolate_one_linear(target));
        }
        SingleSpectrum::new(target_wavelengths.to_vec(), values)
    }

    pub fn cie_interp_linear(
        &self,
        target_wavelengths: &[f64],
        negative_values_allowed: bool,
    ) -> LuxResult<Self> {
        let mut interpolated = self.interpolate_linear(target_wavelengths)?;
        if !negative_values_allowed {
            for value in &mut interpolated.values {
                if *value < 0.0 {
                    *value = 0.0;
                }
            }
        }
        Ok(interpolated)
    }

    fn interpolate_one_linear(&self, target: f64) -> f64 {
        let wavelengths = &self.wavelengths;
        let values = &self.values;

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

    pub fn normalize(
        &self,
        mode: SpectrumNormalization,
        observer: Option<&TristimulusObserver>,
    ) -> LuxResult<Self> {
        let scale = match mode {
            SpectrumNormalization::Max(target_max) => {
                target_max
                    / self
                        .values
                        .iter()
                        .copied()
                        .fold(f64::NEG_INFINITY, f64::max)
            }
            SpectrumNormalization::Area(target_area) => {
                let area: f64 = self
                    .values
                    .iter()
                    .zip(self.spacing()?.iter())
                    .map(|(value, dl)| value * dl)
                    .sum();
                target_area / area
            }
            SpectrumNormalization::Lambda(target_wavelength) => {
                let (index, _) = self
                    .wavelengths
                    .iter()
                    .enumerate()
                    .map(|(index, wavelength)| (index, (wavelength - target_wavelength).abs()))
                    .min_by(|(_, lhs), (_, rhs)| lhs.partial_cmp(rhs).unwrap())
                    .ok_or(LuxError::EmptyInput)?;
                1.0 / self.values[index]
            }
            SpectrumNormalization::Radiometric(target_power) => {
                target_power / spd_to_power(self, PowerType::Radiometric, None)?
            }
            SpectrumNormalization::Photometric(target_power) => {
                target_power / spd_to_power(self, PowerType::Photometric, observer)?
            }
            SpectrumNormalization::Quantal(target_power) => {
                target_power / spd_to_power(self, PowerType::Quantal, None)?
            }
        };

        SingleSpectrum::new(
            self.wavelengths.clone(),
            self.values.iter().map(|value| value * scale).collect(),
        )
    }

    pub fn spd_to_xyz(
        &self,
        observer: &TristimulusObserver,
        relative: bool,
    ) -> LuxResult<[f64; 3]> {
        crate::photometry::spd_to_xyz(self, observer, relative)
    }

    pub fn spd_to_ler(&self, observer: &TristimulusObserver) -> LuxResult<f64> {
        crate::photometry::spd_to_ler(self, observer)
    }

    pub fn spd_to_power(
        &self,
        power_type: PowerType,
        observer: Option<&TristimulusObserver>,
    ) -> LuxResult<f64> {
        spd_to_power(self, power_type, observer)
    }

    pub fn spd_to_ciera(&self) -> LuxResult<f64> {
        spd_to_ciera(self)
    }

    pub fn spd_to_ciera_special(&self) -> LuxResult<Vec<f64>> {
        spd_to_ciera_special(self)
    }

    pub fn spd_to_ciera_result(&self) -> LuxResult<CieRaResult> {
        spd_to_ciera_result(self)
    }

    pub fn spd_to_cierf(&self) -> LuxResult<f64> {
        spd_to_cierf(self)
    }

    pub fn spd_to_iesrf(&self) -> LuxResult<f64> {
        spd_to_iesrf(self)
    }

    pub fn spd_to_cierg(&self) -> LuxResult<f64> {
        spd_to_cierg(self)
    }

    pub fn spd_to_iesrg(&self) -> LuxResult<f64> {
        spd_to_iesrg(self)
    }

    pub fn spd_to_cierf_special(&self) -> LuxResult<Vec<f64>> {
        spd_to_cierf_special(self)
    }

    pub fn spd_to_iesrf_special(&self) -> LuxResult<Vec<f64>> {
        spd_to_iesrf_special(self)
    }

    pub fn spd_to_cierf_result(&self) -> LuxResult<CieRfResult> {
        spd_to_cierf_result(self)
    }

    pub fn spd_to_iesrf_result(&self) -> LuxResult<CieRfResult> {
        spd_to_iesrf_result(self)
    }

    pub fn spd_to_tm30_result(&self) -> LuxResult<Tm30Result> {
        spd_to_tm30_result(self)
    }

    pub fn spd_to_ies_tm30_result(&self) -> LuxResult<Tm30Result> {
        spd_to_ies_tm30_result(self)
    }

    pub fn spectral_mismatch_f1prime(
        &self,
        calibration_illuminant: &Spectrum,
        target_responsivity: &Spectrum,
    ) -> LuxResult<f64> {
        spectral_mismatch_f1prime(self, calibration_illuminant, target_responsivity)
    }

    pub fn spectral_mismatch_correction_factor(
        &self,
        detector: &Spectrum,
        calibration_illuminant: &Spectrum,
        target_responsivity: &Spectrum,
    ) -> LuxResult<f64> {
        spectral_mismatch_correction_factor(
            self,
            detector,
            calibration_illuminant,
            target_responsivity,
        )
    }
}

impl Spectrum {
    pub fn new<T: IntoSpectrumRows>(wavelengths: Vec<f64>, spectra: T) -> LuxResult<Self> {
        let spectra = spectra.into_rows();
        if wavelengths.is_empty() || spectra.is_empty() {
            return Err(LuxError::EmptyInput);
        }
        if wavelengths
            .windows(2)
            .any(|pair| !(pair[1].is_finite() && pair[0].is_finite() && pair[1] > pair[0]))
        {
            return Err(LuxError::NonMonotonicWavelengths);
        }
        for values in &spectra {
            if values.len() != wavelengths.len() {
                return Err(LuxError::MismatchedLengths {
                    wavelengths: wavelengths.len(),
                    values: values.len(),
                });
            }
        }

        Ok(Self {
            wavelengths,
            spectra,
        })
    }

    pub fn wavelengths(&self) -> &[f64] {
        &self.wavelengths
    }

    pub fn values(&self) -> &[f64] {
        self.spectra
            .first()
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn spectra(&self) -> &[Vec<f64>] {
        &self.spectra
    }

    pub fn spectrum_count(&self) -> usize {
        self.spectra.len()
    }

    pub fn wavelength_count(&self) -> usize {
        self.wavelengths.len()
    }

    pub fn spacing(&self) -> LuxResult<Vec<f64>> {
        getwld(&self.wavelengths)
    }

    pub fn into_vec(self) -> Vec<Vec<f64>> {
        self.spectra
    }

    pub fn interpolate_linear(&self, target_wavelengths: &[f64]) -> LuxResult<Self> {
        let mut spectra = Vec::with_capacity(self.spectra.len());
        for values in &self.spectra {
            let spectrum = SingleSpectrum::new(self.wavelengths.clone(), values.clone())?;
            let interpolated = spectrum.interpolate_linear(target_wavelengths)?;
            spectra.push(interpolated.values().to_vec());
        }
        Spectrum::new(target_wavelengths.to_vec(), spectra)
    }

    pub fn cie_interp_linear(
        &self,
        target_wavelengths: &[f64],
        negative_values_allowed: bool,
    ) -> LuxResult<Self> {
        let mut spectra = Vec::with_capacity(self.spectra.len());
        for values in &self.spectra {
            let spectrum = SingleSpectrum::new(self.wavelengths.clone(), values.clone())?;
            let interpolated =
                spectrum.cie_interp_linear(target_wavelengths, negative_values_allowed)?;
            spectra.push(interpolated.values().to_vec());
        }
        Spectrum::new(target_wavelengths.to_vec(), spectra)
    }

    pub fn normalize(
        &self,
        mode: SpectrumNormalization,
        observer: Option<&TristimulusObserver>,
    ) -> LuxResult<Self> {
        self.normalize_each(&[mode], observer)
    }

    pub fn normalize_each(
        &self,
        modes: &[SpectrumNormalization],
        observer: Option<&TristimulusObserver>,
    ) -> LuxResult<Self> {
        if modes.is_empty() {
            return Err(LuxError::EmptyInput);
        }

        let mut spectra = Vec::with_capacity(self.spectra.len());
        for (index, values) in self.spectra.iter().enumerate() {
            let mode = modes.get(index).copied().unwrap_or(modes[0]);
            let spectrum = SingleSpectrum::new(self.wavelengths.clone(), values.clone())?;
            let normalized = spectrum.normalize(mode, observer)?;
            spectra.push(normalized.values().to_vec());
        }

        Spectrum::new(self.wavelengths.clone(), spectra)
    }

    pub fn spd_to_xyz(
        &self,
        observer: &TristimulusObserver,
        relative: bool,
    ) -> LuxResult<Vec<[f64; 3]>> {
        let wavelengths = self.wavelengths();
        let x_bar = observer.x_bar_spectrum()?.interpolate_linear(wavelengths)?;
        let y_bar = observer.vl_spectrum()?.interpolate_linear(wavelengths)?;
        let z_bar = observer.z_bar_spectrum()?.interpolate_linear(wavelengths)?;

        self.spectra
            .iter()
            .map(|values| {
                let spectrum = SingleSpectrum::new(wavelengths.to_vec(), values.clone())?;
                crate::photometry::integrate_xyz(
                    &spectrum,
                    x_bar.values(),
                    y_bar.values(),
                    z_bar.values(),
                    observer.k,
                    relative,
                )
            })
            .collect()
    }

    pub fn spd_to_ler(&self, observer: &TristimulusObserver) -> LuxResult<Vec<f64>> {
        self.spectra
            .iter()
            .map(|values| {
                let spectrum = SingleSpectrum::new(self.wavelengths.to_vec(), values.clone())?;
                spectrum.spd_to_ler(observer)
            })
            .collect()
    }

    pub fn spd_to_ciera(&self) -> LuxResult<Vec<f64>> {
        spds_to_ciera(self)
    }

    pub fn spd_to_ciera_special(&self) -> LuxResult<Vec<Vec<f64>>> {
        spds_to_ciera_special(self)
    }

    pub fn spd_to_ciera_result(&self) -> LuxResult<Vec<CieRaResult>> {
        spds_to_ciera_result(self)
    }

    pub fn spd_to_cierf(&self) -> LuxResult<Vec<f64>> {
        spds_to_cierf(self)
    }

    pub fn spd_to_iesrf(&self) -> LuxResult<Vec<f64>> {
        spds_to_iesrf(self)
    }

    pub fn spd_to_cierg(&self) -> LuxResult<Vec<f64>> {
        spds_to_cierg(self)
    }

    pub fn spd_to_iesrg(&self) -> LuxResult<Vec<f64>> {
        spds_to_iesrg(self)
    }

    pub fn spd_to_cierf_special(&self) -> LuxResult<Vec<Vec<f64>>> {
        spds_to_cierf_special(self)
    }

    pub fn spd_to_iesrf_special(&self) -> LuxResult<Vec<Vec<f64>>> {
        spds_to_iesrf_special(self)
    }

    pub fn spd_to_cierf_result(&self) -> LuxResult<Vec<CieRfResult>> {
        spds_to_cierf_result(self)
    }

    pub fn spd_to_iesrf_result(&self) -> LuxResult<Vec<CieRfResult>> {
        spds_to_iesrf_result(self)
    }

    pub fn spd_to_tm30_result(&self) -> LuxResult<Vec<Tm30Result>> {
        spds_to_tm30_result(self)
    }

    pub fn spd_to_ies_tm30_result(&self) -> LuxResult<Vec<Tm30Result>> {
        spds_to_ies_tm30_result(self)
    }

    pub fn spectral_mismatch_f1primes(
        &self,
        calibration_illuminant: &Spectrum,
        target_responsivity: &Spectrum,
    ) -> LuxResult<Vec<f64>> {
        spectral_mismatch_f1primes(self, calibration_illuminant, target_responsivity)
    }

    pub fn spectral_mismatch_correction_factors(
        &self,
        detectors: &Spectrum,
        calibration_illuminant: &Spectrum,
        target_responsivity: &Spectrum,
    ) -> LuxResult<Vec<Vec<f64>>> {
        spectral_mismatch_correction_factors(
            self,
            detectors,
            calibration_illuminant,
            target_responsivity,
        )
    }
}

fn linear_segment(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    y0 + (y1 - y0) * ((x - x0) / (x1 - x0))
}

pub fn getwlr(grid: WavelengthGrid) -> LuxResult<Vec<f64>> {
    let mut wavelengths = Vec::new();
    let mut current = grid.start;
    let epsilon = grid.step.abs() * 1e-9;

    while current <= grid.end + epsilon {
        wavelengths.push(current);
        current += grid.step;
    }

    if wavelengths.is_empty() {
        return Err(LuxError::InvalidGridSpec);
    }

    Ok(wavelengths)
}

pub fn getwld(wavelengths: &[f64]) -> LuxResult<Vec<f64>> {
    if wavelengths.is_empty() {
        return Err(LuxError::EmptyInput);
    }
    if wavelengths.len() == 1 {
        return Ok(vec![0.0]);
    }
    if wavelengths
        .windows(2)
        .any(|pair| !(pair[1].is_finite() && pair[0].is_finite() && pair[1] > pair[0]))
    {
        return Err(LuxError::NonMonotonicWavelengths);
    }

    let diffs: Vec<f64> = wavelengths
        .windows(2)
        .map(|pair| pair[1] - pair[0])
        .collect();
    let mut spacing = Vec::with_capacity(wavelengths.len());
    spacing.push(diffs[0]);
    for idx in 1..wavelengths.len() - 1 {
        spacing.push((diffs[idx - 1] + diffs[idx]) / 2.0);
    }
    spacing.push(*diffs.last().unwrap());
    Ok(spacing)
}
