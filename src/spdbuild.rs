use crate::color::{xyz_to_yxy, Observer};
use crate::error::{LuxError, LuxResult};
use crate::illuminants::cct_to_xyz;
use crate::spectrum::{getwlr, Spectrum, WavelengthGrid};

pub const DEFAULT_WL_GRID: WavelengthGrid = WavelengthGrid {
    start: 360.0,
    end: 830.0,
    step: 1.0,
};

#[derive(Debug, Clone, PartialEq)]
pub struct RoundedTriangleParams {
    pub peakwl: f64,
    pub fwhm: Option<f64>,
    pub rounding: f64,
    pub min_v: f64,
    pub max_v: f64,
    pub fw: f64,
    pub rw: f64,
}

impl Default for RoundedTriangleParams {
    fn default() -> Self {
        Self {
            peakwl: 530.0,
            fwhm: Some(100.0),
            rounding: 0.5,
            min_v: 0.0,
            max_v: 1.0,
            fw: 100.0,
            rw: 100.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MonoLedParams {
    pub peakwl: f64,
    pub fwhm: f64,
    pub strength_shoulder: f64,
    pub bw_order: f64,
}

impl Default for MonoLedParams {
    fn default() -> Self {
        Self {
            peakwl: 530.0,
            fwhm: 20.0,
            strength_shoulder: 2.0,
            bw_order: -1.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PhosphorLedParams {
    pub peakwl: f64,
    pub fwhm: f64,
    pub bw_order: f64,
    pub strength_shoulder: f64,
    pub strength_ph: Option<f64>,
    pub peakwl_ph1: f64,
    pub fwhm_ph1: f64,
    pub strength_ph1: f64,
    pub peakwl_ph2: f64,
    pub fwhm_ph2: f64,
    pub strength_ph2: Option<f64>,
    pub use_piecewise_fcn: bool,
}

impl Default for PhosphorLedParams {
    fn default() -> Self {
        Self {
            peakwl: 450.0,
            fwhm: 20.0,
            bw_order: -1.0,
            strength_shoulder: 2.0,
            strength_ph: Some(0.0),
            peakwl_ph1: 530.0,
            fwhm_ph1: 80.0,
            strength_ph1: 1.0,
            peakwl_ph2: 560.0,
            fwhm_ph2: 80.0,
            strength_ph2: None,
            use_piecewise_fcn: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PhosphorLedComponents {
    pub spd: Spectrum,
    pub components: Spectrum,
}

/// Generate Gaussian spectrum.
pub fn gaussian_spd(
    peakwls: &[f64],
    fwhms: &[f64],
    grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    if peakwls.is_empty() || fwhms.is_empty() {
        return Err(LuxError::EmptyInput);
    }
    if fwhms.len() != 1 && fwhms.len() != peakwls.len() {
        return Err(LuxError::InvalidGridSpec);
    }

    let grid = grid.unwrap_or(DEFAULT_WL_GRID);
    let wavelengths = getwlr(grid)?;
    let num_peaks = peakwls.len();
    let fwhm_to_sig = 1.0 / (2.0 * (2.0 * 2.0f64.ln()).sqrt());
    let mut spectra = Vec::with_capacity(num_peaks);

    for i in 0..num_peaks {
        let peakwl = peakwls[i];
        let fwhm = if fwhms.len() == 1 { fwhms[0] } else { fwhms[i] };
        let sig = fwhm * fwhm_to_sig;

        let mut values = Vec::with_capacity(wavelengths.len());
        for &wl in &wavelengths {
            let val = (-0.5 * ((wl - peakwl) / sig).powi(2)).exp();
            values.push(val);
        }
        spectra.push(values);
    }

    Spectrum::new(wavelengths, spectra)
}

/// Generate 2nd order Lorentzian spectrum.
pub fn lorentzian2_spd(
    peakwls: &[f64],
    fwhms: &[f64],
    grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    if peakwls.is_empty() || fwhms.is_empty() {
        return Err(LuxError::EmptyInput);
    }
    if fwhms.len() != 1 && fwhms.len() != peakwls.len() {
        return Err(LuxError::InvalidGridSpec);
    }

    let grid = grid.unwrap_or(DEFAULT_WL_GRID);
    let wavelengths = getwlr(grid)?;
    let num_peaks = peakwls.len();
    let n = 2.0 * (2.0f64.sqrt() - 1.0).sqrt();
    let mut spectra = Vec::with_capacity(num_peaks);

    for i in 0..num_peaks {
        let peakwl = peakwls[i];
        let fwhm = if fwhms.len() == 1 { fwhms[0] } else { fwhms[i] };

        let mut values = Vec::with_capacity(wavelengths.len());
        for &wl in &wavelengths {
            let val = (1.0 + (n * (wl - peakwl) / fwhm).powi(2)).powf(-2.0);
            values.push(val);
        }
        spectra.push(values);
    }

    Spectrum::new(wavelengths, spectra)
}

/// Generate Butterworth based spectrum.
pub fn butterworth_spd(
    peakwls: &[f64],
    fwhms: &[f64],
    bw_orders: &[f64],
    grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    if peakwls.is_empty() || fwhms.is_empty() || bw_orders.is_empty() {
        return Err(LuxError::EmptyInput);
    }
    if (fwhms.len() != 1 && fwhms.len() != peakwls.len())
        || (bw_orders.len() != 1 && bw_orders.len() != peakwls.len())
    {
        return Err(LuxError::InvalidGridSpec);
    }

    let grid = grid.unwrap_or(DEFAULT_WL_GRID);
    let wavelengths = getwlr(grid)?;
    let num_peaks = peakwls.len();
    let mut spectra = Vec::with_capacity(num_peaks);

    for i in 0..num_peaks {
        let peakwl = peakwls[i];
        let fwhm = if fwhms.len() == 1 { fwhms[0] } else { fwhms[i] };
        let bw_order = if bw_orders.len() == 1 { bw_orders[0] } else { bw_orders[i] };

        let mut values = Vec::with_capacity(wavelengths.len());
        for &wl in &wavelengths {
            let val = 1.0 / (1.0 + (2.0 * (wl - peakwl) / fwhm).abs().powf(2.0 * bw_order));
            values.push(val);
        }
        spectra.push(values);
    }

    Spectrum::new(wavelengths, spectra)
}

/// Generate rounded triangle spectrum.
pub fn roundedtriangle_spd(
    params: &[RoundedTriangleParams],
    grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    if params.is_empty() {
        return Err(LuxError::EmptyInput);
    }

    let grid = grid.unwrap_or(DEFAULT_WL_GRID);
    let wavelengths = getwlr(grid)?;
    let num_peaks = params.len();
    let mut spectra = Vec::with_capacity(num_peaks);

    for param in params {
        let peakwl = param.peakwl.abs();
        let rounding = param.rounding.abs();
        let min_v = param.min_v.abs();
        let max_v = param.max_v.abs();

        let (fw, rw) = match param.fwhm {
            Some(fwhm) => {
                let width = fwhm.abs() / (rounding / 4.0 + 1.0);
                (width, width)
            }
            None => (param.fw.abs(), param.rw.abs()),
        };

        let eps = 1e-308;
        let r_param = if rounding == 0.0 { eps } else { rounding };

        let mut values = Vec::with_capacity(wavelengths.len());
        for &wl in &wavelengths {
            let wlp = wl - peakwl;
            let x = if wlp < 0.0 { wlp / fw } else { wlp / rw };
            let abs_x = x.abs();

            let rraw = if abs_x < r_param / 2.0 {
                1.0 - r_param / 4.0 - (1.0 / r_param) * x.powi(2)
            } else if abs_x >= r_param / 2.0 && abs_x < 1.0 - r_param / 2.0 {
                1.0 - abs_x
            } else if abs_x >= 1.0 - r_param / 2.0 && abs_x < 1.0 + r_param / 2.0 {
                1.0 / (2.0 * r_param) * (abs_x - (1.0 + r_param / 2.0)).powi(2)
            } else {
                0.0
            };

            let spd_val = min_v + (max_v - min_v) * rraw / (1.0 - rounding / 4.0);
            values.push(spd_val);
        }
        spectra.push(values);
    }

    Spectrum::new(wavelengths, spectra)
}

/// Generate monochromatic LED spectrum.
pub fn mono_led_spd(
    params: &[MonoLedParams],
    grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    if params.is_empty() {
        return Err(LuxError::EmptyInput);
    }

    let grid = grid.unwrap_or(DEFAULT_WL_GRID);
    let wavelengths = getwlr(grid)?;
    let num_peaks = params.len();
    let mut spectra = Vec::with_capacity(num_peaks);

    for param in params {
        let peakwl = param.peakwl;
        let fwhm = param.fwhm;
        let strength_shoulder = param.strength_shoulder;
        let bw_order = param.bw_order;

        let mut values = Vec::with_capacity(wavelengths.len());

        if bw_order == -2.0 {
            let spd = lorentzian2_spd(&[peakwl], &[fwhm], Some(grid))?;
            values.extend_from_slice(spd.spectra()[0].as_slice());
        } else {
            let g_spd = gaussian_spd(&[peakwl], &[fwhm], Some(grid))?;
            let g = g_spd.spectra()[0].as_slice();

            let mut ohno = Vec::with_capacity(wavelengths.len());
            for &g_val in g {
                let val = (g_val + strength_shoulder * g_val.powi(5)) / (1.0 + strength_shoulder);
                ohno.push(val);
            }

            if bw_order == -1.0 || bw_order == 0.0 {
                values.extend_from_slice(&ohno);
            } else if bw_order > 0.0 {
                let bw_spd = butterworth_spd(&[peakwl], &[fwhm], &[bw_order], Some(grid))?;
                values.extend_from_slice(bw_spd.spectra()[0].as_slice());
            } else {
                // Mix case for general negative bw_order values
                let bw_spd = butterworth_spd(&[peakwl], &[fwhm], &[bw_order], Some(grid))?;
                let bw = bw_spd.spectra()[0].as_slice();

                let lz_spd = lorentzian2_spd(&[peakwl], &[fwhm], Some(grid))?;
                let lz = lz_spd.spectra()[0].as_slice();

                for j in 0..wavelengths.len() {
                    let mut val = 0.0;
                    if bw_order >= -1.0 && bw_order <= 0.0 {
                        val += ohno[j];
                    }
                    if bw_order > 0.0 {
                        val += bw[j];
                    }
                    if bw_order >= -2.0 && bw_order < -1.0 {
                        val += lz[j];
                    }
                    values.push(val);
                }
            }
        }

        spectra.push(values);
    }

    Spectrum::new(wavelengths, spectra)
}

/// Generate phosphor LED spectrum with up to 2 phosphors.
pub fn phosphor_led_spd(
    params: &[PhosphorLedParams],
    grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    let res = phosphor_led_spd_with_components(params, grid)?;
    Ok(res.spd)
}

/// Generate phosphor LED spectrum and return both combined spectrum and component spectra.
pub fn phosphor_led_spd_with_components(
    params: &[PhosphorLedParams],
    grid: Option<WavelengthGrid>,
) -> LuxResult<PhosphorLedComponents> {
    if params.is_empty() {
        return Err(LuxError::EmptyInput);
    }

    let grid = grid.unwrap_or(DEFAULT_WL_GRID);
    let wavelengths = getwlr(grid)?;
    let num_mixtures = params.len();

    let mut combined_spectra = Vec::with_capacity(num_mixtures);

    // Determine if we have phosphors (if any params has strength_ph > 0)
    let has_phosphors = params.iter().any(|p| p.strength_ph.unwrap_or(0.0) > 0.0);
    let num_components = if has_phosphors { 3 } else { 1 };

    // component_rows will store:
    // Row 0..N: mono_led
    // Row N..2N: ph1 (if has_phosphors)
    // Row 2N..3N: ph2 (if has_phosphors)
    let mut component_rows = vec![vec![0.0; wavelengths.len()]; num_mixtures * num_components];

    for (i, param) in params.iter().enumerate() {
        // 1. Mono led component
        let mono_params = MonoLedParams {
            peakwl: param.peakwl,
            fwhm: param.fwhm,
            bw_order: param.bw_order,
            strength_shoulder: param.strength_shoulder,
        };
        let mono_spd = mono_led_spd(&[mono_params], Some(grid))?;
        let mono_led = mono_spd.spectra()[0].clone();

        let mut spd = mono_led.clone();

        if let Some(s_ph) = param.strength_ph {
            if s_ph > 0.0 && has_phosphors {
                // 2. Phosphor 1
                let ph1_params = MonoLedParams {
                    peakwl: param.peakwl_ph1,
                    fwhm: param.fwhm_ph1,
                    bw_order: -1.0,
                    strength_shoulder: 1.0,
                };
                let ph1_spd = mono_led_spd(&[ph1_params], Some(grid))?;
                let ph1 = ph1_spd.spectra()[0].clone();

                // 3. Phosphor 2
                let ph2_params = MonoLedParams {
                    peakwl: param.peakwl_ph2,
                    fwhm: param.fwhm_ph2,
                    bw_order: -1.0,
                    strength_shoulder: 1.0,
                };
                let ph2_spd = mono_led_spd(&[ph2_params], Some(grid))?;
                let ph2 = ph2_spd.spectra()[0].clone();

                // Mix phosphors
                let mut phosphors = Vec::with_capacity(wavelengths.len());
                let s_ph1 = param.strength_ph1;

                if let Some(s_ph2) = param.strength_ph2 {
                    let sum = s_ph1 + s_ph2;
                    let denom = if sum == 0.0 { 1e-300 } else { sum };
                    for j in 0..wavelengths.len() {
                        let val = (s_ph1 * ph1[j] + s_ph2 * ph2[j]) / denom + 1e-300;
                        phosphors.push(val);
                    }
                } else {
                    for j in 0..wavelengths.len() {
                        let val = s_ph1 * ph1[j] + (1.0 - s_ph1) * ph2[j] + 1e-300;
                        phosphors.push(val);
                    }
                }

                // Normalize phosphors to max = 1
                let max_ph = phosphors
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                let max_ph_val = if max_ph <= 0.0 { 1.0 } else { max_ph };
                for val in &mut phosphors {
                    *val /= max_ph_val;
                }

                // Combined spd
                for j in 0..wavelengths.len() {
                    spd[j] = mono_led[j] + s_ph * phosphors[j];
                }

                // Store in component_rows
                component_rows[i] = mono_led;
                component_rows[num_mixtures + i] = ph1;
                component_rows[2 * num_mixtures + i] = ph2;
            } else {
                component_rows[i] = mono_led;
            }
        } else {
            component_rows[i] = mono_led;
        }

        // Piecewise function modification
        if param.use_piecewise_fcn {
            // mono_led is the piecewise multiplier for wl < peakwl, else 1
            for j in 0..wavelengths.len() {
                let wl = wavelengths[j];
                let factor = if wl < param.peakwl {
                    let mono_val = component_rows[i][j]; // raw mono_led value
                    mono_val
                } else {
                    1.0
                };
                spd[j] *= factor;

                // Also apply to component spectra
                component_rows[i][j] *= factor;
                if has_phosphors && param.strength_ph.unwrap_or(0.0) > 0.0 {
                    component_rows[num_mixtures + i][j] *= factor;
                    component_rows[2 * num_mixtures + i][j] *= factor;
                }
            }
        }

        // Normalize spd to max = 1
        let max_val = spd
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let max_scale = if max_val <= 0.0 { 1.0 } else { max_val };
        for val in &mut spd {
            *val /= max_scale;
        }

        combined_spectra.push(spd);
    }

    // Normalize each row in component_rows to max = 1
    for row in &mut component_rows {
        let max_val = row
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let max_scale = if max_val <= 0.0 { 1.0 } else { max_val };
        for val in &mut *row {
            *val /= max_scale;
        }
    }

    let spd = Spectrum::new(wavelengths.clone(), combined_spectra)?;
    let components = Spectrum::new(wavelengths, component_rows)?;

    Ok(PhosphorLedComponents { spd, components })
}

fn safe_div(val: f64) -> f64 {
    if val.abs() < 1e-300 {
        if val >= 0.0 {
            1e-300
        } else {
            -1e-300
        }
    } else {
        val
    }
}

/// Calculate fluxes required to obtain target chromaticity additively mixing 3 sources.
pub fn color3mixer(
    yxy_target: [f64; 3],
    yxy1: [f64; 3],
    yxy2: [f64; 3],
    yxy3: [f64; 3],
) -> [f64; 3] {
    let y1 = yxy1[0];
    let x1 = yxy1[1];
    let y1_coord = yxy1[2];

    let y2 = yxy2[0];
    let x2 = yxy2[1];
    let y2_coord = yxy2[2];

    let y3 = yxy3[0];
    let x3 = yxy3[1];
    let y3_coord = yxy3[2];

    let yt = yxy_target[0];
    let xt = yxy_target[1];
    let yt_coord = yxy_target[2];

    let denom = (x3 - x2) * y1_coord + (x2 - x1) * y3_coord + (x1 - x3) * y2_coord;
    let m1 = y1_coord * ((xt - x3) * y2_coord - (yt_coord - y3_coord) * x2 + x3 * yt_coord - xt * y3_coord)
        / safe_div(yt_coord * denom);

    let m2 = -y2_coord * ((xt - x3) * y1_coord - (yt_coord - y3_coord) * x1 + x3 * yt_coord - xt * y3_coord)
        / safe_div(yt_coord * denom);

    let denom3 = (x2 - x1) * y3_coord - (y2_coord - y1_coord) * x3 + x1 * y2_coord - x2 * y1_coord;
    let m3 = y3_coord * ((x2 - x1) * yt_coord - (y2_coord - y1_coord) * xt + x1 * y2_coord - x2 * y1_coord)
        / safe_div(yt_coord * denom3);

    [yt * m1 / safe_div(y1), yt * m2 / safe_div(y2), yt * m3 / safe_div(y3)]
}

/// Helper to solve pseudo-inverse equation Ax = b for 3 x N matrix A (represented by its transpose, N x 3 vector).
fn solve_pseudo_inverse_3xn(a: &[[f64; 3]], b: [f64; 3]) -> Vec<f64> {
    let n = a.len();
    let mut aat = [[0.0; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            let mut sum = 0.0;
            for j in 0..n {
                sum += a[j][i] * a[j][k];
            }
            aat[i][k] = sum;
        }
    }

    // Reuse invert_matrix3 from color module
    let aat_inv = crate::color::invert_matrix3(aat);

    // Compute y = (A A^T)^{-1} b
    let mut y = [0.0; 3];
    for i in 0..3 {
        let mut sum = 0.0;
        for k in 0..3 {
            sum += aat_inv[i][k] * b[k];
        }
        y[i] = sum;
    }

    // Compute x = A^T y
    let mut x = vec![0.0; n];
    for j in 0..n {
        let mut sum = 0.0;
        for i in 0..3 {
            sum += a[j][i] * y[i];
        }
        x[j] = sum;
    }

    x
}

/// Additive color mixer of N primaries using Moore-Penrose pseudo-inverse matrix.
pub fn colormixer_pinv(
    yxy_target: [f64; 3],
    yxy_primaries: &[[f64; 3]],
    input_fmt: &str,
) -> Vec<f64> {
    let n = yxy_primaries.len();
    if input_fmt.to_lowercase() == "xyz" {
        solve_pseudo_inverse_3xn(yxy_primaries, yxy_target)
    } else {
        let yt = yxy_target[0];
        let xt = yxy_target[1];
        let yt_coord = yxy_target[2];

        let mut a_cols = vec![[0.0; 3]; n];
        for j in 0..n {
            let y_i = yxy_primaries[j][0];
            let x_i = yxy_primaries[j][1];
            let y_coord_i = yxy_primaries[j][2];

            let ratio = y_i / y_coord_i.max(1e-300);
            a_cols[j][0] = ratio * (x_i - xt);
            a_cols[j][1] = ratio * (y_coord_i - yt_coord);
            a_cols[j][2] = y_i / yt.max(1e-300);
        }

        solve_pseudo_inverse_3xn(&a_cols, [0.0, 0.0, 1.0])
    }
}

/// Calculate fluxes required to obtain target chromaticity additively mixing N light sources.
pub fn colormixer(
    yxy_target: [f64; 3],
    yxy_primaries: &[[f64; 3]],
    pair_strengths: &[f64],
) -> Vec<f64> {
    let n = yxy_primaries.len();
    if n <= 3 {
        // Fall back to color3mixer
        let p1 = if n > 0 { yxy_primaries[0] } else { [100.0, 1.0/3.0, 1.0/3.0] };
        let p2 = if n > 1 { yxy_primaries[1] } else { [100.0, 1.0/3.0, 1.0/3.0] };
        let p3 = if n > 2 { yxy_primaries[2] } else { [100.0, 1.0/3.0, 1.0/3.0] };
        let m = color3mixer(yxy_target, p1, p2, p3);
        return m.to_vec();
    }

    // mlut stores state:
    // col 0: id
    // col 1..4: Y, x, y
    // col 4, 5: parent_A, parent_B
    // col 6, 7: weight_A, weight_B
    #[derive(Debug, Clone, Copy)]
    struct LRow {
        _id: usize,
        yxy: [f64; 3],
        parent_a: usize,
        parent_b: Option<usize>,
        weight_a: f64,
        weight_b: Option<f64>,
    }

    let mut mlut = Vec::new();
    for i in 0..n {
        mlut.push(LRow {
            _id: i,
            yxy: yxy_primaries[i],
            parent_a: i,
            parent_b: None,
            weight_a: 1.0,
            weight_b: None,
        });
    }

    let mut so: Vec<usize> = (0..n).collect();
    let mut ps = pair_strengths.to_vec();
    if ps.len() < n - 3 {
        // Fill up to n-3 with 0.5 default
        ps.resize(n - 3, 0.5);
    }

    let mut k = 0;
    let mut kk = 0;
    let mut su_k = Vec::new();
    let mut sn_k = Vec::new();

    while so.len() > 3 {
        let pair_strength_ab = ps[kk];
        let p_a = so[2 * k];
        let p_b = so[2 * k + 1];

        // Mix A and B
        let yxy_a = mlut[p_a].yxy;
        let yxy_b = mlut[p_b].yxy;

        let y_a = yxy_a[0];
        let x_a = yxy_a[1];
        let y_coord_a = yxy_a[2];

        let y_b = yxy_b[0];
        let x_b = yxy_b[1];
        let y_coord_b = yxy_b[2];

        let x_val_a = x_a * y_a / y_coord_a.max(1e-300);
        let x_val_b = x_b * y_b / y_coord_b.max(1e-300);
        let z_val_a = (1.0 - x_a - y_coord_a) * y_a / y_coord_a.max(1e-300);
        let z_val_b = (1.0 - x_b - y_coord_b) * y_b / y_coord_b.max(1e-300);

        let xm = pair_strength_ab * x_val_a + (1.0 - pair_strength_ab) * x_val_b;
        let ym = pair_strength_ab * y_a + (1.0 - pair_strength_ab) * y_b;
        let zm = pair_strength_ab * z_val_a + (1.0 - pair_strength_ab) * z_val_b;

        let sum = xm + ym + zm;
        let denom = if sum == 0.0 { 1e-300 } else { sum };
        let xm_coord = xm / denom;
        let ym_coord = ym / denom;

        let new_id = mlut.len();
        mlut.push(LRow {
            _id: new_id,
            yxy: [ym, xm_coord, ym_coord],
            parent_a: p_a,
            parent_b: Some(p_b),
            weight_a: pair_strength_ab,
            weight_b: Some(1.0 - pair_strength_ab),
        });

        su_k.push(p_a);
        su_k.push(p_b);
        sn_k.push(new_id);

        let mut rem_so = Vec::new();
        for &item in &so {
            if !su_k.contains(&item) {
                rem_so.push(item);
            }
        }
        rem_so.extend(&sn_k);

        if rem_so.len() <= 3 {
            so = rem_so;
            break;
        }

        let nn = so.len() / 2;
        if k == nn - 1 {
            so = rem_so;
            su_k.clear();
            sn_k.clear();
            k = 0;
        } else {
            k += 1;
        }
        kk += 1;
    }

    // Solve color3mixer for last 3 sources
    let m3 = color3mixer(yxy_target, mlut[so[0]].yxy, mlut[so[1]].yxy, mlut[so[2]].yxy);
    if m3.iter().any(|&val| val < 0.0 || val.is_nan()) {
        return vec![f64::NAN; n];
    }

    // Backward propagation
    let mut flux_acc = vec![0.0; mlut.len()];
    flux_acc[so[0]] = m3[0];
    flux_acc[so[1]] = m3[1];
    flux_acc[so[2]] = m3[2];

    for i in (n..mlut.len()).rev() {
        let m_i = flux_acc[i];
        let p_a = mlut[i].parent_a;
        let w_a = mlut[i].weight_a;
        flux_acc[p_a] += w_a * m_i;

        if let Some(p_b) = mlut[i].parent_b {
            let w_b = mlut[i].weight_b.unwrap_or(0.0);
            flux_acc[p_b] += w_b * m_i;
        }
    }

    flux_acc[0..n].to_vec()
}

/// Helper function to convert target value to Yxy
fn target_to_yxy(target: &[f64], tar_type: &str, observer: Observer) -> LuxResult<[f64; 3]> {
    match tar_type.to_lowercase().as_str() {
        "cct" => {
            let cct = target[0];
            let xyz = cct_to_xyz(cct, observer)?;
            Ok(xyz_to_yxy(xyz))
        }
        "yxy" => {
            if target.len() < 3 {
                return Err(LuxError::EmptyInput);
            }
            Ok([target[0], target[1], target[2]])
        }
        "xyz" => {
            if target.len() < 3 {
                return Err(LuxError::EmptyInput);
            }
            Ok(xyz_to_yxy([target[0], target[1], target[2]]))
        }
        _ => Err(LuxError::InvalidGridSpec),
    }
}

/// High-level spectrum builder.
pub fn spd_builder(
    flux: Option<&[f64]>,
    component_spds: Option<&Spectrum>,
    params: &PhosphorLedParams,
    pair_strengths: Option<&[f64]>,
    target: Option<&[f64]>,
    tar_type: &str,
    observer: Observer,
    grid: Option<WavelengthGrid>,
) -> LuxResult<Spectrum> {
    let grid = grid.unwrap_or(DEFAULT_WL_GRID);
    let wavelengths = getwlr(grid)?;

    // 1. Get components
    let components = match component_spds {
        Some(s) => s.clone(),
        None => {
            let res = phosphor_led_spd_with_components(&[params.clone()], Some(grid))?;
            res.components
        }
    };

    let n_components = components.spectrum_count();

    // 2. Target optimization
    if let Some(tar) = target {
        if n_components < 3 {
            return Err(LuxError::EmptyInput);
        }

        // Calculate xyz of components
        let xyz_components = components.spd_to_xyz(&observer.standard()?, false)?;
        let mut yxy_components = Vec::with_capacity(n_components);
        for &xyz in &xyz_components {
            yxy_components.push(xyz_to_yxy(xyz));
        }

        // Convert target to Yxy
        let yxy_target = target_to_yxy(tar, tar_type, observer)?;

        // Solve for fluxes
        let m = if n_components == 3 {
            color3mixer(
                yxy_target,
                yxy_components[0],
                yxy_components[1],
                yxy_components[2],
            )
            .to_vec()
        } else {
            let p_strengths = pair_strengths.unwrap_or(&[]);
            colormixer(yxy_target, &yxy_components, p_strengths)
        };

        if m.iter().any(|&val| val.is_nan() || val < 0.0) {
            // Out of gamut: return NaN spectrum as in Python
            let spd_values = vec![f64::NAN; wavelengths.len()];
            return Spectrum::new(wavelengths, vec![spd_values]);
        }

        // Build combined spectrum: sum(M_j * component_j)
        let mut spd_values = vec![0.0; wavelengths.len()];
        for j in 0..n_components {
            let factor = m[j];
            let comp_values = &components.spectra()[j];
            for k in 0..wavelengths.len() {
                spd_values[k] += factor * comp_values[k];
            }
        }

        // Normalize combined spectrum to max = 1
        let max_val = spd_values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let max_scale = if max_val <= 0.0 { 1.0 } else { max_val };
        for val in &mut spd_values {
            *val /= max_scale;
        }

        Spectrum::new(wavelengths, vec![spd_values])
    } else {
        // No target, use flux values
        let flux_vals = flux.unwrap_or(&[]);
        if flux_vals.is_empty() {
            // Just return components
            Ok(components)
        } else {
            // Combine components using flux_vals
            let mut spd_values = vec![0.0; wavelengths.len()];
            let num_to_mix = n_components.min(flux_vals.len());
            for j in 0..num_to_mix {
                let factor = flux_vals[j];
                let comp_values = &components.spectra()[j];
                for k in 0..wavelengths.len() {
                    spd_values[k] += factor * comp_values[k];
                }
            }

            // Normalize
            let max_val = spd_values
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let max_scale = if max_val <= 0.0 { 1.0 } else { max_val };
            for val in &mut spd_values {
                *val /= max_scale;
            }

            Spectrum::new(wavelengths, vec![spd_values])
        }
    }
}

/// Fit single Gaussian peak parameters (peak wavelength, FWHM) to match target chromaticity (x, y).
pub fn fit_gaussian_spd_params(
    target_xy: [f64; 2],
    init_peak: f64,
    init_fwhm: f64,
) -> LuxResult<(f64, f64)> {
    let obj_func = |params: [f64; 2]| -> f64 {
        let peak = params[0];
        let fwhm = params[1];
        if fwhm <= 0.0 || peak < 360.0 || peak > 830.0 {
            return 1e10;
        }
        let grid = DEFAULT_WL_GRID;
        let g = match gaussian_spd(&[peak], &[fwhm], Some(grid)) {
            Ok(s) => s,
            Err(_) => return 1e10,
        };
        let observer = Observer::Cie1931_2;
        let xyz = match g.spd_to_xyz(&observer.standard().unwrap(), false) {
            Ok(v) => v[0],
            Err(_) => return 1e10,
        };
        let yxy = xyz_to_yxy(xyz);
        let x = yxy[1];
        let y = yxy[2];
        (x - target_xy[0]).powi(2) + (y - target_xy[1]).powi(2)
    };

    let opt = nelder_mead_2d(obj_func, [init_peak, init_fwhm], [2.0, 5.0], 1e-12, 1000);
    Ok((opt[0], opt[1]))
}

fn nelder_mead_2d<F>(
    mut obj_func: F,
    init: [f64; 2],
    step: [f64; 2],
    tol: f64,
    max_iter: usize,
) -> [f64; 2]
where
    F: FnMut([f64; 2]) -> f64,
{
    let p0 = init;
    let p1 = [init[0] + step[0], init[1]];
    let p2 = [init[0], init[1] + step[1]];

    let mut points = [
        (p0, obj_func(p0)),
        (p1, obj_func(p1)),
        (p2, obj_func(p2)),
    ];

    for _ in 0..max_iter {
        points.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let diff = (points[2].1 - points[0].1).abs();
        if diff < tol {
            break;
        }

        let centroid = [
            0.5 * (points[0].0[0] + points[1].0[0]),
            0.5 * (points[0].0[1] + points[1].0[1]),
        ];

        let reflected = [
            centroid[0] + 1.0 * (centroid[0] - points[2].0[0]),
            centroid[1] + 1.0 * (centroid[1] - points[2].0[1]),
        ];
        let r_val = obj_func(reflected);

        if r_val < points[1].1 && r_val >= points[0].1 {
            points[2] = (reflected, r_val);
            continue;
        }

        if r_val < points[0].1 {
            let expanded = [
                centroid[0] + 2.0 * (reflected[0] - centroid[0]),
                centroid[1] + 2.0 * (reflected[1] - centroid[1]),
            ];
            let e_val = obj_func(expanded);
            if e_val < r_val {
                points[2] = (expanded, e_val);
            } else {
                points[2] = (reflected, r_val);
            }
            continue;
        }

        let contracted = [
            centroid[0] + 0.5 * (points[2].0[0] - centroid[0]),
            centroid[1] + 0.5 * (points[2].0[1] - centroid[1]),
        ];
        let c_val = obj_func(contracted);
        if c_val < points[2].1 {
            points[2] = (contracted, c_val);
            continue;
        }

        for i in 1..3 {
            points[i].0[0] = points[0].0[0] + 0.5 * (points[i].0[0] - points[0].0[0]);
            points[i].0[1] = points[0].0[1] + 0.5 * (points[i].0[1] - points[0].0[1]);
            points[i].1 = obj_func(points[i].0);
        }
    }

    points[0].0
}

