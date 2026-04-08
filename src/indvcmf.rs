use crate::color::Matrix3;
use crate::error::{LuxError, LuxResult};
use crate::spectrum::Spectrum;

const ASANO_LMS_ABSORBANCE: &str = include_str!("../data/indvcmf/asano_cie2006_Alms.dat");
const ASANO_RELATIVE_MACULAR_DENSITY: &str =
    include_str!("../data/indvcmf/asano_cie2006_RelativeMacularDensity.dat");
const ASANO_OCULAR_DENSITY: &str = include_str!("../data/indvcmf/asano_cie2006_docul.dat");
const WAVELENGTH_START: usize = 390;
const WAVELENGTH_END: usize = 780;
const WAVELENGTH_STEP: usize = 5;
const FIELD_SIZE_MIN: f64 = 2.0;
const FIELD_SIZE_MAX: f64 = 10.0;
const S_CONE_CUTOFF: f64 = 620.0;
const LMS_TO_XYZ_2_DEG: Matrix3 = [
    [0.4151, -0.2424, 0.0425],
    [0.1355, 0.0833, -0.0043],
    [-0.0093, 0.0125, 0.2136],
];
const LMS_TO_XYZ_10_DEG: Matrix3 = [
    [0.4499, -0.2630, 0.0460],
    [0.1617, 0.0726, -0.0011],
    [-0.0036, 0.0054, 0.2291],
];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndividualObserverParameters {
    pub age: f64,
    pub field_size: f64,
    pub lens_density_variation: f64,
    pub macular_density_variation: f64,
    pub cone_density_variation: [f64; 3],
    pub cone_peak_shift: [f64; 3],
    pub allow_negative_xyz_values: bool,
}

impl Default for IndividualObserverParameters {
    fn default() -> Self {
        Self {
            age: 32.0,
            field_size: 10.0,
            lens_density_variation: 0.0,
            macular_density_variation: 0.0,
            cone_density_variation: [0.0, 0.0, 0.0],
            cone_peak_shift: [0.0, 0.0, 0.0],
            allow_negative_xyz_values: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndividualObserverStdDevs {
    pub lens_density: f64,
    pub macular_density: f64,
    pub cone_density: [f64; 3],
    pub cone_peak_shift: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
pub struct IndividualObserverCmf {
    pub lms: Spectrum,
    pub xyz: Spectrum,
    pub lens_transmission: Spectrum,
    pub macular_transmission: Spectrum,
    pub photopigment_sensitivity: Spectrum,
    pub lms_to_xyz_matrix: Matrix3,
}

pub fn individual_observer_default_std_devs() -> IndividualObserverStdDevs {
    IndividualObserverStdDevs {
        lens_density: 19.1,
        macular_density: 37.2,
        cone_density: [17.9, 17.9, 14.7],
        cone_peak_shift: [4.0, 3.0, 2.5],
    }
}

pub fn individual_observer_lms_to_xyz_matrix(field_size: f64) -> Matrix3 {
    let clamped = field_size.clamp(FIELD_SIZE_MIN, FIELD_SIZE_MAX);
    let a = (FIELD_SIZE_MAX - clamped) / (FIELD_SIZE_MAX - FIELD_SIZE_MIN);
    interpolate_matrix3(LMS_TO_XYZ_2_DEG, LMS_TO_XYZ_10_DEG, 1.0 - a)
}

pub fn individual_observer_lms_to_xyz(
    lms: &Spectrum,
    field_size: f64,
    allow_negative_values: bool,
) -> LuxResult<Spectrum> {
    if lms.spectrum_count() != 3 {
        return Err(LuxError::InvalidInput(
            "individual observer LMS input must contain exactly 3 spectra",
        ));
    }

    let matrix = individual_observer_lms_to_xyz_matrix(field_size);
    let wavelengths = lms.wavelengths().to_vec();
    let mut xyz = (0..3)
        .map(|_| Vec::with_capacity(wavelengths.len()))
        .collect::<Vec<_>>();

    for index in 0..wavelengths.len() {
        let lms_sample = [
            lms.spectra()[0][index],
            lms.spectra()[1][index],
            lms.spectra()[2][index],
        ];
        let mut xyz_sample = multiply_matrix3_vector3(matrix, lms_sample);
        if !allow_negative_values {
            for value in &mut xyz_sample {
                if *value < 0.0 {
                    *value = 0.0;
                }
            }
        }
        for axis in 0..3 {
            xyz[axis].push(xyz_sample[axis]);
        }
    }

    Spectrum::new(wavelengths, xyz)
}

pub fn individual_observer_cmf(
    parameters: IndividualObserverParameters,
) -> LuxResult<IndividualObserverCmf> {
    validate_parameters(parameters)?;

    if parameters.cone_peak_shift.iter().any(|shift| *shift != 0.0) {
        return Err(LuxError::InvalidInput(
            "non-zero cone peak shifts are not yet supported in the Rust indvcmf slice",
        ));
    }

    let wavelengths = base_wavelengths();
    let lms_absorbance = parse_columns(ASANO_LMS_ABSORBANCE, 3)?;
    let relative_macular_density = parse_columns(ASANO_RELATIVE_MACULAR_DENSITY, 1)?
        .into_iter()
        .next()
        .ok_or(LuxError::ParseError("missing macular density data"))?;
    let ocular_density = parse_columns(ASANO_OCULAR_DENSITY, 2)?;

    ensure_len(&wavelengths, &relative_macular_density)?;
    ensure_len(&wavelengths, &ocular_density[0])?;
    ensure_len(&wavelengths, &ocular_density[1])?;
    ensure_len(&wavelengths, &lms_absorbance[0])?;
    ensure_len(&wavelengths, &lms_absorbance[1])?;
    ensure_len(&wavelengths, &lms_absorbance[2])?;

    let fs = parameters.field_size;
    let peak_macular_density =
        0.485 * (-fs / 6.132).exp() * (1.0 + parameters.macular_density_variation / 100.0);
    let corrected_macular_density: Vec<f64> = relative_macular_density
        .iter()
        .map(|value| value * peak_macular_density)
        .collect::<Vec<_>>();

    let age_scale = if parameters.age <= 60.0 {
        1.0 + 0.02 * (parameters.age - 32.0)
    } else {
        1.56 + 0.0667 * (parameters.age - 60.0)
    };
    let corrected_ocular_density: Vec<f64> = ocular_density[0]
        .iter()
        .zip(ocular_density[1].iter())
        .map(|(first, second)| {
            (first * age_scale + second) * (1.0 + parameters.lens_density_variation / 100.0)
        })
        .collect::<Vec<_>>();

    let cone_peak_density = [
        (0.38 + 0.54 * (-fs / 1.333).exp()) * (1.0 + parameters.cone_density_variation[0] / 100.0),
        (0.38 + 0.54 * (-fs / 1.333).exp()) * (1.0 + parameters.cone_density_variation[1] / 100.0),
        (0.30 + 0.45 * (-fs / 1.333).exp()) * (1.0 + parameters.cone_density_variation[2] / 100.0),
    ];

    let mut alpha_lms = vec![vec![0.0; wavelengths.len()]; 3];
    for axis in 0..3 {
        for (index, wavelength) in wavelengths.iter().enumerate() {
            alpha_lms[axis][index] = 1.0
                - 10f64.powf(-cone_peak_density[axis] * 10f64.powf(lms_absorbance[axis][index]));
            if axis == 2 && *wavelength >= S_CONE_CUTOFF {
                alpha_lms[axis][index] = 0.0;
            }
        }
    }

    let mut lms = vec![vec![0.0; wavelengths.len()]; 3];
    let mut photopigment_sensitivity = vec![vec![0.0; wavelengths.len()]; 3];
    for axis in 0..3 {
        for (index, wavelength) in wavelengths.iter().enumerate() {
            let lms_quantal = alpha_lms[axis][index]
                * 10f64.powf(-corrected_macular_density[index] - corrected_ocular_density[index]);
            lms[axis][index] = lms_quantal * wavelength;
            photopigment_sensitivity[axis][index] = alpha_lms[axis][index] * wavelength;
        }
        let area: f64 = lms[axis].iter().sum();
        if area == 0.0 {
            return Err(LuxError::InvalidInput(
                "individual observer LMS normalization area must be non-zero",
            ));
        }
        for value in &mut lms[axis] {
            *value = 100.0 * *value / area;
        }
    }

    let lms_matrix = Spectrum::new(wavelengths.clone(), lms)?;
    let xyz_matrix = individual_observer_lms_to_xyz(
        &lms_matrix,
        parameters.field_size,
        parameters.allow_negative_xyz_values,
    )?;
    let lens_transmission = Spectrum::new(
        wavelengths.clone(),
        corrected_ocular_density
            .iter()
            .map(|value| 10f64.powf(-value))
            .collect::<Vec<_>>(),
    )?;
    let macular_transmission = Spectrum::new(
        wavelengths.clone(),
        corrected_macular_density
            .iter()
            .map(|value| 10f64.powf(-value))
            .collect::<Vec<_>>(),
    )?;
    let photopigment_sensitivity = Spectrum::new(wavelengths, photopigment_sensitivity)?;

    Ok(IndividualObserverCmf {
        lms: lms_matrix,
        xyz: xyz_matrix,
        lens_transmission,
        macular_transmission,
        photopigment_sensitivity,
        lms_to_xyz_matrix: individual_observer_lms_to_xyz_matrix(parameters.field_size),
    })
}

fn validate_parameters(parameters: IndividualObserverParameters) -> LuxResult<()> {
    if !parameters.age.is_finite() || parameters.age <= 0.0 {
        return Err(LuxError::InvalidInput("age must be positive and finite"));
    }
    if !parameters.field_size.is_finite()
        || parameters.field_size < FIELD_SIZE_MIN
        || parameters.field_size > FIELD_SIZE_MAX
    {
        return Err(LuxError::InvalidInput(
            "field size must be finite and within 2..=10 degrees",
        ));
    }
    if !parameters.lens_density_variation.is_finite()
        || !parameters.macular_density_variation.is_finite()
        || parameters.lens_density_variation <= -100.0
        || parameters.macular_density_variation <= -100.0
    {
        return Err(LuxError::InvalidInput(
            "lens and macular density variations must be finite and greater than -100%",
        ));
    }
    if parameters
        .cone_density_variation
        .iter()
        .any(|value| !value.is_finite() || *value <= -100.0)
    {
        return Err(LuxError::InvalidInput(
            "cone density variations must be finite and greater than -100%",
        ));
    }
    if parameters
        .cone_peak_shift
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(LuxError::InvalidInput(
            "cone peak shifts must be finite when provided",
        ));
    }
    Ok(())
}

fn base_wavelengths() -> Vec<f64> {
    (WAVELENGTH_START..=WAVELENGTH_END)
        .step_by(WAVELENGTH_STEP)
        .map(|value| value as f64)
        .collect::<Vec<_>>()
}

fn parse_columns(data: &str, expected_columns: usize) -> LuxResult<Vec<Vec<f64>>> {
    let mut columns = vec![Vec::new(); expected_columns];

    for line in data.split(['\n', '\r']) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let values: Vec<f64> = trimmed
            .split(|char: char| char == ',' || char.is_ascii_whitespace())
            .filter(|part| !part.is_empty())
            .map(|part| {
                part.parse::<f64>()
                    .map_err(|_| LuxError::ParseError("invalid indvcmf numeric value"))
            })
            .collect::<LuxResult<Vec<f64>>>()?;

        if values.len() > expected_columns {
            return Err(LuxError::ParseError("unexpected indvcmf column count"));
        }
        let mut padded_values = values;
        while padded_values.len() < expected_columns {
            padded_values.push(0.0);
        }
        for (column, value) in columns.iter_mut().zip(padded_values.into_iter()) {
            column.push(value);
        }
    }

    Ok(columns)
}

fn ensure_len(wavelengths: &[f64], values: &[f64]) -> LuxResult<()> {
    if wavelengths.len() != values.len() {
        return Err(LuxError::MismatchedLengths {
            wavelengths: wavelengths.len(),
            values: values.len(),
        });
    }
    Ok(())
}

fn interpolate_matrix3(lhs: Matrix3, rhs: Matrix3, lhs_weight: f64) -> Matrix3 {
    let rhs_weight = 1.0 - lhs_weight;
    let mut matrix = [[0.0; 3]; 3];
    for row in 0..3 {
        for col in 0..3 {
            matrix[row][col] = lhs[row][col] * lhs_weight + rhs[row][col] * rhs_weight;
        }
    }
    matrix
}

fn multiply_matrix3_vector3(matrix: Matrix3, vector: [f64; 3]) -> [f64; 3] {
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ]
}
