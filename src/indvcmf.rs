use crate::color::Matrix3;
use crate::error::{LuxError, LuxResult};
use crate::spectrum::Spectrum;

const ASANO_LMS_ABSORBANCE: &str = include_str!("../data/indvcmf/asano_cie2006_Alms.dat");
const ASANO_RELATIVE_MACULAR_DENSITY: &str =
    include_str!("../data/indvcmf/asano_cie2006_RelativeMacularDensity.dat");
const ASANO_OCULAR_DENSITY: &str = include_str!("../data/indvcmf/asano_cie2006_docul.dat");
const ASANO_CAT_OBSERVER_FACTORS: &str = include_str!("../data/indvcmf/asano_CatObsPfctr.dat");
const ASANO_US_CENSUS_AGE_DISTRIBUTION: &str =
    include_str!("../data/indvcmf/asano_USCensus2010Population.dat");
const CIETC197_ABSORBANCES: &str = include_str!("../data/indvcmf/cietc197_absorbances0_1nm.dat");
const CIETC197_DOCUL2: &str = include_str!("../data/indvcmf/cietc197_docul2.dat");
const CIE2006_XYZ_2_DEG: &str = include_str!("../data/cmfs/ciexyz_2006_2.dat");
const CIE2006_XYZ_10_DEG: &str = include_str!("../data/cmfs/ciexyz_2006_10.dat");

const WAVELENGTH_START: usize = 390;
const WAVELENGTH_END: usize = 780;
const WAVELENGTH_STEP: usize = 5;
const FIELD_SIZE_MIN: f64 = 2.0;
const FIELD_SIZE_MAX: f64 = 10.0;
const S_CONE_CUTOFF: f64 = 620.0;
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;
const LCG_INCREMENT: u64 = 1;

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
const STOCKMAN2023_LMS_TO_XYZ_2_DEG: Matrix3 = [
    [1.947_354_69, -1.414_451_23, 0.364_763_27],
    [0.689_902_72, 0.348_321_89, 0.0],
    [0.0, 0.0, 1.934_853_43],
];
const STOCKMAN2023_LMS_TO_XYZ_10_DEG: Matrix3 = [
    [1.939_864_43, -1.346_643_59, 0.430_449_35],
    [0.692_839_32, 0.349_675_67, 0.0],
    [0.0, 0.0, 2.146_879_45],
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndividualObserverDataSource {
    Asano,
    CieTc197,
    Stockman2023,
    AicomPlus,
}

impl Default for IndividualObserverDataSource {
    fn default() -> Self {
        Self::Asano
    }
}

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

#[derive(Debug, Clone, PartialEq)]
pub struct IndividualObserverMonteCarloOptions {
    pub n_observers: usize,
    pub field_size: f64,
    pub age_pool: Vec<f64>,
    pub std_devs: IndividualObserverStdDevs,
    pub use_germany_scale_factors: bool,
    pub allow_negative_xyz_values: bool,
    pub data_source: IndividualObserverDataSource,
    pub seed: u64,
}

impl Default for IndividualObserverMonteCarloOptions {
    fn default() -> Self {
        Self {
            n_observers: 1,
            field_size: 10.0,
            age_pool: vec![32.0],
            std_devs: individual_observer_default_std_devs(),
            use_germany_scale_factors: true,
            allow_negative_xyz_values: false,
            data_source: IndividualObserverDataSource::Asano,
            seed: 0xDEC0DED,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct IndividualObserverPopulation {
    pub parameters: Vec<IndividualObserverParameters>,
    pub cmfs: Vec<IndividualObserverCmf>,
}

pub type IndividualObserverModel = IndividualObserverDataSource;

#[derive(Debug, Clone, PartialEq)]
pub struct IndividualObserverRequest {
    pub model: IndividualObserverModel,
    pub parameters: IndividualObserverParameters,
}

impl Default for IndividualObserverRequest {
    fn default() -> Self {
        Self {
            model: IndividualObserverModel::Asano,
            parameters: IndividualObserverParameters::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct IndividualObserverCategoricalOptions {
    pub n_categories: usize,
    pub field_size: f64,
    pub allow_negative_xyz_values: bool,
}

impl Default for IndividualObserverCategoricalOptions {
    fn default() -> Self {
        Self {
            n_categories: 10,
            field_size: 2.0,
            allow_negative_xyz_values: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum IndividualObserverPopulationStrategy {
    MonteCarlo(IndividualObserverMonteCarloOptions),
    Categorical(IndividualObserverCategoricalOptions),
}

#[derive(Debug, Clone, PartialEq)]
pub struct IndividualObserverPopulationRequest {
    pub model: IndividualObserverModel,
    pub strategy: IndividualObserverPopulationStrategy,
}

impl Default for IndividualObserverPopulationRequest {
    fn default() -> Self {
        Self {
            model: IndividualObserverModel::Asano,
            strategy: IndividualObserverPopulationStrategy::MonteCarlo(
                IndividualObserverMonteCarloOptions::default(),
            ),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct ObserverSourceData {
    wavelengths: Vec<f64>,
    lms_absorbance: [Vec<f64>; 3],
    relative_macular_density: Vec<f64>,
    ocular_density: [Vec<f64>; 2],
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

pub fn individual_observer_lms_to_xyz_matrix_stockman2023(field_size: f64) -> Matrix3 {
    let clamped = field_size.clamp(FIELD_SIZE_MIN, FIELD_SIZE_MAX);
    let alpha_10 = (clamped - FIELD_SIZE_MIN) / (FIELD_SIZE_MAX - FIELD_SIZE_MIN);
    let alpha_2 = 1.0 - alpha_10;
    let mut matrix = [[0.0; 3]; 3];
    for row in 0..3 {
        for col in 0..3 {
            matrix[row][col] = STOCKMAN2023_LMS_TO_XYZ_2_DEG[row][col] * alpha_2
                + STOCKMAN2023_LMS_TO_XYZ_10_DEG[row][col] * alpha_10;
        }
    }
    matrix
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
    individual_observer_lms_to_xyz_with_matrix(lms, matrix, allow_negative_values)
}

fn individual_observer_lms_to_xyz_with_matrix(
    lms: &Spectrum,
    matrix: Matrix3,
    allow_negative_values: bool,
) -> LuxResult<Spectrum> {
    if lms.spectrum_count() != 3 {
        return Err(LuxError::InvalidInput(
            "individual observer LMS input must contain exactly 3 spectra",
        ));
    }

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
    individual_observer_cmf_with_source(parameters, IndividualObserverDataSource::Asano)
}

pub fn individual_observer_generate(
    request: IndividualObserverRequest,
) -> LuxResult<IndividualObserverCmf> {
    individual_observer_cmf_with_source(request.parameters, request.model)
}

pub fn individual_observer_cmf_stockman2023(
    parameters: IndividualObserverParameters,
) -> LuxResult<IndividualObserverCmf> {
    individual_observer_cmf_with_source(parameters, IndividualObserverDataSource::Stockman2023)
}

pub fn individual_observer_cmf_aicom_plus(
    parameters: IndividualObserverParameters,
) -> LuxResult<IndividualObserverCmf> {
    individual_observer_cmf_with_source(parameters, IndividualObserverDataSource::AicomPlus)
}

pub fn individual_observer_cmf_with_source(
    parameters: IndividualObserverParameters,
    data_source: IndividualObserverDataSource,
) -> LuxResult<IndividualObserverCmf> {
    validate_parameters(parameters)?;

    let source_data = load_source_data(data_source)?;
    if data_source == IndividualObserverDataSource::Stockman2023 {
        return compute_stockman2023_observer(parameters, source_data);
    }
    if data_source == IndividualObserverDataSource::AicomPlus {
        return compute_aicom_plus_observer(parameters, source_data);
    }
    let wavelengths = source_data.wavelengths;
    let relative_macular_density = source_data.relative_macular_density;
    let ocular_density = source_data.ocular_density;
    let lms_absorbance = source_data.lms_absorbance;

    ensure_len(&wavelengths, &relative_macular_density)?;
    ensure_len(&wavelengths, &ocular_density[0])?;
    ensure_len(&wavelengths, &ocular_density[1])?;
    ensure_len(&wavelengths, &lms_absorbance[0])?;
    ensure_len(&wavelengths, &lms_absorbance[1])?;
    ensure_len(&wavelengths, &lms_absorbance[2])?;

    let shifted_absorbance = (0..3)
        .map(|axis| {
            shift_series_with_linear_extrapolation(
                &wavelengths,
                &lms_absorbance[axis],
                parameters.cone_peak_shift[axis],
            )
        })
        .collect::<LuxResult<Vec<Vec<f64>>>>()?;

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
                - 10f64
                    .powf(-cone_peak_density[axis] * 10f64.powf(shifted_absorbance[axis][index]));
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
    let lms_to_xyz_matrix = match data_source {
        IndividualObserverDataSource::Asano => {
            individual_observer_lms_to_xyz_matrix(parameters.field_size)
        }
        IndividualObserverDataSource::CieTc197 => {
            fit_cietc197_lms_to_xyz_matrix(&lms_matrix, parameters.field_size)?
        }
        IndividualObserverDataSource::Stockman2023 => {
            individual_observer_lms_to_xyz_matrix(parameters.field_size)
        }
        IndividualObserverDataSource::AicomPlus => {
            fit_cietc197_lms_to_xyz_matrix(&lms_matrix, parameters.field_size)?
        }
    };
    let xyz_matrix = individual_observer_lms_to_xyz_with_matrix(
        &lms_matrix,
        lms_to_xyz_matrix,
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
        lms_to_xyz_matrix,
    })
}

pub fn individual_observer_us_census_age_distribution() -> LuxResult<Vec<f64>> {
    let mut ages = Vec::new();

    for line in ASANO_US_CENSUS_AGE_DISTRIBUTION.split(['\n', '\r']) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("Age") {
            continue;
        }
        let fields = split_numeric_tokens(trimmed);
        if fields.len() < 2 {
            return Err(LuxError::ParseError("invalid US census age row"));
        }
        let age = fields[0];
        let population = fields[1];
        if !(10.0..=70.0).contains(&age) {
            continue;
        }
        let repeats = (population / 1000.0).round();
        if repeats <= 0.0 {
            continue;
        }
        for _ in 0..repeats as usize {
            ages.push(age);
        }
    }

    if ages.is_empty() {
        return Err(LuxError::ParseError("empty US census age distribution"));
    }
    Ok(ages)
}

pub fn individual_observer_monte_carlo_parameters(
    options: &IndividualObserverMonteCarloOptions,
) -> LuxResult<Vec<IndividualObserverParameters>> {
    validate_monte_carlo_options(options)?;

    let mut std_devs = options.std_devs;
    if options.use_germany_scale_factors {
        std_devs = scale_std_devs(std_devs);
    }

    let mut rng = LcgRng::new(options.seed);
    let mut parameters = Vec::with_capacity(options.n_observers);

    for _ in 0..options.n_observers {
        let age = draw_age(&options.age_pool, &mut rng)?;
        let lens_density_variation =
            (std_devs.lens_density * rng.next_standard_normal()).max(-100.0);
        let macular_density_variation =
            (std_devs.macular_density * rng.next_standard_normal()).max(-100.0);
        let cone_density_variation = [
            (std_devs.cone_density[0] * rng.next_standard_normal()).max(-100.0),
            (std_devs.cone_density[1] * rng.next_standard_normal()).max(-100.0),
            (std_devs.cone_density[2] * rng.next_standard_normal()).max(-100.0),
        ];
        let cone_peak_shift = [
            std_devs.cone_peak_shift[0] * rng.next_standard_normal(),
            std_devs.cone_peak_shift[1] * rng.next_standard_normal(),
            std_devs.cone_peak_shift[2] * rng.next_standard_normal(),
        ];

        parameters.push(IndividualObserverParameters {
            age,
            field_size: options.field_size,
            lens_density_variation,
            macular_density_variation,
            cone_density_variation,
            cone_peak_shift,
            allow_negative_xyz_values: options.allow_negative_xyz_values,
        });
    }

    Ok(parameters)
}

pub fn individual_observer_monte_carlo(
    options: IndividualObserverMonteCarloOptions,
) -> LuxResult<IndividualObserverPopulation> {
    let parameters = individual_observer_monte_carlo_parameters(&options)?;
    let cmfs = parameters
        .iter()
        .copied()
        .map(|params| individual_observer_cmf_with_source(params, options.data_source))
        .collect::<LuxResult<Vec<_>>>()?;

    Ok(IndividualObserverPopulation { parameters, cmfs })
}

pub fn individual_observer_generate_population(
    request: IndividualObserverPopulationRequest,
) -> LuxResult<IndividualObserverPopulation> {
    match request.strategy {
        IndividualObserverPopulationStrategy::MonteCarlo(mut options) => {
            options.data_source = request.model;
            individual_observer_monte_carlo(options)
        }
        IndividualObserverPopulationStrategy::Categorical(options) => {
            individual_observer_categorical_observers(
                options.n_categories,
                options.field_size,
                request.model,
                options.allow_negative_xyz_values,
            )
        }
    }
}

pub fn individual_observer_categorical_observers(
    n_categories: usize,
    field_size: f64,
    data_source: IndividualObserverDataSource,
    allow_negative_xyz_values: bool,
) -> LuxResult<IndividualObserverPopulation> {
    if n_categories == 0 {
        return Err(LuxError::InvalidInput("category count must be positive"));
    }
    if !field_size.is_finite() || !(FIELD_SIZE_MIN..=FIELD_SIZE_MAX).contains(&field_size) {
        return Err(LuxError::InvalidInput(
            "field size must be finite and within 2..=10 degrees",
        ));
    }

    let (ages, factors) = parse_categorical_observer_factors()?;
    let count = n_categories.min(ages.len());

    let parameters = (0..count)
        .map(|index| IndividualObserverParameters {
            age: ages[index],
            field_size,
            lens_density_variation: factors[0][index],
            macular_density_variation: factors[1][index],
            cone_density_variation: [factors[2][index], factors[3][index], factors[4][index]],
            cone_peak_shift: [factors[5][index], factors[6][index], factors[7][index]],
            allow_negative_xyz_values,
        })
        .collect::<Vec<_>>();

    let cmfs = parameters
        .iter()
        .copied()
        .map(|params| individual_observer_cmf_with_source(params, data_source))
        .collect::<LuxResult<Vec<_>>>()?;

    Ok(IndividualObserverPopulation { parameters, cmfs })
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

fn validate_monte_carlo_options(options: &IndividualObserverMonteCarloOptions) -> LuxResult<()> {
    if options.n_observers == 0 {
        return Err(LuxError::InvalidInput("observer count must be positive"));
    }
    if options.age_pool.is_empty() {
        return Err(LuxError::InvalidInput("age pool cannot be empty"));
    }
    for age in &options.age_pool {
        if !age.is_finite() || *age <= 0.0 {
            return Err(LuxError::InvalidInput(
                "age pool entries must be positive and finite",
            ));
        }
    }

    // Validate field size and variation bounds through the shared parameter validator.
    validate_parameters(IndividualObserverParameters {
        age: options.age_pool[0],
        field_size: options.field_size,
        lens_density_variation: 0.0,
        macular_density_variation: 0.0,
        cone_density_variation: [0.0, 0.0, 0.0],
        cone_peak_shift: [0.0, 0.0, 0.0],
        allow_negative_xyz_values: options.allow_negative_xyz_values,
    })
}

fn scale_std_devs(std_devs: IndividualObserverStdDevs) -> IndividualObserverStdDevs {
    IndividualObserverStdDevs {
        lens_density: std_devs.lens_density * 0.98,
        macular_density: std_devs.macular_density * 0.98,
        cone_density: [
            std_devs.cone_density[0] * 0.5,
            std_devs.cone_density[1] * 0.5,
            std_devs.cone_density[2] * 0.5,
        ],
        cone_peak_shift: [
            std_devs.cone_peak_shift[0] * 0.5,
            std_devs.cone_peak_shift[1] * 0.5,
            std_devs.cone_peak_shift[2] * 0.5,
        ],
    }
}

fn draw_age(age_pool: &[f64], rng: &mut LcgRng) -> LuxResult<f64> {
    if age_pool.is_empty() {
        return Err(LuxError::InvalidInput("age pool cannot be empty"));
    }
    let index = rng.next_usize(age_pool.len());
    Ok(age_pool[index])
}

fn base_wavelengths() -> Vec<f64> {
    (WAVELENGTH_START..=WAVELENGTH_END)
        .step_by(WAVELENGTH_STEP)
        .map(|value| value as f64)
        .collect::<Vec<_>>()
}

fn load_source_data(source: IndividualObserverDataSource) -> LuxResult<ObserverSourceData> {
    match source {
        IndividualObserverDataSource::Asano
        | IndividualObserverDataSource::Stockman2023
        | IndividualObserverDataSource::AicomPlus => {
            let wavelengths = base_wavelengths();
            let lms_absorbance_columns = parse_columns(ASANO_LMS_ABSORBANCE, 3)?;
            let rmd_columns = parse_columns(ASANO_RELATIVE_MACULAR_DENSITY, 1)?;
            let docul_columns = parse_columns(ASANO_OCULAR_DENSITY, 2)?;

            Ok(ObserverSourceData {
                wavelengths,
                lms_absorbance: [
                    lms_absorbance_columns[0].clone(),
                    lms_absorbance_columns[1].clone(),
                    lms_absorbance_columns[2].clone(),
                ],
                relative_macular_density: rmd_columns[0].clone(),
                ocular_density: [docul_columns[0].clone(), docul_columns[1].clone()],
            })
        }
        IndividualObserverDataSource::CieTc197 => load_cietc197_source_data(),
    }
}

fn load_cietc197_source_data() -> LuxResult<ObserverSourceData> {
    let mut wavelengths = Vec::new();
    let mut l_absorbance = Vec::new();
    let mut m_absorbance = Vec::new();
    let mut s_absorbance = Vec::new();
    let mut ocular_sum_32 = Vec::new();
    let mut relative_macular = Vec::new();

    for line in CIETC197_ABSORBANCES.split(['\n', '\r']) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let fields = parse_csv_row_with_empty(trimmed)?;
        if fields.len() < 7 {
            return Err(LuxError::ParseError("invalid cietc197 absorbance row"));
        }

        let wl = fields[0].ok_or(LuxError::ParseError("missing cietc197 wavelength"))?;
        let l = fields[2].ok_or(LuxError::ParseError("missing cietc197 L absorbance"))?;
        let m = fields[3].ok_or(LuxError::ParseError("missing cietc197 M absorbance"))?;
        let s = fields[4].unwrap_or(f64::NEG_INFINITY);
        let ocular_sum = fields[5].ok_or(LuxError::ParseError("missing cietc197 ocular sum"))?;
        let macula = fields[6].ok_or(LuxError::ParseError("missing cietc197 macular density"))?;

        wavelengths.push(wl);
        l_absorbance.push(l);
        m_absorbance.push(m);
        s_absorbance.push(s);
        ocular_sum_32.push(ocular_sum);
        relative_macular.push(macula / 0.35);
    }

    ensure_strictly_increasing(&wavelengths)?;

    if let Some(first_invalid) = s_absorbance.iter().position(|value| !value.is_finite()) {
        for value in &mut s_absorbance[first_invalid..] {
            *value = f64::NEG_INFINITY;
        }
    }

    let (docul2_wavelengths, docul2_values) = parse_two_column_table(CIETC197_DOCUL2)?;
    let interpolated_docul2 = wavelengths
        .iter()
        .map(|wl| interpolate_linear_with_extrapolation(&docul2_wavelengths, &docul2_values, *wl))
        .collect::<Vec<_>>();
    let docul1 = ocular_sum_32
        .iter()
        .zip(interpolated_docul2.iter())
        .map(|(sum_32, second)| sum_32 - second)
        .collect::<Vec<_>>();

    Ok(ObserverSourceData {
        wavelengths,
        lms_absorbance: [l_absorbance, m_absorbance, s_absorbance],
        relative_macular_density: relative_macular,
        ocular_density: [docul1, interpolated_docul2],
    })
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

fn parse_two_column_table(data: &str) -> LuxResult<(Vec<f64>, Vec<f64>)> {
    let mut first = Vec::new();
    let mut second = Vec::new();

    for line in data.split(['\n', '\r']) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let values = split_numeric_tokens(trimmed);
        if values.len() < 2 {
            return Err(LuxError::ParseError("invalid two-column table"));
        }
        first.push(values[0]);
        second.push(values[1]);
    }

    ensure_strictly_increasing(&first)?;
    Ok((first, second))
}

fn parse_three_column_table(data: &str) -> LuxResult<(Vec<f64>, [Vec<f64>; 3])> {
    let mut wavelengths = Vec::new();
    let mut first = Vec::new();
    let mut second = Vec::new();
    let mut third = Vec::new();

    for line in data.split(['\n', '\r']) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let values = split_numeric_tokens(trimmed);
        if values.len() < 4 {
            return Err(LuxError::ParseError("invalid three-column table"));
        }
        wavelengths.push(values[0]);
        first.push(values[1]);
        second.push(values[2]);
        third.push(values[3]);
    }
    ensure_strictly_increasing(&wavelengths)?;
    Ok((wavelengths, [first, second, third]))
}

fn split_numeric_tokens(line: &str) -> Vec<f64> {
    line.split(|char: char| char == ',' || char.is_ascii_whitespace())
        .filter(|part| !part.is_empty())
        .filter_map(|part| part.parse::<f64>().ok())
        .collect::<Vec<_>>()
}

fn parse_csv_row_with_empty(line: &str) -> LuxResult<Vec<Option<f64>>> {
    line.split(',')
        .map(|cell| {
            let trimmed = cell.trim();
            if trimmed.is_empty() {
                Ok(None)
            } else {
                trimmed
                    .parse::<f64>()
                    .map(Some)
                    .map_err(|_| LuxError::ParseError("invalid cietc197 numeric value"))
            }
        })
        .collect::<LuxResult<Vec<Option<f64>>>>()
}

fn parse_categorical_observer_factors() -> LuxResult<(Vec<f64>, [Vec<f64>; 8])> {
    let mut rows = Vec::new();
    for line in ASANO_CAT_OBSERVER_FACTORS.split(['\n', '\r']) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        rows.push(split_numeric_tokens(trimmed));
    }

    if rows.len() < 9 {
        return Err(LuxError::ParseError(
            "invalid categorical observer factor table",
        ));
    }

    let category_count = rows[0].len();
    if category_count == 0 {
        return Err(LuxError::ParseError(
            "empty categorical observer factor table",
        ));
    }
    for row in &rows[1..9] {
        if row.len() != category_count {
            return Err(LuxError::ParseError(
                "categorical observer factor table has inconsistent row lengths",
            ));
        }
    }

    let ages = rows[0].clone();
    let factors = [
        rows[1].clone(),
        rows[2].clone(),
        rows[3].clone(),
        rows[4].clone(),
        rows[5].clone(),
        rows[6].clone(),
        rows[7].clone(),
        rows[8].clone(),
    ];

    Ok((ages, factors))
}

fn compute_stockman2023_observer(
    parameters: IndividualObserverParameters,
    source_data: ObserverSourceData,
) -> LuxResult<IndividualObserverCmf> {
    let wavelengths = source_data.wavelengths;
    let relative_macular_density = source_data.relative_macular_density;
    let ocular_density = source_data.ocular_density;
    let lms_absorbance = source_data.lms_absorbance;

    ensure_len(&wavelengths, &relative_macular_density)?;
    ensure_len(&wavelengths, &ocular_density[0])?;
    ensure_len(&wavelengths, &ocular_density[1])?;
    ensure_len(&wavelengths, &lms_absorbance[0])?;
    ensure_len(&wavelengths, &lms_absorbance[1])?;
    ensure_len(&wavelengths, &lms_absorbance[2])?;

    let shifted_absorbance = (0..3)
        .map(|axis| {
            shift_series_log_wavelength(
                &wavelengths,
                &lms_absorbance[axis],
                parameters.cone_peak_shift[axis],
            )
        })
        .collect::<LuxResult<Vec<Vec<f64>>>>()?;

    let fs = parameters.field_size.clamp(FIELD_SIZE_MIN, FIELD_SIZE_MAX);
    let alpha_10 = (fs - FIELD_SIZE_MIN) / (FIELD_SIZE_MAX - FIELD_SIZE_MIN);
    let alpha_2 = 1.0 - alpha_10;
    let cone_peak_density = [
        (alpha_2 * 0.50 + alpha_10 * 0.38) * (1.0 + parameters.cone_density_variation[0] / 100.0),
        (alpha_2 * 0.50 + alpha_10 * 0.38) * (1.0 + parameters.cone_density_variation[1] / 100.0),
        (alpha_2 * 0.40 + alpha_10 * 0.30) * (1.0 + parameters.cone_density_variation[2] / 100.0),
    ];
    let k_mac =
        (alpha_2 * 1.0 + alpha_10 * 0.271) * (1.0 + parameters.macular_density_variation / 100.0);
    let k_lens = 1.0 * (1.0 + parameters.lens_density_variation / 100.0);

    let corrected_macular_density = relative_macular_density
        .iter()
        .map(|value| value * k_mac)
        .collect::<Vec<_>>();
    let corrected_ocular_density = ocular_density[0]
        .iter()
        .zip(ocular_density[1].iter())
        .map(|(first, second)| (first + second) * k_lens)
        .collect::<Vec<_>>();

    let mut absorptance = vec![vec![0.0; wavelengths.len()]; 3];
    for axis in 0..3 {
        for (index, wavelength) in wavelengths.iter().enumerate() {
            absorptance[axis][index] = 1.0
                - 10f64
                    .powf(-cone_peak_density[axis] * 10f64.powf(shifted_absorbance[axis][index]));
            if axis == 2 && *wavelength >= S_CONE_CUTOFF {
                absorptance[axis][index] = 0.0;
            }
        }
    }

    let mut lms = vec![vec![0.0; wavelengths.len()]; 3];
    let mut photopigment_sensitivity = vec![vec![0.0; wavelengths.len()]; 3];
    for axis in 0..3 {
        for (index, wavelength) in wavelengths.iter().enumerate() {
            let quantal = absorptance[axis][index]
                * 10f64.powf(-(corrected_macular_density[index] + corrected_ocular_density[index]));
            lms[axis][index] = quantal * wavelength;
            photopigment_sensitivity[axis][index] = absorptance[axis][index] * wavelength;
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
    let lms_to_xyz_matrix =
        individual_observer_lms_to_xyz_matrix_stockman2023(parameters.field_size);
    let xyz_matrix = individual_observer_lms_to_xyz_with_matrix(
        &lms_matrix,
        lms_to_xyz_matrix,
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
        lms_to_xyz_matrix,
    })
}

fn compute_aicom_plus_observer(
    parameters: IndividualObserverParameters,
    source_data: ObserverSourceData,
) -> LuxResult<IndividualObserverCmf> {
    let wavelengths = source_data.wavelengths;
    let relative_macular_density = source_data.relative_macular_density;
    let ocular_density = source_data.ocular_density;
    let lms_absorbance = source_data.lms_absorbance;

    ensure_len(&wavelengths, &relative_macular_density)?;
    ensure_len(&wavelengths, &ocular_density[0])?;
    ensure_len(&wavelengths, &ocular_density[1])?;
    ensure_len(&wavelengths, &lms_absorbance[0])?;
    ensure_len(&wavelengths, &lms_absorbance[1])?;
    ensure_len(&wavelengths, &lms_absorbance[2])?;

    let shifted_absorbance = (0..3)
        .map(|axis| {
            shift_series_with_linear_extrapolation(
                &wavelengths,
                &lms_absorbance[axis],
                parameters.cone_peak_shift[axis],
            )
        })
        .collect::<LuxResult<Vec<Vec<f64>>>>()?;

    let fs = parameters.field_size;
    let peak_macular_density =
        0.485 * (-fs / 6.132).exp() * (1.0 + parameters.macular_density_variation / 100.0);
    let corrected_macular_density: Vec<f64> = relative_macular_density
        .iter()
        .map(|value| value * peak_macular_density)
        .collect::<Vec<_>>();

    // AICOM+ adaption: replace the default ocular media template with the CIE 203 style template.
    let cie203_docul = cie203_ocular_density_at(&wavelengths)?;
    let age_scale = if parameters.age <= 60.0 {
        1.0 + 0.02 * (parameters.age - 32.0)
    } else {
        1.56 + 0.0667 * (parameters.age - 60.0)
    };
    let corrected_ocular_density: Vec<f64> = ocular_density[0]
        .iter()
        .zip(cie203_docul.iter())
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
                - 10f64
                    .powf(-cone_peak_density[axis] * 10f64.powf(shifted_absorbance[axis][index]));
            if axis == 2 && *wavelength >= S_CONE_CUTOFF {
                alpha_lms[axis][index] = 0.0;
            }
        }
    }

    let mut lms_raw = vec![vec![0.0; wavelengths.len()]; 3];
    let mut lms_normalized = vec![vec![0.0; wavelengths.len()]; 3];
    let mut photopigment_sensitivity = vec![vec![0.0; wavelengths.len()]; 3];
    for axis in 0..3 {
        for (index, wavelength) in wavelengths.iter().enumerate() {
            let lms_quantal = alpha_lms[axis][index]
                * 10f64.powf(-corrected_macular_density[index] - corrected_ocular_density[index]);
            let energy = lms_quantal * wavelength;
            lms_raw[axis][index] = energy;
            lms_normalized[axis][index] = energy;
            photopigment_sensitivity[axis][index] = alpha_lms[axis][index] * wavelength;
        }
        let area: f64 = lms_normalized[axis].iter().sum();
        if area == 0.0 {
            return Err(LuxError::InvalidInput(
                "individual observer LMS normalization area must be non-zero",
            ));
        }
        for value in &mut lms_normalized[axis] {
            *value = 100.0 * *value / area;
        }
    }

    // AICOM+ adaption: skip LMS normalization before LMS->XYZ conversion.
    let lms_for_xyz = Spectrum::new(wavelengths.clone(), lms_raw)?;
    let lms_output = Spectrum::new(wavelengths.clone(), lms_normalized)?;
    let lms_to_xyz_matrix = fit_cietc197_lms_to_xyz_matrix(&lms_for_xyz, parameters.field_size)?;
    let xyz_matrix = individual_observer_lms_to_xyz_with_matrix(
        &lms_for_xyz,
        lms_to_xyz_matrix,
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
        lms: lms_output,
        xyz: xyz_matrix,
        lens_transmission,
        macular_transmission,
        photopigment_sensitivity,
        lms_to_xyz_matrix,
    })
}

fn cie203_ocular_density_at(wavelengths: &[f64]) -> LuxResult<Vec<f64>> {
    let (docul2_wavelengths, docul2_values) = parse_two_column_table(CIETC197_DOCUL2)?;
    Ok(wavelengths
        .iter()
        .map(|wl| interpolate_linear_with_extrapolation(&docul2_wavelengths, &docul2_values, *wl))
        .collect::<Vec<_>>())
}

fn shift_series_log_wavelength(
    wavelengths: &[f64],
    values: &[f64],
    shift_nm: f64,
) -> LuxResult<Vec<f64>> {
    ensure_len(wavelengths, values)?;
    if wavelengths.is_empty() {
        return Ok(Vec::new());
    }
    if wavelengths.len() == 1 || shift_nm == 0.0 {
        return Ok(values.to_vec());
    }

    let peak_index = values
        .iter()
        .enumerate()
        .filter(|(_, value)| value.is_finite())
        .max_by(|lhs, rhs| {
            lhs.1
                .partial_cmp(rhs.1)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(index, _)| index)
        .ok_or(LuxError::InvalidInput(
            "cannot infer cone absorbance peak for log-shift",
        ))?;
    let lambda_max = wavelengths[peak_index];
    let lambda_max_shifted = lambda_max + shift_nm;
    if lambda_max_shifted <= 0.0 {
        return Err(LuxError::InvalidInput(
            "cone peak shift leads to non-positive shifted lambda_max",
        ));
    }
    let scale = lambda_max / lambda_max_shifted;

    let mut shifted = Vec::with_capacity(wavelengths.len());
    for wavelength in wavelengths {
        let query = wavelength * scale;
        let mut value = interpolate_linear_with_extrapolation(wavelengths, values, query);
        if !value.is_finite() {
            value = f64::NEG_INFINITY;
        }
        shifted.push(value);
    }
    Ok(shifted)
}

fn fit_cietc197_lms_to_xyz_matrix(lms: &Spectrum, field_size: f64) -> LuxResult<Matrix3> {
    let target_xyz = cie2006_xyz_reference(field_size, lms.wavelengths())?;
    fit_lms_to_xyz_matrix(lms, &target_xyz)
}

fn cie2006_xyz_reference(field_size: f64, wavelengths: &[f64]) -> LuxResult<[Vec<f64>; 3]> {
    let (wl2, xyz2) = parse_three_column_table(CIE2006_XYZ_2_DEG)?;
    let (wl10, xyz10) = parse_three_column_table(CIE2006_XYZ_10_DEG)?;

    let fs = field_size.clamp(FIELD_SIZE_MIN, FIELD_SIZE_MAX);
    let alpha_10 = (fs - FIELD_SIZE_MIN) / (FIELD_SIZE_MAX - FIELD_SIZE_MIN);
    let alpha_2 = 1.0 - alpha_10;

    let mut out = [Vec::new(), Vec::new(), Vec::new()];
    for wl in wavelengths {
        let x2 = interpolate_linear_with_extrapolation(&wl2, &xyz2[0], *wl);
        let y2 = interpolate_linear_with_extrapolation(&wl2, &xyz2[1], *wl);
        let z2 = interpolate_linear_with_extrapolation(&wl2, &xyz2[2], *wl);
        let x10 = interpolate_linear_with_extrapolation(&wl10, &xyz10[0], *wl);
        let y10 = interpolate_linear_with_extrapolation(&wl10, &xyz10[1], *wl);
        let z10 = interpolate_linear_with_extrapolation(&wl10, &xyz10[2], *wl);

        out[0].push(alpha_2 * x2 + alpha_10 * x10);
        out[1].push(alpha_2 * y2 + alpha_10 * y10);
        out[2].push(alpha_2 * z2 + alpha_10 * z10);
    }
    Ok(out)
}

fn fit_lms_to_xyz_matrix(lms: &Spectrum, target_xyz: &[Vec<f64>; 3]) -> LuxResult<Matrix3> {
    if lms.spectrum_count() != 3 {
        return Err(LuxError::InvalidInput(
            "individual observer LMS input must contain exactly 3 spectra",
        ));
    }
    for channel in target_xyz {
        ensure_len(lms.wavelengths(), channel)?;
    }

    let l = &lms.spectra()[0];
    let m = &lms.spectra()[1];
    let s = &lms.spectra()[2];

    let mut ata = [[0.0; 3]; 3];
    for index in 0..l.len() {
        let row = [l[index], m[index], s[index]];
        for r in 0..3 {
            for c in 0..3 {
                ata[r][c] += row[r] * row[c];
            }
        }
    }

    let mut atb_rows = [[0.0; 3]; 3];
    for row in 0..3 {
        for index in 0..l.len() {
            let target = target_xyz[row][index];
            atb_rows[row][0] += l[index] * target;
            atb_rows[row][1] += m[index] * target;
            atb_rows[row][2] += s[index] * target;
        }
    }

    let trace = ata[0][0] + ata[1][1] + ata[2][2];
    let base_lambda = (trace / 3.0).max(1e-18);

    for attempt in 0..10 {
        let lambda = base_lambda * 10f64.powi(attempt - 12);
        let regularized = [
            [ata[0][0] + lambda, ata[0][1], ata[0][2]],
            [ata[1][0], ata[1][1] + lambda, ata[1][2]],
            [ata[2][0], ata[2][1], ata[2][2] + lambda],
        ];
        if let Some(ata_inverse) = invert_3x3(regularized) {
            let mut matrix = [[0.0; 3]; 3];
            for row in 0..3 {
                matrix[row] = multiply_matrix3_vector3(ata_inverse, atb_rows[row]);
            }
            if matrix
                .iter()
                .flat_map(|row| row.iter())
                .all(|value| value.is_finite())
            {
                return Ok(matrix);
            }
        }
    }

    Err(LuxError::InvalidInput(
        "unable to solve stable cietc197 lms->xyz matrix fit",
    ))
}

fn invert_3x3(matrix: Matrix3) -> Option<Matrix3> {
    let m = matrix;
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-20 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ])
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

fn ensure_strictly_increasing(values: &[f64]) -> LuxResult<()> {
    for pair in values.windows(2) {
        if pair[1] <= pair[0] {
            return Err(LuxError::NonMonotonicWavelengths);
        }
    }
    Ok(())
}

fn shift_series_with_linear_extrapolation(
    wavelengths: &[f64],
    values: &[f64],
    shift_nm: f64,
) -> LuxResult<Vec<f64>> {
    ensure_len(wavelengths, values)?;

    if wavelengths.is_empty() {
        return Ok(Vec::new());
    }
    if wavelengths.len() == 1 {
        return Ok(vec![values[0]]);
    }

    let mut shifted = Vec::with_capacity(wavelengths.len());
    for wavelength in wavelengths {
        let mut value =
            interpolate_linear_with_extrapolation(wavelengths, values, wavelength - shift_nm);
        if !value.is_finite() {
            value = f64::NEG_INFINITY;
        }
        shifted.push(value);
    }
    Ok(shifted)
}

fn interpolate_linear_with_extrapolation(x: &[f64], y: &[f64], query: f64) -> f64 {
    if query <= x[0] {
        return linear_segment(x[0], y[0], x[1], y[1], query);
    }
    let last = x.len() - 1;
    if query >= x[last] {
        return linear_segment(x[last - 1], y[last - 1], x[last], y[last], query);
    }

    for index in 0..last {
        if query <= x[index + 1] {
            return linear_segment(x[index], y[index], x[index + 1], y[index + 1], query);
        }
    }
    y[last]
}

fn linear_segment(x0: f64, y0: f64, x1: f64, y1: f64, query: f64) -> f64 {
    if x1 == x0 {
        return y0;
    }
    y0 + (query - x0) * (y1 - y0) / (x1 - x0)
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

#[derive(Debug, Clone)]
struct LcgRng {
    state: u64,
}

impl LcgRng {
    fn new(seed: u64) -> Self {
        let state = if seed == 0 { 1 } else { seed };
        Self { state }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(LCG_MULTIPLIER)
            .wrapping_add(LCG_INCREMENT);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        // [0, 1)
        let value = self.next_u64() >> 11;
        (value as f64) / ((1u64 << 53) as f64)
    }

    fn next_usize(&mut self, upper_bound: usize) -> usize {
        if upper_bound <= 1 {
            return 0;
        }
        (self.next_u64() as usize) % upper_bound
    }

    fn next_standard_normal(&mut self) -> f64 {
        let mut u1 = self.next_f64();
        while u1 <= f64::MIN_POSITIVE {
            u1 = self.next_f64();
        }
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        r * theta.cos()
    }
}

// --- Measured individual observer support ---

pub const LMS_TO_XYZ_2DEG_FIXED: Matrix3 = [
    [ 1.94735469, -1.41445123,  0.36476327],
    [ 0.68990272,  0.34832189,  0.00000000],
    [ 0.00000000,  0.00000000,  1.93485343],
];

pub const LMS_TO_XYZ_10DEG_FIXED: Matrix3 = [
    [ 1.93906444, -1.37420684,  0.39960583],
    [ 0.69784280,  0.34538187,  0.00000000],
    [ 0.00000000,  0.00000000,  2.03077344],
];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndividualObserverMeasuredParameters {
    pub lshift: f64,
    pub mshift: f64,
    pub sshift: f64,
    pub lod: f64,
    pub mod_: f64,
    pub sod: f64,
    pub mac: f64,
    pub lens: f64,
    pub field_size: f64,
}

const LSER_COEFFS: [f64; 18] = [
    -42.417608560, -2.656791612, 75.011093607, 56.477062776, 7.509397607, 9.061442173,
    -38.068488495, -20.974610259, -6.642746250, -3.785039126, 9.322071459, 3.134494745,
    1.603799055, 0.439302358, -0.676958684, -0.072988371, -0.078857510, -0.004264105
];

const MCONE_COEFFS: [f64; 18] = [
    -210.6568853069, -0.1458073553, 386.7319763250, 305.4710584670, 5.0218382813, 6.8386224350,
    -208.2062335724, -118.4890200521, -5.7625866330, -3.7973553168, 55.1803460639, 19.9728512548,
    1.8990456325, 0.6913410864, -5.0891806213, -0.7070689492, -0.1419926703, 0.0005894876
];

const SCONE_COEFFS: [f64; 18] = [
    207.3880950935, -6.3065623516, -393.7100478026, -315.6650602846, 19.2917535553, 19.6414743488,
    214.2211570447, 121.8584683485, -15.1820737886, -8.6774057156, -56.7596380441, -20.6318720369,
    3.6934875040, 1.0483022480, 5.3656615075, 0.7898783086, -0.1480357836, 0.0002358232
];

const MACULAR_COEFFS: [f64; 24] = [
    3712.2037792986, 374.1811575175, -7007.6989637831, -5887.2857515364, -633.0475233043,
    -716.0429039473, 4386.8811254914, 2882.1092658881, 638.1347550701, 468.4980700497,
    -1653.7567388120, -817.1240899995, -286.4038978705, -144.7996457395, 340.3364828167,
    115.5652804221, 59.1650826447, 18.6678197694, -30.2344535413, -5.4683753172, -4.1335064207,
    -0.5043959566, 0.5094171266, 1.0050048550
];

const LENS_COEFFS: [f64; 20] = [
    -313.9508632762, -70.3216819666, 585.4719725809, 471.5395862431, 117.3539102044,
    127.0168222865, -324.4700544731, -188.1638078982, -104.5512488013, -68.3078486904,
    89.7815373733, 33.4498264952, 35.2723638870, 13.6524086627, -8.7568168893, -1.2825766708,
    -3.5126531075, -0.4477840959, 0.0428291365, 1.0091871745
];

fn evaluate_fourier_series(x: f64, c: &[f64]) -> f64 {
    let mut sum = c[0];
    for k in 1..=8 {
        let k_f = k as f64;
        sum += c[2 * k - 1] * (k_f * x).cos() + c[2 * k] * (k_f * x).sin();
    }
    sum + c[17]
}

fn evaluate_macular_fourier(x: f64, c: &[f64; 24]) -> f64 {
    let mut sum = c[0];
    for k in 1..=11 {
        let k_f = k as f64;
        sum += c[2 * k - 1] * (k_f * x).cos() + c[2 * k] * (k_f * x).sin();
    }
    sum * c[23]
}

fn evaluate_lens_fourier(x: f64, c: &[f64; 20]) -> f64 {
    let mut sum = c[0];
    for k in 1..=9 {
        let k_f = k as f64;
        sum += c[2 * k - 1] * (k_f * x).cos() + c[2 * k] * (k_f * x).sin();
    }
    sum * c[19]
}

fn lserconelog(nm: f64, lshift: f64) -> f64 {
    let x = (nm.log10() - 2.5563025007672873) / 0.11876664675818423;
    let xshift = (553.1 / (553.1 + lshift)).log10() / 0.11876664675818423;
    evaluate_fourier_series(x + xshift, &LSER_COEFFS)
}

fn mconelog(nm: f64, mshift: f64) -> f64 {
    let x = (nm.log10() - 2.5563025007672873) / 0.11876664675818423;
    let xshift = (529.9 / (529.9 + mshift)).log10() / 0.11876664675818423;
    evaluate_fourier_series(x + xshift, &MCONE_COEFFS)
}

fn sconelog(nm: f64, sshift: f64) -> f64 {
    let x = (nm.log10() - 2.5563025007672873) / 0.11876664675818423;
    let xshift = (416.9 / (416.9 + sshift)).log10() / 0.11876664675818423;
    evaluate_fourier_series(x + xshift, &SCONE_COEFFS)
}

fn macular_density_template(nm: f64) -> f64 {
    let x = (nm - 375.0) / 55.70423008;
    if x >= 0.0 && x <= ((550.0 - 375.0) / 55.70423008) {
        evaluate_macular_fourier(x, &MACULAR_COEFFS)
    } else {
        0.0
    }
}

fn lens_density_template(nm: f64) -> f64 {
    let x = (nm - 360.0) / 95.49296586;
    if x >= 0.0 && x <= ((660.0 - 360.0) / 95.49296586) {
        evaluate_lens_fourier(x, &LENS_COEFFS)
    } else {
        0.0
    }
}

pub fn individual_observer_cmf_from_measured(
    wavelengths: &[f64],
    params: IndividualObserverMeasuredParameters,
    custom_matrix: Option<Matrix3>,
) -> LuxResult<IndividualObserverCmf> {
    if wavelengths.is_empty() {
        return Err(LuxError::EmptyInput);
    }

    let n = wavelengths.len();
    let mut lms_quantal = vec![vec![0.0; n]; 3];
    let mut lens_trans = Vec::with_capacity(n);
    let mut macular_trans = Vec::with_capacity(n);

    let mac_scale = params.mac / 0.35;
    let lens_scale = params.lens / 1.7649;

    for (i, &w) in wavelengths.iter().enumerate() {
        let l_log_abs = lserconelog(w, params.lshift);
        let m_log_abs = mconelog(w, params.mshift);
        let s_log_abs = sconelog(w, params.sshift);

        let l_abs = 10f64.powf(l_log_abs);
        let m_abs = 10f64.powf(m_log_abs);
        let s_abs = 10f64.powf(s_log_abs);

        let l_retina = (1.0 - 10f64.powf(-params.lod * l_abs)) / (1.0 - 10f64.powf(-params.lod));
        let m_retina = (1.0 - 10f64.powf(-params.mod_ * m_abs)) / (1.0 - 10f64.powf(-params.mod_));
        let s_retina = (1.0 - 10f64.powf(-params.sod * s_abs)) / (1.0 - 10f64.powf(-params.sod));

        let mac_temp = macular_density_template(w);
        let lens_temp = lens_density_template(w);

        let mac_d = mac_temp * mac_scale;
        let lens_d = lens_temp * lens_scale;

        let mac_transmission = 10f64.powf(-mac_d);
        let lens_transmission = 10f64.powf(-lens_d);

        lens_trans.push(lens_transmission);
        macular_trans.push(mac_transmission);

        let mac_lens_factor = mac_transmission * lens_transmission;
        lms_quantal[0][i] = l_retina * mac_lens_factor;
        lms_quantal[1][i] = m_retina * mac_lens_factor;
        lms_quantal[2][i] = s_retina * mac_lens_factor;
    }

    for axis in 0..3 {
        let max_val = lms_quantal[axis].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if max_val > 0.0 {
            for val in &mut lms_quantal[axis] {
                *val /= max_val;
            }
        }
    }

    let mut lms_energy = vec![vec![0.0; n]; 3];
    let mut photopigment = vec![vec![0.0; n]; 3];
    for axis in 0..3 {
        for i in 0..n {
            lms_energy[axis][i] = lms_quantal[axis][i] * wavelengths[i];
        }
        let max_val = lms_energy[axis].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if max_val > 0.0 {
            for val in &mut lms_energy[axis] {
                *val /= max_val;
            }
        }

        for i in 0..n {
            let w = wavelengths[i];
            let abs_log = match axis {
                0 => lserconelog(w, params.lshift),
                1 => mconelog(w, params.mshift),
                2 => sconelog(w, params.sshift),
                _ => 0.0,
            };
            let abs_lin = 10f64.powf(abs_log);
            let od = match axis {
                0 => params.lod,
                1 => params.mod_,
                2 => params.sod,
                _ => 0.0,
            };
            let retina = (1.0 - 10f64.powf(-od * abs_lin)) / (1.0 - 10f64.powf(-od));
            photopigment[axis][i] = retina * w;
        }
    }

    let lms_spectrum = Spectrum::new(wavelengths.to_vec(), lms_energy)?;

    let matrix = custom_matrix.unwrap_or_else(|| {
        if params.field_size <= 4.0 {
            LMS_TO_XYZ_2DEG_FIXED
        } else {
            LMS_TO_XYZ_10DEG_FIXED
        }
    });

    let mut xyz_values = vec![vec![0.0; n]; 3];
    for i in 0..n {
        let lms_val = [
            lms_spectrum.spectra()[0][i],
            lms_spectrum.spectra()[1][i],
            lms_spectrum.spectra()[2][i],
        ];
        let xyz_val = multiply_matrix3_vector3(matrix, lms_val);
        for axis in 0..3 {
            xyz_values[axis][i] = xyz_val[axis];
        }
    }

    let xyz_spectrum = Spectrum::new(wavelengths.to_vec(), xyz_values)?;
    let lens_transmission = Spectrum::new(wavelengths.to_vec(), lens_trans)?;
    let macular_transmission = Spectrum::new(wavelengths.to_vec(), macular_trans)?;
    let photopigment_sensitivity = Spectrum::new(wavelengths.to_vec(), photopigment)?;

    Ok(IndividualObserverCmf {
        lms: lms_spectrum,
        xyz: xyz_spectrum,
        lens_transmission,
        macular_transmission,
        photopigment_sensitivity,
        lms_to_xyz_matrix: matrix,
    })
}

