use crate::color::{CatTransform, Matrix3, Observer};
use crate::error::{LuxError, LuxResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CamModel {
    Ciecam02,
    Cam16,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CamSurround {
    Average,
    Dim,
    Dark,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CamSurroundParameters {
    pub c: f64,
    pub nc: f64,
    pub f: f64,
    pub fll: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CamNakaRushtonParameters {
    pub n: f64,
    pub sig: f64,
    pub scaling: f64,
    pub noise: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CamViewingConditions {
    pub model: CamModel,
    pub white_point: [f64; 3],
    pub luminance_factor_white: f64,
    pub adapting_luminance: f64,
    pub background_luminance: f64,
    pub surround: CamSurround,
    pub surround_parameters: CamSurroundParameters,
    pub cat_transform: CatTransform,
    pub degree_of_adaptation: f64,
    pub naka_rushton_parameters: CamNakaRushtonParameters,
    pub k: f64,
    pub fl: f64,
    pub n: f64,
    pub nbb: f64,
    pub ncb: f64,
    pub z: f64,
    pub rgb_white: [f64; 3],
    pub adapted_rgb_white: [f64; 3],
    pub cone_white: [f64; 3],
    pub compressed_cone_white: [f64; 3],
    pub achromatic_response_white: f64,
    pub brightness_white: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CamAppearance {
    pub lightness: f64,
    pub brightness: f64,
    pub chroma: f64,
    pub colorfulness: f64,
    pub saturation: f64,
    pub hue_angle: f64,
    pub a_m: f64,
    pub b_m: f64,
    pub a_c: f64,
    pub b_c: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CamUcsType {
    Ucs,
    Lcd,
    Scd,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CamUcsParameters {
    pub k_l: f64,
    pub c1: f64,
    pub c2: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CamUcsAppearance {
    pub j_prime: f64,
    pub a_prime: f64,
    pub b_prime: f64,
}

impl CamModel {
    pub fn default_cat_transform(self) -> CatTransform {
        match self {
            Self::Ciecam02 => CatTransform::Cat02,
            Self::Cam16 => CatTransform::Cat16,
        }
    }

    pub fn default_naka_rushton_parameters(self) -> CamNakaRushtonParameters {
        match self {
            Self::Ciecam02 | Self::Cam16 => CamNakaRushtonParameters {
                n: 0.42,
                sig: 27.13_f64.powf(1.0 / 0.42),
                scaling: 400.0,
                noise: 0.1,
            },
        }
    }
}

impl CamSurround {
    pub fn parameters(self) -> CamSurroundParameters {
        match self {
            Self::Average => CamSurroundParameters {
                c: 0.69,
                nc: 1.0,
                f: 1.0,
                fll: 1.0,
            },
            Self::Dim => CamSurroundParameters {
                c: 0.59,
                nc: 0.9,
                f: 0.9,
                fll: 1.0,
            },
            Self::Dark => CamSurroundParameters {
                c: 0.525,
                nc: 0.8,
                f: 0.8,
                fll: 1.0,
            },
        }
    }
}

impl CamUcsType {
    pub fn parameters(self) -> CamUcsParameters {
        match self {
            Self::Ucs => CamUcsParameters {
                k_l: 1.0,
                c1: 0.007,
                c2: 0.0228,
            },
            Self::Lcd => CamUcsParameters {
                k_l: 0.77,
                c1: 0.007,
                c2: 0.0053,
            },
            Self::Scd => CamUcsParameters {
                k_l: 1.24,
                c1: 0.007,
                c2: 0.0363,
            },
        }
    }
}

impl CamViewingConditions {
    pub fn new(
        model: CamModel,
        white_point: [f64; 3],
        luminance_factor_white: Option<f64>,
        adapting_luminance: f64,
        background_luminance: f64,
        surround: CamSurround,
        degree_of_adaptation: Option<f64>,
        cat_transform: Option<CatTransform>,
    ) -> LuxResult<Self> {
        validate_xyz_triplet(white_point, "white point values must be finite")?;
        validate_positive_finite(
            adapting_luminance,
            "adapting_luminance must be finite and positive",
        )?;
        validate_positive_finite(
            background_luminance,
            "background_luminance must be finite and positive",
        )?;

        let yw = luminance_factor_white.unwrap_or(white_point[1]);
        validate_positive_finite(yw, "luminance_factor_white must be finite and positive")?;

        let surround_parameters = surround.parameters();
        let d = degree_of_adaptation.unwrap_or_else(|| {
            let raw = surround_parameters.f
                * (1.0 - (1.0 / 3.6) * ((-adapting_luminance - 42.0) / 92.0).exp());
            clamp(raw, 0.0, 1.0)
        });
        validate_degree(d, "degree_of_adaptation must be finite and within 0..=1")?;

        let cat_transform = cat_transform.unwrap_or(model.default_cat_transform());
        let naka_rushton_parameters = model.default_naka_rushton_parameters();
        let normalized_white = normalize_white_point(white_point, yw)?;
        let rgb_white = multiply_matrix3_vector3(cat_transform.matrix(), normalized_white);
        let adapted_rgb_white = apply_von_kries_to_white(rgb_white, yw, d)?;
        let cone_white = match model {
            CamModel::Ciecam02 => {
                let hpe = Observer::Cie1931_2.xyz_to_lms_matrix()?;
                let cone_matrix = multiply_matrix3(hpe, invert_matrix3(cat_transform.matrix()));
                multiply_matrix3_vector3(cone_matrix, adapted_rgb_white)
            }
            CamModel::Cam16 => adapted_rgb_white,
        };

        let k = 1.0 / (5.0 * adapting_luminance + 1.0);
        let fl = 0.2 * k.powi(4) * (5.0 * adapting_luminance)
            + 0.1 * (1.0 - k.powi(4)).powi(2) * (5.0 * adapting_luminance).powf(1.0 / 3.0);
        let n = background_luminance / yw;
        let nbb = 0.725 * (1.0 / n).powf(0.2);
        let ncb = nbb;
        let z = 1.48 + surround_parameters.fll * n.sqrt();
        let compressed_cone_white =
            cam_naka_rushton_scaled(cone_white, fl / 100.0, naka_rushton_parameters, true);
        let achromatic_response_white = (2.0 * compressed_cone_white[0]
            + compressed_cone_white[1]
            + (1.0 / 20.0) * compressed_cone_white[2]
            - 0.305)
            * nbb;
        let brightness_white =
            (4.0 / surround_parameters.c) * (achromatic_response_white + 4.0) * fl.powf(0.25);

        Ok(Self {
            model,
            white_point,
            luminance_factor_white: yw,
            adapting_luminance,
            background_luminance,
            surround,
            surround_parameters,
            cat_transform,
            degree_of_adaptation: d,
            naka_rushton_parameters,
            k,
            fl,
            n,
            nbb,
            ncb,
            z,
            rgb_white,
            adapted_rgb_white,
            cone_white,
            compressed_cone_white,
            achromatic_response_white,
            brightness_white,
        })
    }

    pub fn forward(self, xyz: [f64; 3]) -> LuxResult<CamAppearance> {
        cam_forward(xyz, self)
    }

    pub fn forward_ucs(self, xyz: [f64; 3], ucs_type: CamUcsType) -> LuxResult<CamUcsAppearance> {
        cam_ucs_forward(xyz, self, ucs_type)
    }

    pub fn inverse(self, appearance: CamAppearance) -> LuxResult<[f64; 3]> {
        cam_inverse(appearance, self)
    }

    pub fn inverse_ucs(
        self,
        appearance: CamUcsAppearance,
        ucs_type: CamUcsType,
    ) -> LuxResult<[f64; 3]> {
        cam_ucs_inverse(appearance, self, ucs_type)
    }
}

pub fn cam_naka_rushton(value: f64, parameters: CamNakaRushtonParameters, forward: bool) -> f64 {
    if forward {
        value.signum() * parameters.scaling * value.abs().powf(parameters.n)
            / (value.abs().powf(parameters.n) + parameters.sig.powf(parameters.n))
            + parameters.noise
    } else {
        let z = (value - parameters.noise) / parameters.scaling;
        let sign = z.signum();
        sign * (((z.abs() * parameters.sig.powf(parameters.n)) / (sign - z)).abs())
            .powf(1.0 / parameters.n)
    }
}

pub fn cam16_viewing_conditions(
    white_point: [f64; 3],
    luminance_factor_white: Option<f64>,
    adapting_luminance: f64,
    background_luminance: f64,
    surround: CamSurround,
    degree_of_adaptation: Option<f64>,
    cat_transform: Option<CatTransform>,
) -> LuxResult<CamViewingConditions> {
    CamViewingConditions::new(
        CamModel::Cam16,
        white_point,
        luminance_factor_white,
        adapting_luminance,
        background_luminance,
        surround,
        degree_of_adaptation,
        cat_transform,
    )
}

pub fn ciecam02_viewing_conditions(
    white_point: [f64; 3],
    luminance_factor_white: Option<f64>,
    adapting_luminance: f64,
    background_luminance: f64,
    surround: CamSurround,
    degree_of_adaptation: Option<f64>,
    cat_transform: Option<CatTransform>,
) -> LuxResult<CamViewingConditions> {
    CamViewingConditions::new(
        CamModel::Ciecam02,
        white_point,
        luminance_factor_white,
        adapting_luminance,
        background_luminance,
        surround,
        degree_of_adaptation,
        cat_transform,
    )
}

pub fn cam_forward(xyz: [f64; 3], conditions: CamViewingConditions) -> LuxResult<CamAppearance> {
    validate_xyz_triplet(xyz, "XYZ stimulus values must be finite")?;

    let normalized_xyz = scale_xyz_to_luminance_factor(
        xyz,
        conditions.white_point[1],
        conditions.luminance_factor_white,
    )?;
    let rgb = multiply_matrix3_vector3(conditions.cat_transform.matrix(), normalized_xyz);
    let adapted_rgb = apply_von_kries_to_channels(
        rgb,
        conditions.rgb_white,
        conditions.luminance_factor_white,
        conditions.degree_of_adaptation,
    )?;
    let cone =
        cat_sensor_to_cone_responses(conditions.model, conditions.cat_transform, adapted_rgb)?;
    let compressed = cam_naka_rushton_scaled(
        cone,
        conditions.fl / 100.0,
        conditions.naka_rushton_parameters,
        true,
    );
    let achromatic_response =
        (2.0 * compressed[0] + compressed[1] + compressed[2] / 20.0 - 0.305) * conditions.nbb;

    let a = compressed[0] - 12.0 * compressed[1] / 11.0 + compressed[2] / 11.0;
    let b = (compressed[0] + compressed[1] - 2.0 * compressed[2]) / 9.0;
    let hue_angle = hue_angle_degrees(a, b);
    let eccentricity_factor = (f64::cos(degrees_to_radians(hue_angle) + 2.0) + 3.8) / 4.0;
    let lightness = 100.0
        * (achromatic_response / conditions.achromatic_response_white)
            .powf(conditions.surround_parameters.c * conditions.z);
    let brightness = (4.0 / conditions.surround_parameters.c)
        * (lightness / 100.0).sqrt()
        * (conditions.achromatic_response_white + 4.0)
        * conditions.fl.powf(0.25);
    let t = ((50000.0 / 13.0)
        * conditions.surround_parameters.nc
        * conditions.ncb
        * eccentricity_factor
        * (a * a + b * b).sqrt())
        / (compressed[0] + compressed[1] + (21.0 / 20.0) * compressed[2]);
    let chroma =
        t.powf(0.9) * (lightness / 100.0).sqrt() * (1.64 - 0.29_f64.powf(conditions.n)).powf(0.73);
    let colorfulness = chroma * conditions.fl.powf(0.25);
    let saturation = 100.0 * (colorfulness / brightness).sqrt();
    let hue_radians = degrees_to_radians(hue_angle);

    Ok(CamAppearance {
        lightness,
        brightness,
        chroma,
        colorfulness,
        saturation,
        hue_angle,
        a_m: colorfulness * hue_radians.cos(),
        b_m: colorfulness * hue_radians.sin(),
        a_c: chroma * hue_radians.cos(),
        b_c: chroma * hue_radians.sin(),
    })
}

pub fn cam16_forward(xyz: [f64; 3], conditions: CamViewingConditions) -> LuxResult<CamAppearance> {
    if conditions.model != CamModel::Cam16 {
        return Err(LuxError::InvalidInput(
            "cam16_forward requires CamModel::Cam16 viewing conditions",
        ));
    }
    cam_forward(xyz, conditions)
}

pub fn ciecam02_forward(
    xyz: [f64; 3],
    conditions: CamViewingConditions,
) -> LuxResult<CamAppearance> {
    if conditions.model != CamModel::Ciecam02 {
        return Err(LuxError::InvalidInput(
            "ciecam02_forward requires CamModel::Ciecam02 viewing conditions",
        ));
    }
    cam_forward(xyz, conditions)
}

pub fn cam_inverse(
    appearance: CamAppearance,
    conditions: CamViewingConditions,
) -> LuxResult<[f64; 3]> {
    validate_positive_finite(
        conditions.luminance_factor_white,
        "luminance_factor_white must be finite and positive",
    )?;

    let lightness = appearance.lightness;
    let hue_angle = hue_angle_degrees(appearance.a_m, appearance.b_m);
    let colorfulness = (appearance.a_m * appearance.a_m + appearance.b_m * appearance.b_m).sqrt();
    let chroma = colorfulness / conditions.fl.powf(0.25);

    let (a, b) = if chroma <= 1e-15 || lightness <= 1e-15 {
        (0.0, 0.0)
    } else {
        let t = (chroma
            / ((lightness / 100.0).sqrt() * (1.64 - 0.29_f64.powf(conditions.n)).powf(0.73)))
        .powf(1.0 / 0.9);
        let eccentricity_factor = (f64::cos(degrees_to_radians(hue_angle) + 2.0) + 3.8) / 4.0;
        let achromatic_response = conditions.achromatic_response_white
            * (lightness / 100.0).powf(1.0 / (conditions.surround_parameters.c * conditions.z));
        let at = f64::cos(degrees_to_radians(hue_angle));
        let bt = f64::sin(degrees_to_radians(hue_angle));
        let p1 = (50000.0 / 13.0)
            * conditions.surround_parameters.nc
            * conditions.ncb
            * eccentricity_factor
            / t;
        let p2 = achromatic_response / conditions.nbb + 0.305;
        let p3 = 21.0 / 20.0;

        if bt.abs() >= at.abs() {
            let b = p2 * (2.0 + p3) * (460.0 / 1403.0)
                / (p1 / bt + (2.0 + p3) * (220.0 / 1403.0) * (at / bt) - (27.0 / 1403.0)
                    + p3 * (6300.0 / 1403.0));
            (b * (at / bt), b)
        } else {
            let a = p2 * (2.0 + p3) * (460.0 / 1403.0)
                / (p1 / at + (2.0 + p3) * (220.0 / 1403.0)
                    - ((27.0 / 1403.0) - p3 * (6300.0 / 1403.0)) * (bt / at));
            (a, a * (bt / at))
        }
    };

    let p2 = if lightness <= 1e-15 {
        0.305
    } else {
        let achromatic_response = conditions.achromatic_response_white
            * (lightness / 100.0).powf(1.0 / (conditions.surround_parameters.c * conditions.z));
        achromatic_response / conditions.nbb + 0.305
    };
    let compressed = [
        (460.0 * p2 + 451.0 * a + 288.0 * b) / 1403.0,
        (460.0 * p2 - 891.0 * a - 261.0 * b) / 1403.0,
        (460.0 * p2 - 220.0 * a - 6300.0 * b) / 1403.0,
    ];
    let cone = [
        (100.0 / conditions.fl)
            * cam_naka_rushton(compressed[0], conditions.naka_rushton_parameters, false),
        (100.0 / conditions.fl)
            * cam_naka_rushton(compressed[1], conditions.naka_rushton_parameters, false),
        (100.0 / conditions.fl)
            * cam_naka_rushton(compressed[2], conditions.naka_rushton_parameters, false),
    ];
    let adapted_rgb =
        cone_responses_to_cat_sensor(conditions.model, conditions.cat_transform, cone)?;
    let rgb = invert_von_kries_from_channels(
        adapted_rgb,
        conditions.rgb_white,
        conditions.luminance_factor_white,
        conditions.degree_of_adaptation,
    )?;
    let normalized_xyz =
        multiply_matrix3_vector3(invert_matrix3(conditions.cat_transform.matrix()), rgb);
    scale_xyz_to_luminance_factor(
        normalized_xyz,
        conditions.luminance_factor_white,
        conditions.white_point[1],
    )
}

pub fn cam_ucs_forward(
    xyz: [f64; 3],
    conditions: CamViewingConditions,
    ucs_type: CamUcsType,
) -> LuxResult<CamUcsAppearance> {
    let appearance = cam_forward(xyz, conditions)?;
    Ok(cam_ucs_from_appearance(appearance, ucs_type))
}

pub fn cam16_ucs_forward(
    xyz: [f64; 3],
    conditions: CamViewingConditions,
    ucs_type: CamUcsType,
) -> LuxResult<CamUcsAppearance> {
    if conditions.model != CamModel::Cam16 {
        return Err(LuxError::InvalidInput(
            "cam16_ucs_forward requires CamModel::Cam16 viewing conditions",
        ));
    }
    cam_ucs_forward(xyz, conditions, ucs_type)
}

pub fn ciecam02_ucs_forward(
    xyz: [f64; 3],
    conditions: CamViewingConditions,
    ucs_type: CamUcsType,
) -> LuxResult<CamUcsAppearance> {
    if conditions.model != CamModel::Ciecam02 {
        return Err(LuxError::InvalidInput(
            "ciecam02_ucs_forward requires CamModel::Ciecam02 viewing conditions",
        ));
    }
    cam_ucs_forward(xyz, conditions, ucs_type)
}

pub fn cam_ucs_inverse(
    appearance: CamUcsAppearance,
    conditions: CamViewingConditions,
    ucs_type: CamUcsType,
) -> LuxResult<[f64; 3]> {
    cam_inverse(cam_appearance_from_ucs(appearance, ucs_type), conditions)
}

pub fn cam16_ucs_inverse(
    appearance: CamUcsAppearance,
    conditions: CamViewingConditions,
    ucs_type: CamUcsType,
) -> LuxResult<[f64; 3]> {
    if conditions.model != CamModel::Cam16 {
        return Err(LuxError::InvalidInput(
            "cam16_ucs_inverse requires CamModel::Cam16 viewing conditions",
        ));
    }
    cam_ucs_inverse(appearance, conditions, ucs_type)
}

pub fn ciecam02_ucs_inverse(
    appearance: CamUcsAppearance,
    conditions: CamViewingConditions,
    ucs_type: CamUcsType,
) -> LuxResult<[f64; 3]> {
    if conditions.model != CamModel::Ciecam02 {
        return Err(LuxError::InvalidInput(
            "ciecam02_ucs_inverse requires CamModel::Ciecam02 viewing conditions",
        ));
    }
    cam_ucs_inverse(appearance, conditions, ucs_type)
}

fn cam_naka_rushton_scaled(
    values: [f64; 3],
    scale: f64,
    parameters: CamNakaRushtonParameters,
    forward: bool,
) -> [f64; 3] {
    [
        cam_naka_rushton(values[0] * scale, parameters, forward),
        cam_naka_rushton(values[1] * scale, parameters, forward),
        cam_naka_rushton(values[2] * scale, parameters, forward),
    ]
}

fn cam_ucs_from_appearance(appearance: CamAppearance, ucs_type: CamUcsType) -> CamUcsAppearance {
    let parameters = ucs_type.parameters();
    let j_prime = (1.0 + 100.0 * parameters.c1) * appearance.lightness
        / (1.0 + parameters.c1 * appearance.lightness);
    let m_prime = if parameters.c2 == 0.0 {
        appearance.colorfulness
    } else {
        (1.0 / parameters.c2) * f64::ln(1.0 + parameters.c2 * appearance.colorfulness)
    };
    let hue_radians = degrees_to_radians(appearance.hue_angle);
    CamUcsAppearance {
        j_prime,
        a_prime: m_prime * hue_radians.cos(),
        b_prime: m_prime * hue_radians.sin(),
    }
}

fn cam_appearance_from_ucs(appearance: CamUcsAppearance, ucs_type: CamUcsType) -> CamAppearance {
    let parameters = ucs_type.parameters();
    let hue_angle = hue_angle_degrees(appearance.a_prime, appearance.b_prime);
    let m_prime =
        (appearance.a_prime * appearance.a_prime + appearance.b_prime * appearance.b_prime).sqrt();
    let colorfulness = if parameters.c2 == 0.0 {
        m_prime
    } else {
        (f64::exp(parameters.c2 * m_prime) - 1.0) / parameters.c2
    };
    let lightness = appearance.j_prime / (1.0 + (100.0 - appearance.j_prime) * parameters.c1);
    let chroma = 0.0;
    let saturation = 0.0;
    let hue_radians = degrees_to_radians(hue_angle);
    CamAppearance {
        lightness,
        brightness: 0.0,
        chroma,
        colorfulness,
        saturation,
        hue_angle,
        a_m: colorfulness * hue_radians.cos(),
        b_m: colorfulness * hue_radians.sin(),
        a_c: chroma * hue_radians.cos(),
        b_c: chroma * hue_radians.sin(),
    }
}

fn cat_sensor_to_cone_responses(
    model: CamModel,
    cat_transform: CatTransform,
    adapted_rgb: [f64; 3],
) -> LuxResult<[f64; 3]> {
    match model {
        CamModel::Ciecam02 => {
            let hpe = Observer::Cie1931_2.xyz_to_lms_matrix()?;
            let cone_matrix = multiply_matrix3(hpe, invert_matrix3(cat_transform.matrix()));
            Ok(multiply_matrix3_vector3(cone_matrix, adapted_rgb))
        }
        CamModel::Cam16 => Ok(adapted_rgb),
    }
}

fn cone_responses_to_cat_sensor(
    model: CamModel,
    cat_transform: CatTransform,
    cone: [f64; 3],
) -> LuxResult<[f64; 3]> {
    match model {
        CamModel::Ciecam02 => {
            let hpe = Observer::Cie1931_2.xyz_to_lms_matrix()?;
            let matrix = multiply_matrix3(cat_transform.matrix(), invert_matrix3(hpe));
            Ok(multiply_matrix3_vector3(matrix, cone))
        }
        CamModel::Cam16 => Ok(cone),
    }
}

fn apply_von_kries_to_channels(
    rgb: [f64; 3],
    rgb_white: [f64; 3],
    yw: f64,
    d: f64,
) -> LuxResult<[f64; 3]> {
    let mut out = [0.0; 3];
    for index in 0..3 {
        if rgb_white[index].abs() <= 1e-15 {
            return Err(LuxError::InvalidInput(
                "white point produces zero CAM sensor response",
            ));
        }
        let scale = (d * yw / rgb_white[index]) + (1.0 - d);
        out[index] = scale * rgb[index];
    }
    Ok(out)
}

fn invert_von_kries_from_channels(
    adapted_rgb: [f64; 3],
    rgb_white: [f64; 3],
    yw: f64,
    d: f64,
) -> LuxResult<[f64; 3]> {
    let mut out = [0.0; 3];
    for index in 0..3 {
        if rgb_white[index].abs() <= 1e-15 {
            return Err(LuxError::InvalidInput(
                "white point produces zero CAM sensor response",
            ));
        }
        let scale = (d * yw / rgb_white[index]) + (1.0 - d);
        out[index] = adapted_rgb[index] / scale;
    }
    Ok(out)
}

fn normalize_white_point(white_point: [f64; 3], yw: f64) -> LuxResult<[f64; 3]> {
    if white_point[1].abs() <= 1e-15 {
        return Err(LuxError::InvalidInput(
            "white point Y must be non-zero for CAM normalization",
        ));
    }
    Ok([
        yw * white_point[0] / white_point[1],
        yw,
        yw * white_point[2] / white_point[1],
    ])
}

fn scale_xyz_to_luminance_factor(
    xyz: [f64; 3],
    source_y: f64,
    target_y: f64,
) -> LuxResult<[f64; 3]> {
    if source_y.abs() <= 1e-15 {
        return Err(LuxError::InvalidInput(
            "source Y must be non-zero for CAM normalization",
        ));
    }
    Ok([
        target_y * xyz[0] / source_y,
        target_y * xyz[1] / source_y,
        target_y * xyz[2] / source_y,
    ])
}

fn apply_von_kries_to_white(rgb_white: [f64; 3], yw: f64, d: f64) -> LuxResult<[f64; 3]> {
    let mut out = [0.0; 3];
    for index in 0..3 {
        if rgb_white[index].abs() <= 1e-15 {
            return Err(LuxError::InvalidInput(
                "white point produces zero CAM sensor response",
            ));
        }
        out[index] = ((d * yw / rgb_white[index]) + (1.0 - d)) * rgb_white[index];
    }
    Ok(out)
}

fn multiply_matrix3_vector3(matrix: Matrix3, vector: [f64; 3]) -> [f64; 3] {
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ]
}

fn multiply_matrix3(left: Matrix3, right: Matrix3) -> Matrix3 {
    let mut out = [[0.0; 3]; 3];
    for row in 0..3 {
        for col in 0..3 {
            out[row][col] = left[row][0] * right[0][col]
                + left[row][1] * right[1][col]
                + left[row][2] * right[2][col];
        }
    }
    out
}

fn invert_matrix3(matrix: Matrix3) -> Matrix3 {
    let a = matrix[0][0];
    let b = matrix[0][1];
    let c = matrix[0][2];
    let d = matrix[1][0];
    let e = matrix[1][1];
    let f = matrix[1][2];
    let g = matrix[2][0];
    let h = matrix[2][1];
    let i = matrix[2][2];

    let cofactor00 = e * i - f * h;
    let cofactor01 = -(d * i - f * g);
    let cofactor02 = d * h - e * g;
    let cofactor10 = -(b * i - c * h);
    let cofactor11 = a * i - c * g;
    let cofactor12 = -(a * h - b * g);
    let cofactor20 = b * f - c * e;
    let cofactor21 = -(a * f - c * d);
    let cofactor22 = a * e - b * d;

    let determinant = a * cofactor00 + b * cofactor01 + c * cofactor02;
    let inv_det = 1.0 / determinant;

    [
        [
            cofactor00 * inv_det,
            cofactor10 * inv_det,
            cofactor20 * inv_det,
        ],
        [
            cofactor01 * inv_det,
            cofactor11 * inv_det,
            cofactor21 * inv_det,
        ],
        [
            cofactor02 * inv_det,
            cofactor12 * inv_det,
            cofactor22 * inv_det,
        ],
    ]
}

fn clamp(value: f64, min: f64, max: f64) -> f64 {
    value.max(min).min(max)
}

fn degrees_to_radians(value: f64) -> f64 {
    value * std::f64::consts::PI / 180.0
}

fn hue_angle_degrees(a: f64, b: f64) -> f64 {
    let angle = b.atan2(a).to_degrees();
    if angle < 0.0 {
        angle + 360.0
    } else {
        angle
    }
}

fn validate_xyz_triplet(xyz: [f64; 3], label: &'static str) -> LuxResult<()> {
    if xyz.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(LuxError::InvalidInput(label))
    }
}

fn validate_degree(value: f64, label: &'static str) -> LuxResult<()> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        Err(LuxError::InvalidInput(label))
    } else {
        Ok(())
    }
}

fn validate_positive_finite(value: f64, label: &'static str) -> LuxResult<()> {
    if !value.is_finite() || value <= 0.0 {
        Err(LuxError::InvalidInput(label))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{
        cam16_forward, cam16_ucs_forward, cam16_viewing_conditions, cam_inverse, cam_naka_rushton,
        ciecam02_forward, ciecam02_ucs_forward, ciecam02_ucs_inverse, ciecam02_viewing_conditions,
        CamModel, CamNakaRushtonParameters, CamSurround, CamUcsAppearance, CamUcsType,
        CamViewingConditions,
    };
    use crate::color::CatTransform;

    #[test]
    fn computes_cam16_viewing_conditions() {
        let vc = cam16_viewing_conditions(
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        assert_eq!(vc.model, CamModel::Cam16);
        assert_eq!(vc.cat_transform, CatTransform::Cat16);
        assert!((vc.fl - 0.793_700_527_546_167_3).abs() < 1e-12);
        assert!((vc.n - 0.2).abs() < 1e-12);
        assert!((vc.nbb - 1.000_304_004_559_380_7).abs() < 1e-12);
        assert!((vc.z - 1.927_213_595_499_957_9).abs() < 1e-12);
    }

    #[test]
    fn computes_ciecam02_viewing_conditions() {
        let vc = ciecam02_viewing_conditions(
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        assert_eq!(vc.model, CamModel::Ciecam02);
        assert_eq!(vc.cat_transform, CatTransform::Cat02);
        assert!((vc.fl - 0.793_700_527_546_167_3).abs() < 1e-12);
        assert!((vc.achromatic_response_white - 39.501_102_128_317_285).abs() < 1e-9);
        assert!((vc.brightness_white - 238.026_509_521_529_2).abs() < 1e-9);
    }

    #[test]
    fn computes_cam16_white_responses() {
        let vc = CamViewingConditions::new(
            CamModel::Cam16,
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        assert!((vc.achromatic_response_white - 39.500_996_860_311_794).abs() < 1e-9);
        assert!((vc.brightness_white - 238.025_933_522_883_52).abs() < 1e-9);
    }

    #[test]
    fn applies_naka_rushton_forward() {
        let parameters = CamNakaRushtonParameters {
            n: 0.42,
            sig: 27.13_f64.powf(1.0 / 0.42),
            scaling: 400.0,
            noise: 0.1,
        };
        let result = cam_naka_rushton(1.0, parameters, true);
        assert!((result - 14.319_694_276_573_056).abs() < 1e-12);
    }

    #[test]
    fn computes_cam16_forward_correlates() {
        let vc = cam16_viewing_conditions(
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        let appearance = cam16_forward([19.01, 20.0, 21.78], vc).unwrap();
        assert!((appearance.lightness - 41.579_587_651_004_54).abs() < 1e-11);
        assert!((appearance.brightness - 153.484_444_297_693).abs() < 1e-11);
        assert!((appearance.chroma - 0.017_728_932_726_109_602).abs() < 1e-12);
        assert!((appearance.colorfulness - 0.016_733_884_199_670_097).abs() < 1e-12);
        assert!((appearance.saturation - 1.044_157_943_071_202_5).abs() < 1e-11);
        assert!((appearance.hue_angle - 296.915_507_076_488_44).abs() < 1e-11);
        assert!((appearance.a_m - 0.007_575_028_723_342_962).abs() < 1e-12);
        assert!((appearance.b_m + 0.014_921_186_958_432_55).abs() < 1e-12);
    }

    #[test]
    fn computes_ciecam02_forward_correlates() {
        let vc = ciecam02_viewing_conditions(
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        let appearance = ciecam02_forward([19.01, 20.0, 21.78], vc).unwrap();
        assert!((appearance.lightness - 41.579_506_441_924_01).abs() < 1e-11);
        assert!((appearance.brightness - 153.484_665_828_851_66).abs() < 1e-11);
        assert!((appearance.chroma - 0.018_926_054_653_553_997).abs() < 1e-11);
        assert!((appearance.colorfulness - 0.017_863_816_836_688_81).abs() < 1e-11);
        assert!((appearance.saturation - 1.078_834_017_553_237).abs() < 1e-11);
        assert!((appearance.hue_angle - 302.725_637_816_980_2).abs() < 1e-10);
        assert!((appearance.a_m - 0.009_657_479_716_078_998).abs() < 1e-11);
        assert!((appearance.b_m + 0.015_028_274_601_839_338).abs() < 1e-11);
    }

    #[test]
    fn computes_cam16_ucs_forward_correlates() {
        let vc = cam16_viewing_conditions(
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        let appearance = cam16_ucs_forward([19.01, 20.0, 21.78], vc, CamUcsType::Ucs).unwrap();
        assert!((appearance.j_prime - 54.749_939_614_956_66).abs() < 1e-10);
        assert!((appearance.a_prime - 0.007_573_584_030_746_738).abs() < 1e-10);
        assert!((appearance.b_prime + 0.014_918_341_222_909_555).abs() < 1e-10);
    }

    #[test]
    fn computes_ciecam02_ucs_forward_correlates() {
        let vc = ciecam02_viewing_conditions(
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        let appearance = ciecam02_ucs_forward([19.01, 20.0, 21.78], vc, CamUcsType::Ucs).unwrap();
        assert!((appearance.j_prime - 54.749_856_789_698_91).abs() < 1e-10);
        assert!((appearance.a_prime - 0.009_655_513_528_223_68).abs() < 1e-10);
        assert!((appearance.b_prime + 0.015_025_214_961_863_152).abs() < 1e-10);
    }

    #[test]
    fn inverts_cam16_jabm_to_xyz() {
        let vc = cam16_viewing_conditions(
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        let xyz = cam_inverse(
            super::CamAppearance {
                lightness: 41.579_587_651_004_54,
                brightness: 0.0,
                chroma: 0.0,
                colorfulness: 0.016_733_884_199_670_097,
                saturation: 0.0,
                hue_angle: 0.0,
                a_m: 0.007_575_028_723_342_962,
                b_m: -0.014_921_186_958_432_55,
                a_c: 0.0,
                b_c: 0.0,
            },
            vc,
        )
        .unwrap();
        assert!((xyz[0] - 19.01).abs() < 1e-10);
        assert!((xyz[1] - 20.0).abs() < 1e-10);
        assert!((xyz[2] - 21.78).abs() < 1e-10);
    }

    #[test]
    fn inverts_ciecam02_ucs_to_xyz() {
        let vc = ciecam02_viewing_conditions(
            [95.047, 100.0, 108.883],
            None,
            100.0,
            20.0,
            CamSurround::Average,
            Some(1.0),
            None,
        )
        .unwrap();
        let xyz = ciecam02_ucs_inverse(
            CamUcsAppearance {
                j_prime: 54.749_856_789_698_91,
                a_prime: 0.009_655_513_528_223_68,
                b_prime: -0.015_025_214_961_863_152,
            },
            vc,
            CamUcsType::Ucs,
        )
        .unwrap();
        assert!((xyz[0] - 19.01).abs() < 1e-10);
        assert!((xyz[1] - 20.0).abs() < 1e-10);
        assert!((xyz[2] - 21.78).abs() < 1e-10);
    }
}
