use crate::cam::{
    cam16_forward, cam16_ucs_forward, cam16_ucs_inverse, cam_forward, cam_inverse, cam_ucs_forward,
    ciecam02_forward, ciecam02_ucs_forward, ciecam02_ucs_inverse, CamAppearance, CamUcsAppearance,
    CamUcsType, CamViewingConditions as ModelCamViewingConditions,
};
use crate::error::{LuxError, LuxResult};
use crate::spectrum::{SpectralMatrix, Spectrum};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Observer {
    Cie1931_2,
    Cie1964_10,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TristimulusObserver {
    pub wavelengths: Vec<f64>,
    pub x_bar: Vec<f64>,
    pub y_bar: Vec<f64>,
    pub z_bar: Vec<f64>,
    pub k: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MesopicLuminousEfficiency {
    pub curves: SpectralMatrix,
    pub k_mesopic: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tristimulus {
    values: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
pub struct TristimulusSet {
    values: Vec<[f64; 3]>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeltaEFormula {
    Cie76,
    Ciede2000,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CatTransform {
    Bradford,
    Cat02,
    Cat16,
    Sharp,
    Bianco,
    Cmc,
    Kries,
    Judd1945,
    Judd1945Cie016,
    Judd1935,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CatSurround {
    Average,
    Dim,
    Dark,
    Display,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CatMode {
    OneStep,
    SourceToBaseline,
    BaselineToTarget,
    TwoStep,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CatViewingConditions {
    pub surround: CatSurround,
    pub adapting_luminance: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CatContext {
    pub source_white: [f64; 3],
    pub target_white: [f64; 3],
    pub baseline_white: Option<[f64; 3]>,
    pub transform: CatTransform,
    pub mode: CatMode,
    pub source_conditions: CatViewingConditions,
    pub target_conditions: CatViewingConditions,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CatAdapter {
    matrix: Matrix3,
}

pub type Matrix3 = [[f64; 3]; 3];

const EPSILON: f64 = 1e-15;
const LAB_LINEAR_THRESHOLD: f64 = (24.0 / 116.0) * (24.0 / 116.0) * (24.0 / 116.0);
const LAB_LINEAR_SCALE: f64 = 841.0 / 108.0;
const LAB_INVERSE_LINEAR_SCALE: f64 = 108.0 / 841.0;
const LUV_LINEAR_THRESHOLD: f64 = (6.0 / 29.0) * (6.0 / 29.0) * (6.0 / 29.0);
const LUV_LINEAR_SCALE: f64 = (29.0 / 3.0) * (29.0 / 3.0) * (29.0 / 3.0);
const SRGB_XYZ_TO_RGB: Matrix3 = [
    [3.2404542, -1.5371385, -0.4985314],
    [-0.9692660, 1.8760108, 0.0415560],
    [0.0556434, -0.2040259, 1.0572252],
];
const SRGB_RGB_TO_XYZ: Matrix3 = [
    [0.4124564, 0.3575761, 0.1804375],
    [0.2126729, 0.7151522, 0.0721750],
    [0.0193339, 0.1191920, 0.9503041],
];

impl Observer {
    pub fn standard(self) -> LuxResult<TristimulusObserver> {
        match self {
            Self::Cie1931_2 => TristimulusObserver::from_csv(
                include_str!("../data/cmfs/ciexyz_1931_2.dat"),
                683.002,
            ),
            Self::Cie1964_10 => TristimulusObserver::from_csv(
                include_str!("../data/cmfs/ciexyz_1964_10.dat"),
                683.599,
            ),
        }
    }

    pub fn xyzbar(self) -> LuxResult<SpectralMatrix> {
        self.standard()?.xyz_spectra()
    }

    pub fn xyzbar_linear(self, target_wavelengths: &[f64]) -> LuxResult<SpectralMatrix> {
        self.xyzbar()?.cie_interp_linear(target_wavelengths, false)
    }

    pub fn vlbar(self) -> LuxResult<(Spectrum, f64)> {
        let observer = self.standard()?;
        Ok((observer.vl_spectrum()?, observer.k))
    }

    pub fn vlbar_linear(self, target_wavelengths: &[f64]) -> LuxResult<(Spectrum, f64)> {
        let (vl, k) = self.vlbar()?;
        Ok((vl.cie_interp_linear(target_wavelengths, false)?, k))
    }

    pub fn xyz_to_lms_matrix(self) -> LuxResult<Matrix3> {
        match self {
            Self::Cie1931_2 => Ok([
                [0.38971, 0.68898, -0.07868],
                [-0.22981, 1.1834, 0.04641],
                [0.0, 0.0, 1.0],
            ]),
            Self::Cie1964_10 => Ok([
                [
                    0.217_010_449_691_388_16,
                    0.835_733_670_117_584_4,
                    -0.043_510_597_212_556_935,
                ],
                [
                    -0.429_979_507_573_619_8,
                    1.203_889_456_462_98,
                    0.086_210_895_329_211_28,
                ],
                [0.0, 0.0, 0.465_792_338_736_113],
            ]),
        }
    }
}

impl Tristimulus {
    pub fn new(values: [f64; 3]) -> Self {
        Self { values }
    }

    pub fn values(self) -> [f64; 3] {
        self.values
    }

    pub fn xyz_to_yxy(self) -> Self {
        Self::new(xyz_to_yxy(self.values))
    }

    pub fn yxy_to_xyz(self) -> Self {
        Self::new(yxy_to_xyz(self.values))
    }

    pub fn xyz_to_yuv(self) -> Self {
        Self::new(xyz_to_yuv(self.values))
    }

    pub fn yuv_to_xyz(self) -> Self {
        Self::new(yuv_to_xyz(self.values))
    }

    pub fn xyz_to_lab(self, white_point: [f64; 3]) -> Self {
        Self::new(xyz_to_lab(self.values, white_point))
    }

    pub fn lab_to_xyz(self, white_point: [f64; 3]) -> Self {
        Self::new(lab_to_xyz(self.values, white_point))
    }

    pub fn xyz_to_luv(self, white_point: [f64; 3]) -> Self {
        Self::new(xyz_to_luv(self.values, white_point))
    }

    pub fn luv_to_xyz(self, white_point: [f64; 3]) -> Self {
        Self::new(luv_to_xyz(self.values, white_point))
    }

    pub fn xyz_to_lms(self, observer: Observer) -> LuxResult<Self> {
        Ok(Self::new(xyz_to_lms(self.values, observer)?))
    }

    pub fn lms_to_xyz(self, observer: Observer) -> LuxResult<Self> {
        Ok(Self::new(lms_to_xyz(self.values, observer)?))
    }

    pub fn xyz_to_srgb(self, gamma: f64, offset: f64, use_linear_part: bool) -> Self {
        Self::new(xyz_to_srgb(self.values, gamma, offset, use_linear_part))
    }

    pub fn srgb_to_xyz(self, gamma: f64, offset: f64, use_linear_part: bool) -> Self {
        Self::new(srgb_to_xyz(self.values, gamma, offset, use_linear_part))
    }

    pub fn cat_apply(
        self,
        source_white: [f64; 3],
        target_white: [f64; 3],
        transform: CatTransform,
        degree_of_adaptation: f64,
    ) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_degree(
            source_white,
            target_white,
            transform,
            degree_of_adaptation,
        )?)
    }

    pub fn cat_apply_mode(
        self,
        source_white: [f64; 3],
        target_white: [f64; 3],
        baseline_white: Option<[f64; 3]>,
        transform: CatTransform,
        mode: CatMode,
        degrees_of_adaptation: [f64; 2],
    ) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_mode(
            source_white,
            target_white,
            baseline_white,
            transform,
            mode,
            degrees_of_adaptation,
        )?)
    }

    pub fn cat_apply_with_conditions(
        self,
        source_white: [f64; 3],
        target_white: [f64; 3],
        transform: CatTransform,
        surround: CatSurround,
        adapting_luminance: f64,
    ) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_degree(
            source_white,
            target_white,
            transform,
            cat_degree_of_adaptation(surround, adapting_luminance)?,
        )?)
    }

    pub fn cat_apply_mode_with_conditions(
        self,
        source_white: [f64; 3],
        target_white: [f64; 3],
        baseline_white: Option<[f64; 3]>,
        transform: CatTransform,
        mode: CatMode,
        source_conditions: CatViewingConditions,
        target_conditions: CatViewingConditions,
    ) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_conditions(
            source_white,
            target_white,
            baseline_white,
            transform,
            mode,
            source_conditions,
            target_conditions,
        )?)
    }

    pub fn cat_apply_context(self, context: CatContext) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_context(context)?)
    }

    pub fn cat_apply_adapter(self, adapter: CatAdapter) -> LuxResult<Self> {
        Ok(Self::new(adapter.apply(self.values)?))
    }

    pub fn delta_e(self, other: Self, white_point: [f64; 3], formula: DeltaEFormula) -> f64 {
        delta_e(self.values, other.values, white_point, formula)
    }

    pub fn cam_forward(self, conditions: ModelCamViewingConditions) -> LuxResult<CamAppearance> {
        cam_forward(self.values, conditions)
    }

    pub fn cam16_forward(self, conditions: ModelCamViewingConditions) -> LuxResult<CamAppearance> {
        cam16_forward(self.values, conditions)
    }

    pub fn ciecam02_forward(
        self,
        conditions: ModelCamViewingConditions,
    ) -> LuxResult<CamAppearance> {
        ciecam02_forward(self.values, conditions)
    }

    pub fn cam_ucs_forward(
        self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<CamUcsAppearance> {
        cam_ucs_forward(self.values, conditions, ucs_type)
    }

    pub fn cam16_ucs_forward(
        self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<CamUcsAppearance> {
        cam16_ucs_forward(self.values, conditions, ucs_type)
    }

    pub fn ciecam02_ucs_forward(
        self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<CamUcsAppearance> {
        ciecam02_ucs_forward(self.values, conditions, ucs_type)
    }

    pub fn cam_inverse(self, conditions: ModelCamViewingConditions) -> LuxResult<Self> {
        Ok(Self::new(cam_inverse(
            CamAppearance {
                lightness: self.values[0],
                brightness: 0.0,
                chroma: 0.0,
                colorfulness: (self.values[1] * self.values[1] + self.values[2] * self.values[2])
                    .sqrt(),
                saturation: 0.0,
                hue_angle: 0.0,
                a_m: self.values[1],
                b_m: self.values[2],
                a_c: 0.0,
                b_c: 0.0,
            },
            conditions,
        )?))
    }

    pub fn cam16_ucs_inverse(
        self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<Self> {
        Ok(Self::new(cam16_ucs_inverse(
            CamUcsAppearance {
                j_prime: self.values[0],
                a_prime: self.values[1],
                b_prime: self.values[2],
            },
            conditions,
            ucs_type,
        )?))
    }

    pub fn ciecam02_ucs_inverse(
        self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<Self> {
        Ok(Self::new(ciecam02_ucs_inverse(
            CamUcsAppearance {
                j_prime: self.values[0],
                a_prime: self.values[1],
                b_prime: self.values[2],
            },
            conditions,
            ucs_type,
        )?))
    }
}

impl From<[f64; 3]> for Tristimulus {
    fn from(values: [f64; 3]) -> Self {
        Self::new(values)
    }
}

impl From<Tristimulus> for [f64; 3] {
    fn from(value: Tristimulus) -> Self {
        value.values
    }
}

impl TristimulusSet {
    pub fn new(values: Vec<[f64; 3]>) -> Self {
        Self { values }
    }

    pub fn from_single(value: Tristimulus) -> Self {
        Self::new(vec![value.values()])
    }

    pub fn values(&self) -> &[[f64; 3]] {
        &self.values
    }

    pub fn iter(&self) -> impl Iterator<Item = Tristimulus> + '_ {
        self.values.iter().copied().map(Tristimulus::new)
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn into_vec(self) -> Vec<[f64; 3]> {
        self.values
    }

    pub fn xyz_to_yxy(&self) -> Self {
        Self::new(self.values.iter().copied().map(xyz_to_yxy).collect())
    }

    pub fn yxy_to_xyz(&self) -> Self {
        Self::new(self.values.iter().copied().map(yxy_to_xyz).collect())
    }

    pub fn xyz_to_yuv(&self) -> Self {
        Self::new(self.values.iter().copied().map(xyz_to_yuv).collect())
    }

    pub fn yuv_to_xyz(&self) -> Self {
        Self::new(self.values.iter().copied().map(yuv_to_xyz).collect())
    }

    pub fn xyz_to_lab(&self, white_point: [f64; 3]) -> Self {
        Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| xyz_to_lab(value, white_point))
                .collect(),
        )
    }

    pub fn lab_to_xyz(&self, white_point: [f64; 3]) -> Self {
        Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| lab_to_xyz(value, white_point))
                .collect(),
        )
    }

    pub fn xyz_to_luv(&self, white_point: [f64; 3]) -> Self {
        Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| xyz_to_luv(value, white_point))
                .collect(),
        )
    }

    pub fn luv_to_xyz(&self, white_point: [f64; 3]) -> Self {
        Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| luv_to_xyz(value, white_point))
                .collect(),
        )
    }

    pub fn xyz_to_lms(&self, observer: Observer) -> LuxResult<Self> {
        Ok(Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| xyz_to_lms(value, observer))
                .collect::<LuxResult<Vec<_>>>()?,
        ))
    }

    pub fn lms_to_xyz(&self, observer: Observer) -> LuxResult<Self> {
        Ok(Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| lms_to_xyz(value, observer))
                .collect::<LuxResult<Vec<_>>>()?,
        ))
    }

    pub fn xyz_to_srgb(&self, gamma: f64, offset: f64, use_linear_part: bool) -> Self {
        Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| xyz_to_srgb(value, gamma, offset, use_linear_part))
                .collect(),
        )
    }

    pub fn srgb_to_xyz(&self, gamma: f64, offset: f64, use_linear_part: bool) -> Self {
        Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| srgb_to_xyz(value, gamma, offset, use_linear_part))
                .collect(),
        )
    }

    pub fn cat_apply(
        &self,
        source_white: [f64; 3],
        target_white: [f64; 3],
        transform: CatTransform,
        degree_of_adaptation: f64,
    ) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_degree(
            source_white,
            target_white,
            transform,
            degree_of_adaptation,
        )?)
    }

    pub fn cat_apply_mode(
        &self,
        source_white: [f64; 3],
        target_white: [f64; 3],
        baseline_white: Option<[f64; 3]>,
        transform: CatTransform,
        mode: CatMode,
        degrees_of_adaptation: [f64; 2],
    ) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_mode(
            source_white,
            target_white,
            baseline_white,
            transform,
            mode,
            degrees_of_adaptation,
        )?)
    }

    pub fn cat_apply_with_conditions(
        &self,
        source_white: [f64; 3],
        target_white: [f64; 3],
        transform: CatTransform,
        surround: CatSurround,
        adapting_luminance: f64,
    ) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_degree(
            source_white,
            target_white,
            transform,
            cat_degree_of_adaptation(surround, adapting_luminance)?,
        )?)
    }

    pub fn cat_apply_mode_with_conditions(
        &self,
        source_white: [f64; 3],
        target_white: [f64; 3],
        baseline_white: Option<[f64; 3]>,
        transform: CatTransform,
        mode: CatMode,
        source_conditions: CatViewingConditions,
        target_conditions: CatViewingConditions,
    ) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_conditions(
            source_white,
            target_white,
            baseline_white,
            transform,
            mode,
            source_conditions,
            target_conditions,
        )?)
    }

    pub fn cat_apply_context(&self, context: CatContext) -> LuxResult<Self> {
        self.cat_apply_adapter(CatAdapter::from_context(context)?)
    }

    pub fn cat_apply_adapter(&self, adapter: CatAdapter) -> LuxResult<Self> {
        Ok(Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| adapter.apply(value))
                .collect::<LuxResult<Vec<_>>>()?,
        ))
    }

    pub fn delta_e(
        &self,
        other: &Self,
        white_point: [f64; 3],
        formula: DeltaEFormula,
    ) -> LuxResult<Vec<f64>> {
        if self.len() != other.len() {
            return Err(LuxError::MismatchedLengths {
                wavelengths: self.len(),
                values: other.len(),
            });
        }

        Ok(self
            .values
            .iter()
            .copied()
            .zip(other.values.iter().copied())
            .map(|(left, right)| delta_e(left, right, white_point, formula))
            .collect())
    }

    pub fn cam_forward(
        &self,
        conditions: ModelCamViewingConditions,
    ) -> LuxResult<Vec<CamAppearance>> {
        self.values
            .iter()
            .copied()
            .map(|value| cam_forward(value, conditions))
            .collect()
    }

    pub fn cam16_forward(
        &self,
        conditions: ModelCamViewingConditions,
    ) -> LuxResult<Vec<CamAppearance>> {
        self.values
            .iter()
            .copied()
            .map(|value| cam16_forward(value, conditions))
            .collect()
    }

    pub fn ciecam02_forward(
        &self,
        conditions: ModelCamViewingConditions,
    ) -> LuxResult<Vec<CamAppearance>> {
        self.values
            .iter()
            .copied()
            .map(|value| ciecam02_forward(value, conditions))
            .collect()
    }

    pub fn cam_ucs_forward(
        &self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<Vec<CamUcsAppearance>> {
        self.values
            .iter()
            .copied()
            .map(|value| cam_ucs_forward(value, conditions, ucs_type))
            .collect()
    }

    pub fn cam16_ucs_forward(
        &self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<Vec<CamUcsAppearance>> {
        self.values
            .iter()
            .copied()
            .map(|value| cam16_ucs_forward(value, conditions, ucs_type))
            .collect()
    }

    pub fn ciecam02_ucs_forward(
        &self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<Vec<CamUcsAppearance>> {
        self.values
            .iter()
            .copied()
            .map(|value| ciecam02_ucs_forward(value, conditions, ucs_type))
            .collect()
    }

    pub fn cam_inverse(&self, conditions: ModelCamViewingConditions) -> LuxResult<Self> {
        Ok(Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| {
                    cam_inverse(
                        CamAppearance {
                            lightness: value[0],
                            brightness: 0.0,
                            chroma: 0.0,
                            colorfulness: (value[1] * value[1] + value[2] * value[2]).sqrt(),
                            saturation: 0.0,
                            hue_angle: 0.0,
                            a_m: value[1],
                            b_m: value[2],
                            a_c: 0.0,
                            b_c: 0.0,
                        },
                        conditions,
                    )
                })
                .collect::<LuxResult<Vec<_>>>()?,
        ))
    }

    pub fn cam16_ucs_inverse(
        &self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<Self> {
        Ok(Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| {
                    cam16_ucs_inverse(
                        CamUcsAppearance {
                            j_prime: value[0],
                            a_prime: value[1],
                            b_prime: value[2],
                        },
                        conditions,
                        ucs_type,
                    )
                })
                .collect::<LuxResult<Vec<_>>>()?,
        ))
    }

    pub fn ciecam02_ucs_inverse(
        &self,
        conditions: ModelCamViewingConditions,
        ucs_type: CamUcsType,
    ) -> LuxResult<Self> {
        Ok(Self::new(
            self.values
                .iter()
                .copied()
                .map(|value| {
                    ciecam02_ucs_inverse(
                        CamUcsAppearance {
                            j_prime: value[0],
                            a_prime: value[1],
                            b_prime: value[2],
                        },
                        conditions,
                        ucs_type,
                    )
                })
                .collect::<LuxResult<Vec<_>>>()?,
        ))
    }
}

impl From<Vec<[f64; 3]>> for TristimulusSet {
    fn from(values: Vec<[f64; 3]>) -> Self {
        Self::new(values)
    }
}

impl From<Tristimulus> for TristimulusSet {
    fn from(value: Tristimulus) -> Self {
        Self::from_single(value)
    }
}

impl FromIterator<Tristimulus> for TristimulusSet {
    fn from_iter<T: IntoIterator<Item = Tristimulus>>(iter: T) -> Self {
        Self::new(iter.into_iter().map(Tristimulus::values).collect())
    }
}

impl CatTransform {
    pub fn matrix(self) -> Matrix3 {
        match self {
            Self::Bradford => [
                [0.8951, 0.2664, -0.1614],
                [-0.7502, 1.7135, 0.0367],
                [0.0389, -0.0685, 1.0296],
            ],
            Self::Cat02 => [
                [0.7328, 0.4296, -0.1624],
                [-0.7036, 1.6975, 0.0061],
                [0.0030, 0.0136, 0.9834],
            ],
            Self::Cat16 => [
                [0.401288, 0.650173, -0.051461],
                [-0.250268, 1.204414, 0.045854],
                [-0.002079, 0.048952, 0.953127],
            ],
            Self::Sharp => [
                [1.2694, -0.0988, -0.1706],
                [-0.8364, 1.8006, 0.0357],
                [0.0297, -0.0315, 1.0018],
            ],
            Self::Bianco => [
                [0.8752, 0.2787, -0.1539],
                [-0.8904, 1.8709, 0.0195],
                [-0.0061, 0.0162, 0.9899],
            ],
            Self::Cmc => [
                [0.7982, 0.3389, -0.1371],
                [-0.5918, 1.5512, 0.0406],
                [0.0008, 0.0239, 0.9753],
            ],
            Self::Kries => [
                [0.40024, 0.70760, -0.08081],
                [-0.22630, 1.16532, 0.04570],
                [0.0, 0.0, 0.91822],
            ],
            Self::Judd1945 => [[0.0, 1.0, 0.0], [-0.460, 1.359, 0.101], [0.0, 0.0, 1.0]],
            Self::Judd1945Cie016 => [[0.0, 1.0, 0.0], [-0.460, 1.360, 0.102], [0.0, 0.0, 1.0]],
            Self::Judd1935 => [
                [3.1956, 2.4478, -0.1434],
                [-2.5455, 7.0942, 0.9963],
                [0.0, 0.0, 1.0],
            ],
        }
    }
}

impl CatSurround {
    pub fn factor(self) -> f64 {
        match self {
            Self::Average => 1.0,
            Self::Dim => 0.9,
            Self::Dark => 0.8,
            Self::Display => 0.0,
        }
    }
}

impl CatMode {
    pub fn default_baseline_white(self) -> [f64; 3] {
        match self {
            Self::SourceToBaseline | Self::BaselineToTarget | Self::TwoStep => {
                [100.0, 100.0, 100.0]
            }
            Self::OneStep => [100.0, 100.0, 100.0],
        }
    }
}

impl CatViewingConditions {
    pub fn new(surround: CatSurround, adapting_luminance: f64) -> LuxResult<Self> {
        validate_adapting_luminance(adapting_luminance)?;
        Ok(Self {
            surround,
            adapting_luminance,
        })
    }

    pub fn degree_of_adaptation(self) -> LuxResult<f64> {
        cat_degree_of_adaptation(self.surround, self.adapting_luminance)
    }
}

impl CatContext {
    pub fn new(
        source_white: [f64; 3],
        target_white: [f64; 3],
        baseline_white: Option<[f64; 3]>,
        transform: CatTransform,
        mode: CatMode,
        source_conditions: CatViewingConditions,
        target_conditions: CatViewingConditions,
    ) -> Self {
        Self {
            source_white,
            target_white,
            baseline_white,
            transform,
            mode,
            source_conditions,
            target_conditions,
        }
    }

    pub fn baseline_white_or_default(self) -> [f64; 3] {
        self.baseline_white
            .unwrap_or(self.mode.default_baseline_white())
    }

    pub fn degrees_of_adaptation(self) -> LuxResult<[f64; 2]> {
        cat_mode_degrees_from_conditions(self.mode, self.source_conditions, self.target_conditions)
    }
}

impl CatAdapter {
    pub fn new(matrix: Matrix3) -> Self {
        Self { matrix }
    }

    pub fn matrix(self) -> Matrix3 {
        self.matrix
    }

    pub fn apply(self, xyz: [f64; 3]) -> LuxResult<[f64; 3]> {
        validate_xyz_triplet(xyz, "xyz values must be finite")?;
        Ok(multiply_matrix3_vector3(self.matrix, xyz))
    }

    pub fn from_degree(
        source_white: [f64; 3],
        target_white: [f64; 3],
        transform: CatTransform,
        degree_of_adaptation: f64,
    ) -> LuxResult<Self> {
        validate_xyz_triplet(source_white, "source white values must be finite")?;
        validate_xyz_triplet(target_white, "target white values must be finite")?;
        validate_degree(
            degree_of_adaptation,
            "degree_of_adaptation must be finite and within 0..=1",
        )?;

        let sensor_matrix = transform.matrix();
        let inverse = invert_matrix3(sensor_matrix);
        let rgbw_source = multiply_matrix3_vector3(sensor_matrix, source_white);
        let rgbw_target = multiply_matrix3_vector3(sensor_matrix, target_white);
        let mut diagonal = [[0.0; 3]; 3];

        for index in 0..3 {
            if rgbw_source[index].abs() <= EPSILON {
                return Err(LuxError::InvalidInput(
                    "source white produces zero CAT sensor response",
                ));
            }
            let ratio = rgbw_target[index] / rgbw_source[index];
            diagonal[index][index] = degree_of_adaptation * ratio + (1.0 - degree_of_adaptation);
        }

        Ok(Self::new(multiply_matrix3(
            inverse,
            multiply_matrix3(diagonal, sensor_matrix),
        )))
    }

    pub fn from_mode(
        source_white: [f64; 3],
        target_white: [f64; 3],
        baseline_white: Option<[f64; 3]>,
        transform: CatTransform,
        mode: CatMode,
        degrees_of_adaptation: [f64; 2],
    ) -> LuxResult<Self> {
        validate_xyz_triplet(source_white, "source white values must be finite")?;
        validate_xyz_triplet(target_white, "target white values must be finite")?;
        validate_degree(
            degrees_of_adaptation[0],
            "degrees_of_adaptation[0] must be finite and within 0..=1",
        )?;
        validate_degree(
            degrees_of_adaptation[1],
            "degrees_of_adaptation[1] must be finite and within 0..=1",
        )?;

        let baseline_white = baseline_white.unwrap_or(mode.default_baseline_white());
        validate_xyz_triplet(baseline_white, "baseline white values must be finite")?;

        match mode {
            CatMode::OneStep => Self::from_degree(
                source_white,
                target_white,
                transform,
                degrees_of_adaptation[0],
            ),
            CatMode::SourceToBaseline => Self::from_degree(
                source_white,
                baseline_white,
                transform,
                degrees_of_adaptation[0],
            ),
            CatMode::BaselineToTarget => Self::from_degree(
                baseline_white,
                target_white,
                transform,
                degrees_of_adaptation[0],
            ),
            CatMode::TwoStep => {
                let sensor_matrix = transform.matrix();
                let inverse = invert_matrix3(sensor_matrix);
                let rgbw1 = multiply_matrix3_vector3(sensor_matrix, source_white);
                let rgbw2 = multiply_matrix3_vector3(sensor_matrix, target_white);
                let rgbw0 = multiply_matrix3_vector3(sensor_matrix, baseline_white);
                let mut diagonal = [[0.0; 3]; 3];

                for index in 0..3 {
                    if rgbw1[index].abs() <= EPSILON {
                        return Err(LuxError::InvalidInput(
                            "source white produces zero CAT sensor response",
                        ));
                    }
                    if rgbw2[index].abs() <= EPSILON {
                        return Err(LuxError::InvalidInput(
                            "target white produces zero CAT sensor response",
                        ));
                    }
                    let scale10 = degrees_of_adaptation[0] * (rgbw0[index] / rgbw1[index])
                        + (1.0 - degrees_of_adaptation[0]);
                    let scale20 = degrees_of_adaptation[1] * (rgbw0[index] / rgbw2[index])
                        + (1.0 - degrees_of_adaptation[1]);
                    diagonal[index][index] = scale10 / scale20;
                }

                Ok(Self::new(multiply_matrix3(
                    inverse,
                    multiply_matrix3(diagonal, sensor_matrix),
                )))
            }
        }
    }

    pub fn from_conditions(
        source_white: [f64; 3],
        target_white: [f64; 3],
        baseline_white: Option<[f64; 3]>,
        transform: CatTransform,
        mode: CatMode,
        source_conditions: CatViewingConditions,
        target_conditions: CatViewingConditions,
    ) -> LuxResult<Self> {
        let degrees = cat_mode_degrees_from_conditions(mode, source_conditions, target_conditions)?;
        Self::from_mode(
            source_white,
            target_white,
            baseline_white,
            transform,
            mode,
            degrees,
        )
    }

    pub fn from_context(context: CatContext) -> LuxResult<Self> {
        Self::from_conditions(
            context.source_white,
            context.target_white,
            context.baseline_white,
            context.transform,
            context.mode,
            context.source_conditions,
            context.target_conditions,
        )
    }
}

impl TristimulusObserver {
    pub fn from_csv(csv: &str, k: f64) -> LuxResult<Self> {
        let mut wavelengths = Vec::new();
        let mut x_bar = Vec::new();
        let mut y_bar = Vec::new();
        let mut z_bar = Vec::new();

        for line in csv.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            let mut parts = trimmed.split(',');
            let wl = parts
                .next()
                .ok_or(LuxError::ParseError("missing wavelength"))?
                .trim()
                .parse::<f64>()
                .map_err(|_| LuxError::ParseError("invalid wavelength"))?;
            let x = parts
                .next()
                .ok_or(LuxError::ParseError("missing x_bar"))?
                .trim()
                .parse::<f64>()
                .map_err(|_| LuxError::ParseError("invalid x_bar"))?;
            let y = parts
                .next()
                .ok_or(LuxError::ParseError("missing y_bar"))?
                .trim()
                .parse::<f64>()
                .map_err(|_| LuxError::ParseError("invalid y_bar"))?;
            let z = parts
                .next()
                .ok_or(LuxError::ParseError("missing z_bar"))?
                .trim()
                .parse::<f64>()
                .map_err(|_| LuxError::ParseError("invalid z_bar"))?;

            wavelengths.push(wl);
            x_bar.push(x);
            y_bar.push(y);
            z_bar.push(z);
        }

        if wavelengths.is_empty() {
            return Err(LuxError::EmptyInput);
        }

        Ok(Self {
            wavelengths,
            x_bar,
            y_bar,
            z_bar,
            k,
        })
    }

    pub fn vl_spectrum(&self) -> LuxResult<Spectrum> {
        Spectrum::new(self.wavelengths.clone(), self.y_bar.clone())
    }

    pub fn xyz_spectra(&self) -> LuxResult<SpectralMatrix> {
        SpectralMatrix::new(
            self.wavelengths.clone(),
            vec![self.x_bar.clone(), self.y_bar.clone(), self.z_bar.clone()],
        )
    }

    pub fn x_bar_spectrum(&self) -> LuxResult<Spectrum> {
        Spectrum::new(self.wavelengths.clone(), self.x_bar.clone())
    }

    pub fn z_bar_spectrum(&self) -> LuxResult<Spectrum> {
        Spectrum::new(self.wavelengths.clone(), self.z_bar.clone())
    }
}

fn nonzero(value: f64) -> f64 {
    if value == 0.0 {
        EPSILON
    } else {
        value
    }
}

fn lab_response_curve(value: f64, white: f64) -> f64 {
    let ratio = value / white;
    if ratio <= LAB_LINEAR_THRESHOLD {
        LAB_LINEAR_SCALE * ratio + 16.0 / 116.0
    } else {
        ratio.cbrt()
    }
}

fn lab_inverse_response_curve(response: f64, white: f64) -> f64 {
    if response <= 24.0 / 116.0 {
        white * ((response - 16.0 / 116.0) * LAB_INVERSE_LINEAR_SCALE)
    } else {
        white * response.powi(3)
    }
}

fn cie_lightness_from_ratio(y_ratio: f64) -> f64 {
    if y_ratio <= LUV_LINEAR_THRESHOLD {
        LUV_LINEAR_SCALE * y_ratio
    } else {
        116.0 * y_ratio.cbrt() - 16.0
    }
}

fn cie_y_ratio_from_lightness(lightness: f64) -> f64 {
    let y_ratio = ((lightness + 16.0) / 116.0).powi(3);
    if y_ratio < LUV_LINEAR_THRESHOLD {
        lightness / LUV_LINEAR_SCALE
    } else {
        y_ratio
    }
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

fn validate_adapting_luminance(adapting_luminance: f64) -> LuxResult<()> {
    if !adapting_luminance.is_finite() || adapting_luminance < 0.0 {
        Err(LuxError::InvalidInput(
            "adapting_luminance must be finite and non-negative",
        ))
    } else {
        Ok(())
    }
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

// CIE chromaticity transforms.

pub fn xyz_to_yxy(xyz: [f64; 3]) -> [f64; 3] {
    let sum = xyz[0] + xyz[1] + xyz[2];
    let denominator = nonzero(sum);
    [xyz[1], xyz[0] / denominator, xyz[1] / denominator]
}

pub fn yxy_to_xyz(yxy: [f64; 3]) -> [f64; 3] {
    let y = nonzero(yxy[2]);
    [
        yxy[0] * yxy[1] / y,
        yxy[0],
        yxy[0] * (1.0 - yxy[1] - yxy[2]) / y,
    ]
}

pub fn xyz_to_yuv(xyz: [f64; 3]) -> [f64; 3] {
    let denominator = xyz[0] + 15.0 * xyz[1] + 3.0 * xyz[2];
    let denominator = nonzero(denominator);
    [
        xyz[1],
        4.0 * xyz[0] / denominator,
        9.0 * xyz[1] / denominator,
    ]
}

pub fn yuv_to_xyz(yuv: [f64; 3]) -> [f64; 3] {
    let v = nonzero(yuv[2]);
    [
        yuv[0] * (9.0 * yuv[1]) / (4.0 * v),
        yuv[0],
        yuv[0] * (12.0 - 3.0 * yuv[1] - 20.0 * yuv[2]) / (4.0 * v),
    ]
}

// LMS transforms.

pub fn xyz_to_lms_with_matrix(xyz: [f64; 3], matrix: Matrix3) -> [f64; 3] {
    multiply_matrix3_vector3(matrix, xyz)
}

pub fn lms_to_xyz_with_matrix(lms: [f64; 3], matrix: Matrix3) -> [f64; 3] {
    multiply_matrix3_vector3(invert_matrix3(matrix), lms)
}

pub fn xyz_to_lms(xyz: [f64; 3], observer: Observer) -> LuxResult<[f64; 3]> {
    Ok(xyz_to_lms_with_matrix(xyz, observer.xyz_to_lms_matrix()?))
}

pub fn lms_to_xyz(lms: [f64; 3], observer: Observer) -> LuxResult<[f64; 3]> {
    Ok(lms_to_xyz_with_matrix(lms, observer.xyz_to_lms_matrix()?))
}

// Chromatic adaptation.

pub fn cat_apply(
    xyz: [f64; 3],
    source_white: [f64; 3],
    target_white: [f64; 3],
    transform: CatTransform,
    degree_of_adaptation: f64,
) -> LuxResult<[f64; 3]> {
    cat_compile(source_white, target_white, transform, degree_of_adaptation)?.apply(xyz)
}

pub fn cat_apply_mode(
    xyz: [f64; 3],
    source_white: [f64; 3],
    target_white: [f64; 3],
    baseline_white: Option<[f64; 3]>,
    transform: CatTransform,
    mode: CatMode,
    degrees_of_adaptation: [f64; 2],
) -> LuxResult<[f64; 3]> {
    cat_compile_mode(
        source_white,
        target_white,
        baseline_white,
        transform,
        mode,
        degrees_of_adaptation,
    )?
    .apply(xyz)
}

pub fn cat_compile(
    source_white: [f64; 3],
    target_white: [f64; 3],
    transform: CatTransform,
    degree_of_adaptation: f64,
) -> LuxResult<CatAdapter> {
    CatAdapter::from_degree(source_white, target_white, transform, degree_of_adaptation)
}

pub fn cat_compile_mode(
    source_white: [f64; 3],
    target_white: [f64; 3],
    baseline_white: Option<[f64; 3]>,
    transform: CatTransform,
    mode: CatMode,
    degrees_of_adaptation: [f64; 2],
) -> LuxResult<CatAdapter> {
    CatAdapter::from_mode(
        source_white,
        target_white,
        baseline_white,
        transform,
        mode,
        degrees_of_adaptation,
    )
}

pub fn cat_degree_of_adaptation(surround: CatSurround, adapting_luminance: f64) -> LuxResult<f64> {
    validate_adapting_luminance(adapting_luminance)?;

    let factor = surround.factor();
    let degree = factor * (1.0 - (1.0 / 3.6) * ((-adapting_luminance - 42.0) / 92.0).exp());
    Ok(clamp(degree, 0.0, 1.0))
}

pub fn cat_mode_degrees_from_conditions(
    mode: CatMode,
    source_conditions: CatViewingConditions,
    target_conditions: CatViewingConditions,
) -> LuxResult<[f64; 2]> {
    let source_degree = source_conditions.degree_of_adaptation()?;
    let target_degree = target_conditions.degree_of_adaptation()?;

    Ok(match mode {
        CatMode::OneStep | CatMode::SourceToBaseline => [source_degree, source_degree],
        CatMode::BaselineToTarget => [target_degree, target_degree],
        CatMode::TwoStep => [source_degree, target_degree],
    })
}

pub fn cat_apply_with_conditions(
    xyz: [f64; 3],
    source_white: [f64; 3],
    target_white: [f64; 3],
    transform: CatTransform,
    surround: CatSurround,
    adapting_luminance: f64,
) -> LuxResult<[f64; 3]> {
    cat_compile_with_conditions(
        source_white,
        target_white,
        transform,
        surround,
        adapting_luminance,
    )?
    .apply(xyz)
}

pub fn cat_apply_mode_with_conditions(
    xyz: [f64; 3],
    source_white: [f64; 3],
    target_white: [f64; 3],
    baseline_white: Option<[f64; 3]>,
    transform: CatTransform,
    mode: CatMode,
    source_conditions: CatViewingConditions,
    target_conditions: CatViewingConditions,
) -> LuxResult<[f64; 3]> {
    cat_compile_mode_with_conditions(
        source_white,
        target_white,
        baseline_white,
        transform,
        mode,
        source_conditions,
        target_conditions,
    )
    .and_then(|adapter| adapter.apply(xyz))
}

pub fn cat_apply_context(xyz: [f64; 3], context: CatContext) -> LuxResult<[f64; 3]> {
    cat_compile_context(context)?.apply(xyz)
}

pub fn cat_compile_with_conditions(
    source_white: [f64; 3],
    target_white: [f64; 3],
    transform: CatTransform,
    surround: CatSurround,
    adapting_luminance: f64,
) -> LuxResult<CatAdapter> {
    let degree = cat_degree_of_adaptation(surround, adapting_luminance)?;
    cat_compile(source_white, target_white, transform, degree)
}

pub fn cat_compile_mode_with_conditions(
    source_white: [f64; 3],
    target_white: [f64; 3],
    baseline_white: Option<[f64; 3]>,
    transform: CatTransform,
    mode: CatMode,
    source_conditions: CatViewingConditions,
    target_conditions: CatViewingConditions,
) -> LuxResult<CatAdapter> {
    CatAdapter::from_conditions(
        source_white,
        target_white,
        baseline_white,
        transform,
        mode,
        source_conditions,
        target_conditions,
    )
}

pub fn cat_compile_context(context: CatContext) -> LuxResult<CatAdapter> {
    CatAdapter::from_context(context)
}

// sRGB transforms.

pub fn xyz_to_srgb(xyz: [f64; 3], gamma: f64, offset: f64, use_linear_part: bool) -> [f64; 3] {
    let linear = multiply_matrix3_vector3(
        SRGB_XYZ_TO_RGB,
        [xyz[0] / 100.0, xyz[1] / 100.0, xyz[2] / 100.0],
    );

    let mut rgb = [0.0; 3];
    for (index, linear_value) in linear.iter().enumerate() {
        let srgb = clamp(*linear_value, 0.0, 1.0);
        let mut encoded = ((1.0 - offset) * srgb.powf(1.0 / gamma) + offset) * 255.0;
        if use_linear_part && srgb <= 0.0031308 {
            encoded = srgb * 12.92 * 255.0;
        }
        rgb[index] = clamp(encoded, 0.0, 255.0);
    }

    rgb
}

pub fn srgb_to_xyz(rgb: [f64; 3], gamma: f64, offset: f64, use_linear_part: bool) -> [f64; 3] {
    let scaled = [rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0];
    let mut linear = [0.0; 3];

    for (index, encoded) in scaled.iter().enumerate() {
        let mut value = ((*encoded - offset) / (1.0 - offset)).powf(gamma);
        if use_linear_part && value < 0.0031308 {
            value = *encoded / 12.92;
        }
        linear[index] = value;
    }

    let xyz = multiply_matrix3_vector3(SRGB_RGB_TO_XYZ, linear);
    [xyz[0] * 100.0, xyz[1] * 100.0, xyz[2] * 100.0]
}

// CIE perceptual color spaces with explicit white point input.

pub fn xyz_to_lab(xyz: [f64; 3], white_point: [f64; 3]) -> [f64; 3] {
    let fx = lab_response_curve(xyz[0], white_point[0]);
    let fy = lab_response_curve(xyz[1], white_point[1]);
    let fz = lab_response_curve(xyz[2], white_point[2]);
    let l = cie_lightness_from_ratio(xyz[1] / white_point[1]);

    [l, 500.0 * (fx - fy), 200.0 * (fy - fz)]
}

pub fn lab_to_xyz(lab: [f64; 3], white_point: [f64; 3]) -> [f64; 3] {
    let fy = (lab[0] + 16.0) / 116.0;
    let fx = lab[1] / 500.0 + fy;
    let fz = fy - lab[2] / 200.0;

    [
        lab_inverse_response_curve(fx, white_point[0]),
        lab_inverse_response_curve(fy, white_point[1]),
        lab_inverse_response_curve(fz, white_point[2]),
    ]
}

pub fn xyz_to_luv(xyz: [f64; 3], white_point: [f64; 3]) -> [f64; 3] {
    let yuv = xyz_to_yuv(xyz);
    let white_yuv = xyz_to_yuv(white_point);
    let y_ratio = yuv[0] / white_yuv[0];
    let l = cie_lightness_from_ratio(y_ratio);

    [
        l,
        13.0 * l * (yuv[1] - white_yuv[1]),
        13.0 * l * (yuv[2] - white_yuv[2]),
    ]
}

pub fn luv_to_xyz(luv: [f64; 3], white_point: [f64; 3]) -> [f64; 3] {
    let white_yuv = xyz_to_yuv(white_point);
    let mut yuv = [0.0; 3];
    if luv[0] == 0.0 {
        yuv[1] = 0.0;
        yuv[2] = 0.0;
    } else {
        yuv[1] = luv[1] / (13.0 * luv[0]) + white_yuv[1];
        yuv[2] = luv[2] / (13.0 * luv[0]) + white_yuv[2];
    }

    yuv[0] = white_yuv[0] * cie_y_ratio_from_lightness(luv[0]);

    yuv_to_xyz(yuv)
}

// Color difference.

fn delta_e_lab(lab1: [f64; 3], lab2: [f64; 3], formula: DeltaEFormula) -> f64 {
    match formula {
        DeltaEFormula::Cie76 => delta_e_cie76_lab(lab1, lab2),
        DeltaEFormula::Ciede2000 => delta_e_ciede2000_lab(lab1, lab2),
    }
}

pub fn delta_e(
    xyz1: [f64; 3],
    xyz2: [f64; 3],
    white_point: [f64; 3],
    formula: DeltaEFormula,
) -> f64 {
    delta_e_lab(
        xyz_to_lab(xyz1, white_point),
        xyz_to_lab(xyz2, white_point),
        formula,
    )
}

pub fn delta_e_cie76(xyz1: [f64; 3], xyz2: [f64; 3], white_point: [f64; 3]) -> f64 {
    delta_e(xyz1, xyz2, white_point, DeltaEFormula::Cie76)
}

pub fn delta_e_ciede2000(xyz1: [f64; 3], xyz2: [f64; 3], white_point: [f64; 3]) -> f64 {
    delta_e(xyz1, xyz2, white_point, DeltaEFormula::Ciede2000)
}

fn delta_e_cie76_lab(lab1: [f64; 3], lab2: [f64; 3]) -> f64 {
    let dl = lab1[0] - lab2[0];
    let da = lab1[1] - lab2[1];
    let db = lab1[2] - lab2[2];
    (dl * dl + da * da + db * db).sqrt()
}

fn delta_e_ciede2000_lab(lab1: [f64; 3], lab2: [f64; 3]) -> f64 {
    let (l1, a1, b1) = (lab1[0], lab1[1], lab1[2]);
    let (l2, a2, b2) = (lab2[0], lab2[1], lab2[2]);

    let c1 = (a1 * a1 + b1 * b1).sqrt();
    let c2 = (a2 * a2 + b2 * b2).sqrt();
    let c_bar = (c1 + c2) / 2.0;
    let c_bar7 = c_bar.powi(7);
    let g = 0.5 * (1.0 - (c_bar7 / (c_bar7 + 25_f64.powi(7))).sqrt());

    let a1_prime = (1.0 + g) * a1;
    let a2_prime = (1.0 + g) * a2;
    let c1_prime = (a1_prime * a1_prime + b1 * b1).sqrt();
    let c2_prime = (a2_prime * a2_prime + b2 * b2).sqrt();
    let h1_prime = if c1_prime == 0.0 {
        0.0
    } else {
        hue_angle_degrees(a1_prime, b1)
    };
    let h2_prime = if c2_prime == 0.0 {
        0.0
    } else {
        hue_angle_degrees(a2_prime, b2)
    };

    let delta_l_prime = l2 - l1;
    let delta_c_prime = c2_prime - c1_prime;

    let delta_h_prime = if c1_prime == 0.0 || c2_prime == 0.0 {
        0.0
    } else {
        let mut delta = h2_prime - h1_prime;
        if delta > 180.0 {
            delta -= 360.0;
        } else if delta < -180.0 {
            delta += 360.0;
        }
        delta
    };
    let delta_big_h_prime =
        2.0 * (c1_prime * c2_prime).sqrt() * degrees_to_radians(delta_h_prime / 2.0).sin();

    let l_bar_prime = (l1 + l2) / 2.0;
    let c_bar_prime = (c1_prime + c2_prime) / 2.0;
    let h_bar_prime = if c1_prime == 0.0 || c2_prime == 0.0 {
        h1_prime + h2_prime
    } else if (h1_prime - h2_prime).abs() > 180.0 {
        if h1_prime + h2_prime < 360.0 {
            (h1_prime + h2_prime + 360.0) / 2.0
        } else {
            (h1_prime + h2_prime - 360.0) / 2.0
        }
    } else {
        (h1_prime + h2_prime) / 2.0
    };

    let t = 1.0 - 0.17 * degrees_to_radians(h_bar_prime - 30.0).cos()
        + 0.24 * degrees_to_radians(2.0 * h_bar_prime).cos()
        + 0.32 * degrees_to_radians(3.0 * h_bar_prime + 6.0).cos()
        - 0.20 * degrees_to_radians(4.0 * h_bar_prime - 63.0).cos();
    let delta_theta = 30.0 * (-(((h_bar_prime - 275.0) / 25.0).powi(2))).exp();
    let c_bar_prime7 = c_bar_prime.powi(7);
    let r_c = 2.0 * (c_bar_prime7 / (c_bar_prime7 + 25_f64.powi(7))).sqrt();
    let s_l =
        1.0 + (0.015 * (l_bar_prime - 50.0).powi(2)) / (20.0 + (l_bar_prime - 50.0).powi(2)).sqrt();
    let s_c = 1.0 + 0.045 * c_bar_prime;
    let s_h = 1.0 + 0.015 * c_bar_prime * t;
    let r_t = -degrees_to_radians(2.0 * delta_theta).sin() * r_c;

    let l_term = delta_l_prime / s_l;
    let c_term = delta_c_prime / s_c;
    let h_term = delta_big_h_prime / s_h;

    (l_term * l_term + c_term * c_term + h_term * h_term + r_t * c_term * h_term).sqrt()
}

pub fn get_cie_mesopic_adaptation(
    photopic_luminance: &[f64],
    scotopic_luminance: Option<&[f64]>,
    s_p_ratio: Option<&[f64]>,
) -> LuxResult<(Vec<f64>, Vec<f64>)> {
    if photopic_luminance.is_empty() {
        return Err(LuxError::EmptyInput);
    }
    if scotopic_luminance.is_some() == s_p_ratio.is_some() {
        return Err(LuxError::InvalidInput(
            "provide exactly one of scotopic_luminance or s_p_ratio",
        ));
    }

    let len = photopic_luminance.len();
    if let Some(ls) = scotopic_luminance {
        if ls.len() != len {
            return Err(LuxError::MismatchedLengths {
                wavelengths: len,
                values: ls.len(),
            });
        }
    }
    if let Some(sp) = s_p_ratio {
        if sp.len() != len {
            return Err(LuxError::MismatchedLengths {
                wavelengths: len,
                values: sp.len(),
            });
        }
    }

    let mut lmes = Vec::with_capacity(len);
    let mut m_values = Vec::with_capacity(len);

    for index in 0..len {
        let lp = photopic_luminance[index];
        if !lp.is_finite() || lp <= 0.0 {
            return Err(LuxError::InvalidInput(
                "photopic luminance values must be finite and positive",
            ));
        }

        let sp = if let Some(ls) = scotopic_luminance {
            let scotopic = ls[index];
            if !scotopic.is_finite() || scotopic < 0.0 {
                return Err(LuxError::InvalidInput(
                    "scotopic luminance values must be finite and non-negative",
                ));
            }
            scotopic / lp
        } else {
            let ratio = s_p_ratio.unwrap()[index];
            if !ratio.is_finite() || ratio < 0.0 {
                return Err(LuxError::InvalidInput(
                    "S/P ratio values must be finite and non-negative",
                ));
            }
            ratio
        };

        let f_lmes = |m: f64| {
            ((m * lp) + (1.0 - m) * sp * 683.0 / 1699.0) / (m + (1.0 - m) * 683.0 / 1699.0)
        };
        let f_m = |m: f64| 0.767 + 0.3334 * f_lmes(m).log10();

        let mut previous = 0.5;
        let mut current = f_m(previous);
        let mut iterations = 0;
        while (current - previous).abs() > 1e-12 && iterations < 100 {
            previous = current;
            current = f_m(previous);
            iterations += 1;
        }

        lmes.push(f_lmes(current));
        m_values.push(current.clamp(0.0, 1.0));
    }

    Ok((lmes, m_values))
}

pub fn vlbar_cie_mesopic(
    m_values: &[f64],
    target_wavelengths: Option<&[f64]>,
) -> LuxResult<MesopicLuminousEfficiency> {
    if m_values.is_empty() {
        return Err(LuxError::EmptyInput);
    }

    let photopic = Observer::Cie1931_2.vlbar()?.0;
    let wavelengths = photopic.wavelengths().to_vec();
    let scotopic = load_scotopic_vlbar_on(&wavelengths)?;
    let peak_index = wavelengths
        .iter()
        .position(|&wavelength| (wavelength - 555.0).abs() < 1e-12)
        .ok_or(LuxError::ParseError(
            "missing 555 nm in mesopic source data",
        ))?;

    let mut curves = Vec::with_capacity(m_values.len());
    let mut k_mesopic = Vec::with_capacity(m_values.len());

    for &m in m_values {
        let m = m.clamp(0.0, 1.0);
        let values: Vec<f64> = photopic
            .values()
            .iter()
            .zip(scotopic.values().iter())
            .map(|(vp, vs)| m * vp + (1.0 - m) * vs)
            .collect();

        let k = 683.0 / values[peak_index];
        curves.push(values);
        k_mesopic.push(k);
    }

    let curves = SpectralMatrix::new(wavelengths, curves)?;
    let curves = if let Some(target_wavelengths) = target_wavelengths {
        curves.cie_interp_linear(target_wavelengths, false)?
    } else {
        curves
    };
    let normalization =
        vec![crate::spectrum::SpectrumNormalization::Max(1.0); curves.spectrum_count()];
    let curves = curves.normalize_each(&normalization, None)?;

    Ok(MesopicLuminousEfficiency { curves, k_mesopic })
}

fn load_scotopic_vlbar_on(target_wavelengths: &[f64]) -> LuxResult<Spectrum> {
    let base = TristimulusObserver::from_csv(
        include_str!("../data/cmfs/ciexyz_1951_20_scotopic.dat"),
        1699.0,
    )?
    .vl_spectrum()?;

    let source_wavelengths = base.wavelengths().to_vec();
    let interpolated = base.cie_interp_linear(target_wavelengths, false)?;
    let clipped = target_wavelengths
        .iter()
        .zip(interpolated.values().iter())
        .map(|(&wavelength, &value)| {
            if wavelength < source_wavelengths[0]
                || wavelength > source_wavelengths[source_wavelengths.len() - 1]
            {
                0.0
            } else if value.is_sign_negative() {
                0.0
            } else {
                value
            }
        })
        .collect();

    Spectrum::new(target_wavelengths.to_vec(), clipped)
}

#[cfg(test)]
mod tests {
    use crate::cam::{
        cam16_viewing_conditions, ciecam02_viewing_conditions, CamSurround as ModelCamSurround,
        CamUcsType,
    };

    use super::{
        cat_apply, cat_apply_context, cat_apply_mode, cat_apply_mode_with_conditions,
        cat_apply_with_conditions, cat_compile, cat_compile_context, cat_compile_mode,
        cat_compile_mode_with_conditions, cat_compile_with_conditions, cat_degree_of_adaptation,
        cat_mode_degrees_from_conditions, delta_e, delta_e_cie76, delta_e_cie76_lab,
        delta_e_ciede2000, delta_e_ciede2000_lab, get_cie_mesopic_adaptation, lab_to_xyz,
        lms_to_xyz, luv_to_xyz, srgb_to_xyz, vlbar_cie_mesopic, xyz_to_lab, xyz_to_lms, xyz_to_luv,
        xyz_to_srgb, xyz_to_yuv, xyz_to_yxy, yuv_to_xyz, yxy_to_xyz, CatAdapter, CatContext,
        CatMode, CatSurround, CatTransform, CatViewingConditions, DeltaEFormula, Observer,
        Tristimulus, TristimulusSet,
    };

    #[test]
    fn loads_standard_observer() {
        let observer = Observer::Cie1931_2.standard().unwrap();
        assert_eq!(observer.wavelengths.first().copied(), Some(360.0));
        assert_eq!(observer.wavelengths.last().copied(), Some(830.0));
        assert_eq!(observer.wavelengths.len(), 471);
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
        let yxy = xyz_to_yxy([0.25, 0.5, 0.25]);
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
        let lab = xyz_to_lab([0.25, 0.5, 0.25], [0.5, 0.5, 0.5]);
        assert!((lab[0] - 100.0).abs() < 1e-12);
        assert!((lab[1] + 103.149_737_007_950_17).abs() < 1e-9);
        assert!((lab[2] - 41.259_894_803_180_07).abs() < 1e-9);
    }

    #[test]
    fn converts_lab_to_xyz() {
        let xyz = lab_to_xyz(
            [100.0, -103.149_737_007_950_17, 41.259_894_803_180_07],
            [0.5, 0.5, 0.5],
        );
        assert!((xyz[0] - 0.25).abs() < 1e-9);
        assert!((xyz[1] - 0.5).abs() < 1e-12);
        assert!((xyz[2] - 0.25).abs() < 1e-9);
    }

    #[test]
    fn converts_xyz_to_luv() {
        let luv = xyz_to_luv([0.25, 0.5, 0.25], [0.5, 0.5, 0.5]);
        assert!((luv[0] - 100.0).abs() < 1e-12);
        assert!((luv[1] + 120.743_034_055_727_58).abs() < 1e-9);
        assert!((luv[2] - 72.445_820_433_436_54).abs() < 1e-9);
    }

    #[test]
    fn converts_luv_to_xyz() {
        let xyz = luv_to_xyz(
            [100.0, -120.743_034_055_727_58, 72.445_820_433_436_54],
            [0.5, 0.5, 0.5],
        );
        assert!((xyz[0] - 0.25).abs() < 1e-9);
        assert!((xyz[1] - 0.5).abs() < 1e-12);
        assert!((xyz[2] - 0.25).abs() < 1e-9);
    }

    #[test]
    fn computes_delta_e_cie76() {
        let white = [95.047, 100.0, 108.883];
        let xyz1 = lab_to_xyz([50.0, 2.5, -80.0], white);
        let xyz2 = lab_to_xyz([50.0, 0.0, -82.5], white);
        let delta = delta_e_cie76(xyz1, xyz2, white);
        assert!((delta - 3.535_533_905_932_737_8).abs() < 1e-12);
    }

    #[test]
    fn computes_delta_e_ciede2000() {
        let white = [95.047, 100.0, 108.883];
        let xyz1 = lab_to_xyz([50.0, 2.6772, -79.7751], white);
        let xyz2 = lab_to_xyz([50.0, 0.0, -82.7485], white);
        let delta = delta_e_ciede2000(xyz1, xyz2, white);
        assert!((delta - 2.042_459_680_156_574).abs() < 1e-12);
    }

    #[test]
    fn computes_delta_e_with_formula_dispatch() {
        let white = [95.047, 100.0, 108.883];
        let xyz1 = lab_to_xyz([50.0, 2.6772, -79.7751], white);
        let xyz2 = lab_to_xyz([50.0, 0.0, -82.7485], white);
        let delta = delta_e(xyz1, xyz2, white, DeltaEFormula::Ciede2000);
        assert!((delta - 2.042_459_680_156_574).abs() < 1e-12);
    }

    #[test]
    fn computes_delta_e_from_xyz() {
        let white = [95.047, 100.0, 108.883];
        let lab1 = [50.0, 2.5, -80.0];
        let lab2 = [50.0, 0.0, -82.5];
        let xyz1 = lab_to_xyz(lab1, white);
        let xyz2 = lab_to_xyz(lab2, white);
        let delta = delta_e(xyz1, xyz2, white, DeltaEFormula::Cie76);
        assert!((delta - delta_e_cie76_lab(lab1, lab2)).abs() < 1e-9);
    }

    #[test]
    fn internal_lab_paths_match_xyz_paths() {
        let white = [95.047, 100.0, 108.883];
        let lab1 = [50.0, 2.6772, -79.7751];
        let lab2 = [50.0, 0.0, -82.7485];
        let xyz1 = lab_to_xyz(lab1, white);
        let xyz2 = lab_to_xyz(lab2, white);
        assert!((delta_e_cie76(xyz1, xyz2, white) - delta_e_cie76_lab(lab1, lab2)).abs() < 1e-12);
        assert!(
            (delta_e_ciede2000(xyz1, xyz2, white) - delta_e_ciede2000_lab(lab1, lab2)).abs()
                < 1e-12
        );
    }

    #[test]
    fn applies_bradford_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_203_102_921_214).abs() < 1e-12);
        assert!((adapted[1] - 19.999_901_615_516_674).abs() < 1e-12);
        assert!((adapted[2] - 7.118_055_791_689_174).abs() < 1e-12);
    }

    #[test]
    fn applies_cat02_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Cat02,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_153_635_389_728).abs() < 1e-12);
        assert!((adapted[1] - 19.999_847_882_170_943).abs() < 1e-12);
        assert!((adapted[2] - 7.118_149_458_933_564).abs() < 1e-12);
    }

    #[test]
    fn applies_cat16_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Cat16,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_301_223_531_525).abs() < 1e-12);
        assert!((adapted[1] - 20.000_021_021_038_33).abs() < 1e-12);
        assert!((adapted[2] - 7.118_208_448_159_319).abs() < 1e-12);
    }

    #[test]
    fn applies_sharp_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Sharp,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_337_526_821_627).abs() < 1e-12);
        assert!((adapted[1] - 19.999_953_773_116_744).abs() < 1e-12);
        assert!((adapted[2] - 7.118_112_433_644_779).abs() < 1e-12);
    }

    #[test]
    fn applies_bianco_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bianco,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_237_997_793_32).abs() < 1e-12);
        assert!((adapted[1] - 19.999_886_291_132_23).abs() < 1e-12);
        assert!((adapted[2] - 7.118_132_338_899_736).abs() < 1e-12);
    }

    #[test]
    fn applies_cmc_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Cmc,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_245_614_232_79).abs() < 1e-12);
        assert!((adapted[1] - 19.999_939_194_170_24).abs() < 1e-12);
        assert!((adapted[2] - 7.118_164_955_975_901).abs() < 1e-12);
    }

    #[test]
    fn applies_kries_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Kries,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_131_711_895_99).abs() < 1e-12);
        assert!((adapted[1] - 19.999_997_693_482_22).abs() < 1e-12);
        assert!((adapted[2] - 7.118_111_183_564_011).abs() < 1e-12);
    }

    #[test]
    fn applies_judd1945_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Judd1945,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_117_638_953_994).abs() < 1e-12);
        assert!((adapted[1] - 20.0).abs() < 1e-12);
        assert!((adapted[2] - 7.118_111_183_564_01).abs() < 1e-12);
    }

    #[test]
    fn applies_judd1945_cie016_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Judd1945Cie016,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_113_747_706_712).abs() < 1e-12);
        assert!((adapted[1] - 20.0).abs() < 1e-12);
        assert!((adapted[2] - 7.118_111_183_564_01).abs() < 1e-12);
    }

    #[test]
    fn applies_judd1935_chromatic_adaptation() {
        let adapted = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Judd1935,
            1.0,
        )
        .unwrap();
        assert!((adapted[0] - 21.970_394_658_200_817).abs() < 1e-12);
        assert!((adapted[1] - 20.000_197_359_403_444).abs() < 1e-12);
        assert!((adapted[2] - 7.118_111_183_564_01).abs() < 1e-12);
    }

    #[test]
    fn rejects_invalid_adaptation_degree() {
        let err = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            1.5,
        )
        .unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid input: degree_of_adaptation must be finite and within 0..=1"
        );
    }

    #[test]
    fn computes_degree_of_adaptation_for_average_surround() {
        let degree = cat_degree_of_adaptation(CatSurround::Average, 318.31).unwrap();
        assert!((degree - 0.994_468_780_088_437_4).abs() < 1e-12);
    }

    #[test]
    fn computes_degree_of_adaptation_for_dim_surround() {
        let degree = cat_degree_of_adaptation(CatSurround::Dim, 20.0).unwrap();
        assert!((degree - 0.772_572_461_903_455_1).abs() < 1e-12);
    }

    #[test]
    fn computes_degree_of_adaptation_for_dark_surround() {
        let degree = cat_degree_of_adaptation(CatSurround::Dark, 0.0).unwrap();
        assert!((degree - 0.659_225_947_140_2).abs() < 1e-12);
    }

    #[test]
    fn computes_degree_of_adaptation_from_viewing_conditions() {
        let conditions = CatViewingConditions::new(CatSurround::Average, 318.31).unwrap();
        let degree = conditions.degree_of_adaptation().unwrap();
        assert!((degree - 0.994_468_780_088_437_4).abs() < 1e-12);
    }

    #[test]
    fn applies_chromatic_adaptation_with_conditions() {
        let adapted = cat_apply_with_conditions(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            CatSurround::Average,
            318.31,
        )
        .unwrap();
        assert!((adapted[0] - 21.953_829_568_576_072).abs() < 1e-12);
        assert!((adapted[1] - 19.999_902_159_702_89).abs() < 1e-12);
        assert!((adapted[2] - 7.199_154_229_436_402).abs() < 1e-12);
    }

    #[test]
    fn resolves_mode_degrees_from_viewing_conditions() {
        let source = CatViewingConditions::new(CatSurround::Average, 318.31).unwrap();
        let target = CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap();
        let degrees = cat_mode_degrees_from_conditions(CatMode::TwoStep, source, target).unwrap();
        assert!((degrees[0] - 0.994_468_780_088_437_4).abs() < 1e-12);
        assert!((degrees[1] - 0.772_572_461_903_455_1).abs() < 1e-12);

        let target_only =
            cat_mode_degrees_from_conditions(CatMode::BaselineToTarget, source, target).unwrap();
        assert!((target_only[0] - 0.772_572_461_903_455_1).abs() < 1e-12);
    }

    #[test]
    fn applies_mode_chromatic_adaptation_with_conditions() {
        let source = CatViewingConditions::new(CatSurround::Average, 318.31).unwrap();
        let target = CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap();
        let adapted = cat_apply_mode_with_conditions(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            source,
            target,
        )
        .unwrap();
        let manual = cat_apply_mode(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            [
                source.degree_of_adaptation().unwrap(),
                target.degree_of_adaptation().unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(adapted, manual);
    }

    #[test]
    fn baseline_to_target_mode_uses_target_conditions_degree() {
        let source = CatViewingConditions::new(CatSurround::Average, 318.31).unwrap();
        let target = CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap();
        let adapted = cat_apply_mode_with_conditions(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::BaselineToTarget,
            source,
            target,
        )
        .unwrap();
        let manual = cat_apply_mode(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::BaselineToTarget,
            [
                target.degree_of_adaptation().unwrap(),
                target.degree_of_adaptation().unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(adapted, manual);
    }

    #[test]
    fn context_exposes_default_baseline_and_mode_degrees() {
        let context = CatContext::new(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            None,
            CatTransform::Bradford,
            CatMode::TwoStep,
            CatViewingConditions::new(CatSurround::Average, 318.31).unwrap(),
            CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap(),
        );
        assert_eq!(context.baseline_white_or_default(), [100.0, 100.0, 100.0]);
        let degrees = context.degrees_of_adaptation().unwrap();
        assert!((degrees[0] - 0.994_468_780_088_437_4).abs() < 1e-12);
        assert!((degrees[1] - 0.772_572_461_903_455_1).abs() < 1e-12);
    }

    #[test]
    fn applies_chromatic_adaptation_from_context() {
        let context = CatContext::new(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            CatViewingConditions::new(CatSurround::Average, 318.31).unwrap(),
            CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap(),
        );
        let adapted = cat_apply_context([19.01, 20.0, 21.78], context).unwrap();
        let manual = cat_apply_mode_with_conditions(
            [19.01, 20.0, 21.78],
            context.source_white,
            context.target_white,
            context.baseline_white,
            context.transform,
            context.mode,
            context.source_conditions,
            context.target_conditions,
        )
        .unwrap();
        assert_eq!(adapted, manual);
    }

    #[test]
    fn compiled_adapter_matches_single_step_helper() {
        let adapter = CatAdapter::from_degree(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            1.0,
        )
        .unwrap();
        let adapted = adapter.apply([19.01, 20.0, 21.78]).unwrap();
        let helper = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            1.0,
        )
        .unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn cat_compile_matches_single_step_helper() {
        let adapter = cat_compile(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            1.0,
        )
        .unwrap();
        let adapted = adapter.apply([19.01, 20.0, 21.78]).unwrap();
        let helper = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            1.0,
        )
        .unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn compiled_mode_adapter_matches_mode_helper() {
        let adapter = CatAdapter::from_mode(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            [0.8, 0.6],
        )
        .unwrap();
        let adapted = adapter.apply([19.01, 20.0, 21.78]).unwrap();
        let helper = cat_apply_mode(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            [0.8, 0.6],
        )
        .unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn cat_compile_mode_matches_mode_helper() {
        let adapter = cat_compile_mode(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            [0.8, 0.6],
        )
        .unwrap();
        let adapted = adapter.apply([19.01, 20.0, 21.78]).unwrap();
        let helper = cat_apply_mode(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            [0.8, 0.6],
        )
        .unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn compiled_context_adapter_matches_context_helper() {
        let context = CatContext::new(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            CatViewingConditions::new(CatSurround::Average, 318.31).unwrap(),
            CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap(),
        );
        let adapter = CatAdapter::from_context(context).unwrap();
        let adapted = adapter.apply([19.01, 20.0, 21.78]).unwrap();
        let helper = cat_apply_context([19.01, 20.0, 21.78], context).unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn cat_compile_with_conditions_matches_helper() {
        let adapter = cat_compile_with_conditions(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            CatSurround::Average,
            318.31,
        )
        .unwrap();
        let adapted = adapter.apply([19.01, 20.0, 21.78]).unwrap();
        let helper = cat_apply_with_conditions(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            CatSurround::Average,
            318.31,
        )
        .unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn cat_compile_mode_with_conditions_matches_helper() {
        let source = CatViewingConditions::new(CatSurround::Average, 318.31).unwrap();
        let target = CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap();
        let adapter = cat_compile_mode_with_conditions(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            source,
            target,
        )
        .unwrap();
        let adapted = adapter.apply([19.01, 20.0, 21.78]).unwrap();
        let helper = cat_apply_mode_with_conditions(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            source,
            target,
        )
        .unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn cat_compile_context_matches_helper() {
        let context = CatContext::new(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            CatViewingConditions::new(CatSurround::Average, 318.31).unwrap(),
            CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap(),
        );
        let adapter = cat_compile_context(context).unwrap();
        let adapted = adapter.apply([19.01, 20.0, 21.78]).unwrap();
        let helper = cat_apply_context([19.01, 20.0, 21.78], context).unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn tristimulus_adapter_wrapper_matches_adapter() {
        let adapter = CatAdapter::from_degree(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            CatTransform::Bradford,
            1.0,
        )
        .unwrap();
        let adapted = Tristimulus::new([19.01, 20.0, 21.78])
            .cat_apply_adapter(adapter)
            .unwrap()
            .values();
        assert_eq!(adapted, adapter.apply([19.01, 20.0, 21.78]).unwrap());
    }

    #[test]
    fn applies_two_step_bradford_chromatic_adaptation() {
        let adapted = cat_apply_mode(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            [0.8, 0.6],
        )
        .unwrap();
        assert!((adapted[0] - 20.321_183_547_718_547).abs() < 1e-12);
        assert!((adapted[1] - 19.738_985_345_802_1).abs() < 1e-12);
        assert!((adapted[2] - 9.694_619_002_109_818).abs() < 1e-12);
    }

    #[test]
    fn applies_two_step_cat16_chromatic_adaptation() {
        let adapted = cat_apply_mode(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Cat16,
            CatMode::TwoStep,
            [0.8, 0.6],
        )
        .unwrap();
        assert!((adapted[0] - 20.564_644_514_788_387).abs() < 1e-12);
        assert!((adapted[1] - 20.001_575_611_836_01).abs() < 1e-12);
        assert!((adapted[2] - 9.934_161_728_245_801).abs() < 1e-12);
    }

    #[test]
    fn source_to_baseline_mode_matches_one_step_helper() {
        let adapted = cat_apply_mode(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::SourceToBaseline,
            [1.0, 1.0],
        )
        .unwrap();
        let helper = cat_apply(
            [19.01, 20.0, 21.78],
            [95.047, 100.0, 108.883],
            [100.0, 100.0, 100.0],
            CatTransform::Bradford,
            1.0,
        )
        .unwrap();
        assert_eq!(adapted, helper);
    }

    #[test]
    fn batch_chromaticity_transforms_match_scalar_versions() {
        let xyz = [[0.25, 0.5, 0.25], [0.2, 0.3, 0.4]];
        assert_eq!(
            TristimulusSet::new(xyz.to_vec()).xyz_to_yxy().into_vec(),
            vec![xyz_to_yxy(xyz[0]), xyz_to_yxy(xyz[1])]
        );
        let yxy = [
            [0.5, 0.25, 0.5],
            [0.3, 0.222_222_222_222_222_2, 0.333_333_333_333_333_3],
        ];
        assert_eq!(
            TristimulusSet::new(yxy.to_vec()).yxy_to_xyz().into_vec(),
            vec![yxy_to_xyz(yxy[0]), yxy_to_xyz(yxy[1])]
        );
        assert_eq!(
            TristimulusSet::new(xyz.to_vec()).xyz_to_yuv().into_vec(),
            vec![xyz_to_yuv(xyz[0]), xyz_to_yuv(xyz[1])]
        );
        let yuv = [
            [0.5, 0.117_647_058_823_529_41, 0.529_411_764_705_882_4],
            [0.3, 0.129_032_258_064_516_13, 0.435_483_870_967_741_94],
        ];
        assert_eq!(
            TristimulusSet::new(yuv.to_vec()).yuv_to_xyz().into_vec(),
            vec![yuv_to_xyz(yuv[0]), yuv_to_xyz(yuv[1])]
        );
    }

    #[test]
    fn batch_color_space_transforms_match_scalar_versions() {
        let xyz = [[0.25, 0.5, 0.25], [0.2, 0.3, 0.4]];
        let white = [0.5, 0.5, 0.5];
        let xyz_set = TristimulusSet::new(xyz.to_vec());
        let lab = xyz_set.xyz_to_lab(white).into_vec();
        assert_eq!(
            lab,
            vec![xyz_to_lab(xyz[0], white), xyz_to_lab(xyz[1], white)]
        );
        assert_eq!(
            TristimulusSet::new(lab.clone())
                .lab_to_xyz(white)
                .into_vec(),
            vec![lab_to_xyz(lab[0], white), lab_to_xyz(lab[1], white)]
        );
        let luv = xyz_set.xyz_to_luv(white).into_vec();
        assert_eq!(
            luv,
            vec![xyz_to_luv(xyz[0], white), xyz_to_luv(xyz[1], white)]
        );
        assert_eq!(
            TristimulusSet::new(luv.clone())
                .luv_to_xyz(white)
                .into_vec(),
            vec![luv_to_xyz(luv[0], white), luv_to_xyz(luv[1], white)]
        );
    }

    #[test]
    fn batch_lms_and_srgb_transforms_match_scalar_versions() {
        let xyz = [[0.25, 0.5, 0.25], [20.0, 21.0, 22.0]];
        let lms_many = TristimulusSet::new(vec![xyz[0], xyz[0]])
            .xyz_to_lms(Observer::Cie1931_2)
            .unwrap()
            .into_vec();
        assert_eq!(
            lms_many,
            vec![xyz_to_lms(xyz[0], Observer::Cie1931_2).unwrap(); 2]
        );
        let lms_input = [lms_many[0], lms_many[1]];
        assert_eq!(
            TristimulusSet::new(lms_input.to_vec())
                .lms_to_xyz(Observer::Cie1931_2)
                .unwrap()
                .into_vec(),
            vec![
                lms_to_xyz(lms_input[0], Observer::Cie1931_2).unwrap(),
                lms_to_xyz(lms_input[1], Observer::Cie1931_2).unwrap()
            ]
        );
        assert_eq!(
            TristimulusSet::new(vec![xyz[1], xyz[1]])
                .xyz_to_srgb(2.4, -0.055, true)
                .into_vec(),
            vec![xyz_to_srgb(xyz[1], 2.4, -0.055, true); 2]
        );
        let rgb = [[64.0, 128.0, 192.0], [32.0, 64.0, 96.0]];
        assert_eq!(
            TristimulusSet::new(rgb.to_vec())
                .srgb_to_xyz(2.4, -0.055, true)
                .into_vec(),
            vec![
                srgb_to_xyz(rgb[0], 2.4, -0.055, true),
                srgb_to_xyz(rgb[1], 2.4, -0.055, true)
            ]
        );
    }

    #[test]
    fn batch_cat_transforms_match_scalar_versions() {
        let xyz = [[19.01, 20.0, 21.78], [20.0, 21.0, 22.0]];
        let source_conditions = CatViewingConditions::new(CatSurround::Average, 318.31).unwrap();
        let target_conditions = CatViewingConditions::new(CatSurround::Dim, 20.0).unwrap();
        let context = CatContext::new(
            [95.047, 100.0, 108.883],
            [109.85, 100.0, 35.585],
            Some([100.0, 100.0, 100.0]),
            CatTransform::Bradford,
            CatMode::TwoStep,
            source_conditions,
            target_conditions,
        );
        let many = TristimulusSet::new(xyz.to_vec())
            .cat_apply(
                [95.047, 100.0, 108.883],
                [109.85, 100.0, 35.585],
                CatTransform::Bradford,
                1.0,
            )
            .unwrap()
            .into_vec();
        assert_eq!(
            many,
            vec![
                cat_apply(
                    xyz[0],
                    [95.047, 100.0, 108.883],
                    [109.85, 100.0, 35.585],
                    CatTransform::Bradford,
                    1.0
                )
                .unwrap(),
                cat_apply(
                    xyz[1],
                    [95.047, 100.0, 108.883],
                    [109.85, 100.0, 35.585],
                    CatTransform::Bradford,
                    1.0
                )
                .unwrap()
            ]
        );
        let context_many = TristimulusSet::new(xyz.to_vec())
            .cat_apply_context(context)
            .unwrap()
            .into_vec();
        assert_eq!(context_many[0], cat_apply_context(xyz[0], context).unwrap());
        let adapter = CatAdapter::from_context(context).unwrap();
        let adapter_many = TristimulusSet::new(xyz.to_vec())
            .cat_apply_adapter(adapter)
            .unwrap()
            .into_vec();
        assert_eq!(adapter_many[0], adapter.apply(xyz[0]).unwrap());
        assert_eq!(adapter_many[1], adapter.apply(xyz[1]).unwrap());
        let conditioned = TristimulusSet::new(xyz.to_vec())
            .cat_apply_with_conditions(
                [95.047, 100.0, 108.883],
                [109.85, 100.0, 35.585],
                CatTransform::Bradford,
                CatSurround::Average,
                318.31,
            )
            .unwrap()
            .into_vec();
        assert_eq!(
            conditioned[0],
            cat_apply_with_conditions(
                xyz[0],
                [95.047, 100.0, 108.883],
                [109.85, 100.0, 35.585],
                CatTransform::Bradford,
                CatSurround::Average,
                318.31
            )
            .unwrap()
        );
        let mode_many = TristimulusSet::new(xyz.to_vec())
            .cat_apply_mode(
                [95.047, 100.0, 108.883],
                [109.85, 100.0, 35.585],
                Some([100.0, 100.0, 100.0]),
                CatTransform::Bradford,
                CatMode::TwoStep,
                [0.8, 0.6],
            )
            .unwrap()
            .into_vec();
        assert_eq!(
            mode_many[0],
            cat_apply_mode(
                xyz[0],
                [95.047, 100.0, 108.883],
                [109.85, 100.0, 35.585],
                Some([100.0, 100.0, 100.0]),
                CatTransform::Bradford,
                CatMode::TwoStep,
                [0.8, 0.6]
            )
            .unwrap()
        );
        let mode_conditioned = TristimulusSet::new(xyz.to_vec())
            .cat_apply_mode_with_conditions(
                [95.047, 100.0, 108.883],
                [109.85, 100.0, 35.585],
                Some([100.0, 100.0, 100.0]),
                CatTransform::Bradford,
                CatMode::TwoStep,
                source_conditions,
                target_conditions,
            )
            .unwrap()
            .into_vec();
        assert_eq!(
            mode_conditioned[0],
            cat_apply_mode_with_conditions(
                xyz[0],
                [95.047, 100.0, 108.883],
                [109.85, 100.0, 35.585],
                Some([100.0, 100.0, 100.0]),
                CatTransform::Bradford,
                CatMode::TwoStep,
                source_conditions,
                target_conditions
            )
            .unwrap()
        );
    }

    #[test]
    fn tristimulus_wrapper_matches_scalar_transforms() {
        let xyz = Tristimulus::new([0.25, 0.5, 0.25]);
        assert_eq!(xyz.xyz_to_yxy().values(), xyz_to_yxy([0.25, 0.5, 0.25]));
        assert_eq!(
            xyz.xyz_to_lab([0.5, 0.5, 0.5]).values(),
            xyz_to_lab([0.25, 0.5, 0.25], [0.5, 0.5, 0.5])
        );
        assert_eq!(
            xyz.xyz_to_lms(Observer::Cie1931_2).unwrap().values(),
            xyz_to_lms([0.25, 0.5, 0.25], Observer::Cie1931_2).unwrap()
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
        let cam16_scalar = crate::cam::cam16_forward([0.25, 0.5, 0.25], cam16_conditions).unwrap();
        assert_eq!(cam16, cam16_scalar);
        let cam16_ucs = xyz
            .cam16_ucs_forward(cam16_conditions, CamUcsType::Ucs)
            .unwrap();
        let cam16_ucs_scalar =
            crate::cam::cam16_ucs_forward([0.25, 0.5, 0.25], cam16_conditions, CamUcsType::Ucs)
                .unwrap();
        assert_eq!(cam16_ucs, cam16_ucs_scalar);
        let cam16_back =
            Tristimulus::new([cam16_ucs.j_prime, cam16_ucs.a_prime, cam16_ucs.b_prime])
                .cam16_ucs_inverse(cam16_conditions, CamUcsType::Ucs)
                .unwrap();
        assert!((cam16_back.values()[0] - 0.25).abs() < 1e-10);
        assert!((cam16_back.values()[1] - 0.5).abs() < 1e-10);
        assert!((cam16_back.values()[2] - 0.25).abs() < 1e-10);
    }

    #[test]
    fn tristimulus_set_wrapper_matches_batch_transforms() {
        let xyz = TristimulusSet::new(vec![[0.25, 0.5, 0.25], [0.2, 0.3, 0.4]]);
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
            .map(|value| crate::cam::ciecam02_forward(value, ciecam02_conditions).unwrap())
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
                crate::cam::ciecam02_ucs_forward(value, ciecam02_conditions, CamUcsType::Ucs)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        assert_eq!(cam_ucs_many, cam_ucs_scalar);
        let ucs_triplets = TristimulusSet::new(
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
    fn tristimulus_set_delta_e_matches_pairwise_scalar_computation() {
        let left = TristimulusSet::new(vec![
            lab_to_xyz([50.0, 2.5, -80.0], [95.047, 100.0, 108.883]),
            lab_to_xyz([50.0, 2.6772, -79.7751], [95.047, 100.0, 108.883]),
        ]);
        let right = TristimulusSet::new(vec![
            lab_to_xyz([50.0, 0.0, -82.5], [95.047, 100.0, 108.883]),
            lab_to_xyz([50.0, 0.0, -82.7485], [95.047, 100.0, 108.883]),
        ]);
        let result = left
            .delta_e(&right, [95.047, 100.0, 108.883], DeltaEFormula::Ciede2000)
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
}
