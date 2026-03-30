pub mod cam;
pub mod color;
pub mod error;
pub mod illuminants;
pub mod photometry;
pub mod spectrum;

#[allow(deprecated)]
pub use cam::{
    cam16_forward, cam16_ucs_forward, cam16_ucs_inverse, cam16_viewing_conditions, cam_forward,
    cam_inverse, cam_naka_rushton, cam_ucs_forward, cam_ucs_inverse, ciecam02_forward,
    ciecam02_ucs_forward, ciecam02_ucs_inverse, ciecam02_viewing_conditions, CamAppearance,
    CamModel, CamNakaRushtonParameters, CamSurround, CamSurroundParameters, CamUcsAppearance,
    CamUcsParameters, CamUcsType, CamViewingConditions,
};
#[allow(deprecated)]
pub use color::{
    cat_apply, cat_apply_context, cat_apply_mode, cat_apply_mode_with_conditions,
    cat_apply_with_conditions, cat_compile, cat_compile_context, cat_compile_mode,
    cat_compile_mode_with_conditions, cat_compile_with_conditions, cat_degree_of_adaptation,
    cat_mode_degrees_from_conditions, delta_e, delta_e_cie76, delta_e_ciede2000,
    get_cie_mesopic_adaptation, lab_to_xyz, lms_to_xyz, lms_to_xyz_with_matrix, luv_to_xyz,
    srgb_to_xyz, vlbar_cie_mesopic, xyz_to_lab, xyz_to_lms, xyz_to_lms_with_matrix, xyz_to_luv,
    xyz_to_srgb, xyz_to_yuv, xyz_to_yxy, yuv_to_xyz, yxy_to_xyz, CatAdapter, CatContext, CatMode,
    CatSurround, CatTransform, CatViewingConditions, DeltaEFormula, Matrix3,
    MesopicLuminousEfficiency, Observer, Tristimulus, TristimulusObserver, TristimulusSet,
};
pub use error::{LuxError, LuxResult};
pub use illuminants::{
    blackbody, cct_to_xyz, cri_ref, daylightlocus, daylightphase, standard_illuminant,
    standard_illuminant_names, xyz_to_cct,
};
pub use photometry::{spd_to_ler, spd_to_power, spd_to_xyz, PowerType};
pub use spectrum::{
    getwld, getwlr, SpectralMatrix, Spectrum, SpectrumNormalization, WavelengthGrid,
};
