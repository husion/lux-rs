from pathlib import Path

import luxpy as lx
import numpy as np

from baseline_common import scalar_line, vec_line


def generate_source_baselines(root: Path) -> list[tuple[str, str]]:
    d65 = np.loadtxt(root / "data" / "spds" / "CIE_D65.csv", delimiter=",").T
    f_series = np.loadtxt(root / "data" / "spds" / "CIE_F_1to12_1nm.csv", delimiter=",")
    f4 = np.vstack((f_series[:, 0], f_series[:, 4]))
    x_d, y_d = lx.daylightlocus(6500.0)
    cct_sample, duv_sample = lx.xyz_to_cct(np.array([[100.0, 100.0, 100.0]]), out="cct,duv")

    return [
        vec_line(
            "blackbody_relative_6500",
            lx.blackbody(6500.0, wl3=[360.0, 365.0, 1.0], relative=True).ravel(),
        ),
        vec_line(
            "blackbody_absolute_6500_560",
            lx.blackbody(6500.0, wl3=[560.0, 560.0, 1.0], relative=False).ravel(),
        ),
        vec_line("daylightlocus_6500", [x_d, y_d]),
        vec_line(
            "daylightphase_6500",
            lx.daylightphase(6500.0, wl3=[360.0, 365.0, 1.0]).ravel(),
        ),
        vec_line(
            "daylightphase_3500",
            lx.daylightphase(3500.0, wl3=[360.0, 365.0, 1.0]).ravel(),
        ),
        scalar_line("ciera_d65", lx.cri.spd_to_ciera(d65).ravel()[0]),
        vec_line("ciera_d65_ri", lx.cri.spd_to_ciera(d65, out="Rfi").ravel()),
        scalar_line("ciera_f4", lx.cri.spd_to_ciera(f4).ravel()[0]),
        vec_line("ciera_f4_ri", lx.cri.spd_to_ciera(f4, out="Rfi").ravel()),
        scalar_line("cierf_d65", lx.cri.spd_to_cierf(d65).ravel()[0]),
        vec_line("cierf_d65_rfi", lx.cri.spd_to_cierf(d65, out="Rfi").ravel()),
        scalar_line("cierg_d65", lx.cri.spd_to_cri(d65, cri_type="cierf", out="Rg").ravel()[0]),
        scalar_line("cierf_f4", lx.cri.spd_to_cierf(f4).ravel()[0]),
        vec_line("cierf_f4_rfi", lx.cri.spd_to_cierf(f4, out="Rfi").ravel()),
        scalar_line("cierg_f4", lx.cri.spd_to_cri(f4, cri_type="cierf", out="Rg").ravel()[0]),
        vec_line(
            "cri_ref_3000_6500",
            lx.cri_ref([3000.0, 6500.0], wl3=[360.0, 365.0, 1.0]).ravel(),
        ),
        vec_line(
            "xyz_to_cct_sample",
            [np.ravel(cct_sample)[0], np.ravel(duv_sample)[0]],
        ),
        vec_line("cct_to_xyz_6500", lx.cct_to_xyz(6500.0).ravel()),
        vec_line("illuminant_A", lx._CIE_ILLUMINANTS["A"][:, 0:6].ravel()),
        vec_line("illuminant_D65", lx._CIE_ILLUMINANTS["D65"][:, 0:6].ravel()),
        vec_line("illuminant_F4", lx._CIE_ILLUMINANTS["F4"][:, 0:6].ravel()),
        vec_line("illuminant_LED_B1", lx._CIE_ILLUMINANTS["LED_B1"][:, 0:6].ravel()),
        vec_line(
            "illuminant_D50",
            lx.daylightphase(5000.0, wl3=[360.0, 365.0, 1.0], cct_is_nominal=True).ravel(),
        ),
        vec_line(
            "spdbuild_gaussian",
            lx.toolboxes.spdbuild.gaussian_spd(530, 20, wl=[380, 780, 5]).ravel(),
        ),
        vec_line(
            "spdbuild_lorentzian2",
            lx.toolboxes.spdbuild.lorentzian2_spd(530, 20, wl=[380, 780, 5]).ravel(),
        ),
        vec_line(
            "spdbuild_butterworth",
            lx.toolboxes.spdbuild.butterworth_spd(530, 20, 2, wl=[380, 780, 5]).ravel(),
        ),
        vec_line(
            "spdbuild_roundedtriangle",
            lx.toolboxes.spdbuild.roundedtriangle_spd(530, 100, 0.5, wl=[380, 780, 5]).ravel(),
        ),
        vec_line(
            "spdbuild_mono_led",
            lx.toolboxes.spdbuild.mono_led_spd(530, 20, wl=[380, 780, 5], strength_shoulder=2, bw_order=-1).ravel(),
        ),
        vec_line(
            "spdbuild_phosphor_led",
            lx.toolboxes.spdbuild.phosphor_led_spd(450, 20, wl=[380, 780, 5], strength_ph=0.5, peakwl_ph1=530, fwhm_ph1=80, strength_ph1=0.8, peakwl_ph2=560, fwhm_ph2=80, use_piecewise_fcn=True, with_wl=False).ravel(),
        ),
        vec_line(
            "spdbuild_color3mixer",
            lx.toolboxes.spdbuild.color3mixer(np.array([[100.0, 0.4, 0.4]]), np.array([[100.0, 0.2, 0.2]]), np.array([[100.0, 0.6, 0.2]]), np.array([[100.0, 0.3, 0.7]])).ravel(),
        ),
        vec_line(
            "spdbuild_colormixer_pinv",
            lx.toolboxes.spdbuild.colormixer_pinv(np.array([[100.0, 0.4, 0.4]]), np.array([[100.0, 0.2, 0.2], [100.0, 0.6, 0.2], [100.0, 0.3, 0.7]]), input_fmt='Yxy').ravel(),
        ),
        vec_line(
            "spdbuild_colormixer",
            lx.toolboxes.spdbuild.colormixer(np.array([[100.0, 0.4, 0.4]]), np.array([[100.0, 0.2, 0.2], [100.0, 0.6, 0.2], [100.0, 0.3, 0.7], [100.0, 0.5, 0.5]]), pair_strengths=np.array([0.5])).ravel(),
        ),
        vec_line(
            "spdbuild_spd_builder",
            lx.toolboxes.spdbuild.spd_builder(wl=[380, 780, 5], peakwl=450, fwhm=20, strength_ph=0.5, peakwl_ph1=530, fwhm_ph1=80, strength_ph1=0.8, peakwl_ph2=560, fwhm_ph2=80, target=np.array([[100.0, 0.4, 0.4]]), tar_type='Yxy', with_wl=False).ravel(),
        ),
    ]
