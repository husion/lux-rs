"""Color-space conversion and delta-E LuxPy baseline values."""

import luxpy as lx
import numpy as np
from luxpy.color.deltaE import DE2000

from baseline_helpers import fmt_scalar, fmt_vec


def iter_color_baselines():
    xyz_sample = np.array([[0.25, 0.5, 0.25]])
    white_sample = np.array([[0.5, 0.5, 0.5]])
    yxy_sample = lx.xyz_to_Yxy(xyz_sample)
    yuv_sample = lx.xyz_to_Yuv(xyz_sample)
    lab_sample = lx.xyz_to_lab(xyz_sample, xyzw=white_sample)
    luv_sample = lx.xyz_to_luv(xyz_sample, xyzw=white_sample)
    lms_1931 = lx.xyz_to_lms(xyz_sample, cieobs="1931_2")
    lms_1964 = lx.xyz_to_lms(xyz_sample, cieobs="1964_10")

    yield "xyz_to_yxy", fmt_vec(yxy_sample.ravel())
    yield "yxy_to_xyz", fmt_vec(lx.Yxy_to_xyz(yxy_sample).ravel())
    yield "xyz_to_yuv", fmt_vec(yuv_sample.ravel())
    yield "yuv_to_xyz", fmt_vec(lx.Yuv_to_xyz(yuv_sample).ravel())
    yield "xyz_to_lab", fmt_vec(lab_sample.ravel())
    yield "lab_to_xyz", fmt_vec(lx.lab_to_xyz(lab_sample, xyzw=white_sample).ravel())
    yield "xyz_to_luv", fmt_vec(luv_sample.ravel())
    yield "luv_to_xyz", fmt_vec(lx.luv_to_xyz(luv_sample, xyzw=white_sample).ravel())

    white_d65 = np.array([[95.047, 100.0, 108.883]])
    delta_xyz1_cie76 = lx.lab_to_xyz(np.array([[50.0, 2.5, -80.0]]), xyzw=white_d65)
    delta_xyz2_cie76 = lx.lab_to_xyz(np.array([[50.0, 0.0, -82.5]]), xyzw=white_d65)
    delta_xyz1_ciede2000 = lx.lab_to_xyz(np.array([[50.0, 2.6772, -79.7751]]), xyzw=white_d65)
    delta_xyz2_ciede2000 = lx.lab_to_xyz(np.array([[50.0, 0.0, -82.7485]]), xyzw=white_d65)

    yield "delta_e_cie76", fmt_scalar(
        np.linalg.norm(
            lx.xyz_to_lab(delta_xyz1_cie76, xyzw=white_d65).ravel()
            - lx.xyz_to_lab(delta_xyz2_cie76, xyzw=white_d65).ravel()
        )
    )
    yield "delta_e_ciede2000", fmt_scalar(
        DE2000(
            delta_xyz1_ciede2000,
            delta_xyz2_ciede2000,
            dtype="xyz",
            xyzwt=white_d65,
            xyzwr=white_d65,
        ).ravel()[0]
    )
    yield "xyz_to_lms_1931", fmt_vec(lms_1931.ravel())
    yield "lms_to_xyz_1931", fmt_vec(lx.lms_to_xyz(lms_1931, cieobs="1931_2").ravel())
    yield "xyz_to_lms_1964", fmt_vec(lms_1964.ravel())


def iter_display_color_baselines():
    yield "xyz_to_srgb", fmt_vec(lx.xyz_to_srgb(np.array([[20.0, 21.0, 22.0]])).ravel())
    yield "srgb_to_xyz", fmt_vec(lx.srgb_to_xyz(np.array([[64.0, 128.0, 192.0]])).ravel())
