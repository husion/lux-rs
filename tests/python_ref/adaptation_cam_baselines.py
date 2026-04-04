"""Chromatic adaptation and CAM LuxPy baseline values."""

import luxpy as lx
import numpy as np
from luxpy.color import cat
from luxpy.color.cam.cam02ucs import run as cam02ucs
from luxpy.color.cam.cam16ucs import run as cam16ucs
from luxpy.color.cam.ciecam02 import run as ciecam02
from luxpy.color.cam.ciecam16 import run as ciecam16

from baseline_helpers import fmt_scalar, fmt_vec


def iter_adaptation_cam_baselines():
    cat_xyz = np.array([[19.01, 20.0, 21.78]])
    cat_w1 = np.array([[95.047, 100.0, 108.883]])
    cat_w2 = np.array([[109.85, 100.0, 35.585]])
    cat_d_avg = cat.get_degree_of_adaptation(Dtype="cat02", F="avg", La=318.31)[0]
    cat_d_dim = cat.get_degree_of_adaptation(Dtype="cat16", F="dim", La=20.0)[0]
    cat_d_dark = cat.get_degree_of_adaptation(Dtype="cat16", F="dark", La=0.0)[0]

    yield "cat_d_avg", fmt_scalar(cat_d_avg)
    yield "cat_d_dim", fmt_scalar(cat_d_dim)
    yield "cat_d_dark", fmt_scalar(cat_d_dark)
    yield "cat_bradford", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="bfd").ravel())
    yield "cat_cat02", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="cat02").ravel())
    yield "cat_cat16", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="cat16").ravel())
    yield "cat_sharp", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="sharp").ravel())
    yield "cat_bianco", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="bianco").ravel())
    yield "cat_cmc", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="cmc").ravel())
    yield "cat_kries", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="kries").ravel())
    yield "cat_judd1945", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="judd-1945").ravel())
    yield "cat_judd1945_cie016", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="judd-1945-CIE016").ravel())
    yield "cat_judd1935", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="judd-1935").ravel())
    yield "cat_bradford_avg", fmt_vec(cat.apply_vonkries1(cat_xyz, xyzw1=cat_w1, xyzw2=cat_w2, mcat="bfd", D=cat_d_avg).ravel())
    yield "cat_two_step_bradford", fmt_vec(cat.apply_vonkries(cat_xyz, cat_w1, cat_w2, xyzw0=np.array([[100.0, 100.0, 100.0]]), D=[0.8, 0.6], mcat="bfd", catmode="1>0>2").ravel())
    yield "cat_two_step_cat16", fmt_vec(cat.apply_vonkries(cat_xyz, cat_w1, cat_w2, xyzw0=np.array([[100.0, 100.0, 100.0]]), D=[0.8, 0.6], mcat="cat16", catmode="1>0>2").ravel())

    cam_conditions = {"La": 100.0, "Yb": 20.0, "surround": "avg", "D": 1.0, "Dtype": None}
    yield "cam16_forward", fmt_vec(ciecam16(cat_xyz, xyzw=cat_w1, conditions=cam_conditions, outin="J,Q,C,M,s,h,aM,bM,aC,bC").ravel())
    yield "ciecam02_forward", fmt_vec(ciecam02(cat_xyz, xyzw=cat_w1, conditions=cam_conditions, outin="J,Q,C,M,s,h,aM,bM,aC,bC").ravel())
    yield "cam16_ucs", fmt_vec(cam16ucs(cat_xyz, xyzw=cat_w1, conditions=cam_conditions).ravel())
    yield "cam16_lcd", fmt_vec(cam16ucs(cat_xyz, xyzw=cat_w1, conditions=cam_conditions, ucstype="lcd").ravel())
    yield "cam16_scd", fmt_vec(cam16ucs(cat_xyz, xyzw=cat_w1, conditions=cam_conditions, ucstype="scd").ravel())
    yield "cam02_ucs", fmt_vec(cam02ucs(cat_xyz, xyzw=cat_w1, conditions=cam_conditions).ravel())
    yield "cam02_lcd", fmt_vec(cam02ucs(cat_xyz, xyzw=cat_w1, conditions=cam_conditions, ucstype="lcd").ravel())
    yield "cam02_scd", fmt_vec(cam02ucs(cat_xyz, xyzw=cat_w1, conditions=cam_conditions, ucstype="scd").ravel())
    yield "cam16_inverse", fmt_vec(np.array([[19.01, 20.0, 21.78]]).ravel())
    yield "cam02_inverse", fmt_vec(np.array([[19.01, 20.0, 21.78]]).ravel())
    yield "cam16ucs_inverse", fmt_vec(np.array([[19.01, 20.0, 21.78]]).ravel())
    yield "cam02ucs_inverse", fmt_vec(np.array([[19.01, 20.0, 21.78]]).ravel())
