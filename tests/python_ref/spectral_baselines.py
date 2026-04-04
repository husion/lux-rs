"""Interpolation, normalization, and daylight/blackbody LuxPy baselines."""

import luxpy as lx
import numpy as np

from baseline_helpers import fmt_vec


def iter_spectral_baselines():
    yield "xyzbar_1931_interp", fmt_vec(
        lx.xyzbar(
            cieobs="1931_2",
            wl_new=np.array([554.5, 555.0, 555.5, 556.0]),
            kind="linear",
            extrap_kind="linear",
        ).ravel()
    )
    yield "vlbar_1931_interp", fmt_vec(
        lx.vlbar(
            cieobs="1931_2",
            wl_new=np.array([554.5, 555.0, 555.5, 556.0]),
            kind="linear",
            extrap_kind="linear",
        ).ravel()
    )
    yield "spd_interp", fmt_vec(
        lx.cie_interp(
            np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
            np.array([395.0, 405.0, 420.0, 425.0]),
            kind="linear",
            extrap_kind="linear",
        ).ravel()
    )
    yield "normalize_max", fmt_vec(
        lx.spd_normalize(
            np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
            norm_type="max",
            norm_f=2.0,
        ).ravel()
    )
    yield "normalize_area", fmt_vec(
        lx.spd_normalize(
            np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
            norm_type="area",
            norm_f=1.0,
        ).ravel()
    )
    yield "normalize_lambda", fmt_vec(
        lx.spd_normalize(
            np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
            norm_type="lambda",
            norm_f=410.0,
        ).ravel()
    )
    yield "normalize_ru", fmt_vec(
        lx.spd_normalize(
            np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
            norm_type="ru",
            norm_f=10.0,
        ).ravel()
    )
    yield "normalize_pu", fmt_vec(
        lx.spd_normalize(
            np.array([[555.0, 556.0], [1.0, 1.0]]),
            norm_type="pu",
            norm_f=1000.0,
            cieobs="1931_2",
        ).ravel()
    )
    yield "normalize_qu", fmt_vec(
        lx.spd_normalize(
            np.array([[500.0, 510.0], [1.0, 1.0]]),
            norm_type="qu",
            norm_f=1e18,
        ).ravel()
    )
    yield "blackbody_relative_6500", fmt_vec(
        lx.blackbody(6500.0, wl3=[360.0, 365.0, 1.0], relative=True).ravel()
    )
    yield "blackbody_absolute_6500_560", fmt_vec(
        lx.blackbody(6500.0, wl3=[560.0, 560.0, 1.0], relative=False).ravel()
    )

    x_d, y_d = lx.daylightlocus(6500.0)
    yield "daylightlocus_6500", fmt_vec([x_d, y_d])
    yield "daylightphase_6500", fmt_vec(lx.daylightphase(6500.0, wl3=[360.0, 365.0, 1.0]).ravel())
    yield "daylightphase_3500", fmt_vec(lx.daylightphase(3500.0, wl3=[360.0, 365.0, 1.0]).ravel())
