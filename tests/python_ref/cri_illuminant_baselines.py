"""CRI, illuminant, and SPD-derived LuxPy baseline values."""

from pathlib import Path

import luxpy as lx
import numpy as np

from baseline_helpers import fmt_scalar, fmt_vec


def iter_cri_illuminant_baselines(root: Path):
    d65 = np.loadtxt(root / "data" / "spds" / "CIE_D65.csv", delimiter=",").T
    f_series = np.loadtxt(root / "data" / "spds" / "CIE_F_1to12_1nm.csv", delimiter=",")
    f4 = np.vstack((f_series[:, 0], f_series[:, 4]))

    yield "ciera_d65", fmt_scalar(lx.cri.spd_to_ciera(d65).ravel()[0])
    yield "ciera_d65_ri", fmt_vec(lx.cri.spd_to_ciera(d65, out="Rfi").ravel())
    yield "ciera_f4", fmt_scalar(lx.cri.spd_to_ciera(f4).ravel()[0])
    yield "ciera_f4_ri", fmt_vec(lx.cri.spd_to_ciera(f4, out="Rfi").ravel())
    yield "cierf_d65", fmt_scalar(lx.cri.spd_to_cierf(d65).ravel()[0])
    yield "cierf_d65_rfi", fmt_vec(lx.cri.spd_to_cierf(d65, out="Rfi").ravel())
    yield "cierg_d65", fmt_scalar(lx.cri.spd_to_cri(d65, cri_type="cierf", out="Rg").ravel()[0])
    yield "cierf_f4", fmt_scalar(lx.cri.spd_to_cierf(f4).ravel()[0])
    yield "cierf_f4_rfi", fmt_vec(lx.cri.spd_to_cierf(f4, out="Rfi").ravel())
    yield "cierg_f4", fmt_scalar(lx.cri.spd_to_cri(f4, cri_type="cierf", out="Rg").ravel()[0])
    yield "cri_ref_3000_6500", fmt_vec(lx.cri_ref([3000.0, 6500.0], wl3=[360.0, 365.0, 1.0]).ravel())

    xyz_to_cct_sample = np.array([[100.0, 100.0, 100.0]])
    cct_sample, duv_sample = lx.xyz_to_cct(xyz_to_cct_sample, out="cct,duv")
    yield "xyz_to_cct_sample", fmt_vec([np.ravel(cct_sample)[0], np.ravel(duv_sample)[0]])
    yield "cct_to_xyz_6500", fmt_vec(lx.cct_to_xyz(6500.0).ravel())
    yield "illuminant_A", fmt_vec(lx._CIE_ILLUMINANTS["A"][:, 0:6].ravel())
    yield "illuminant_D65", fmt_vec(lx._CIE_ILLUMINANTS["D65"][:, 0:6].ravel())
    yield "illuminant_F4", fmt_vec(lx._CIE_ILLUMINANTS["F4"][:, 0:6].ravel())
    yield "illuminant_LED_B1", fmt_vec(lx._CIE_ILLUMINANTS["LED_B1"][:, 0:6].ravel())
    yield "illuminant_D50", fmt_vec(
        lx.daylightphase(5000.0, wl3=[360.0, 365.0, 1.0], cct_is_nominal=True).ravel()
    )

    yield "getwlr", fmt_vec(lx.getwlr([360, 365, 1]))
    yield "getwld_equal_scalar", fmt_scalar(lx.getwld(np.array([400.0, 410.0, 420.0])))
    yield "getwld_unequal", fmt_vec(lx.getwld(np.array([400.0, 410.0, 430.0])))
    yield "power_ru", fmt_scalar(
        lx.spd_to_power(np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]), "ru")[0, 0]
    )
    yield "power_pu", fmt_scalar(
        lx.spd_to_power(
            np.array([[555.0, 556.0], [1.0, 1.0]]),
            "pu",
            cieobs="1931_2",
        )[0, 0]
    )
    yield "power_qu", fmt_scalar(
        lx.spd_to_power(np.array([[500.0, 510.0], [1.0, 1.0]]), "qu")[0, 0]
    )
    yield "ler_1931", fmt_scalar(
        lx.spd_to_ler(np.array([[555.0, 556.0], [1.0, 1.0]]), cieobs="1931_2")[0, 0]
    )
    yield "ler_1964", fmt_scalar(
        lx.spd_to_ler(np.array([[555.0, 556.0], [1.0, 1.0]]), cieobs="1964_10")[0, 0]
    )
    yield "ler_many_1931", fmt_vec(
        lx.spd_to_ler(np.array([[555.0, 556.0], [1.0, 1.0], [2.0, 2.0]]), cieobs="1931_2").ravel()
    )
    yield "xyz_relative", fmt_vec(
        lx.spd_to_xyz(np.array([[555.0, 556.0], [1.0, 1.0]]), cieobs="1931_2")[0]
    )
    yield "xyz_absolute", fmt_vec(
        lx.spd_to_xyz(
            np.array([[555.0, 556.0], [1.0, 1.0]]),
            cieobs="1931_2",
            relative=False,
        )[0]
    )
    yield "xyz_relative_1964_10", fmt_vec(
        lx.spd_to_xyz(np.array([[555.0, 556.0], [1.0, 1.0]]), cieobs="1964_10")[0]
    )
    yield "xyz_relative_many", fmt_vec(
        lx.spd_to_xyz(np.array([[555.0, 556.0], [1.0, 1.0], [2.0, 2.0]]), cieobs="1931_2").ravel()
    )
    yield "xyz_absolute_many", fmt_vec(
        lx.spd_to_xyz(
            np.array([[555.0, 556.0], [1.0, 1.0], [2.0, 2.0]]),
            cieobs="1931_2",
            relative=False,
        ).ravel()
    )
