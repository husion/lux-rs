"""Observer and mesopic LuxPy baseline values."""

import luxpy as lx
import numpy as np

from baseline_helpers import fmt_scalar, fmt_vec


def iter_observer_baselines():
    xyzbar_1931 = lx.xyzbar(cieobs="1931_2")
    vlbar_1931, k_1931 = lx.vlbar(cieobs="1931_2", out=2)
    xyzbar_1964 = lx.xyzbar(cieobs="1964_10")
    vlbar_1964, k_1964 = lx.vlbar(cieobs="1964_10", out=2)

    yield "xyzbar_1931_shape", f"{xyzbar_1931.shape[0]},{xyzbar_1931.shape[1]}"
    yield "xyzbar_1931_555", fmt_vec(xyzbar_1931[:, xyzbar_1931[0] == 555.0].ravel())
    yield "vlbar_1931_shape", f"{vlbar_1931.shape[0]},{vlbar_1931.shape[1]}"
    yield "vlbar_1931_555", fmt_vec(vlbar_1931[:, vlbar_1931[0] == 555.0].ravel())
    yield "vlbar_1931_k", fmt_scalar(k_1931)
    yield "xyzbar_1964_555", fmt_vec(xyzbar_1964[:, xyzbar_1964[0] == 555.0].ravel())
    yield "vlbar_1964_555", fmt_vec(vlbar_1964[:, vlbar_1964[0] == 555.0].ravel())
    yield "vlbar_1964_k", fmt_scalar(k_1964)

    lmes_1, m_1 = lx.get_cie_mesopic_adaptation(1.0, SP=1.0)
    vlbar_mesopic, k_mesopic = lx.vlbar_cie_mesopic(
        m=[0.5, 1.0],
        wl_new=np.array([555.0]),
        out=2,
    )
    yield "mesopic_lmes_sp_1", fmt_scalar(np.ravel(lmes_1)[0])
    yield "mesopic_m_sp_1", fmt_scalar(np.ravel(m_1)[0])
    yield "mesopic_vlbar_shape", f"{vlbar_mesopic.shape[0]},{vlbar_mesopic.shape[1]}"
    yield "mesopic_vlbar_555", fmt_vec(vlbar_mesopic.ravel())
    yield "mesopic_k", fmt_vec(np.ravel(k_mesopic))
