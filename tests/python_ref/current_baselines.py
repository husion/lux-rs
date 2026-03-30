import luxpy as lx
import numpy as np


def fmt_scalar(value: float) -> str:
    return repr(float(value))


def fmt_vec(values) -> str:
    return ",".join(repr(float(value)) for value in values)


def main() -> None:
    xyzbar_1931 = lx.xyzbar(cieobs="1931_2")
    vlbar_1931, k_1931 = lx.vlbar(cieobs="1931_2", out=2)
    xyzbar_1964 = lx.xyzbar(cieobs="1964_10")
    vlbar_1964, k_1964 = lx.vlbar(cieobs="1964_10", out=2)

    print(f"xyzbar_1931_shape={xyzbar_1931.shape[0]},{xyzbar_1931.shape[1]}")
    print(f"xyzbar_1931_555={fmt_vec(xyzbar_1931[:, xyzbar_1931[0] == 555.0].ravel())}")
    print(f"vlbar_1931_shape={vlbar_1931.shape[0]},{vlbar_1931.shape[1]}")
    print(f"vlbar_1931_555={fmt_vec(vlbar_1931[:, vlbar_1931[0] == 555.0].ravel())}")
    print(f"vlbar_1931_k={fmt_scalar(k_1931)}")
    print(f"xyzbar_1964_555={fmt_vec(xyzbar_1964[:, xyzbar_1964[0] == 555.0].ravel())}")
    print(f"vlbar_1964_555={fmt_vec(vlbar_1964[:, vlbar_1964[0] == 555.0].ravel())}")
    print(f"vlbar_1964_k={fmt_scalar(k_1964)}")
    lmes_1, m_1 = lx.get_cie_mesopic_adaptation(1.0, SP=1.0)
    vlbar_mesopic, k_mesopic = lx.vlbar_cie_mesopic(
        m=[0.5, 1.0],
        wl_new=np.array([555.0]),
        out=2,
    )
    print(f"mesopic_lmes_sp_1={fmt_scalar(np.ravel(lmes_1)[0])}")
    print(f"mesopic_m_sp_1={fmt_scalar(np.ravel(m_1)[0])}")
    print(f"mesopic_vlbar_shape={vlbar_mesopic.shape[0]},{vlbar_mesopic.shape[1]}")
    print(f"mesopic_vlbar_555={fmt_vec(vlbar_mesopic.ravel())}")
    print(f"mesopic_k={fmt_vec(np.ravel(k_mesopic))}")
    xyz_sample = np.array([[0.25, 0.5, 0.25]])
    white_sample = np.array([[0.5, 0.5, 0.5]])
    yxy_sample = lx.xyz_to_Yxy(xyz_sample)
    yuv_sample = lx.xyz_to_Yuv(xyz_sample)
    lab_sample = lx.xyz_to_lab(xyz_sample, xyzw=white_sample)
    luv_sample = lx.xyz_to_luv(xyz_sample, xyzw=white_sample)
    lms_1931 = lx.xyz_to_lms(xyz_sample, cieobs="1931_2")
    lms_1964 = lx.xyz_to_lms(xyz_sample, cieobs="1964_10")
    srgb_sample = lx.xyz_to_srgb(np.array([[20.0, 21.0, 22.0]]))
    print(f"xyz_to_yxy={fmt_vec(yxy_sample.ravel())}")
    print(f"yxy_to_xyz={fmt_vec(lx.Yxy_to_xyz(yxy_sample).ravel())}")
    print(f"xyz_to_yuv={fmt_vec(yuv_sample.ravel())}")
    print(f"yuv_to_xyz={fmt_vec(lx.Yuv_to_xyz(yuv_sample).ravel())}")
    print(f"xyz_to_lab={fmt_vec(lab_sample.ravel())}")
    print(f"lab_to_xyz={fmt_vec(lx.lab_to_xyz(lab_sample, xyzw=white_sample).ravel())}")
    print(f"xyz_to_luv={fmt_vec(luv_sample.ravel())}")
    print(f"luv_to_xyz={fmt_vec(lx.luv_to_xyz(luv_sample, xyzw=white_sample).ravel())}")
    print(f"xyz_to_lms_1931={fmt_vec(lms_1931.ravel())}")
    print(f"lms_to_xyz_1931={fmt_vec(lx.lms_to_xyz(lms_1931, cieobs='1931_2').ravel())}")
    print(f"xyz_to_lms_1964={fmt_vec(lms_1964.ravel())}")
    print(f"xyz_to_srgb={fmt_vec(srgb_sample.ravel())}")
    print(
        f"srgb_to_xyz={fmt_vec(lx.srgb_to_xyz(np.array([[64.0, 128.0, 192.0]])).ravel())}"
    )
    print(
        "xyzbar_1931_interp="
        + fmt_vec(
            lx.xyzbar(
                cieobs="1931_2",
                wl_new=np.array([554.5, 555.0, 555.5, 556.0]),
                kind="linear",
                extrap_kind="linear",
            ).ravel()
        )
    )
    print(
        "vlbar_1931_interp="
        + fmt_vec(
            lx.vlbar(
                cieobs="1931_2",
                wl_new=np.array([554.5, 555.0, 555.5, 556.0]),
                kind="linear",
                extrap_kind="linear",
            ).ravel()
        )
    )
    print(
        "spd_interp="
        + fmt_vec(
            lx.cie_interp(
                np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
                np.array([395.0, 405.0, 420.0, 425.0]),
                kind="linear",
                extrap_kind="linear",
            ).ravel()
        )
    )
    print(
        "normalize_max="
        + fmt_vec(
            lx.spd_normalize(
                np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
                norm_type="max",
                norm_f=2.0,
            ).ravel()
        )
    )
    print(
        "normalize_area="
        + fmt_vec(
            lx.spd_normalize(
                np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
                norm_type="area",
                norm_f=1.0,
            ).ravel()
        )
    )
    print(
        "normalize_lambda="
        + fmt_vec(
            lx.spd_normalize(
                np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
                norm_type="lambda",
                norm_f=410.0,
            ).ravel()
        )
    )
    print(
        "normalize_ru="
        + fmt_vec(
            lx.spd_normalize(
                np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]),
                norm_type="ru",
                norm_f=10.0,
            ).ravel()
        )
    )
    print(
        "normalize_pu="
        + fmt_vec(
            lx.spd_normalize(
                np.array([[555.0, 556.0], [1.0, 1.0]]),
                norm_type="pu",
                norm_f=1000.0,
                cieobs="1931_2",
            ).ravel()
        )
    )
    print(
        "normalize_qu="
        + fmt_vec(
            lx.spd_normalize(
                np.array([[500.0, 510.0], [1.0, 1.0]]),
                norm_type="qu",
                norm_f=1e18,
            ).ravel()
        )
    )
    print(
        "blackbody_relative_6500="
        + fmt_vec(lx.blackbody(6500.0, wl3=[360.0, 365.0, 1.0], relative=True).ravel())
    )
    print(
        "blackbody_absolute_6500_560="
        + fmt_vec(lx.blackbody(6500.0, wl3=[560.0, 560.0, 1.0], relative=False).ravel())
    )
    x_d, y_d = lx.daylightlocus(6500.0)
    print(f"daylightlocus_6500={fmt_vec([x_d, y_d])}")
    print(
        "daylightphase_6500="
        + fmt_vec(lx.daylightphase(6500.0, wl3=[360.0, 365.0, 1.0]).ravel())
    )
    print(
        "daylightphase_3500="
        + fmt_vec(lx.daylightphase(3500.0, wl3=[360.0, 365.0, 1.0]).ravel())
    )
    print(
        "cri_ref_3000_6500="
        + fmt_vec(lx.cri_ref([3000.0, 6500.0], wl3=[360.0, 365.0, 1.0]).ravel())
    )
    xyz_to_cct_sample = np.array([[100.0, 100.0, 100.0]])
    cct_sample, duv_sample = lx.xyz_to_cct(xyz_to_cct_sample, out="cct,duv")
    print(
        "xyz_to_cct_sample="
        + fmt_vec([np.ravel(cct_sample)[0], np.ravel(duv_sample)[0]])
    )
    print("cct_to_xyz_6500=" + fmt_vec(lx.cct_to_xyz(6500.0).ravel()))
    print("illuminant_A=" + fmt_vec(lx._CIE_ILLUMINANTS["A"][:, 0:6].ravel()))
    print("illuminant_D65=" + fmt_vec(lx._CIE_ILLUMINANTS["D65"][:, 0:6].ravel()))
    print("illuminant_F4=" + fmt_vec(lx._CIE_ILLUMINANTS["F4"][:, 0:6].ravel()))
    print("illuminant_LED_B1=" + fmt_vec(lx._CIE_ILLUMINANTS["LED_B1"][:, 0:6].ravel()))
    print(
        "illuminant_D50="
        + fmt_vec(lx.daylightphase(5000.0, wl3=[360.0, 365.0, 1.0], cct_is_nominal=True).ravel())
    )

    print(f"getwlr={fmt_vec(lx.getwlr([360, 365, 1]))}")
    print(f"getwld_equal_scalar={fmt_scalar(lx.getwld(np.array([400.0, 410.0, 420.0])))}")
    print(f"getwld_unequal={fmt_vec(lx.getwld(np.array([400.0, 410.0, 430.0])))}")
    print(
        "power_ru="
        + fmt_scalar(lx.spd_to_power(np.array([[400.0, 410.0, 420.0], [1.0, 2.0, 3.0]]), "ru")[0, 0])
    )
    print(
        "power_pu="
        + fmt_scalar(
            lx.spd_to_power(
                np.array([[555.0, 556.0], [1.0, 1.0]]),
                "pu",
                cieobs="1931_2",
            )[0, 0]
        )
    )
    print(
        "power_qu="
        + fmt_scalar(lx.spd_to_power(np.array([[500.0, 510.0], [1.0, 1.0]]), "qu")[0, 0])
    )
    print(
        "ler_1931="
        + fmt_scalar(lx.spd_to_ler(np.array([[555.0, 556.0], [1.0, 1.0]]), cieobs="1931_2")[0, 0])
    )
    print(
        "ler_1964="
        + fmt_scalar(lx.spd_to_ler(np.array([[555.0, 556.0], [1.0, 1.0]]), cieobs="1964_10")[0, 0])
    )
    print(
        "ler_many_1931="
        + fmt_vec(
            lx.spd_to_ler(np.array([[555.0, 556.0], [1.0, 1.0], [2.0, 2.0]]), cieobs="1931_2").ravel()
        )
    )
    print(
        "xyz_relative="
        + fmt_vec(lx.spd_to_xyz(np.array([[555.0, 556.0], [1.0, 1.0]]), cieobs="1931_2")[0])
    )
    print(
        "xyz_absolute="
        + fmt_vec(
            lx.spd_to_xyz(
                np.array([[555.0, 556.0], [1.0, 1.0]]),
                cieobs="1931_2",
                relative=False,
            )[0]
        )
    )
    print(
        "xyz_relative_1964_10="
        + fmt_vec(lx.spd_to_xyz(np.array([[555.0, 556.0], [1.0, 1.0]]), cieobs="1964_10")[0])
    )
    print(
        "xyz_relative_many="
        + fmt_vec(
            lx.spd_to_xyz(np.array([[555.0, 556.0], [1.0, 1.0], [2.0, 2.0]]), cieobs="1931_2").ravel()
        )
    )
    print(
        "xyz_absolute_many="
        + fmt_vec(
            lx.spd_to_xyz(
                np.array([[555.0, 556.0], [1.0, 1.0], [2.0, 2.0]]),
                cieobs="1931_2",
                relative=False,
            ).ravel()
        )
    )


if __name__ == "__main__":
    main()
