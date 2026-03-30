# lux

Pure Rust lighting and color science library, built as a rewrite path for the computational core of [`luxpy`](https://github.com/ksmet1977/luxpy).

## Overview

`lux` focuses on a native Rust implementation of core lighting and color science workflows that are widely used in `luxpy`, without calling Python at runtime.

Current repository status, aligned with [`TODO_REFACTOR.md`](./TODO_REFACTOR.md):

- `P0` base spectral kernel: completed
- `P1` reference-source and CCT path: completed
- `P1.5` first standard-illuminant registry: completed
- next priority: `CAT / deltaE` (`P2`)

At the moment the crate already covers:

- wavelength grid generation compatible with `luxpy.getwlr`
- wavelength spacing calculation compatible with `luxpy.getwld`
- `Spectrum` and `SpectralMatrix` data models
- linear interpolation with linear extrapolation
- spectrum normalization:
  - `max`
  - `area`
  - `lambda`
  - radiometric / photometric / quantal targets
- embedded standard observers:
  - `1931_2`
  - `1964_10`
- photometric and tristimulus pipeline:
  - `spd_to_power`
  - `spd_to_ler`
  - `spd_to_xyz`
- CIE 191:2010 mesopic support:
  - `get_cie_mesopic_adaptation`
  - `vlbar_cie_mesopic`
- reference-source models:
  - `blackbody`
  - `daylightlocus`
  - `daylightphase`
  - `cri_ref`
- standard illuminant registry:
  - `standard_illuminant(name, wl_grid)`
  - `A`
  - `D50 / D55 / D65 / D75`
  - `F1..F12`
  - CIE LED series
- CCT / Duv:
  - `xyz_to_cct`
  - `cct_to_xyz`
- color transforms:
  - `XYZ <-> Yxy`
  - `XYZ <-> Yuv`
  - `XYZ <-> Lab`
  - `XYZ <-> Luv`
  - `XYZ <-> LMS`
  - `XYZ <-> sRGB`

## Why This Repo Exists

This repository is not a Rust binding around Python. It is a direct Rust implementation route for the parts of `luxpy` that matter most to numerical core workflows:

- predictable native deployment
- easier integration into Rust systems
- clearer data ownership and API design
- parity testing against an existing scientific reference implementation

A local `luxpy/` reference directory may be present in the development workspace for parity work, while the upstream project lives at <https://github.com/ksmet1977/luxpy>.

## Verification Status

The current crate state is backed by:

- `56` Rust unit tests
- `1` Python parity integration test
- local parity baselines against `luxpy` for:
  - spectral grid helpers
  - interpolation and normalization
  - observer / CMF access
  - `spd_to_power`, `spd_to_ler`, `spd_to_xyz`
  - `blackbody`, `daylightphase`, `cri_ref`
  - `xyz_to_cct`, `cct_to_xyz`
  - first-batch standard illuminants

Parity checks currently run through [`tests/python_parity.rs`](./tests/python_parity.rs) and [`tests/python_ref/current_baselines.py`](./tests/python_ref/current_baselines.py).

## Quick Example

```rust
use lux::{spd_to_power, Observer, PowerType, Spectrum};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::Cie1931_2.standard()?;
    let spectrum = Spectrum::new(vec![555.0, 556.0], vec![1.0, 1.0])?;
    let lumens = spd_to_power(&spectrum, PowerType::Photometric, Some(&observer))?;

    println!("{lumens}");
    Ok(())
}
```

## Public API Snapshot

The root crate currently re-exports these main entry points:

- spectrum:
  - `getwlr`
  - `getwld`
  - `Spectrum`
  - `SpectralMatrix`
  - `SpectrumNormalization`
  - `WavelengthGrid`
- photometry:
  - `spd_to_power`
  - `spd_to_ler`
  - `spd_to_ler_many`
  - `spd_to_xyz`
  - `spd_to_xyz_many`
- illuminants:
  - `blackbody`
  - `daylightlocus`
  - `daylightphase`
  - `cri_ref`
  - `standard_illuminant`
  - `standard_illuminant_names`
  - `xyz_to_cct`
  - `cct_to_xyz`
- color:
  - `xyz_to_yxy`, `yxy_to_xyz`
  - `xyz_to_yuv`, `yuv_to_xyz`
  - `xyz_to_lab`, `lab_to_xyz`
  - `xyz_to_luv`, `luv_to_xyz`
  - `xyz_to_lms`, `lms_to_xyz`
  - `xyz_to_srgb`, `srgb_to_xyz`
  - `vlbar_cie_mesopic`
  - `get_cie_mesopic_adaptation`

## Roadmap

Near-term work, following [`TODO_REFACTOR.md`](./TODO_REFACTOR.md):

1. finish illuminant naming cleanup and alias normalization
2. implement `CAT` main paths
3. implement `deltaE`
4. expand observer sets and interpolation semantics toward broader `luxpy` coverage

Longer-term items such as CAM, CRI/TM-30, photobiological metrics, individual observers, and hyperspectral tooling remain intentionally deferred until the core kernel stays stable.

## Relationship To LuxPy

[`luxpy`](https://github.com/ksmet1977/luxpy) is a comprehensive Python toolbox for lighting and color science. This repository currently targets the computational subset that forms a solid native core for future Rust-first use.

Scope difference, in short:

- `luxpy`: broad toolbox including CAM, CAT, CRI/TM-30, photobiology, hyperspectral imaging, instrument/toolbox integrations, and more
- `lux`: currently focused on the spectral kernel, observers, integration, reference illuminants, CCT, and common color transforms

That means `lux` is currently narrower than `luxpy`, but already suitable for a meaningful subset of core calculations and for parity-driven incremental expansion.

## Citing LuxPy

If this repository or its design work benefits from `luxpy`, please cite the original `luxpy` project and tutorial paper.

Recommended citation from the upstream [`luxpy` README](https://github.com/ksmet1977/luxpy/blob/master/README.md):

> Smet, K. A. G. (2020). Tutorial: The LuxPy Python Toolbox for Lighting and Color Science. LEUKOS, 1-23. https://doi.org/10.1080/15502724.2018.1518717

Useful upstream references:

- LuxPy repository: <https://github.com/ksmet1977/luxpy>
- LuxPy tutorial paper: <https://www.tandfonline.com/doi/full/10.1080/15502724.2018.1518717>
- LuxPy Zenodo DOI: <https://doi.org/10.5281/zenodo.1298963>

## License

This crate is licensed under `GPL-3.0-only`. See [`Cargo.toml`](./Cargo.toml) and the repository license terms for details.
