# lux

Pure Rust lighting and color science library, built as a rewrite path for the computational core of [`luxpy`](https://github.com/ksmet1977/luxpy).

## Overview

`lux` focuses on a native Rust implementation of core lighting and color science workflows that are widely used in `luxpy`, without calling Python at runtime.

Current repository status, aligned with [`TODO_REFACTOR.md`](./TODO_REFACTOR.md):

- `P0` base spectral kernel: completed
- `P1` reference-source and CCT path: completed
- `P1.5` first standard-illuminant registry: completed
- `P2` status: color transforms, `deltaE`, and CAT utilities completed
- `P3` status: first CAM / CAM-UCS forward and inverse paths completed
- `P3` status: first CRI path completed for `CIE Ra`, `CIE Rf / Rg`, and `TM-30` result objects
- next priority: move into photobiological metrics

At the moment the crate already covers:

- spectral foundations: wavelength grids, spacing helpers, interpolation, normalization, and single/batch spectrum models
- observers and photometry: embedded standard observers, tristimulus integration, radiometric / photometric / quantal power, and mesopic support
- illuminants and reference sources: blackbody, daylight family, CRI reference sources, and a registry for common CIE illuminants and LED series
- color kernels: CCT, common XYZ-derived transforms, color difference, and chromatic adaptation including viewing-condition and compiled-adapter workflows
- appearance models: first-pass `CIECAM02`, `CAM16`, `CAM02-UCS`, and `CAM16-UCS` forward / inverse paths plus wrapper APIs on top of the color data models
- color quality metrics: `CIE Ra`, `CIE Rf / Rg`, and structured `TM-30` result objects for single and batch spectral workflows

## Why This Repo Exists

This repository is not a Rust binding around Python. It is a direct Rust implementation route for the parts of `luxpy` that matter most to numerical core workflows:

- predictable native deployment
- easier integration into Rust systems
- clearer data ownership and API design
- parity testing against an existing scientific reference implementation

A local `luxpy/` reference directory may be present in the development workspace for parity work, while the upstream project lives at <https://github.com/ksmet1977/luxpy>.

## Verification Status

The current crate state is backed by:

- `130` Rust unit tests
- `1` Python parity integration test
- local parity baselines against `luxpy` for:
  - spectral grid helpers
  - interpolation and normalization
  - observer / CMF access
  - `spd_to_power`, `spd_to_ler`, `spd_to_xyz`
  - `blackbody`, `daylightphase`, `cri_ref`
  - `xyz_to_cct`, `cct_to_xyz`
  - first-batch standard illuminants
  - one-step and two-step `CAT`
  - CAT adaptation degree and viewing-condition paths
  - CAM viewing-condition kernels for `CIECAM02` and `CAM16`
  - first CAM forward correlates for `CIECAM02` and `CAM16`
  - CAM-UCS forward and inverse paths for `CAM02-UCS` and `CAM16-UCS`
  - `CIE Ra`
  - `CIE Rf / Rg`
  - `TM-30` result objects

Parity checks currently run through [`tests/python_parity.rs`](./tests/python_parity.rs) and [`tests/python_ref/current_baselines.py`](./tests/python_ref/current_baselines.py).

## Quick Example

```rust
use lux::{CatTransform, Tristimulus, TristimulusSet};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sample = Tristimulus::new([19.01, 20.0, 21.78]);
    let adapted = sample.cat_apply(
        [95.047, 100.0, 108.883],
        [109.85, 100.0, 35.585],
        CatTransform::Bradford,
        1.0,
    )?;

    let batch = TristimulusSet::new(vec![[0.25, 0.5, 0.25], [0.2, 0.3, 0.4]]);
    let lab = batch.xyz_to_lab([0.5, 0.5, 0.5]);

    println!("{:?}", adapted.values());
    println!("{:?}", lab.values());
    Ok(())
}
```

For color calculations, use `Tristimulus` and `TristimulusSet` as the primary API.
For spectral calculations, use `Spectrum` and `SpectralMatrix` as the primary API.

## Roadmap

Near-term work, following [`TODO_REFACTOR.md`](./TODO_REFACTOR.md):

1. finish illuminant naming cleanup and alias normalization
2. start `photobiochem` base metrics
3. return to broader observer coverage and remaining result-layer polish

Longer-term items such as broader CAM families, TM-30 graphics, photobiological metrics, individual observers, and hyperspectral tooling remain intentionally deferred until the current core stays stable.

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
