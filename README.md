# lux-rs

Pure Rust lighting and color science library for spectral, photometric, and colorimetric workflows.

## Overview

`lux-rs` provides a native Rust API for core lighting and color science calculations without requiring Python at runtime.

The crate currently includes:

- `P0` base spectral kernel: completed
- `P1` reference-source and CCT path: completed
- `P1.5` first standard-illuminant registry: completed
- `P2` status: color transforms, `deltaE`, and CAT utilities completed
- `P3` status: first CAM / CAM-UCS forward and inverse paths completed
- `P3` status: first CRI path completed for `CIE Ra`, `CIE Rf / Rg`, and `TM-30` result objects
- next priority: move into photobiological metrics

## Design Goals

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

## Verification

The crate is validated with Rust tests and parity checks against [`luxpy`](https://github.com/ksmet1977/luxpy). Covered paths include:

- spectral grid helpers
- interpolation and normalization
- observer and CMF access
- `spd_to_power`, `spd_to_ler`, `spd_to_xyz`
- `blackbody`, `daylightphase`, `cri_ref`
- `xyz_to_cct`, `cct_to_xyz`
- standard illuminants
- one-step and two-step `CAT`
- `CIECAM02`, `CAM16`, `CAM02-UCS`, and `CAM16-UCS`
- `CIE Ra`
- `CIE Rf / Rg`
- `TM-30` result objects

## Install

```bash
cargo add lux-rs
```

## Quick Example

```rust
use lux_rs::{spd_to_ler, spd_to_xyz, standard_illuminant, xyz_to_cct, Observer};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let observer = Observer::Cie1931_2.standard()?;
    let d65 = standard_illuminant("D65", None)?;

    let xyz = spd_to_xyz(&d65, &observer, true)?;
    let ler = spd_to_ler(&d65, &observer)?;
    let (cct, duv) = xyz_to_cct(xyz, Observer::Cie1931_2)?;

    println!("XYZ: {:?}", xyz);
    println!("LER: {:.3} lm/W", ler);
    println!("CCT: {:.1} K, Duv: {:.6}", cct, duv);
    Ok(())
}
```

For color calculations, use `Tristimulus` and `TristimulusSet` as the primary API.
For spectral calculations, use `Spectrum` and `SpectralMatrix` as the primary API.

## API Shape Conventions

Phase 1 convergence keeps the current numerical behavior, but makes the intended public API shape explicit:

- use `Spectrum` for one spectral power distribution and `SpectralMatrix` for aligned batch spectral workflows
- use `Tristimulus` for one XYZ-like color value and `TristimulusSet` for aligned batch color workflows
- prefer single-item public entry points first; add batch variants only when the same operation needs to run over aligned multi-row data
- when both forms exist, keep singular/plural naming aligned, for example `spd_to_*` for `Spectrum` and `spds_to_*` for `SpectralMatrix`
- source constructors that produce one SPD should return `Spectrum`; constructors that naturally produce multiple SPDs should return `SpectralMatrix`
- fixed-size leaf kernels may still use raw `[f64; 3]` values internally, but public wrappers should prefer `Tristimulus` / `TristimulusSet` when they represent semantic color results
- wrapper APIs should stay thin: the single-item and batch forms should share one core implementation rather than fork behavior

## Roadmap

Near-term work, following [`TODO_REFACTOR.md`](./TODO_REFACTOR.md):

1. finish illuminant naming cleanup and alias normalization
2. start `photobiochem` base metrics
3. return to broader observer coverage and remaining result-layer polish

Longer-term items such as broader CAM families, TM-30 graphics, photobiological metrics, individual observers, and hyperspectral tooling remain intentionally deferred until the current core stays stable.

## Relationship To LuxPy

[`luxpy`](https://github.com/ksmet1977/luxpy) is a comprehensive Python toolbox for lighting and color science. `lux-rs` is not a binding layer around Python; it is a native Rust implementation that draws on the same problem domain and uses LuxPy for parity-oriented validation during development.

Scope difference, in short:

- `luxpy`: broad toolbox including CAM, CAT, CRI/TM-30, photobiology, hyperspectral imaging, instrument/toolbox integrations, and more
- `lux-rs`: focused on spectral kernels, observers, integration, reference illuminants, photometry, CCT, color transforms, CAM, and CRI/TM-30 core workflows

That means `lux-rs` is narrower in scope than `luxpy`, while still being suitable for a meaningful subset of core numerical workflows.

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
