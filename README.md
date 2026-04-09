# lux-rs

Pure Rust lighting and color science library for spectral, photometric, and colorimetric workflows.

[![crates.io](https://img.shields.io/crates/v/lux-rs.svg)](https://crates.io/crates/lux-rs)
[![docs.rs](https://docs.rs/lux-rs/badge.svg)](https://docs.rs/lux-rs)
[![license](https://img.shields.io/crates/l/lux-rs.svg)](https://github.com/husion/lux-rs/blob/main/LICENSE)

## Quick Links

- [Overview](#overview)
- [Install](#install)
- [Quick Example](#quick-example)
- [What's New in 0.1.1](#whats-new-in-011)
- [API Shape Conventions](#api-shape-conventions)
- [Roadmap](#roadmap)

## Overview

`lux-rs` provides a native Rust API for core lighting and color science calculations without requiring Python at runtime.

Current modules and responsibilities:

- `spectrum`: wavelength grids, spacing helpers, interpolation, and spectrum normalization for single and batch SPD workflows.
- `illuminants`: standard illuminant registry, blackbody/daylight/reference source generation, and CCT-XYZ conversion helpers.
- `photometry`: SPD integration to tristimulus values plus radiometric/photometric/quantal power and luminous efficacy.
- `color`: observer access, XYZ/Lab/Luv/Yuv/Yxy/sRGB/LMS transforms, color-difference metrics, and chromatic adaptation utilities.
- `cam`: `CIECAM02`, `CAM16`, and CAM-UCS forward/inverse appearance-model pipelines with viewing-condition helpers.
- `cri`: color rendering metrics including `CIE Ra`, `CIE Rf / Rg`, and structured `TM-30` results for single and batch spectra.
- `spectral_mismatch`: detector spectral mismatch metrics (`f1-prime`) and correction-factor computation utilities.
- `indvcmf`: deterministic individual-observer CMF construction (Asano-style slice), including LMS-to-XYZ conversion support.
- `error`: shared `LuxError` / `LuxResult` types used across modules.

## What's New in 0.1.1

- Added individual observer CMF support via the new `indvcmf` module (`individual_observer_cmf`, LMS-to-XYZ helpers, and defaults).
- Added detector spectral mismatch utilities via the new `spectral_mismatch` module (`f1′` and correction factor helpers).
- Expanded regression and parity coverage with new API-level tests for color, illuminants, photometry, spectral mismatch, spectra, and individual observer paths.
- Tightened API consistency and warning-free quality across modules (public API alignment and clippy cleanup).

## Design Goals

- spectral foundations: wavelength grids, spacing helpers, interpolation, normalization, and unified single/batch `Spectrum` workflows
- observers and photometry: embedded standard observers, tristimulus integration, radiometric / photometric / quantal power, and mesopic support
- illuminants and reference sources: blackbody, daylight family, CRI reference sources, and a registry for common CIE illuminants and LED series
- color kernels: CCT, common XYZ-derived transforms, color difference, and chromatic adaptation including viewing-condition and compiled-adapter workflows
- appearance models: first-pass `CIECAM02`, `CAM16`, `CAM02-UCS`, and `CAM16-UCS` forward / inverse paths plus wrapper APIs on top of the color data models
- color quality metrics: `CIE Ra`, `CIE Rf / Rg`, and structured `TM-30` result objects for single and batch spectral workflows
- advanced detector utilities: first-pass spectral mismatch (`f1′`) and correction-factor workflows on top of `Spectrum`
- advanced observer utilities: first-pass `indvcmf` deterministic LMS / XYZ CMF generation for one observer profile

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

For color calculations, use `Tristimulus` as the primary batch API (`1` row represents a single XYZ-like sample).
For spectral calculations, use `Spectrum` as the primary API (`1` row represents a single SPD).

## API Shape Conventions

Phase 1 convergence keeps the current numerical behavior, but makes the intended public API shape explicit:

- use `Spectrum` for both one SPD and aligned multi-SPD workflows; represent a single SPD as a one-row batch
- use `Tristimulus` for aligned XYZ-like color workflows; represent a single item as a one-row batch
- keep scalar and batch paths numerically aligned, for example free `spd_to_*` helpers and row-wise `Spectrum::spd_to_*` batch methods
- prefer constructors that keep row alignment explicit (`Spectrum::new(...)` and `Tristimulus::new(...)`), with `Tristimulus::from_single(...)` as convenience for one item
- fixed-size leaf kernels may still use raw `[f64; 3]` values internally, but public wrappers should prefer `Tristimulus` when they represent semantic color results
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
