# lux

`lux` is a standalone Rust library for lighting and color science calculations.

This crate is being built as a pure Rust rewrite path for LuxPy's computational
core. It does not call Python code at runtime.

Current scope:

- wavelength grid generation compatible with LuxPy `getwlr`
- wavelength spacing calculation compatible with LuxPy `getwld`
- `Spectrum` data model
- `SpectralMatrix` batch spectral data model
- linear interpolation with linear extrapolation
- spectrum normalization modes:
  - `max`
  - `area`
  - `lambda`
  - radiometric / photometric / quantal targets
- SPD power integration for:
  - radiometric units
  - photometric units using embedded standard observers
  - quantal units
- embedded CIE standard observers:
  - `1931_2`
  - `1964_10`
- CIE 191:2010 mesopic support:
  - `get_cie_mesopic_adaptation`
  - `vlbar_cie_mesopic`
- reference source models:
  - `blackbody`
  - `daylightlocus`
  - `daylightphase`
  - `cri_ref` (default CIE Ra path)
- standard illuminant registry:
  - `standard_illuminant(name, wl_grid)`
  - `A`
  - `D50/D55/D65/D75`
  - `F1..F12`
  - CIE LED series
- CCT / Duv:
  - `xyz_to_cct`
  - `cct_to_xyz`
- color transforms:
  - `XYZ <-> Yxy`
  - `XYZ <-> Yuv`
  - `XYZ <-> Lab` (explicit white point)
  - `XYZ <-> Luv` (explicit white point)
  - `XYZ <-> LMS` (observer/default matrix or explicit matrix)
  - `XYZ <-> sRGB`
- Python parity baselines for the current P0 core plus implemented color transforms

Short example:

```rust
use lux::{spd_to_power, Observer, PowerType, Spectrum};

let observer = Observer::Cie1931_2.standard()?;
let spectrum = Spectrum::new(vec![555.0, 556.0], vec![1.0, 1.0])?;
let lumens = spd_to_power(&spectrum, PowerType::Photometric, Some(&observer))?;
```

Planned next:

- reference illuminant mixing for non-default CRI / TM-30 paths
- wider standard illuminant coverage and alias cleanup:
  - additional CIE / IEC fixed datasets
  - alias normalization and lookup ergonomics
- CLI and optional Python bindings as thin wrappers on top of Rust
