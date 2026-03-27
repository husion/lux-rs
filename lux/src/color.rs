use crate::error::{LuxError, LuxResult};
use crate::spectrum::Spectrum;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Observer {
    Cie1931_2,
    Cie1964_10,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TristimulusObserver {
    pub wavelengths: Vec<f64>,
    pub x_bar: Vec<f64>,
    pub y_bar: Vec<f64>,
    pub z_bar: Vec<f64>,
    pub k: f64,
}

impl Observer {
    pub fn standard(self) -> LuxResult<TristimulusObserver> {
        match self {
            Self::Cie1931_2 => TristimulusObserver::from_csv(
                include_str!("../data/cmfs/ciexyz_1931_2.dat"),
                683.002,
            ),
            Self::Cie1964_10 => TristimulusObserver::from_csv(
                include_str!("../data/cmfs/ciexyz_1964_10.dat"),
                683.599,
            ),
        }
    }
}

impl TristimulusObserver {
    pub fn from_csv(csv: &str, k: f64) -> LuxResult<Self> {
        let mut wavelengths = Vec::new();
        let mut x_bar = Vec::new();
        let mut y_bar = Vec::new();
        let mut z_bar = Vec::new();

        for line in csv.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            let mut parts = trimmed.split(',');
            let wl = parts
                .next()
                .ok_or(LuxError::ParseError("missing wavelength"))?
                .parse::<f64>()
                .map_err(|_| LuxError::ParseError("invalid wavelength"))?;
            let x = parts
                .next()
                .ok_or(LuxError::ParseError("missing x_bar"))?
                .parse::<f64>()
                .map_err(|_| LuxError::ParseError("invalid x_bar"))?;
            let y = parts
                .next()
                .ok_or(LuxError::ParseError("missing y_bar"))?
                .parse::<f64>()
                .map_err(|_| LuxError::ParseError("invalid y_bar"))?;
            let z = parts
                .next()
                .ok_or(LuxError::ParseError("missing z_bar"))?
                .parse::<f64>()
                .map_err(|_| LuxError::ParseError("invalid z_bar"))?;

            wavelengths.push(wl);
            x_bar.push(x);
            y_bar.push(y);
            z_bar.push(z);
        }

        if wavelengths.is_empty() {
            return Err(LuxError::EmptyInput);
        }

        Ok(Self {
            wavelengths,
            x_bar,
            y_bar,
            z_bar,
            k,
        })
    }

    pub fn vl_spectrum(&self) -> LuxResult<Spectrum> {
        Spectrum::new(self.wavelengths.clone(), self.y_bar.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::Observer;

    #[test]
    fn loads_standard_observer() {
        let observer = Observer::Cie1931_2.standard().unwrap();
        assert_eq!(observer.wavelengths.first().copied(), Some(360.0));
        assert_eq!(observer.wavelengths.last().copied(), Some(830.0));
        assert_eq!(observer.wavelengths.len(), 471);
    }
}
