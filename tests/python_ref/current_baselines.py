"""Top-level smoke script that preserves the existing parity invocation contract."""

from itertools import chain
from pathlib import Path

from adaptation_cam_baselines import iter_adaptation_cam_baselines
from baseline_helpers import print_pairs
from color_baselines import iter_color_baselines, iter_display_color_baselines
from cri_illuminant_baselines import iter_cri_illuminant_baselines
from observer_baselines import iter_observer_baselines
from spectral_baselines import iter_spectral_baselines


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    print_pairs(
        chain(
            iter_observer_baselines(),
            iter_color_baselines(),
            iter_adaptation_cam_baselines(),
            iter_display_color_baselines(),
            iter_spectral_baselines(),
            iter_cri_illuminant_baselines(root),
        )
    )


if __name__ == "__main__":
    main()
