"""Shared formatting and emission helpers for Python parity baselines."""

def fmt_scalar(value: float) -> str:
    return repr(float(value))


def fmt_vec(values) -> str:
    return ",".join(repr(float(value)) for value in values)


def print_pairs(pairs) -> None:
    for key, value in pairs:
        print(f"{key}={value}")
