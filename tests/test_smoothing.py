from __future__ import annotations

import numpy as np

import eic_smoothing_simulator as sim


def test_ppm_window() -> None:
    low, high = sim.ppm_window(500.0, 10.0)
    assert np.isclose(low, 499.995)
    assert np.isclose(high, 500.005)


def test_constant_signal_is_preserved() -> None:
    x = np.arange(21, dtype=float)
    y = np.full(x.size, 7.5)
    for method in ("SMA", "LWMA", "TBLWMA", "MSDIAL_SG", "BF", "LOESS"):
        output, _ = sim.apply_smoothing(
            y,
            method,
            msdial_level=2,
            mzmine_sg_width=5,
            mzmine_loess_width=5,
            x=x,
        )
        assert output.shape == y.shape
        assert np.all(np.isfinite(output))
        assert np.allclose(output, y, atol=1e-9)


def test_reference_impulse_responses() -> None:
    x = np.arange(21, dtype=float)
    y = np.zeros(x.size)
    center = x.size // 2
    y[center] = 1.0

    lwma, _ = sim.apply_smoothing(y, "LWMA", msdial_level=2, x=x)
    tblwma, _ = sim.apply_smoothing(y, "TBLWMA", msdial_level=2, x=x)
    binomial, _ = sim.apply_smoothing(y, "BF", msdial_level=2, x=x)

    expected_linear = np.asarray([1, 2, 3, 2, 1], dtype=float) / 9.0
    expected_binomial = np.asarray([1, 4, 6, 4, 1], dtype=float) / 16.0
    assert np.allclose(lwma[center - 2 : center + 3], expected_linear)
    assert np.allclose(tblwma[center - 2 : center + 3], expected_linear)
    assert np.allclose(binomial[center - 2 : center + 3], expected_binomial)


def test_msdial_sg_level_two_coefficients() -> None:
    expected = np.asarray([-3, 12, 17, 12, -3], dtype=float) / 35.0
    assert np.allclose(sim._msdial_sg_coefficients(2), expected)


def test_mzmine_zero_positions_remain_zero() -> None:
    x = np.arange(21, dtype=float)
    y = np.full(x.size, 10.0)
    center = x.size // 2
    y[center] = 0.0

    sg, _ = sim.apply_smoothing(y, "MZMINE_SG", mzmine_sg_width=5, x=x)
    loess, _ = sim.apply_smoothing(y, "LOESS", mzmine_loess_width=5, x=x)

    assert sg[center] == 0.0
    assert loess[center] == 0.0


def test_loess_preserves_linear_signal() -> None:
    x = np.arange(21, dtype=float)
    y = 2.0 * x + 1.0
    output, _ = sim.apply_smoothing(y, "LOESS", mzmine_loess_width=5, x=x)
    assert np.allclose(output, y, atol=1e-9)
