"""Tests for the linearized L2 halo-orbit seed."""

import pytest

from halo_insertion import constants as c
from halo_insertion.equilibria import find_L2
from halo_insertion.halo_seed import linear_halo_seed

MU = c.MU


def test_seed_symmetric_ic_form():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    x0, y0, z0, xdot0, ydot0, zdot0 = seed.state0
    assert y0 == 0.0
    assert xdot0 == 0.0
    assert zdot0 == 0.0
    assert z0 == pytest.approx(0.03)
    assert x0 != l2.x  # offset from L2
    assert ydot0 != 0.0


def test_seed_near_l2():
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    x0 = seed.state0[0]
    # The seed's in-plane offset should be small relative to the
    # Earth-Moon distance scale (a few percent of 1 DU), i.e. genuinely
    # "near" L2, not some unrelated part of the system.
    assert abs(x0 - l2.x) < 0.1


def test_seed_z_sign_convention():
    l2 = find_L2(MU)
    seed_pos = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    seed_neg = linear_halo_seed(l2.x, MU, 0.03, z_sign=-1.0)
    assert seed_pos.state0[2] == pytest.approx(0.03)
    assert seed_neg.state0[2] == pytest.approx(-0.03)


def test_seed_rejects_nonpositive_amplitude():
    l2 = find_L2(MU)
    with pytest.raises(ValueError):
        linear_halo_seed(l2.x, MU, 0.0, z_sign=1.0)
    with pytest.raises(ValueError):
        linear_halo_seed(l2.x, MU, -0.01, z_sign=1.0)


def test_seed_frequencies_are_close_but_not_forced_equal():
    # Documents the key M3 finding: the linear in-plane and out-of-plane
    # frequencies at L2 differ by a few percent (this is exactly why a
    # linear seed needs nonlinear differential correction to close into
    # a genuine periodic 3D orbit).
    l2 = find_L2(MU)
    seed = linear_halo_seed(l2.x, MU, 0.03, z_sign=1.0)
    ratio = seed.wz / seed.wxy
    assert 0.9 < ratio < 1.0  # close to 1, but not equal
