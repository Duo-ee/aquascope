"""Flood it (#658): the local-inertial solver against still water, the volume budget, Manning's normal depth and
Ritter's dam break; the PNG and Terrarium decoding; the channel; the agreement measure."""

from __future__ import annotations

import math
import struct
import zlib

import numpy as np
import pytest

from aquascope import floodsim as fs


def test_still_water_stays_still():
    """Well-balanced: a lake at rest over a bumpy bed, with dry islands, does not move (closed or open edges)."""
    rng = np.random.default_rng(1)
    z = rng.uniform(0, 3, size=(40, 50))
    z[10:14, 10:14] = 5.0          # an island above the water
    level = 2.5
    h0 = np.maximum(level - z, 0.0)
    # Raise the edges above the lake so the open-edge case has nothing to drain.
    z[0, :] = z[-1, :] = z[:, 0] = z[:, -1] = 6.0
    h0[0, :] = h0[-1, :] = h0[:, 0] = h0[:, -1] = 0.0
    for closed in (True, False):
        r = fs.simulate(z, 10.0, duration_s=600.0, h0=h0, n=0.03, closed=closed)
        assert r["steps"] > 50
        assert np.abs(r["h"] - h0).max() < 1e-12
        assert r["volume"]["outflow_m3"] == 0.0


def test_mass_conservation_closed_and_open():
    rng = np.random.default_rng(2)
    z = np.cumsum(rng.normal(0, 0.2, size=(60, 60)), axis=0) + rng.uniform(0, 1, size=(60, 60))
    # Closed: everything poured in stays in.
    r = fs.simulate(z, 20.0, duration_s=3600.0, n=0.05, sources=[(30, 30, 50.0)], closed=True)
    v = r["volume"]
    assert v["inflow_m3"] == pytest.approx(50.0 * 3600.0)
    assert v["outflow_m3"] == 0.0
    assert abs(v["error_rel"]) < 1e-12
    assert r["h"].min() >= 0.0
    # Open: what leaves at the edges is counted, and the budget still closes.
    r = fs.simulate(z, 20.0, duration_s=6 * 3600.0, n=0.05, sources=[(30, 30, 50.0)])
    v = r["volume"]
    assert v["outflow_m3"] > 0.0
    assert abs(v["error_rel"]) < 1e-12


def test_normal_depth_on_a_slope():
    """Steady uniform flow down a plane reaches Manning's normal depth, h = (n q / sqrt(S))^(3/5)."""
    s, n, q, dx, nx = 1e-3, 0.05, 2.0, 10.0, 400
    z = -(np.arange(nx) * dx * s)[None, :]
    r = fs.simulate(z, dx, duration_s=6 * 3600.0, n=n, sources=[(0, 5, q * dx)], closed="nsw")
    expected = (n * q / math.sqrt(s)) ** 0.6
    mid = r["h"][0, 150:250]
    assert np.abs(mid / expected - 1).max() < 0.02


def test_ritter_dam_break_within_the_stated_tolerance():
    """Ritter (1892): dry, flat, frictionless bed. The local-inertial scheme has no advection term, so it cannot match
    the supercritical front (the literature says so); the stated tolerance is a relative L1 depth error under 0.20
    at 30 s, which holds as the grid is refined (0.177 to 0.185 for dx = 4 m to 0.5 m), with the volume exact."""
    h0, dx, t = 1.0, 1.0, 30.0
    nx = 600
    x = (np.arange(nx) + 0.5) * dx - 200.0
    z = np.zeros((1, nx))
    h = np.where(x < 0, h0, 0.0)[None, :]
    r = fs.simulate(z, dx, duration_s=t, h0=h, n=0.0, closed=True)
    hm, ha = r["h"][0], fs.ritter(x, r["t_s"], h0)
    l1 = np.abs(hm - ha).sum() / ha.sum()
    assert l1 < 0.20
    # Upstream of the rarefaction (beyond a 20 m numerical spread) the water has not moved yet.
    c0 = math.sqrt(fs.G * h0)
    calm = x < -c0 * t - 20
    assert np.abs(hm[calm] - h0).max() < 1e-6
    # The front lags Ritter's 2 c0 t (the known under-estimate of the scheme) but has moved.
    front = x[hm > 0.01].max()
    assert 0.25 * 2 * c0 * t < front < 2 * c0 * t
    assert abs(r["volume"]["error_rel"]) < 1e-12


def test_ritter_formula():
    c0 = math.sqrt(fs.G * 2.0)
    assert fs.ritter([-100.0], 1.0, 2.0)[0] == 2.0
    assert fs.ritter([0.0], 10.0, 2.0)[0] == pytest.approx(4 / 9 * 2.0)
    assert fs.ritter([2 * c0 * 10 + 1], 10.0, 2.0)[0] == 0.0


def _png(rows: list[list[tuple[int, int, int]]], filters: list[int]) -> bytes:
    """A tiny RGB PNG with the given per-row filters (0 none, 1 sub, 2 up)."""
    h, w = len(rows), len(rows[0])
    raw = b""
    prev = bytes(3 * w)
    for row, f in zip(rows, filters):
        line = bytes(v for px in row for v in px)
        if f == 1:
            enc = bytes((line[i] - (line[i - 3] if i >= 3 else 0)) & 255 for i in range(len(line)))
        elif f == 2:
            enc = bytes((line[i] - prev[i]) & 255 for i in range(len(line)))
        else:
            enc = line
        raw += bytes([f]) + enc
        prev = line

    def chunk(kind: bytes, body: bytes) -> bytes:
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body))

    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


def test_png_and_terrarium_decoding():
    rows = [[(128, 0, 0), (128, 10, 128)], [(127, 255, 0), (129, 2, 64)], [(128, 100, 200), (128, 1, 1)]]
    rgb = fs.decode_png_rgb(_png(rows, [0, 1, 2]))
    assert rgb.tolist() == [[list(p) for p in r] for r in rows]
    z = fs.terrarium_elevation(rgb)
    assert z[0, 0] == 0.0                      # 128*256 - 32768
    assert z[0, 1] == pytest.approx(10.5)
    assert z[1, 0] == pytest.approx(-1.0)
    assert z[1, 1] == pytest.approx(258.25)


def test_channel_depth_and_burning():
    # d = (n Q / (W sqrt S))^(3/5)
    assert fs.channel_depth(100.0, 25.0, 1e-3) == pytest.approx((0.035 * 100 / (25 * math.sqrt(1e-3))) ** 0.6)
    # The slope is kept within the range.
    assert fs.channel_depth(100.0, 25.0, 0.0) == fs.channel_depth(100.0, 25.0, fs.SLOPE_RANGE[0])
    z = np.array([[10.0, 9.0, 9.5, 8.0, 8.2]])
    out = fs.burn_channel(z, [(0, i) for i in range(5)], 2.0)
    # The bank is a running minimum, so the bed never rises downstream.
    assert out[0].tolist() == [8.0, 7.0, 7.0, 6.0, 6.0]


def test_pixel_size_and_cells():
    assert fs.pixel_size_m(0.0, 0) == pytest.approx(fs.EARTH_CIRCUMFERENCE_M / 256)
    grid = {"zoom": 12, "x0": 0, "y0": 0, "size": 512}
    x, y = fs._world_px(-2.2, 52.06, 12)
    grid["x0"], grid["y0"] = int(x) - 256, int(y) - 256
    assert fs.cell_of(grid, -2.2, 52.06) == (256, 256)
    lon, lat = fs.cell_lonlat(grid, 256, 256)
    assert abs(float(lon) + 2.2) < 3e-4 and abs(float(lat) - 52.06) < 3e-4


def test_compare_depths():
    m = np.array([[0.0, 1.0, 2.0], [0.5, 0.0, 0.0]])
    ref = np.array([[0.0, 2.0, 0.0], [1.0, 1.0, np.nan]])
    c = fs.compare_depths(m, ref)
    assert (c["wet_both"], c["model_only"], c["reference_only"]) == (2, 1, 1)
    assert c["csi"] == 0.5
    assert c["depth_diff_m"]["mean"] == pytest.approx(-0.75)


def test_matches_the_javascript_twin_fixture():
    """explorer/tests/floodsim-core.test.mjs runs the same grid in JS against this fixture; here it is regenerated
    and checked so the two cannot drift apart silently."""
    import json
    from pathlib import Path

    fx = json.loads((Path(__file__).parent / "fixtures" / "floodsim_twin.json").read_text())
    r = fs.simulate(np.array(fx["z"]), fx["dx"], duration_s=fx["duration_s"], n=np.array(fx["n"]),
                    sources=[tuple(fx["source"])])
    assert r["steps"] == fx["steps"]
    assert np.abs(r["h"] - np.array(fx["h"])).max() < 1e-9
