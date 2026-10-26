"""Flood it (#658): a 2D screening flood from a river reach, with the local-inertial scheme. Experimental (a spike).

The scheme is the local-inertial approximation of the shallow water equations used by LISFLOOD-FP (Bates,
Horritt and Fewtrell, 2010, J. Hydrol. 387, 33-45, doi:10.1016/j.jhydrol.2010.03.027), in the q-centred form of
de Almeida, Bates, Freer and Souvignet (2012, Water Resour. Res. 48, W05528, doi:10.1029/2011WR011570). On a
regular grid of cells of side ``dx``, depth ``h`` and bed ``z`` at cell centres and unit discharge ``q`` (m2/s) on
the faces, one step is:

* the flow depth at a face: ``hf = max(eta_l, eta_r) - max(z_l, z_r)`` with ``eta = z + h``;
* the face discharge (Manning friction, semi-implicit)::

      q_new = [theta q + (1 - theta)/2 (q_left + q_right) - g hf dt (eta_r - eta_l)/dx]
              / [1 + g dt n^2 |q| / hf^(7/3)]

  where ``q_left`` and ``q_right`` are the old discharges at the two neighbouring faces in the same direction
  (``theta = 1`` gives Bates et al. 2010); a face with ``hf`` at or below :data:`H_DRY` carries nothing;
* the depth: ``h_new = h + dt/dx (q_in - q_out)`` summed over the four faces, plus any inflow ``Q dt / dx^2``;
* the time step: ``dt = alpha dx / sqrt(g h_max)`` (the CFL condition of Bates et al. 2010, eq. 12).

The forms were checked against the equations quoted in Shaw et al. (2021, GMD 14, 3577, eqs. 18-19, the
LISFLOOD-FP 8.0 paper) and Cozzolino et al. (2019, Adv. Water Resour. 133, 103422, eqs. 46-48), which also give
the coefficients used here: ``alpha = 0.7`` (de Almeida et al. 2012) and ``theta = 0.7`` (de Almeida and Bates
2013).

Wet and dry: a face carries water only where ``hf > H_DRY``. Depth never goes negative: when a cell's outgoing
faces would take more than it holds in a step, they are scaled down so that it empties exactly (a flux limiter
that keeps the volume exact, rather than clipping negative depths, which would create water). The edges of the
grid are a free outfall: a ghost cell beyond the edge has the edge cell's bed and no water, and the edge faces only
let water out. The volume that leaves is counted, so the budget closes: final = initial + inflow - outflow.

What the scheme is not good at: supercritical flow. Without the advection term it under-estimates fast fronts
(de Almeida and Bates 2013; Cozzolino et al. 2019 show it fails for supercritical flow), so the Ritter dam break
below is a check with a stated tolerance, not an exact match. Floods on flood plains are slow, which is what the
scheme is for.

Terrain: the Terrarium tiles the Explorer's hillshade uses (``elevation-tiles-prod``, Mapzen/Tilezen joerd, a
mix of SRTM, GMTED, ETOPO1 and national lidar; decode ``h = R*256 + G + B/256 - 32768``). No bathymetry: the river
is burnt in (:func:`burn_channel`). This is a screening flood: terrain errors of metres, no levees unless the
terrain shows them, no rain, a constant inflow.

numpy and the standard library only (the PNG tiles are decoded with zlib), so it runs in Pyodide too.
"""

from __future__ import annotations

import math
import struct
import time
import zlib
from collections.abc import Callable
from typing import Any

import numpy as np

__all__ = [
    "ALPHA", "G", "H_DRY", "N_CHANNEL", "N_FLOODPLAIN", "THETA", "TERRARIUM", "burn_channel", "channel_depth",
    "compare_depths", "decode_png_rgb", "flood_reach", "jrc_on_grid", "pixel_size_m", "reach_grid", "ritter",
    "simulate", "terrain_grid", "terrarium_elevation",
]

G = 9.81
#: The CFL coefficient (de Almeida et al. 2012, as used by Cozzolino et al. 2019).
ALPHA = 0.7
#: The q-centred weight (de Almeida and Bates 2013, as used by Cozzolino et al. 2019); 1 is Bates et al. 2010.
THETA = 0.7
#: A face whose flow depth is at or below this (m) carries no water.
H_DRY = 1e-3
#: Manning's n for the burnt channel and the flood plain: screening defaults, not calibrated. Chow (1959) gives
#: about 0.025-0.04 for clean natural channels and 0.03-0.08 for pasture to light brush; JRC's maps are made with
#: their own values.
N_CHANNEL = 0.035
N_FLOODPLAIN = 0.06
#: The slope used to size the channel is kept in this range (the terrain's slope over a few km is noisy).
SLOPE_RANGE = (1e-4, 1e-2)

TERRARIUM = "https://elevation-tiles-prod.s3.amazonaws.com/terrarium/{z}/{x}/{y}.png"
TERRAIN_CREDIT = ("Terrain: Terrarium tiles (Mapzen/Tilezen joerd on AWS Open Data: SRTM, GMTED2010, ETOPO1, "
                  "national lidar and others; see the joerd attribution list)")
EARTH_CIRCUMFERENCE_M = 40075016.686


# ── the solver ──────────────────────────────────────────────────────────────


def _faces(z: np.ndarray, eta: np.ndarray, n: np.ndarray, q: np.ndarray, dt: float, dx: float, theta: float,
           axis: int) -> np.ndarray:
    """New discharges on every face along ``axis`` (1: x, between columns; 0: y, between rows), edges included.

    ``q`` has one more entry than the cells along ``axis``. Ghost cells beyond the edges have the edge cell's
    bed and no water, so the edge faces drain the grid; they are clipped to outflow only.
    """
    zl = np.concatenate([np.take(z, [0], axis=axis), z], axis=axis)
    zr = np.concatenate([z, np.take(z, [-1], axis=axis)], axis=axis)
    el = np.concatenate([np.take(z, [0], axis=axis), eta], axis=axis)
    er = np.concatenate([eta, np.take(z, [-1], axis=axis)], axis=axis)
    nl = np.concatenate([np.take(n, [0], axis=axis), n], axis=axis)
    nr = np.concatenate([n, np.take(n, [-1], axis=axis)], axis=axis)
    hf = np.maximum(el, er) - np.maximum(zl, zr)
    # The q-centred average: the old discharge here and at the two faces either side (the edge faces use their own).
    lo = np.concatenate([np.take(q, [0], axis=axis), np.take(q, range(q.shape[axis] - 1), axis=axis)], axis=axis)
    hi = np.concatenate([np.take(q, range(1, q.shape[axis]), axis=axis), np.take(q, [-1], axis=axis)], axis=axis)
    qc = theta * q + 0.5 * (1.0 - theta) * (lo + hi)
    nf = 0.5 * (nl + nr)
    wet = hf > H_DRY
    hs = np.where(wet, hf, 1.0)
    out = (qc - G * hs * dt * (er - el) / dx) / (1.0 + G * dt * nf * nf * np.abs(q) / hs ** (7.0 / 3.0))
    out = np.where(wet, out, 0.0)
    # Edge faces only let water out of the grid.
    first = [slice(None)] * 2
    last = [slice(None)] * 2
    first[axis] = 0
    last[axis] = -1
    out[tuple(first)] = np.minimum(out[tuple(first)], 0.0)
    out[tuple(last)] = np.maximum(out[tuple(last)], 0.0)
    return out


def _limit(h: np.ndarray, qx: np.ndarray, qy: np.ndarray, dt: float, dx: float) -> None:
    """Scale each cell's outgoing faces so that the cell cannot go below zero in this step (in place)."""
    out = (np.maximum(qx[:, 1:], 0) + np.maximum(-qx[:, :-1], 0)
           + np.maximum(qy[1:, :], 0) + np.maximum(-qy[:-1, :], 0)) * dt / dx
    f = np.where(out > h, h / np.where(out > 0, out, 1.0), 1.0)
    # Each face takes the factor of the cell its water leaves (its donor).
    fx_l = np.concatenate([f[:, :1], f], axis=1)     # the cell left of each x face (the edge: the edge cell)
    fx_r = np.concatenate([f, f[:, -1:]], axis=1)
    qx *= np.where(qx > 0, fx_l, fx_r)
    fy_t = np.concatenate([f[:1, :], f], axis=0)
    fy_b = np.concatenate([f, f[-1:, :]], axis=0)
    qy *= np.where(qy > 0, fy_t, fy_b)


def simulate(z: Any, dx: float, *, duration_s: float, h0: Any = None, n: Any = N_FLOODPLAIN,
             sources: list[tuple[int, int, float]] | None = None, frames_every_s: float | None = None,
             theta: float = THETA, alpha: float = ALPHA, dt_max: float = 10.0, closed: bool | str = False,
             max_steps: int = 5_000_000) -> dict[str, Any]:
    """Run the local-inertial scheme on bed ``z`` (m, rows north to south) with cells of side ``dx`` (m).

    ``sources`` are ``(row, col, flow m3/s)`` inflows, constant over the run. ``closed=True`` walls every edge
    in, and a string walls some of them (``"ns"``: the north and south edges); the others are a free outfall.
    Returns the final depth ``h``, the deepest depth of the run ``h_max``, the time each cell first got wetter
    than 0.1 m (``arrival_s``, NaN if never), the depth ``frames`` every ``frames_every_s`` with their
    ``times_s``, the step count, the wall time and the volume budget.
    """
    z = np.asarray(z, dtype="f8")
    if z.ndim != 2:
        raise ValueError("z must be a 2D grid")
    ny, nx = z.shape
    h = np.zeros_like(z) if h0 is None else np.array(h0, dtype="f8")
    nn = np.broadcast_to(np.asarray(n, dtype="f8"), z.shape)
    qx = np.zeros((ny, nx + 1))
    qy = np.zeros((ny + 1, nx))
    src = [(int(r), int(c), float(q)) for r, c, q in (sources or [])]
    walls = "nsew" if closed is True else (closed or "")
    h_max = h.copy()
    arrival = np.full(z.shape, np.nan)
    v0 = float(h.sum()) * dx * dx
    v_in = v_out = 0.0
    frames: list[np.ndarray] = []
    times: list[float] = []
    next_frame = 0.0 if frames_every_s else math.inf
    t = 0.0
    steps = 0
    started = time.perf_counter()
    while t < duration_s - 1e-9 and steps < max_steps:
        if t >= next_frame - 1e-9:
            frames.append(h.astype("f4"))
            times.append(t)
            next_frame += frames_every_s or math.inf
        hm = float(h.max())
        dt = min(dt_max, duration_s - t, alpha * dx / math.sqrt(G * hm) if hm > H_DRY else dt_max)
        if frames_every_s:
            dt = min(dt, max(next_frame - t, 1e-6))
        eta = z + h
        qx = _faces(z, eta, nn, qx, dt, dx, theta, axis=1)
        qy = _faces(z, eta, nn, qy, dt, dx, theta, axis=0)
        if "w" in walls:
            qx[:, 0] = 0.0
        if "e" in walls:
            qx[:, -1] = 0.0
        if "n" in walls:
            qy[0, :] = 0.0
        if "s" in walls:
            qy[-1, :] = 0.0
        _limit(h, qx, qy, dt, dx)
        v_out += float((qx[:, -1].sum() - qx[:, 0].sum() + qy[-1, :].sum() - qy[0, :].sum()) * dx * dt)
        h = h + dt / dx * (qx[:, :-1] - qx[:, 1:] + qy[:-1, :] - qy[1:, :])
        h = np.maximum(h, 0.0)   # only round-off can go below zero after the limiter
        for r, c, q in src:
            h[r, c] += q * dt / (dx * dx)
            v_in += q * dt
        t += dt
        steps += 1
        np.maximum(h_max, h, out=h_max)
        newly = np.isnan(arrival) & (h > 0.1)
        arrival[newly] = t
    if frames_every_s:
        frames.append(h.astype("f4"))
        times.append(t)
    wall = time.perf_counter() - started
    v1 = float(h.sum()) * dx * dx
    budget_err = v1 - (v0 + v_in - v_out)
    scale = max(v0 + v_in, 1e-12)
    return {
        "h": h, "h_max": h_max, "arrival_s": arrival, "frames": frames, "times_s": times, "t_s": t, "steps": steps,
        "wall_s": wall, "steps_per_s": steps / wall if wall > 0 else None,
        "volume": {"initial_m3": v0, "inflow_m3": v_in, "outflow_m3": v_out, "final_m3": v1,
                   "error_m3": budget_err, "error_rel": budget_err / scale},
    }


# ── the analytic check ──────────────────────────────────────────────────────


def ritter(x: Any, t: float, h0: float, x0: float = 0.0) -> np.ndarray:
    """Ritter's (1892) dam break on a dry, flat, frictionless bed: the depth at ``x`` a time ``t`` after a dam at
    ``x0`` holding ``h0`` upstream (x < x0) vanishes. ``c0 = sqrt(g h0)``; the depth is ``h0`` for
    ``x - x0 <= -c0 t``, ``(2 c0 - (x - x0)/t)^2 / (9 g)`` up to the front at ``x - x0 = 2 c0 t``, and zero beyond."""
    x = np.asarray(x, dtype="f8")
    c0 = math.sqrt(G * h0)
    xi = (x - x0) / t
    h = (2 * c0 - xi) ** 2 / (9 * G)
    return np.where(xi <= -c0, h0, np.where(xi >= 2 * c0, 0.0, h))


# ── the terrain ─────────────────────────────────────────────────────────────


def decode_png_rgb(data: bytes) -> np.ndarray:
    """An 8-bit RGB or RGBA, non-interlaced PNG (what the Terrarium tiles are) as a (rows, cols, 3) uint8 array."""
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")
    i = 8
    idat = b""
    width = height = channels = 0
    while i < len(data):
        (length,) = struct.unpack(">I", data[i:i + 4])
        kind = data[i + 4:i + 8]
        body = data[i + 8:i + 8 + length]
        if kind == b"IHDR":
            width, height, depth, colour, _, _, interlace = struct.unpack(">IIBBBBB", body)
            if depth != 8 or colour not in (2, 6) or interlace:
                raise ValueError("only 8-bit RGB(A) non-interlaced PNGs are read")
            channels = 3 if colour == 2 else 4
        elif kind == b"IDAT":
            idat += body
        elif kind == b"IEND":
            break
        i += 12 + length
    raw = np.frombuffer(zlib.decompress(idat), dtype="u1")
    stride = width * channels
    rows = raw.reshape(height, stride + 1)
    out = np.zeros((height, stride), dtype="i4")
    prev = np.zeros(stride, dtype="i4")
    for r in range(height):
        ftype = rows[r, 0]
        line = rows[r, 1:].astype("i4")
        if ftype == 0:
            cur = line
        elif ftype == 2:
            cur = (line + prev) & 255
        else:
            cur = np.zeros(stride, dtype="i4")
            for c in range(stride):
                a = cur[c - channels] if c >= channels else 0
                b = prev[c]
                if ftype == 1:
                    p = a
                elif ftype == 3:
                    p = (a + b) >> 1
                else:   # 4, Paeth
                    cc = prev[c - channels] if c >= channels else 0
                    pa, pb, pc = abs(b - cc), abs(a - cc), abs(a + b - 2 * cc)
                    p = a if pa <= pb and pa <= pc else (b if pb <= pc else cc)
                cur[c] = (line[c] + p) & 255
        out[r] = cur
        prev = cur
    return out.reshape(height, width, channels)[:, :, :3].astype("u1")


def terrarium_elevation(rgb: np.ndarray) -> np.ndarray:
    """Terrarium encoding to metres: ``R*256 + G + B/256 - 32768``."""
    rgb = np.asarray(rgb, dtype="f8")
    return rgb[..., 0] * 256.0 + rgb[..., 1] + rgb[..., 2] / 256.0 - 32768.0


def pixel_size_m(lat: float, zoom: int) -> float:
    """Ground size (m) of a 256-pixel Web Mercator tile's pixel at ``lat`` and ``zoom``."""
    return EARTH_CIRCUMFERENCE_M * math.cos(math.radians(lat)) / (256 * 2 ** zoom)


def _world_px(lon: float, lat: float, zoom: int) -> tuple[float, float]:
    s = 256 * 2 ** zoom
    x = (lon + 180.0) / 360.0 * s
    y = (1 - math.log(math.tan(math.radians(lat)) + 1 / math.cos(math.radians(lat))) / math.pi) / 2 * s
    return x, y


def _px_lonlat(x: float, y: float, zoom: int) -> tuple[float, float]:
    s = 256 * 2 ** zoom
    lon = x / s * 360.0 - 180.0
    lat = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / s))))
    return lon, lat


def _http_get(url: str) -> bytes:
    import httpx

    r = httpx.get(url, timeout=60, follow_redirects=True)
    r.raise_for_status()
    return r.content


def terrain_grid(lon: float, lat: float, *, zoom: int = 12, size: int = 512,
                 fetch: Callable[[str], bytes] | None = None) -> dict[str, Any]:
    """A ``size`` x ``size`` window of Terrarium terrain centred on (lon, lat), at ``zoom``'s pixel size.

    Returns ``z`` (m, rows north to south), ``dx`` (m, the pixel's ground size at the centre), the window's
    world-pixel origin ``x0``/``y0`` (Web Mercator pixels at ``zoom``) and its ``bounds`` [W, S, E, N]."""
    get = fetch or _http_get
    cx, cy = _world_px(lon, lat, zoom)
    x0, y0 = int(round(cx - size / 2)), int(round(cy - size / 2))
    tx0, ty0 = x0 // 256, y0 // 256
    tx1, ty1 = (x0 + size - 1) // 256, (y0 + size - 1) // 256
    mosaic = np.zeros(((ty1 - ty0 + 1) * 256, (tx1 - tx0 + 1) * 256))
    n = 2 ** zoom
    for ty in range(ty0, ty1 + 1):
        for tx in range(tx0, tx1 + 1):
            rgb = decode_png_rgb(get(TERRARIUM.format(z=zoom, x=tx % n, y=ty)))
            mosaic[(ty - ty0) * 256:(ty - ty0 + 1) * 256, (tx - tx0) * 256:(tx - tx0 + 1) * 256] = \
                terrarium_elevation(rgb)
    ox, oy = x0 - tx0 * 256, y0 - ty0 * 256
    z = mosaic[oy:oy + size, ox:ox + size]
    w, n_ = _px_lonlat(x0, y0, zoom)
    e, s = _px_lonlat(x0 + size, y0 + size, zoom)
    return {"z": z, "dx": pixel_size_m(lat, zoom), "x0": x0, "y0": y0, "zoom": zoom, "size": size,
            "bounds": [w, s, e, n_]}


def cell_of(grid: dict[str, Any], lon: float, lat: float) -> tuple[int, int]:
    """The (row, col) of the grid cell holding (lon, lat)."""
    x, y = _world_px(lon, lat, grid["zoom"])
    return int(math.floor(y - grid["y0"])), int(math.floor(x - grid["x0"]))


def cell_lonlat(grid: dict[str, Any], row: Any, col: Any) -> tuple[np.ndarray, np.ndarray]:
    """Cell-centre longitudes and latitudes for arrays of rows and columns."""
    s = 256 * 2 ** grid["zoom"]
    x = grid["x0"] + np.asarray(col, dtype="f8") + 0.5
    y = grid["y0"] + np.asarray(row, dtype="f8") + 0.5
    lon = x / s * 360.0 - 180.0
    lat = np.degrees(np.arctan(np.sinh(np.pi * (1 - 2 * y / s))))
    return lon, lat


# ── the channel ─────────────────────────────────────────────────────────────


def channel_depth(q_bankfull: float, width_m: float, slope: float, n: float = N_CHANNEL) -> float:
    """The depth (m) at which a wide rectangular channel of ``width_m`` carries ``q_bankfull`` (m3/s) at normal
    depth by Manning's equation, ``Q = (1/n) W d^(5/3) S^(1/2)``, so ``d = (n Q / (W S^(1/2)))^(3/5)``. The slope
    is kept within :data:`SLOPE_RANGE`."""
    s = min(max(float(slope), SLOPE_RANGE[0]), SLOPE_RANGE[1])
    return (n * float(q_bankfull) / (float(width_m) * math.sqrt(s))) ** 0.6


def _line_cells(grid: dict[str, Any], line: list[tuple[float, float]]) -> list[tuple[int, int]]:
    """The cells a (lon, lat) polyline passes through, in order, without repeats (sampled every half cell)."""
    zoom = grid["zoom"]
    pts = [_world_px(lo, la, zoom) for lo, la in line]
    out: list[tuple[int, int]] = []
    size = grid["size"]
    for (ax, ay), (bx, by) in zip(pts, pts[1:]):
        k = max(1, int(math.ceil(2 * math.hypot(bx - ax, by - ay))))
        for i in range(k + 1):
            x = ax + (bx - ax) * i / k - grid["x0"]
            y = ay + (by - ay) * i / k - grid["y0"]
            r, c = int(math.floor(y)), int(math.floor(x))
            if 0 <= r < size and 0 <= c < size and (not out or out[-1] != (r, c)):
                if (r, c) not in out[-3:]:
                    out.append((r, c))
    return out


def burn_channel(z: np.ndarray, cells: list[tuple[int, int]], depth_m: float) -> np.ndarray:
    """A copy of ``z`` with a channel ``depth_m`` deep along ``cells`` (ordered upstream to downstream).

    The bank line is first made non-increasing downstream (a running minimum), so the burnt bed has no pits or
    bumps from the terrain's noise; the bed is that line minus ``depth_m``."""
    out = np.array(z, dtype="f8")
    if not cells:
        return out
    bank = np.minimum.accumulate(np.array([z[r, c] for r, c in cells], dtype="f8"))
    for (r, c), b in zip(cells, bank):
        out[r, c] = min(out[r, c], b - depth_m)
    return out


def _orient(z: np.ndarray, cells: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Cells ordered so that the higher end comes first (upstream)."""
    if len(cells) > 1 and z[cells[0]] < z[cells[-1]]:
        return cells[::-1]
    return cells


# ── comparing with a reference map ──────────────────────────────────────────


def compare_depths(model: Any, reference: Any, *, wet_m: float = 0.05, mask: Any = None) -> dict[str, Any]:
    """Agreement between two depth grids on the same cells: the cells wet in both, in one only, the critical
    success index ``hits / (hits + misses + false alarms)`` (the share of the cells wet in either that are wet in
    both), the hit rate (share of the reference's wet cells the model wets), and the depth difference (model -
    reference) where both are wet."""
    m = np.asarray(model, dtype="f8")
    ref = np.asarray(reference, dtype="f8")
    ok = np.isfinite(ref) & np.isfinite(m)
    if mask is not None:
        ok &= np.asarray(mask, dtype=bool)
    mw = (m > wet_m) & ok
    rw = (ref > wet_m) & ok
    hits = int((mw & rw).sum())
    misses = int((~mw & rw).sum())
    false = int((mw & ~rw).sum())
    either = hits + misses + false
    d = (m - ref)[mw & rw]
    return {
        "cells": int(ok.sum()), "wet_model": int(mw.sum()), "wet_reference": int(rw.sum()), "wet_both": hits,
        "reference_only": misses, "model_only": false,
        "csi": round(hits / either, 3) if either else None,
        "hit_rate": round(hits / (hits + misses), 3) if hits + misses else None,
        "false_alarm_ratio": round(false / (hits + false), 3) if hits + false else None,
        "depth_diff_m": None if not d.size else {
            "mean": round(float(d.mean()), 2), "median": round(float(np.median(d)), 2),
            "mean_abs": round(float(np.abs(d).mean()), 2), "p10": round(float(np.percentile(d, 10)), 2),
            "p90": round(float(np.percentile(d, 90)), 2)},
        "wet_m": wet_m,
    }


def jrc_on_grid(grid: dict[str, Any], return_period: int = 100) -> np.ndarray:
    """The JRC CEMS-GloFAS flood depth (m, 3 arc-seconds) of ``return_period`` at every grid cell's centre (nearest
    pixel; NaN outside the map, 0 where JRC has no flood)."""
    from aquascope.flood_depth import NODATA, depth_url, tile_for
    from aquascope.utils.cog import COGNotFound, open_cog

    size = grid["size"]
    rr, cc = np.mgrid[0:size, 0:size]
    lon, lat = cell_lonlat(grid, rr, cc)
    out = np.full((size, size), np.nan)
    names = {tile_for(float(la), float(lo)) for la, lo in ((lat.min(), lon.min()), (lat.min(), lon.max()),
                                                           (lat.max(), lon.min()), (lat.max(), lon.max()))}
    for name in sorted(x for x in names if x):
        try:
            cog = open_cog(depth_url(name, return_period))
        except COGNotFound:
            continue
        img = cog.images[0]
        x0, dxp, _, y0, _, dyp = img.transform
        col = np.floor((lon - x0) / dxp).astype(int)
        row = np.floor((lat - y0) / dyp).astype(int)
        inside = (col >= 0) & (col < img.width) & (row >= 0) & (row < img.height)
        if not inside.any():
            continue
        c0, c1 = int(col[inside].min()), int(col[inside].max())
        r0, r1 = int(row[inside].min()), int(row[inside].max())
        win = np.asarray(cog.read_window(c0, r0, c1 + 1, r1 + 1), dtype="f8")
        vals = win[row[inside] - r0, col[inside] - c0]
        vals = np.where((vals == NODATA) | ~np.isfinite(vals), 0.0, vals)
        out[inside] = vals
    return out


# ── the whole thing, from a reach ───────────────────────────────────────────


def _join(segments: list[tuple[float, float, float, float]]) -> list[list[tuple[float, float]]]:
    """(lon1, lat1, lon2, lat2) segments chained into polylines where one ends where the next starts."""
    lines: list[list[tuple[float, float]]] = []
    for x1, y1, x2, y2 in segments:
        if lines and lines[-1][-1] == (x1, y1):
            lines[-1].append((x2, y2))
        else:
            lines.append([(x1, y1), (x2, y2)])
    return lines


def reach_grid(lat: float, lon: float, *, river_id: int, bankfull_cms: float, order: int | None = None,
               zoom: int = 12, size: int = 512, lines: list[list[tuple[float, float]]] | None = None,
               grid: dict[str, Any] | None = None) -> dict[str, Any]:
    """The grid a reach floods on: Terrarium terrain around (lon, lat), the GEOGLOWS river lines (the reach and the
    rivers at least its Strahler order) burnt in to the depth that carries ``bankfull_cms`` (use the 2-year flow)
    in a one-cell-wide channel (:func:`channel_depth`), Manning's n per cell, and the burnt cell nearest the click.
    ``grid`` and ``lines`` skip the downloads."""
    grid = grid or terrain_grid(lon, lat, zoom=zoom, size=size)
    z = grid["z"]
    dx = grid["dx"]
    if lines is None:
        from aquascope.flood_depth import river_lines

        half_km = grid["size"] * dx / 2000 * math.sqrt(2)
        got = river_lines(lon, lat, half_km, int(river_id), order)
        lines = _join(got[0]) + _join(got[1]) if got else []
    paths = [_orient(z, _line_cells(grid, ln)) for ln in lines]
    cells = [c for p in paths for c in p]
    # The slope that sizes the channel: the network's drop (5th to 95th percentile of its banks) over its extent.
    slope = SLOPE_RANGE[0]
    if len(cells) > 2:
        zs = np.array([z[c] for c in cells])
        rr = np.array([c[0] for c in cells], dtype="f8")
        cc = np.array([c[1] for c in cells], dtype="f8")
        span = math.hypot(rr.max() - rr.min(), cc.max() - cc.min()) * dx
        slope = float((np.percentile(zs, 95) - np.percentile(zs, 5)) / max(span, dx))
    depth = channel_depth(bankfull_cms, dx, slope)
    zb = np.array(z)
    for p in paths:
        zb = burn_channel(zb, p, depth)
    n = np.full(z.shape, N_FLOODPLAIN)
    for r, c in cells:
        n[r, c] = N_CHANNEL
    r0, c0 = cell_of(grid, lon, lat)
    if cells:   # the nearest burnt cell to the click
        r0, c0 = min(cells, key=lambda rc: (rc[0] - r0) ** 2 + (rc[1] - c0) ** 2)
    return {"grid": grid, "z_burnt": zb, "n": n, "click_cell": [r0, c0],
            "channel": {"depth_m": round(depth, 2), "width_m": round(dx, 1), "slope": round(slope, 5),
                        "bankfull_cms": bankfull_cms, "cells": len(cells)}}


def flood_reach(lat: float, lon: float, *, river_id: int, flow_cms: float, bankfull_cms: float,
                order: int | None = None, hours: float = 6.0, zoom: int = 12, size: int = 512,
                frames_every_s: float | None = 1800.0, lines: list[list[tuple[float, float]]] | None = None,
                grid: dict[str, Any] | None = None, inflow_cell: tuple[int, int] | None = None,
                walls: str = "") -> dict[str, Any]:
    """A screening flood: ``flow_cms`` (m3/s) poured into the reach's cell nearest (lon, lat) for ``hours``, on the
    Terrarium terrain around it, with the GEOGLOWS river lines (the reach and the rivers at least its size) burnt
    in to the depth that carries ``bankfull_cms`` (use the 2-year flow) in a one-cell-wide channel.

    ``inflow_cell`` (row, col) puts the inflow elsewhere, for example where the river enters the grid, so the
    whole reach in view floods as in a hazard map; ``walls`` ("n", "s", "e", "w") closes those edges, for example the
    upstream one, which is otherwise a free outfall like the others. ``grid`` and ``lines`` skip the downloads."""
    g = reach_grid(lat, lon, river_id=river_id, bankfull_cms=bankfull_cms, order=order, zoom=zoom, size=size,
                   lines=lines, grid=grid)
    r0, c0 = g["click_cell"] if inflow_cell is None else (int(inflow_cell[0]), int(inflow_cell[1]))
    res = simulate(g["z_burnt"], g["grid"]["dx"], duration_s=hours * 3600.0, n=g["n"],
                   sources=[(r0, c0, float(flow_cms))], frames_every_s=frames_every_s, closed=walls)
    res.update({"grid": g["grid"], "z_burnt": g["z_burnt"], "source_cell": [r0, c0], "flow_cms": flow_cms,
                "hours": hours, "channel": g["channel"]})
    return res
