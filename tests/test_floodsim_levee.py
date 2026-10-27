"""Levees and breaches in Flood it (#663): the supercover line, a levee that holds (straight and diagonal), one that
overtops with the volume still balanced, a breach against a broad-crested weir, the twin fixture with a levee in
place, and the three faces."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from aquascope import floodsim as fs

ZOOM = 12


def _grid(z: np.ndarray, lon: float = -2.2, lat: float = 52.05) -> dict:
    """A prepare_grid-like answer around (lon, lat) at zoom 12 (cells of about 23.5 m), walled all round."""
    x0, y0 = fs._world_px(lon, lat, ZOOM)
    ny, nx = z.shape
    win = fs.window(int(x0), int(y0), nx, ny, ZOOM)
    return {"window": win, "nx": nx, "ny": ny, "dx": win["dx"], "z": z.astype("f8"), "z_burnt": z.astype("f8"),
            "n": np.full(z.shape, fs.N_FLOODPLAIN), "edges": fs.edge_codes(ny, nx, closed=True)}


def _ll(g: dict, r: float, c: float) -> list[float]:
    """The (lon, lat) of a point at (row, col) in cell units (cell centres at .5)."""
    w = g["window"]
    return list(fs._px_lonlat(w["x0"] + c, w["y0"] + r, w["zoom"]))


def _run(g: dict, sources, hours: float, **kw) -> dict:
    return fs.simulate(g["z_burnt"], g["dx"], duration_s=hours * 3600.0, n=g["n"], sources=sources,
                       edges=g["edges"], **kw)


# ── the line ────────────────────────────────────────────────────────────────


def test_a_supercover_line_is_four_connected_and_takes_both_cells_at_a_corner():
    g = _grid(np.zeros((40, 40)))
    # A diagonal through cell corners exactly: both side cells at every corner.
    cells = fs._supercover(2.0, 2.0, 6.0, 6.0)
    assert cells[:4] == [(2, 2), (2, 3), (3, 2), (3, 3)]
    # Any drawn line: consecutive cells share a face (or both sides of a corner were taken), no repeats.
    line = [_ll(g, 3.3, 2.7), _ll(g, 30.6, 21.2), _ll(g, 12.1, 37.9)]
    got = fs.supercover_cells(g, line)
    assert len(got) == len(set(got))
    pos = {rc: i for i, rc in enumerate(got)}
    for (r, c) in got[1:]:
        assert any(abs(r - rr) + abs(c - cc) == 1 and pos[(rr, cc)] < pos[(r, c)] for rr, cc in got), (r, c)
    # The same line, the same cells; the old half-cell sampling (the channel's) is eight-connected only.
    assert fs.supercover_cells(g, line) == got
    sampled = fs._line_cells(g["window"], line)
    assert set(sampled) <= set(got)
    # Cells off the grid are left out; a single point is not a line.
    assert all(0 <= r < 40 and 0 <= c < 40 for r, c in fs.supercover_cells(g, [_ll(g, -5, 5), _ll(g, 45, 30)]))
    with pytest.raises(ValueError):
        fs.supercover_cells(g, [_ll(g, 3, 3)])


def test_levee_and_breach_edit_a_copy_and_say_what_they_changed():
    z = np.full((20, 30), 10.0)
    z[:, 12] = 9.0                                      # a burnt channel down column 12
    g = _grid(z)
    g["z_burnt"] = g["z"].copy()
    g["z_burnt"][:, 12] = 7.0
    line = [_ll(g, 10.5, 2.5), _ll(g, 10.5, 25.5)]
    out = fs.apply_levee(g, line, 12.0)
    assert out["cells"] == [[10, c] for c in range(2, 26)]
    assert out["terrain"][10, 2:26].min() == 12.0 and out["ground"][10, 2:26].min() == 12.0
    assert out["raised"] == 24 and out["channel_cells"] == 1          # it dams the channel, and says so
    assert out["height_m"] == [2.0, 3.0] and out["crest_m"] == 12.0
    assert (g["z_burnt"][10] == np.where(np.arange(30) == 12, 7.0, 10.0)).all()   # the input is unchanged
    # A crest below the ground leaves the ground alone.
    assert fs.apply_levee(g, line, 5.0)["raised"] == 0
    # The edits in order: a levee, then a breach without an invert, which puts the natural ground back.
    e = fs.edit_grid(g, [{"kind": "levee", "line": line, "crest": 12.0},
                         {"kind": "breach", "line": [_ll(g, 10.5, 15.5), _ll(g, 10.5, 18.5)]}])
    assert e["z"][10, 15:19].tolist() == [10.0] * 4 and e["z"][10, 14] == 12.0 and e["z"][10, 19] == 12.0
    assert e["edits"][1]["lowered"] == 4 and e["edits"][1]["invert_m"] is None
    assert (e["natural"]["z_burnt"] == g["z_burnt"]).all() and e["z_burnt"][10, 12] == 12.0
    # The same edits give the same grid; an invert lowers below the natural ground.
    e2 = fs.edit_grid(g, [{"kind": "levee", "line": line, "crest": 12.0},
                          {"kind": "breach", "line": [_ll(g, 10.5, 15.5), _ll(g, 10.5, 18.5)]}])
    assert (e2["z_burnt"] == e["z_burnt"]).all()
    lo = fs.apply_breach(e, [_ll(g, 10.5, 15.5), _ll(g, 10.5, 16.5)], invert=8.5)
    assert lo["terrain"][10, 15] == 8.5 and lo["ground"][10, 16] == 8.5
    with pytest.raises(ValueError):
        fs.edit_grid(g, [{"kind": "levee", "line": line}])
    with pytest.raises(ValueError):
        fs.edit_grid(g, [{"kind": "wall", "line": line, "crest": 3}])


def test_reading_lines_from_geojson():
    ls = {"type": "LineString", "coordinates": [[-2.2, 52.0], [-2.1, 52.1]]}
    assert fs.read_line(ls) == [[[-2.2, 52.0], [-2.1, 52.1]]]
    fc = {"type": "FeatureCollection", "features": [{"type": "Feature", "geometry": ls, "properties": {}},
                                                    {"type": "Feature", "geometry": {
                                                        "type": "MultiLineString",
                                                        "coordinates": [[[0, 0], [1, 1]], [[2, 2], [3, 3]]]}}]}
    assert len(fs.read_line(json.dumps(fc))) == 3
    assert fs.read_line([[0, 0], [1, 1]]) == [[[0.0, 0.0], [1.0, 1.0]]]
    with pytest.raises(ValueError):
        fs.read_line({"type": "Polygon", "coordinates": []})


# ── the analytic checks ─────────────────────────────────────────────────────


def _basin(ny: int = 60, nx: int = 60) -> np.ndarray:
    """A flat plain at 10 m, slightly tilted down to the south-east so water spreads everywhere it can reach."""
    rows, cols = np.mgrid[0:ny, 0:nx]
    return 10.0 - 0.0005 * (rows + cols)


@pytest.mark.parametrize("diagonal", [False, True])
def test_a_levee_higher_than_the_flood_keeps_every_cell_behind_it_dry(diagonal):
    g = _grid(_basin())
    if diagonal:   # a levee at an angle across the plain, through no cell corner on purpose
        line = [_ll(g, -1.0, 59.3), _ll(g, 61.0, 7.9)]
    else:
        line = [_ll(g, -1.0, 40.5), _ll(g, 61.0, 40.5)]
    e = fs.edit_grid(g, [{"kind": "levee", "line": line, "crest": 12.0}])
    levee = np.zeros(g["z"].shape, bool)
    for r, c in e["edits"][0]["cells"]:
        levee[r, c] = True
    # The protected side: the cells beyond the levee from the inflow at the north-west corner.
    from scipy import ndimage

    lab, _ = ndimage.label(~levee, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    behind = (lab != lab[2, 2]) & ~levee
    assert behind.sum() > 300
    res = _run(e, [(2, 2, 40.0)], 6.0)
    # The flood stands about 1.5 m deep against a 2 m levee, and not a drop gets behind it.
    assert res["h_max"][~behind & ~levee].max() > 1.0
    assert float((e["z_burnt"] + res["h_max"])[~behind & ~levee].max()) < 12.0
    assert res["h_max"][behind].max() == 0.0
    assert abs(res["volume"]["error_rel"]) < 1e-12


def test_a_levee_lower_than_the_water_overtops_and_the_volume_still_balances():
    g = _grid(_basin())
    line = [_ll(g, -1.0, 59.3), _ll(g, 61.0, 7.9)]
    e = fs.edit_grid(g, [{"kind": "levee", "line": line, "crest": 10.6}])
    res = _run(e, [(2, 2, 40.0)], 6.0)
    levee = np.zeros(g["z"].shape, bool)
    for r, c in e["edits"][0]["cells"]:
        levee[r, c] = True
    from scipy import ndimage

    lab, _ = ndimage.label(~levee, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    behind = (lab != lab[2, 2]) & ~levee
    assert float((e["z_burnt"] + res["h_max"])[2, 2]) > 10.6          # the water rises over the crest
    assert (res["h_max"][behind] > 0.05).sum() > 100                  # and runs down the other side
    v = res["volume"]
    assert v["outflow_m3"] == 0.0 and abs(v["error_rel"]) < 1e-12
    assert v["final_m3"] == pytest.approx(40.0 * 6 * 3600.0, rel=1e-12)


def _weir_grid(q_unit: float, b_cells: int = 8) -> tuple[dict, float, int, int]:
    """A pond on a flat floor at 0 m held by a levee with a breach of ``b_cells`` cells down to the floor; beyond
    it the ground falls 2 m to a free outfall, so the gap is never drowned. The inflow is ``q_unit`` m2/s per
    metre of gap."""
    ny, nxr, nxp = 40, 30, 20
    z = np.zeros((ny, nxr + 1 + nxp))
    z[:, nxr + 1:] = -2.0
    g = _grid(z)
    g["edges"]["e"][:] = fs.FREE
    levee = [_ll(g, -1.0, nxr + 0.5), _ll(g, ny + 1.0, nxr + 0.5)]
    r0 = ny // 2 - b_cells // 2
    breach = [_ll(g, r0 + 0.5, nxr + 0.5), _ll(g, r0 + b_cells - 0.5, nxr + 0.5)]
    e = fs.edit_grid(g, [{"kind": "levee", "line": levee, "crest": 5.0},
                         {"kind": "breach", "line": breach, "invert": 0.0}])
    assert e["edits"][1]["lowered"] == b_cells
    return e, q_unit * b_cells * e["dx"], r0, nxr


@pytest.mark.parametrize("q_unit", [1.0, 2.0, 4.0])
def test_a_breach_floods_the_protected_side_like_a_broad_crested_weir(q_unit):
    """At a steady state the flow through the gap is the inflow; the head over the gap's floor, read just upstream,
    gives a broad-crested weir's flow ``C b H^1.5``. The local-inertial scheme has no critical-flow control at a
    crest: the flow over a crest cell is set by the water-surface slope and friction across it, so the match
    depends on the head over the cell size. At Flood it's cells (about 24 m) and n = 0.06, with heads of 0.8 to
    1.7 m, the flow is within 30% of the weir's (docs/floodsim.md has the sweep)."""
    e, q, r0, col = _weir_grid(q_unit)
    b = 8 * e["dx"]
    res = _run(e, [(20, 2, q)], 24.0, steady_tol=1e-4, steady_every_s=900.0, steady_min_s=3600.0)
    assert res["steady_at_s"] is not None
    assert res["h"][:, col + 1:].max() > 0.05                       # the protected side floods
    head = float(np.mean(res["h"][r0:r0 + 8, col - 3:col - 1]))
    ratio = q / fs.weir_flow(b, head)
    assert 0.7 < ratio < 1.3, (head, ratio)
    assert abs(res["volume"]["error_rel"]) < 1e-12


def test_the_weir_law():
    assert fs.WEIR_C == pytest.approx(1.7049, abs=1e-4)
    assert fs.weir_flow(10.0, 1.0) == pytest.approx(17.049, abs=1e-3)
    assert fs.weir_flow(10.0, -1.0) == 0.0


# ── the twin, the profile, the before and after ─────────────────────────────

FIXTURE = Path(__file__).parent / "fixtures" / "floodsim_twin_levee.json"


def _levee_fixture_grid() -> tuple[dict, list]:
    g = _grid(_basin(24, 30))
    g["edges"]["e"][10:14] = 0.002                      # a normal-depth outfall on the east side
    line = [_ll(g, -1.0, 27.2), _ll(g, 25.0, 3.9)]
    e = fs.edit_grid(g, [{"kind": "levee", "line": line, "crest": 10.7},
                         {"kind": "breach", "line": [_ll(g, 11.2, 15.0), _ll(g, 12.8, 13.4)], "invert": 9.9}])
    return e, [[2, 2, 6.0], [3, 5, [0.0, 4.0, 2.0]]]


def test_the_javascript_twin_matches_numpy_with_a_levee_in_place():
    """explorer/tests/floodsim-core.test.mjs runs this fixture's grid (a diagonal levee with a breach) in JS; the WGSL
    was checked on it in headless Chrome. Here it is regenerated and checked so the three cannot drift apart."""
    e, srcs = _levee_fixture_grid()
    fx = json.loads(FIXTURE.read_text())
    assert np.array(fx["z"]) == pytest.approx(e["z_burnt"], abs=1e-9)
    edges = {k: np.array([np.nan if x is None else x for x in v]) for k, v in fx["edges"].items()}
    r = fs.simulate(np.array(fx["z"]), fx["dx"], duration_s=fx["duration_s"], n=np.array(fx["n"]),
                    sources=[tuple(s) for s in fx["sources"]], edges=edges, table_dt_s=fx["table_dt_s"])
    assert r["steps"] == fx["steps"]
    assert np.abs(r["h"] - np.array(fx["h"])).max() < 1e-6
    assert abs(r["volume"]["error_rel"]) < 1e-12


def _write_fixture() -> None:   # pragma: no cover - run by hand when the fixture changes
    e, srcs = _levee_fixture_grid()
    dur = 2 * 3600.0
    r = fs.simulate(e["z_burnt"], e["dx"], duration_s=dur, n=e["n"], sources=[tuple(s) for s in srcs],
                    edges=e["edges"], table_dt_s=1800.0)
    FIXTURE.write_text(json.dumps({
        "z": e["z_burnt"].round(9).tolist(), "n": e["n"].tolist(), "dx": e["dx"], "duration_s": dur,
        "table_dt_s": 1800.0, "sources": srcs,
        "edges": {k: [None if math.isnan(x) else float(x) for x in v] for k, v in e["edges"].items()},
        "steps": r["steps"], "h": r["h"].tolist(), "volume": r["volume"],
        "wet": int((r["h"] > 0.01).sum())}))


def test_the_profile_suggests_a_crest_above_the_flood():
    g = _grid(_basin())
    line = [_ll(g, 5.5, 30.5), _ll(g, 5.5, 50.5)]
    p = fs.levee_profile(g, line)
    assert p["highest_water_m"] is None and p["suggested_crest_m"] == pytest.approx(math.ceil((p["highest_ground_m"]
                                                                                               + 1) * 10) / 10)
    h = np.zeros(g["z"].shape)
    h[5, 35:40] = 0.8
    p = fs.levee_profile(g, line, h)
    top = max(g["z"][5, c] + 0.8 for c in range(35, 40))
    assert p["highest_water_m"] == pytest.approx(top, abs=0.01)
    assert p["suggested_crest_m"] == pytest.approx(math.ceil((round(top, 2) + fs.FREEBOARD_M) * 10) / 10)
    assert p["water_m"][0] is None and len(p["distance_m"]) == len(p["ground_m"]) == 21
    assert p["distance_m"][-1] == pytest.approx(20 * g["dx"], abs=0.06)


def test_before_and_after_on_the_same_grid():
    g = _grid(_basin())
    line = [_ll(g, -1.0, 40.5), _ll(g, 61.0, 40.5)]
    e = fs.edit_grid(g, [{"kind": "levee", "line": line, "crest": 12.0}])
    r0 = _run(g, [(2, 2, 40.0)], 6.0)
    r1 = _run(e, [(2, 2, 40.0)], 6.0)
    c = fs.compare_runs(g, r0["h_max"], e, r1["h_max"])
    assert c["kept_dry_km2"] > 0 and c["wet_km2_after"] < c["wet_km2_before"]
    assert c["wet_km2_change"] == pytest.approx(c["wet_km2_after"] - c["wet_km2_before"], abs=2e-3)
    assert c["depth_change_m"]["max"] > 0                       # the water stands deeper in front of the levee
