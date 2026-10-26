// Flood it (#658), the spike: from a river reach's card, pour its 100-year flow (or the 2- or 10-year) into the
// reach for six hours and watch the water spread over the terrain, on the visitor's GPU (WebGPU) or, without
// one, in a worker on a coarser grid. A screening flood: no bathymetry, terrain errors of metres.
//
// The pieces: the Terrarium terrain the hillshade already uses (512 x 512 cells at zoom 12, about 20 m to 38 m),
// the GEOGLOWS river lines in view (the reach and every river at least its order) burnt in to the depth that
// carries the 2-year flow (Manning, floodsim-core.js channelDepth), the flow from the reach's own return periods
// (the worker's "now" forecast: Log-Pearson III on GEOGLOWS's simulated record since 1940), the solver
// (floodsim-gpu.js or floodsim-worker.js), and the depth drawn as a raster with the time bar on the model's
// own clock. aquascope.floodsim is the same scheme in numpy for the checks; every run goes in the undo log.

import { actions, state } from "./core.js?v=__BUILD__";
import { map } from "./map.js?v=__BUILD__";
import { openCard, updateCard } from "./map-card.js?v=__BUILD__";
import { releaseModelClock, useModelClock } from "./time-ui.js?v=__BUILD__";
import { call } from "./worker-client.js?v=__BUILD__";
import { logAction, registerSlotWriter, slotValue } from "./map-actions.js?v=__BUILD__";
import { announce } from "./a11y.js?v=__BUILD__";
import {
  FLOODSIM_CREDIT, LABEL, N_CHANNEL, N_FLOODPLAIN, burnChannel, channelDepth, coarsen, depthImage, lineCells,
  networkSlope, orient, summary, terrariumElevation, windowTiles, worldPx,
} from "./floodsim-core.js?v=__BUILD__";
import { createGpuSim, gpuAvailable } from "./floodsim-gpu.js?v=__BUILD__";

const SRC = "floodsim-depth";
const ZOOM = 12;
const SIZE = 512;          // the GPU grid; the worker runs SIZE / 2
const HOURS = 6;
const FRAME_S = 900;       // a frame every 15 minutes of model time
const RIVER_SOURCE = "river-net";

let shown = null;          // { target, rp, flow, frames: [{ t, url }], bounds, stats } on the map now
let running = 0;

const hm = (s) => {
  const m = Math.round(s / 60);
  return m < 60 ? `${m} min` : `${Math.floor(m / 60)} h${m % 60 ? ` ${m % 60} min` : ""}`;
};
const fmt = (v, d = 0) => Number(v).toLocaleString("en-GB", { maximumFractionDigits: d, minimumFractionDigits: d });

// ── inputs: terrain, lines, flows ───────────────────────────────────────────

async function readTerrain(win) {
  const canvas = new OffscreenCanvas(win.cols * 256, win.rows * 256);
  const ctx = canvas.getContext("2d", { willReadFrequently: true });
  await Promise.all(win.tiles.map(async (t) => {
    const r = await fetch(t.url);
    if (!r.ok) throw new Error(`terrain tile ${t.x}/${t.y}: HTTP ${r.status}`);
    const bmp = await createImageBitmap(await r.blob(), { colorSpaceConversion: "none", premultiplyAlpha: "none" });
    ctx.drawImage(bmp, t.col * 256, t.row * 256);
  }));
  return terrariumElevation(ctx.getImageData(win.ox, win.oy, win.size, win.size, { colorSpace: "srgb" }).data, 4);
}

function riverLines(t) {
  if (!state.mapOk || !map.getSource(RIVER_SOURCE)) return [];
  const feats = map.querySourceFeatures(RIVER_SOURCE, { sourceLayer: "streams" });
  const own = feats.find((f) => Number(f.properties.riverId) === Number(t.river_id));
  const min = Number(t.order || (own && own.properties.strahlerOrder) || 0);
  const lines = [];
  const seen = new Set();
  for (const f of feats) {
    const p = f.properties;
    if (Number(p.riverId) !== Number(t.river_id) && !(Number(p.strahlerOrder) >= min)) continue;
    const g = f.geometry;
    const parts = g.type === "LineString" ? [g.coordinates] : g.type === "MultiLineString" ? g.coordinates : [];
    for (const part of parts) {
      const key = `${p.riverId}:${part[0]}:${part.at(-1)}`;
      if (part.length < 2 || seen.has(key)) continue;
      seen.add(key);
      lines.push(part);
    }
  }
  return lines;
}

async function readFlows(t) {
  const res = await call("now", { op: "forecast", args: { lat: t.lat, lon: t.lon, river_id: t.river_id, glofas: false } });
  const thr = res && res.thresholds;
  if (!thr || !Array.isArray(thr.q)) throw new Error((thr && thr.error) || "this reach has no return-period flows");
  const at = (T) => thr.q[thr.return_periods.indexOf(T)];
  return { q2: at(2), q10: at(10), q100: at(100), method: thr.method };
}

/** The grid one run floods on: the terrain with the rivers burnt in, Manning's n, the inflow cell. */
export function prepareGrid(z, win, lines, qBankfull, click) {
  const zb = Float64Array.from(z);
  const n = new Float64Array(z.length).fill(N_FLOODPLAIN);
  const paths = lines.map((ln) => orient(z, win.size, lineCells(win, ln)));
  const cells = paths.flat();
  const slope = networkSlope(z, win.size, cells, win.dx);
  const depth = channelDepth(qBankfull, win.dx, slope);
  for (const p of paths) burnChannel(zb, win.size, p, depth);
  for (const [r, c] of cells) n[r * win.size + c] = N_CHANNEL;
  let src = click;
  if (cells.length) {
    let best = Infinity;
    for (const [r, c] of cells) {
      const d = (r - click[0]) ** 2 + (c - click[1]) ** 2;
      if (d < best) { best = d; src = [r, c]; }
    }
  }
  return { z: zb, n, source: src, channel: { depth_m: depth, slope, cells: cells.length } };
}

// ── the run ─────────────────────────────────────────────────────────────────

async function frameUrl(h, size) {
  const c = new OffscreenCanvas(size, size);
  c.getContext("2d").putImageData(new ImageData(depthImage(h, size, size), size, size), 0, 0);
  return URL.createObjectURL(await c.convertToBlob({ type: "image/png" }));
}

async function runGpu(grid, src, q, onFrame) {
  const sim = await createGpuSim(grid, { source: [...src, q] });
  if (!sim) return null;
  try {
    const started = performance.now();
    let steps = 0;
    await onFrame(0, await sim.readDepth());
    for (let t = FRAME_S; t <= HOURS * 3600; t += FRAME_S) {
      const r = await sim.advance(t, 400);
      steps = r.steps;
      await onFrame(r.t, await sim.readDepth());
    }
    return { engine: "WebGPU", steps, wall_s: (performance.now() - started) / 1000, hMax: await sim.readMaxDepth() };
  } finally {
    sim.destroy();
  }
}

function runWorker(grid, src, q, onFrame) {
  return new Promise((resolve, reject) => {
    const w = new Worker(new URL("./floodsim-worker.js?v=__BUILD__", import.meta.url), { type: "module" });
    let chain = Promise.resolve();
    w.onmessage = (e) => {
      const m = e.data;
      if (m.frame) chain = chain.then(() => onFrame(m.frame.t, m.frame.h));
      else if (m.done) { chain.then(() => { w.terminate(); resolve({ engine: "worker (no WebGPU)", ...m.done }); }); }
      else if (m.error) { w.terminate(); reject(new Error(m.error)); }
    };
    w.onerror = (e) => { w.terminate(); reject(new Error(e.message || "the worker failed")); };
    w.postMessage({ ...grid, source: [...src, q], hours: HOURS, frameS: FRAME_S });
  });
}

/** Flood a reach: t = { river_id, lat, lon, order }; rp = 2, 10 or 100. */
export async function floodIt(t, { rp = 100, by = "rules" } = {}) {
  const my = ++running;
  const id = `flood:${t.river_id}`;
  const base = { id, lngLat: [t.lon, t.lat], lift: 14, what: "Flood it", title: `River reach ${t.river_id}`,
    sub: `${LABEL}`, credit: FLOODSIM_CREDIT, buttons: rpButtons(t, rp, true) };
  openCard({ ...base, status: { pending: "Reading this reach's flood flows…" } });
  try {
    const flows = await readFlows(t);
    const q = rp === 2 ? flows.q2 : rp === 10 ? flows.q10 : flows.q100;
    if (!(q > 0)) throw new Error(`no ${rp}-year flow here`);
    if (my !== running) return null;
    updateCard(id, { status: { pending: `Reading the terrain for ${fmt(q)} m³/s…` } });
    const win = windowTiles(t.lon, t.lat, ZOOM, SIZE);
    const z = await readTerrain(win);
    const lines = riverLines(t);
    // navigator.gpu can be there with no adapter behind it (headless, blocklisted GPUs): ask for one.
    const gpu = gpuAvailable() && Boolean(await navigator.gpu.requestAdapter().catch(() => null));
    // Without WebGPU the worker runs a grid of half the size: the terrain averaged 2 x 2, the rivers burnt again.
    const cwin = gpu ? win : { ...win, size: SIZE / 2, zoom: ZOOM - 1, x0: win.x0 / 2, y0: win.y0 / 2, dx: win.dx * 2 };
    const zz = gpu ? z : coarsen(z, SIZE);
    const [px, py] = worldPx(t.lon, t.lat, cwin.zoom);
    const clickCell = [Math.floor(py - cwin.y0), Math.floor(px - cwin.x0)];
    const prep = prepareGrid(zz, cwin, lines, flows.q2, clickCell);
    const grid = { z: prep.z, n: prep.n, nx: cwin.size, ny: cwin.size, dx: cwin.dx };
    const frames = [];
    updateCard(id, { status: { pending: `Flooding on ${gpu ? "your GPU" : "your CPU"}…` } });
    const onFrame = async (time, h) => {
      if (my !== running) return;
      frames.push({ t: time, url: await frameUrl(h, cwin.size), stats: summary(h, cwin.dx) });
      draw({ frames, bounds: win.bounds }, frames.length - 1);
      updateCard(id, { status: { pending: `Flooding… ${hm(time)} of ${HOURS} h` } });
    };
    let res = gpu ? await runGpu(grid, prep.source, q, onFrame).catch(() => null) : null;
    if (!res) res = await runWorker(grid, prep.source, q, onFrame);
    if (my !== running) return null;
    const last = frames.at(-1).stats;
    const value = { target: t, rp, flow: q, frames, bounds: win.bounds, engine: res.engine, steps: res.steps,
      wall_s: res.wall_s, grid: cwin.size, dx: cwin.dx, channel: prep.channel, lines: lines.length, stats: last,
      max: summary(res.hMax, cwin.dx) };
    const before = slotValue("floodsim");
    show(value);
    if (state.mapOk) {   // the whole run in view, the card's anchor still on the reach
      const [w, s, e, n] = win.bounds;
      map.fitBounds([[w, s], [e, n]], { padding: 40, duration: 600 });
    }
    logAction({ action: { type: "flood_reach", river_id: t.river_id, return_period: rp, flow_cms: q }, before,
      after: value, label: `Flood it: ${rp}-year flow`, by, detail: `${fmt(q)} m³/s for ${HOURS} h` });
    announce(`Flood it: ${fmt(value.max.wet_km2, 1)} square kilometres wet after ${HOURS} hours`);
    return value;
  } catch (err) {
    if (my === running) updateCard(id, { status: { text: `Could not flood this reach: ${err.message}`, muted: true }, buttons: rpButtons(t, rp, false) });
    return null;
  }
}

function rpButtons(t, rp, busy) {
  return [2, 10, 100].map((T) => ({
    id: `rp${T}`, label: `${T}-year`, pressed: T === rp, disabled: busy,
    title: `Pour this reach's ${T}-year flow in for ${HOURS} hours`, onClick: () => floodIt(t, { rp: T }),
  })).concat([{ id: "clear", label: "Clear", hidden: busy, title: "Take the flood off the map", onClick: clearByUndo }]);
}

// ── drawing ─────────────────────────────────────────────────────────────────

function draw(v, i) {
  if (!state.mapOk || !v.frames[i]) return;
  const [w, s, e, n] = v.bounds;
  const url = v.frames[i].url;
  const src = map.getSource(SRC);
  if (src && src.updateImage) {
    if (src.url !== url) src.updateImage({ url, coordinates: [[w, n], [e, n], [e, s], [w, s]] });
    return;
  }
  map.addSource(SRC, { type: "image", url, coordinates: [[w, n], [e, n], [e, s], [w, s]] });
  // Under the rivers and the gauges, as the other depth rasters: at the top of the stack on the globe, an image
  // source darkened the whole map in Chrome 155.
  map.addLayer({ id: SRC, type: "raster", source: SRC,
    paint: { "raster-opacity": 0.92, "raster-fade-duration": 0, "raster-resampling": "linear" } }, beforeId());
}

function beforeId() {
  const layers = (map.getStyle() && map.getStyle().layers) || [];
  const hit = layers.find((l) => /^river-/.test(l.id) || ["catchment-fill", "gauge-heat", "clusters", "points"].includes(l.id));
  return hit ? hit.id : undefined;
}

function clearMap() {
  if (!state.mapOk) return;
  if (map.getLayer(SRC)) map.removeLayer(SRC);
  if (map.getSource(SRC)) map.removeSource(SRC);
}

function card(v, i) {
  const f = v.frames[i];
  const t = v.target;
  updateCard(`flood:${t.river_id}`, {
    status: { text: `${fmt(v.flow)} m³/s, the ${v.rp}-year flow, poured in for ${HOURS} h. ${hm(f.t)} in: ` +
      `${fmt(f.stats.wet_km2, 2)} km² under water, deepest ${fmt(f.stats.max_m, 1)} m.` },
    figure: { value: fmt(v.max.wet_km2, 2), unit: " km²", label: `wet at some time in ${HOURS} h` },
    note: `${v.engine}: ${fmt(v.steps)} steps on ${v.grid} x ${v.grid} cells of ${fmt(v.dx)} m in ` +
      `${fmt(v.wall_s, 1)} s. Channel burnt ${fmt(v.channel.depth_m, 1)} m deep along ${v.lines} river lines ` +
      "(Manning, the 2-year flow).",
    buttons: rpButtons(t, v.rp, false),
  });
}

function show(v) {
  shown = v;
  const last = v.frames.length - 1;
  draw(v, last);
  card(v, last);
  useModelClock({
    times: v.frames.map((f) => f.t), index: last, label: (i) => `+${hm(v.frames[i].t)}`,
    summary: `Flood it: ${HOURS} h of the ${v.rp}-year flow, model time.`,
    onStep: (i) => { if (shown === v) { draw(v, i); card(v, i); } },
    onExit: () => { if (shown === v) { draw(v, last); card(v, last); } },
  });
}

function paint(v) {
  if (v) {
    if (state.mapOk) openCard({ id: `flood:${v.target.river_id}`, lngLat: [v.target.lon, v.target.lat], lift: 14,
      what: "Flood it", title: `River reach ${v.target.river_id}`, sub: LABEL, credit: FLOODSIM_CREDIT });
    show(v);
  } else {
    shown = null;
    running++;
    clearMap();
    releaseModelClock();
  }
}

function clearByUndo() {
  const before = slotValue("floodsim");
  if (!before) return;
  paint(null);
  logAction({ action: { type: "flood_reach", river_id: before.target.river_id, return_period: null }, before,
    after: null, label: "Flood it cleared" });
}

export function initFloodIt() {
  actions.openFloodIt = (t) => floodIt(t);
  registerSlotWriter("floodsim", (value) => paint(value || null));
  globalThis.__aq = globalThis.__aq || {};
  globalThis.__aq.floodIt = { run: floodIt, shown: () => shown, prepareGrid };
}
