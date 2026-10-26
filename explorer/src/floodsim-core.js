// Flood it (#658), the pure part: no DOM, no map, so node can test it and the fallback worker can import it.
//
// The local-inertial scheme of Bates, Horritt and Fewtrell (2010) in the q-centred form of de Almeida et al.
// (2012), the same steps as aquascope/floodsim.py (a node test runs both on one grid and compares them):
// the face flow depth, the semi-implicit Manning friction, a flux limiter that keeps depths non-negative and the
// volume exact, a free outfall at the edges, and the CFL step dt = alpha dx / sqrt(g h_max). The WGSL in
// FLOODSIM_WGSL is the same scheme for the GPU (floodsim-gpu.js), in f32.
//
// Also here: the Terrarium decode, the channel burnt along the river lines, its depth by Manning's equation,
// and the depth colours.

export const G = 9.81;
export const ALPHA = 0.7;
export const THETA = 0.7;
export const H_DRY = 1e-3;
export const N_CHANNEL = 0.035;
export const N_FLOODPLAIN = 0.06;
export const SLOPE_RANGE = [1e-4, 1e-2];
export const TERRARIUM = "https://elevation-tiles-prod.s3.amazonaws.com/terrarium/{z}/{x}/{y}.png";
export const EARTH_CIRCUMFERENCE_M = 40075016.686;
export const LABEL = "screening flood: no bathymetry, terrain errors of metres, a constant inflow";
export const FLOODSIM_CREDIT = "Local-inertial scheme (Bates et al. 2010; de Almeida et al. 2012) run on your device; " +
  "terrain: Terrarium tiles (Mapzen/Tilezen joerd: SRTM, GMTED2010, ETOPO1 and others); rivers: GEOGLOWS v2 " +
  "(TDX-Hydro, CC BY-SA 4.0); flows: GEOGLOWS v2 retrospective, CC BY 4.0. A screening flood, not a design flood.";

// ── terrain ─────────────────────────────────────────────────────────────────

/** Terrarium RGB(A) bytes -> metres: R*256 + G + B/256 - 32768. `rgba` has `channels` bytes a pixel. */
export function terrariumElevation(rgba, channels = 4) {
  const n = Math.floor(rgba.length / channels);
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const k = i * channels;
    out[i] = rgba[k] * 256 + rgba[k + 1] + rgba[k + 2] / 256 - 32768;
  }
  return out;
}

export const pixelSizeM = (lat, zoom) => EARTH_CIRCUMFERENCE_M * Math.cos((lat * Math.PI) / 180) / (256 * 2 ** zoom);

export function worldPx(lon, lat, zoom) {
  const s = 256 * 2 ** zoom;
  const r = (lat * Math.PI) / 180;
  return [((lon + 180) / 360) * s, ((1 - Math.log(Math.tan(r) + 1 / Math.cos(r)) / Math.PI) / 2) * s];
}

export function pxLonLat(x, y, zoom) {
  const s = 256 * 2 ** zoom;
  return [(x / s) * 360 - 180, (Math.atan(Math.sinh(Math.PI * (1 - (2 * y) / s))) * 180) / Math.PI];
}

/** The tiles a size x size window centred on (lon, lat) needs, and where the window sits in their mosaic. */
export function windowTiles(lon, lat, zoom, size) {
  const [cx, cy] = worldPx(lon, lat, zoom);
  const x0 = Math.round(cx - size / 2), y0 = Math.round(cy - size / 2);
  const tx0 = Math.floor(x0 / 256), ty0 = Math.floor(y0 / 256);
  const tx1 = Math.floor((x0 + size - 1) / 256), ty1 = Math.floor((y0 + size - 1) / 256);
  const n = 2 ** zoom;
  const tiles = [];
  for (let ty = ty0; ty <= ty1; ty++) {
    for (let tx = tx0; tx <= tx1; tx++) {
      tiles.push({ x: ((tx % n) + n) % n, y: ty, col: tx - tx0, row: ty - ty0,
        url: TERRARIUM.replace("{z}", zoom).replace("{x}", ((tx % n) + n) % n).replace("{y}", ty) });
    }
  }
  const [w, nn] = pxLonLat(x0, y0, zoom);
  const [e, s] = pxLonLat(x0 + size, y0 + size, zoom);
  return { x0, y0, zoom, size, tiles, cols: tx1 - tx0 + 1, rows: ty1 - ty0 + 1, ox: x0 - tx0 * 256, oy: y0 - ty0 * 256,
    bounds: [w, s, e, nn], dx: pixelSizeM(lat, zoom) };
}

/** Halve a grid (mean of 2x2 blocks): the fallback runs on a coarser grid. */
export function coarsen(z, size) {
  const m = size >> 1;
  const out = new Float64Array(m * m);
  for (let r = 0; r < m; r++) {
    for (let c = 0; c < m; c++) {
      const i = 2 * r * size + 2 * c;
      out[r * m + c] = (z[i] + z[i + 1] + z[i + size] + z[i + size + 1]) / 4;
    }
  }
  return out;
}

// ── the channel ─────────────────────────────────────────────────────────────

/** Manning normal depth of a wide rectangular channel: d = (n Q / (W sqrt(S)))^(3/5), S kept in SLOPE_RANGE. */
export function channelDepth(qBankfull, widthM, slope, n = N_CHANNEL) {
  const s = Math.min(Math.max(slope, SLOPE_RANGE[0]), SLOPE_RANGE[1]);
  return ((n * qBankfull) / (widthM * Math.sqrt(s))) ** 0.6;
}

/** The cells (row, col) a [lon, lat] polyline crosses, in order, sampled every half cell (as floodsim._line_cells). */
export function lineCells(win, line) {
  const pts = line.map(([lo, la]) => worldPx(lo, la, win.zoom));
  const out = [];
  const seen = (r, c) => out.slice(-3).some(([a, b]) => a === r && b === c);
  for (let i = 0; i + 1 < pts.length; i++) {
    const [ax, ay] = pts[i], [bx, by] = pts[i + 1];
    const k = Math.max(1, Math.ceil(2 * Math.hypot(bx - ax, by - ay)));
    for (let j = 0; j <= k; j++) {
      const r = Math.floor(ay + ((by - ay) * j) / k - win.y0);
      const c = Math.floor(ax + ((bx - ax) * j) / k - win.x0);
      if (r >= 0 && r < win.size && c >= 0 && c < win.size && !seen(r, c)) out.push([r, c]);
    }
  }
  return out;
}

/** Cells ordered with the higher end first (upstream). */
export function orient(z, size, cells) {
  if (cells.length > 1 && z[cells[0][0] * size + cells[0][1]] < z[cells.at(-1)[0] * size + cells.at(-1)[1]]) return cells.slice().reverse();
  return cells;
}

/** Burn a channel depthM deep along ordered cells (in place): the bank made non-increasing downstream, minus depthM. */
export function burnChannel(z, size, cells, depthM) {
  let bank = Infinity;
  for (const [r, c] of cells) {
    const i = r * size + c;
    bank = Math.min(bank, z[i]);
    z[i] = Math.min(z[i], bank - depthM);
  }
  return z;
}

/** The slope that sizes the channel: the network's bank drop (5th to 95th percentile) over its extent. */
export function networkSlope(z, size, cells, dx) {
  if (cells.length < 3) return SLOPE_RANGE[0];
  const zs = cells.map(([r, c]) => z[r * size + c]).sort((a, b) => a - b);
  const pct = (p) => {
    const k = (zs.length - 1) * p, lo = Math.floor(k), hi = Math.ceil(k);
    return zs[lo] + (zs[hi] - zs[lo]) * (k - lo);
  };
  const rs = cells.map((c) => c[0]), cs = cells.map((c) => c[1]);
  const span = Math.hypot(Math.max(...rs) - Math.min(...rs), Math.max(...cs) - Math.min(...cs)) * dx;
  return (pct(0.95) - pct(0.05)) / Math.max(span, dx);
}

// ── the solver (CPU) ────────────────────────────────────────────────────────

/**
 * A run on the CPU (the fallback, and the twin the node test checks against Python).
 * grid: { z (Float64Array ny*nx, rows north to south), nx, ny, dx, n (array or number) }
 * opts: { durationS, sources: [[row, col, q]], framesEveryS, theta, alpha, dtMax, onFrame(t, h), now }
 */
export function createCpuSim(grid, opts = {}) {
  const { nx, ny, dx } = grid;
  const z = Float64Array.from(grid.z);
  const N = nx * ny;
  const n = typeof grid.n === "number" || grid.n === undefined
    ? new Float64Array(N).fill(grid.n === undefined ? N_FLOODPLAIN : grid.n) : Float64Array.from(grid.n);
  const theta = opts.theta ?? THETA, alpha = opts.alpha ?? ALPHA, dtMax = opts.dtMax ?? 10;
  const h = grid.h0 ? Float64Array.from(grid.h0) : new Float64Array(N);
  let qx = new Float64Array(ny * (nx + 1)), qy = new Float64Array((ny + 1) * nx);
  let qxn = new Float64Array(qx.length), qyn = new Float64Array(qy.length);
  const f = new Float64Array(N);
  const hMax = Float64Array.from(h);
  const sources = (opts.sources || []).map(([r, c, q]) => [r * nx + c, q]);
  const st = { t: 0, steps: 0, vIn: 0, vOut: 0, v0: h.reduce((a, b) => a + b, 0) * dx * dx };

  function face(zl, zr, el, er, nl, nr, qOld, qLo, qHi, dt) {
    const hf = Math.max(el, er) - Math.max(zl, zr);
    if (!(hf > H_DRY)) return 0;
    const qc = theta * qOld + 0.5 * (1 - theta) * (qLo + qHi);
    const nf = 0.5 * (nl + nr);
    return (qc - (G * hf * dt * (er - el)) / dx) / (1 + (G * dt * nf * nf * Math.abs(qOld)) / hf ** (7 / 3));
  }

  function step(dt) {
    // x faces: row r, face k (0..nx): between cell k-1 and k; ghosts at the edges have the edge bed, no water.
    for (let r = 0; r < ny; r++) {
      const row = r * nx, fr = r * (nx + 1);
      for (let k = 0; k <= nx; k++) {
        const il = row + Math.max(k - 1, 0), ir = row + Math.min(k, nx - 1);
        const zl = z[il], zr = z[ir];
        const el = k === 0 ? zl : zl + h[il];
        const er = k === nx ? zr : zr + h[ir];
        const qo = qx[fr + k];
        const lo = k === 0 ? qo : qx[fr + k - 1], hi = k === nx ? qo : qx[fr + k + 1];
        let q = face(zl, zr, el, er, n[il], n[ir], qo, lo, hi, dt);
        if (k === 0) q = Math.min(q, 0); else if (k === nx) q = Math.max(q, 0);
        qxn[fr + k] = q;
      }
    }
    // y faces: face k (0..ny) between row k-1 and k, column c.
    for (let k = 0; k <= ny; k++) {
      for (let c = 0; c < nx; c++) {
        const it = Math.max(k - 1, 0) * nx + c, ib = Math.min(k, ny - 1) * nx + c;
        const zt = z[it], zb = z[ib];
        const et = k === 0 ? zt : zt + h[it];
        const eb = k === ny ? zb : zb + h[ib];
        const qo = qy[k * nx + c];
        const lo = k === 0 ? qo : qy[(k - 1) * nx + c], hi = k === ny ? qo : qy[(k + 1) * nx + c];
        let q = face(zt, zb, et, eb, n[it], n[ib], qo, lo, hi, dt);
        if (k === 0) q = Math.min(q, 0); else if (k === ny) q = Math.max(q, 0);
        qyn[k * nx + c] = q;
      }
    }
    // The limiter: each cell's outgoing faces scaled so that it cannot go below zero.
    for (let r = 0; r < ny; r++) {
      for (let c = 0; c < nx; c++) {
        const i = r * nx + c, fr = r * (nx + 1);
        const out = (Math.max(qxn[fr + c + 1], 0) + Math.max(-qxn[fr + c], 0)
          + Math.max(qyn[(r + 1) * nx + c], 0) + Math.max(-qyn[r * nx + c], 0)) * dt / dx;
        f[i] = out > h[i] ? h[i] / out : 1;
      }
    }
    for (let r = 0; r < ny; r++) {
      const fr = r * (nx + 1);
      for (let k = 0; k <= nx; k++) {
        const q = qxn[fr + k];
        qxn[fr + k] = q * (q > 0 ? f[r * nx + Math.max(k - 1, 0)] : f[r * nx + Math.min(k, nx - 1)]);
      }
    }
    for (let k = 0; k <= ny; k++) {
      for (let c = 0; c < nx; c++) {
        const q = qyn[k * nx + c];
        qyn[k * nx + c] = q * (q > 0 ? f[Math.max(k - 1, 0) * nx + c] : f[Math.min(k, ny - 1) * nx + c]);
      }
    }
    let out = 0;
    for (let r = 0; r < ny; r++) out += qxn[r * (nx + 1) + nx] - qxn[r * (nx + 1)];
    for (let c = 0; c < nx; c++) out += qyn[ny * nx + c] - qyn[c];
    st.vOut += out * dx * dt;
    for (let r = 0; r < ny; r++) {
      for (let c = 0; c < nx; c++) {
        const i = r * nx + c, fr = r * (nx + 1);
        const v = h[i] + (dt / dx) * (qxn[fr + c] - qxn[fr + c + 1] + qyn[r * nx + c] - qyn[(r + 1) * nx + c]);
        h[i] = v > 0 ? v : 0;
      }
    }
    for (const [i, q] of sources) { h[i] += (q * dt) / (dx * dx); st.vIn += q * dt; }
    for (let i = 0; i < N; i++) if (h[i] > hMax[i]) hMax[i] = h[i];
    [qx, qxn] = [qxn, qx];
    [qy, qyn] = [qyn, qy];
  }

  function dtNow(remaining) {
    let hm = 0;
    for (let i = 0; i < N; i++) if (h[i] > hm) hm = h[i];
    return Math.min(dtMax, remaining, hm > H_DRY ? (alpha * dx) / Math.sqrt(G * hm) : dtMax);
  }

  /** Advance to time `until` (s), or by at most `maxSteps`; returns the state. */
  function advance(until, maxSteps = Infinity) {
    let k = 0;
    while (st.t < until - 1e-9 && k < maxSteps) {
      const dt = dtNow(until - st.t);
      step(dt);
      st.t += dt;
      st.steps++;
      k++;
    }
    return st;
  }

  function volume() {
    let v = 0;
    for (let i = 0; i < N; i++) v += h[i];
    const final = v * dx * dx;
    return { initial_m3: st.v0, inflow_m3: st.vIn, outflow_m3: st.vOut, final_m3: final,
      error_m3: final - (st.v0 + st.vIn - st.vOut) };
  }

  return { h, hMax, state: st, advance, volume, nx, ny, dx };
}

// ── the drawing ─────────────────────────────────────────────────────────────

// Depth (m) -> colour, the same stops as the JRC depth map (flood-depth-core.js RAMP), so the two read alike.
export const RAMP = [
  { depth_m: 0.05, rgb: [191, 224, 247], alpha: 0.55 },
  { depth_m: 1.0, rgb: [106, 174, 232], alpha: 0.74 },
  { depth_m: 3.0, rgb: [47, 124, 204], alpha: 0.86 },
  { depth_m: 10.0, rgb: [23, 63, 153], alpha: 0.93 },
];

/** RGBA bytes for one depth (0 alpha when shallower than the first stop). */
export function depthRgba(d) {
  if (!(d >= RAMP[0].depth_m)) return [0, 0, 0, 0];
  for (let i = 1; i < RAMP.length; i++) {
    if (d <= RAMP[i].depth_m) {
      const a = RAMP[i - 1], b = RAMP[i];
      const t = (d - a.depth_m) / (b.depth_m - a.depth_m);
      const mix = (u, v) => Math.round(u + (v - u) * t);
      return [mix(a.rgb[0], b.rgb[0]), mix(a.rgb[1], b.rgb[1]), mix(a.rgb[2], b.rgb[2]), Math.round(255 * (a.alpha + (b.alpha - a.alpha) * t))];
    }
  }
  const l = RAMP.at(-1);
  return [...l.rgb, Math.round(255 * l.alpha)];
}

/** A depth grid as RGBA pixels (ImageData's layout). */
export function depthImage(h, nx, ny) {
  const px = new Uint8ClampedArray(nx * ny * 4);
  for (let i = 0; i < nx * ny; i++) {
    const [r, g, b, a] = depthRgba(h[i]);
    px[4 * i] = r; px[4 * i + 1] = g; px[4 * i + 2] = b; px[4 * i + 3] = a;
  }
  return px;
}

/** The wet area (km2) and the deepest depth of a depth grid. */
export function summary(h, dx, wet = 0.05) {
  let cells = 0, deepest = 0;
  for (const v of h) { if (v > wet) cells++; if (v > deepest) deepest = v; }
  return { wet_km2: (cells * dx * dx) / 1e6, max_m: deepest, wet_cells: cells };
}

// ── the GPU scheme (WGSL) ───────────────────────────────────────────────────
// Seven storage buffers (the default limit is eight a stage), f32 unless said: zn (the bed, then Manning's n),
// h and hmax (cells), lim (the limiter, cells), qO and qN (the x faces, ny*(nx+1), then the y faces, (ny+1)*nx: the
// old and the new discharges, swapped every step through two bind groups), and st: { dt, t, tEnd, hmaxBits (an
// atomic u32, the largest depth's bits: positive floats order like their bits), n (the steps taken) }. Per step: stepDt (one thread),
// faces (both directions), limit, scale, depth. Ghost cells, the limiter and the edge rule as in createCpuSim.

export const FLOODSIM_WGSL = /* wgsl */ `
struct P { nx: u32, ny: u32, src: u32, pad: u32, dx: f32, theta: f32, alpha: f32, dtMax: f32, srcQ: f32, g: f32, hDry: f32, pad2: f32 };
struct S { dt: f32, t: f32, tEnd: f32, hmaxBits: atomic<u32>, n: u32 };
@group(0) @binding(0) var<uniform> p: P;
@group(0) @binding(1) var<storage, read_write> st: S;
@group(0) @binding(2) var<storage, read> zn: array<f32>;
@group(0) @binding(3) var<storage, read_write> h: array<f32>;
@group(0) @binding(4) var<storage, read_write> hmax: array<f32>;
@group(0) @binding(5) var<storage, read_write> lim: array<f32>;
@group(0) @binding(6) var<storage, read> qO: array<f32>;
@group(0) @binding(7) var<storage, read_write> qN: array<f32>;

fn ncells() -> u32 { return p.nx * p.ny; }
fn nqx() -> u32 { return (p.nx + 1u) * p.ny; }

@compute @workgroup_size(1)
fn stepDt() {
  let hm = bitcast<f32>(atomicLoad(&st.hmaxBits));
  var dt = p.dtMax;
  if (hm > p.hDry) { dt = min(dt, p.alpha * p.dx / sqrt(p.g * hm)); }
  dt = max(min(dt, st.tEnd - st.t), 0.0);
  st.dt = dt;
  st.t = st.t + dt;
  if (dt > 0.0) { st.n = st.n + 1u; }
  atomicStore(&st.hmaxBits, 0u);
}

fn face(zl: f32, zr: f32, el: f32, er: f32, nl: f32, nr: f32, qo: f32, lo: f32, hi: f32, dt: f32) -> f32 {
  let hf = max(el, er) - max(zl, zr);
  if (!(hf > p.hDry)) { return 0.0; }
  let qc = p.theta * qo + 0.5 * (1.0 - p.theta) * (lo + hi);
  let nf = 0.5 * (nl + nr);
  return (qc - p.g * hf * dt * (er - el) / p.dx) / (1.0 + p.g * dt * nf * nf * abs(qo) / pow(hf, 7.0 / 3.0));
}

@compute @workgroup_size(16, 16)
fn faces(@builtin(global_invocation_id) id: vec3<u32>) {
  let nx = p.nx; let ny = p.ny; let dt = st.dt; let N = ncells(); let X = nqx();
  let k = id.x; let r = id.y;
  if (r < ny && k <= nx) {            // x face k (0..nx) of row r
    let il = r * nx + select(k - 1u, 0u, k == 0u);
    let ir = r * nx + min(k, nx - 1u);
    let zl = zn[il]; let zr = zn[ir];
    let el = select(zl + h[il], zl, k == 0u);
    let er = select(zr + h[ir], zr, k == nx);
    let fi = r * (nx + 1u) + k;
    let qo = qO[fi];
    let lo = select(qO[max(fi, 1u) - 1u], qo, k == 0u);
    let hi = select(qO[min(fi + 1u, X - 1u)], qo, k == nx);
    var q = face(zl, zr, el, er, zn[N + il], zn[N + ir], qo, lo, hi, dt);
    if (k == 0u) { q = min(q, 0.0); } else if (k == nx) { q = max(q, 0.0); }
    qN[fi] = q;
  }
  let c = k; let kk = r;
  if (c < nx && kk <= ny) {           // y face kk (0..ny) of column c
    let it = select(kk - 1u, 0u, kk == 0u) * nx + c;
    let ib = min(kk, ny - 1u) * nx + c;
    let zt = zn[it]; let zb = zn[ib];
    let et = select(zt + h[it], zt, kk == 0u);
    let eb = select(zb + h[ib], zb, kk == ny);
    let fi = kk * nx + c;
    let qo = qO[X + fi];
    let lo = select(qO[X + max(fi, nx) - nx], qo, kk == 0u);
    let hi = select(qO[X + min(fi + nx, (ny + 1u) * nx - 1u)], qo, kk == ny);
    var q = face(zt, zb, et, eb, zn[N + it], zn[N + ib], qo, lo, hi, dt);
    if (kk == 0u) { q = min(q, 0.0); } else if (kk == ny) { q = max(q, 0.0); }
    qN[X + fi] = q;
  }
}

@compute @workgroup_size(16, 16)
fn limit(@builtin(global_invocation_id) id: vec3<u32>) {
  let c = id.x; let r = id.y; let nx = p.nx; let X = nqx();
  if (c >= nx || r >= p.ny) { return; }
  let i = r * nx + c; let fr = r * (nx + 1u);
  let o = (max(qN[fr + c + 1u], 0.0) + max(-qN[fr + c], 0.0) + max(qN[X + (r + 1u) * nx + c], 0.0) + max(-qN[X + r * nx + c], 0.0)) * st.dt / p.dx;
  lim[i] = select(1.0, h[i] / o, o > h[i]);
}

@compute @workgroup_size(16, 16)
fn scale(@builtin(global_invocation_id) id: vec3<u32>) {
  let nx = p.nx; let ny = p.ny; let k = id.x; let r = id.y; let X = nqx();
  if (r < ny && k <= nx) {
    let fi = r * (nx + 1u) + k;
    let q = qN[fi];
    let donor = select(r * nx + min(k, nx - 1u), r * nx + select(k - 1u, 0u, k == 0u), q > 0.0);
    qN[fi] = q * lim[donor];
  }
  let c = k; let kk = r;
  if (c < nx && kk <= ny) {
    let fi = X + kk * nx + c;
    let q = qN[fi];
    let donor = select(min(kk, ny - 1u) * nx + c, select(kk - 1u, 0u, kk == 0u) * nx + c, q > 0.0);
    qN[fi] = q * lim[donor];
  }
}

var<workgroup> wmax: atomic<u32>;

@compute @workgroup_size(16, 16)
fn depth(@builtin(global_invocation_id) id: vec3<u32>, @builtin(local_invocation_index) li: u32) {
  if (li == 0u) { atomicStore(&wmax, 0u); }
  workgroupBarrier();
  let c = id.x; let r = id.y; let nx = p.nx; let X = nqx();
  if (c < nx && r < p.ny) {
    let i = r * nx + c; let fr = r * (nx + 1u); let dt = st.dt;
    var v = h[i] + dt / p.dx * (qN[fr + c] - qN[fr + c + 1u] + qN[X + r * nx + c] - qN[X + (r + 1u) * nx + c]);
    v = max(v, 0.0);
    if (i == p.src) { v = v + p.srcQ * dt / (p.dx * p.dx); }
    h[i] = v;
    hmax[i] = max(hmax[i], v);
    atomicMax(&wmax, bitcast<u32>(v));
  }
  workgroupBarrier();
  if (li == 0u) { atomicMax(&st.hmaxBits, atomicLoad(&wmax)); }
}
`;
