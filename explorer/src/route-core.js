// Route it (#656), the pure part: no DOM, no map, so node can test it.
// Every number comes from Python (aquascope.routing.route: Muskingum down GEOGLOWS's own network, with its k and x
// per reach). The frames arrive as whole-number shares of the peak that entered (0 to 100) for each reach and
// hour; this module only decides the colours, the widths, the words and which frame is the peak.

// The wave's colour by its share of the peak that entered: pale where little arrives, deep where it is near
// full strength (a sequential blue-purple ramp that reads on the light and dark basemaps alike).
export const RAMP = [
  [0, "#c6dbef"],
  [10, "#9ecae1"],
  [25, "#6baed6"],
  [45, "#4292c6"],
  [65, "#6a51a3"],
  [85, "#54278f"],
  [100, "#3f007d"],
];

const finite = (v) => v !== null && v !== undefined && v !== "" && Number.isFinite(Number(v));

/** The line colour expression: the reach's feature-state "rw" (its share, 0 to 100) on RAMP. */
export function waveColorExpr() {
  const stops = RAMP.flatMap(([v, c]) => [v, c]);
  return ["interpolate", ["linear"], ["coalesce", ["feature-state", "rw"], 0], ...stops];
}

/** The line width: wider as more of the wave is there, and wider as you zoom in. */
export function waveWidthExpr() {
  const share = ["coalesce", ["feature-state", "rw"], 0];
  const at = (k) => ["interpolate", ["linear"], share, 0, 0.8 * k, 20, 2.4 * k, 100, 7 * k];
  return ["interpolate", ["linear"], ["zoom"], 3, at(0.6), 7, at(1), 12, at(1.6)];
}

/** Reaches with no water yet are drawn faintly, so the path is there before the wave is. */
export function waveOpacityExpr() {
  return ["case", [">", ["coalesce", ["feature-state", "rw"], 0], 0], 0.95, 0.35];
}

/** The filter that keeps the path's reaches: a match on riverId, so the tiles' own lines are coloured. */
export function pathFilterExpr(ids) {
  const list = [...new Set((ids || []).map(Number).filter(Number.isFinite))];
  if (!list.length) return ["boolean", false];
  return ["match", ["to-number", ["get", "riverId"]], list, true, false];
}

/** "6 h", "1 d 4 h", "12 days": the time since the water started to enter. */
export function hoursWords(h) {
  if (!finite(h)) return "n/a";
  const n = Number(h);
  if (n < 48) return `${Math.round(n)} h`;
  const d = Math.floor(n / 24), r = Math.round(n - 24 * d);
  if (n >= 240) return `${(n / 24).toFixed(1).replace(/\.0$/, "")} days`;
  return r ? `${d} d ${r} h` : `${d} days`;
}

/** The clock's label for a frame: "+14 h", "+2 d 6 h". */
export function frameLabel(tH) {
  return `+${hoursWords(tH)}`;
}

/** A flow in m³/s, rounded as the card says it: "1,240", "38", "4.2". */
export function flowWords(v) {
  if (!finite(v)) return "n/a";
  const n = Number(v);
  if (n >= 100) return Math.round(n).toLocaleString("en-US");
  if (n >= 10) return n.toFixed(0);
  return n.toFixed(1).replace(/\.0$/, "");
}

/** "−45%" (lower than where it started), "0%". */
export function flatterWords(p) {
  if (!finite(p)) return "";
  const n = Math.round(Number(p));
  return n === 0 ? "0%" : `−${Math.abs(n)}%`;
}

/** The frame where the most of the wave is on the map at once (the most total share), for the still picture. */
export function busiestFrame(frames) {
  const share = (frames && frames.share) || [];
  let best = 0, top = -1;
  share.forEach((row, i) => {
    const s = row.reduce((a, b) => a + b, 0);
    if (s > top) { top = s; best = i; }
  });
  return best;
}

/** Each reach's highest share over the run: the peak everywhere, the picture of how the wave flattens. */
export function peakShares(frames) {
  const ids = (frames && frames.ids) || [];
  const out = new Array(ids.length).fill(0);
  for (const row of (frames && frames.share) || []) row.forEach((v, j) => { if (v > out[j]) out[j] = v; });
  return out;
}

/** The feature-states for one frame: [[river_id, share], ...]; `peak: true` gives each reach's peak instead. */
export function frameStates(frames, i, { peak = false } = {}) {
  const ids = (frames && frames.ids) || [];
  const row = peak ? peakShares(frames) : ((frames && frames.share && frames.share[i]) || []);
  return ids.map((id, j) => [id, row[j] || 0]);
}

/** The places' labels on the map: one line each, "Cologne · peak 3 d 6 h · 1,240 m³/s". */
export function placeFeatures(places) {
  const feats = [];
  for (const [n, p] of (places || []).entries()) {
    if (!finite(p.lat) || !finite(p.lon)) continue;
    const reached = p.reached && finite(p.peak_h);
    const text = reached ? `${p.name} · peak ${hoursWords(p.peak_h)} · ${flowWords(p.peak_cms)} m³/s`
      : `${p.name} · ${p.reached ? "still rising" : "too small to matter"}`;
    feats.push({ type: "Feature", geometry: { type: "Point", coordinates: [Number(p.lon), Number(p.lat)] },
      properties: { n, kind: p.kind === "gauge" ? "gauge" : "town", text, rank: p.kind === "town" ? 0 : 1 } });
  }
  return { type: "FeatureCollection", features: feats };
}

/** The places the card lists: the towns first (people), then the gauges, each in the order the water reaches them,
 * at most `max`. */
export function cardPlaces(places, max = 8) {
  const all = (places || []).filter((p) => p.reached && finite(p.peak_h));
  const towns = all.filter((p) => p.kind !== "gauge");
  const gauges = all.filter((p) => p.kind === "gauge");
  const pick = [...towns.slice(0, max), ...gauges.slice(0, Math.max(0, max - towns.length))].slice(0, max);
  return pick.sort((a, b) => (a.along_km || 0) - (b.along_km || 0));
}

/** The result's first line: what entered and how far its peak got. */
export function headline(res) {
  if (!res || !res.ok) return "";
  const i = res.inflow || {};
  const f = res.furthest || {};
  const what = i.shape === "block" ? `${flowWords(i.peak_cms)} m³/s released for ${hoursWords(i.duration_h)}`
    : `A flood peaking at ${flowWords(i.peak_cms)} m³/s over ${hoursWords(i.duration_h)}`;
  if (!f.index) return `${what}: it does not get past the first reach in ${hoursWords(res.hours)}.`;
  const km = finite(f.along_km) ? `${Math.round(f.along_km).toLocaleString("en-US")} km down` : "downstream";
  return `${what}. ${hoursWords(f.peak_h)} later its peak is ${km} at ${flowWords(f.peak_cms)} m³/s, ` +
    `${flatterWords(f.attenuation_pct)} on the start.`;
}

/** The action log's line for a run. */
export function logLabel(res) {
  const i = (res && res.inflow) || {};
  if (finite(i.return_period)) return `Route the ${Number(i.return_period)}-year flood`;
  return `Route ${flowWords(i.peak_cms)} m³/s${i.shape === "block" ? " released" : ""}`;
}

/** The [[w, s], [e, n]] box of the path's line, or null. */
export function pathBounds(res) {
  const c = (res && res.geometry && res.geometry.coordinates) || [];
  if (c.length < 2) return null;
  let w = Infinity, s = Infinity, e = -Infinity, n = -Infinity;
  for (const [x, y] of c) { w = Math.min(w, x); e = Math.max(e, x); s = Math.min(s, y); n = Math.max(n, y); }
  return [[w, s], [e, n]];
}

/** The part of the path the wave's peak reached in the hours routed, for framing (the rest is the way on). */
export function reachedBounds(res) {
  const c = (res && res.geometry && res.geometry.coordinates) || [];
  const rows = (res && res.reaches) || [];
  const f = res && res.furthest;
  const end = f && finite(f.index) ? rows[f.index] : null;
  if (!end || !finite(end.along_km) || c.length < 2) return pathBounds(res);
  let km = 0;
  const part = [c[0]];
  for (let k = 1; k < c.length; k++) {
    const [x0, y0] = c[k - 1], [x1, y1] = c[k];
    const dx = (x1 - x0) * 111.32 * Math.cos(((y0 + y1) / 2) * Math.PI / 180), dy = (y1 - y0) * 110.57;
    km += Math.hypot(dx, dy);
    part.push(c[k]);
    if (km >= end.along_km) break;
  }
  return pathBounds({ geometry: { coordinates: part } });
}

/** The inputs the panel sends, checked: { peak_cms } or { return_period }, duration_h, shape. Throws on a bad one. */
export function routeArgs({ kind, peak, rp, duration, block }) {
  const d = Number(duration);
  if (!(d > 0 && d <= 720)) throw new Error("the duration must be between 0 and 720 hours");
  const out = { duration_h: d, shape: block ? "block" : "triangle" };
  if (kind === "peak") {
    const q = Number(peak);
    if (!(q > 0 && q <= 1e6)) throw new Error("give the peak flow in m³/s");
    out.peak_cms = q;
  } else {
    const t = Number(rp);
    if (!(t >= 1.01 && t <= 1000)) throw new Error("the return period must be between 1.01 and 1000 years");
    out.return_period = t;
  }
  return out;
}

/** The towns in the basemap's own label tiles (OpenMapTiles "place"): [{ name, lat, lon, kind: "town" }]. */
export function townsFrom(features) {
  const seen = new Set();
  const out = [];
  for (const f of features || []) {
    const p = f.properties || {};
    if (!["city", "town"].includes(p.class)) continue;
    const name = p["name:latin"] || p.name_en || p.name;
    const g = f.geometry;
    if (!name || !g || g.type !== "Point") continue;
    if (seen.has(name)) continue;
    seen.add(name);
    out.push({ name: String(name), lat: g.coordinates[1], lon: g.coordinates[0], kind: "town" });
  }
  return out;
}
