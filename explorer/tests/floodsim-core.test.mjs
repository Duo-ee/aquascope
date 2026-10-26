// Flood it (#658): the JS twin of aquascope.floodsim runs the same grid as the Python fixture and must match it,
// and the WGSL keeps the same scheme's pieces.
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";

import {
  FLOODSIM_WGSL, burnChannel, channelDepth, coarsen, createCpuSim, depthRgba, lineCells, pixelSizeM,
  terrariumElevation, windowTiles, worldPx,
} from "../src/floodsim-core.js";

const fx = JSON.parse(readFileSync(new URL("../../tests/fixtures/floodsim_twin.json", import.meta.url), "utf8"));

test("the JS solver matches the numpy solver step for step", () => {
  const ny = fx.z.length, nx = fx.z[0].length;
  const sim = createCpuSim({ z: fx.z.flat(), n: fx.n.flat(), nx, ny, dx: fx.dx }, { sources: [fx.source] });
  sim.advance(fx.duration_s);
  assert.equal(sim.state.steps, fx.steps);
  const want = fx.h.flat();
  let worst = 0;
  for (let i = 0; i < want.length; i++) worst = Math.max(worst, Math.abs(sim.h[i] - want[i]));
  // Round-off differs (pow, the order of sums), and a wet/dry face can amplify it a little: within a micrometre.
  assert.ok(worst < 1e-6, `largest difference ${worst}`);
  const v = sim.volume();
  assert.ok(Math.abs(v.error_m3) / v.inflow_m3 < 1e-12);
  assert.ok(Math.abs(v.outflow_m3 - fx.volume.outflow_m3) / fx.volume.outflow_m3 < 1e-8);
});

test("still water stays still in JS too", () => {
  const nx = 12, ny = 10, z = [], h0 = [];
  for (let i = 0; i < nx * ny; i++) { const zi = (i * 7919) % 13 / 10; z.push(zi); h0.push(Math.max(2 - zi, 0)); }
  for (let c = 0; c < nx; c++) { z[c] = 5; h0[c] = 0; z[(ny - 1) * nx + c] = 5; h0[(ny - 1) * nx + c] = 0; }
  for (let r = 0; r < ny; r++) { z[r * nx] = 5; h0[r * nx] = 0; z[r * nx + nx - 1] = 5; h0[r * nx + nx - 1] = 0; }
  const sim = createCpuSim({ z, n: 0.03, nx, ny, dx: 10, h0 });
  sim.advance(300);
  for (let i = 0; i < h0.length; i++) assert.ok(Math.abs(sim.h[i] - h0[i]) < 1e-12);
});

test("terrarium, tiles and pixel size", () => {
  assert.deepEqual([...terrariumElevation([128, 0, 0, 255, 128, 10, 128, 255])], [0, 10.5]);
  assert.ok(Math.abs(pixelSizeM(0, 0) - 40075016.686 / 256) < 1e-6);
  const w = windowTiles(-2.2, 52.06, 12, 512);
  assert.ok(w.tiles.length >= 4 && w.tiles.length <= 9);
  assert.match(w.tiles[0].url, /^https:\/\/elevation-tiles-prod\.s3\.amazonaws\.com\/terrarium\/12\/\d+\/\d+\.png$/);
  const [x, y] = worldPx(-2.2, 52.06, 12);
  assert.equal(Math.floor(y - w.y0), 256);
  assert.equal(Math.floor(x - w.x0), 256);
  assert.ok(w.bounds[0] < -2.2 && w.bounds[2] > -2.2 && w.bounds[1] < 52.06 && w.bounds[3] > 52.06);
});

test("channel depth, burning, cells and coarsening", () => {
  assert.ok(Math.abs(channelDepth(100, 25, 1e-3) - ((0.035 * 100) / (25 * Math.sqrt(1e-3))) ** 0.6) < 1e-12);
  assert.equal(channelDepth(100, 25, 0), channelDepth(100, 25, 1e-4));
  const z = [10, 9, 9.5, 8, 8.2];
  assert.deepEqual(burnChannel(z, 5, [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4]], 2), [8, 7, 7, 6, 6]);
  const win = { zoom: 12, size: 512, ...(() => { const [a, b] = worldPx(-2.2, 52.06, 12); return { x0: Math.floor(a) - 256, y0: Math.floor(b) - 256 }; })() };
  const cells = lineCells(win, [[-2.2, 52.06], [-2.19, 52.06]]);
  assert.equal(cells[0][0], 256);
  assert.ok(cells.length > 10 && cells.every(([r]) => r === 256));
  assert.deepEqual([...coarsen([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], 4)], [3.5, 5.5, 11.5, 13.5]);
});

test("depth colours follow the JRC ramp's stops", () => {
  assert.deepEqual(depthRgba(0.01), [0, 0, 0, 0]);
  assert.deepEqual(depthRgba(1), [106, 174, 232, Math.round(255 * 0.74)]);
  assert.deepEqual(depthRgba(50), [23, 63, 153, Math.round(255 * 0.93)]);
});

test("the WGSL has every pass and the same coefficients' slots", () => {
  for (const fn of ["fn stepDt", "fn faces", "fn limit", "fn scale", "fn depth", "pow(hf, 7.0 / 3.0)", "atomicMax"]) {
    assert.ok(FLOODSIM_WGSL.includes(fn), fn);
  }
  // At most eight storage buffers a stage (the WebGPU default limit).
  assert.ok((FLOODSIM_WGSL.match(/var<storage/g) || []).length <= 8);
});
