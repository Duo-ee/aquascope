// Flood it (#658), the fallback without WebGPU: the same local-inertial scheme (floodsim-core.js createCpuSim)
// in a module worker, so the page stays responsive. It gets { z, n, nx, ny, dx, source: [row, col, q], hours,
// frameS } and posts { frame: { t, h } } every frameS of simulated time, then { done: { steps, wall_s, volume,
// hMax } }, or { error }.

import { createCpuSim } from "./floodsim-core.js?v=__BUILD__";

self.onmessage = (e) => {
  const job = e.data;
  try {
    const sim = createCpuSim({ z: job.z, n: job.n, nx: job.nx, ny: job.ny, dx: job.dx }, { sources: [job.source] });
    const end = job.hours * 3600;
    const started = performance.now();
    self.postMessage({ frame: { t: 0, h: new Float32Array(sim.h) } });
    for (let t = job.frameS; t <= end + 1e-6; t += job.frameS) {
      sim.advance(Math.min(t, end));
      self.postMessage({ frame: { t: sim.state.t, h: new Float32Array(sim.h) } });
    }
    self.postMessage({ done: { steps: sim.state.steps, wall_s: (performance.now() - started) / 1000,
      volume: sim.volume(), hMax: new Float32Array(sim.hMax) } });
  } catch (err) {
    self.postMessage({ error: String((err && err.message) || err) });
  }
};
