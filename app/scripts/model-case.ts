/**
 * Разбор варианта «материал 2, слабое предсжатие, удар — рывок — удар» гибридным методом:
 * журнал точных задач (что встретилось, какие волны родились), траектории фронтов и поля
 * для статьи docs/model/model.tex.
 *
 * Запуск (из app/):  node scripts/model-case.ts ../docs/model/figures/case.json
 */
import { writeFileSync } from 'node:fs';
import { normalize } from '../src/core/material.ts';
import { MATERIAL_2 } from '../src/core/presets.ts';
import { type HybridResult, hybridStateAt, solveHybrid } from '../src/core/hybrid.ts';
import { traction } from '../src/core/material1d.ts';

const m = normalize(MATERIAL_2);
const F0: [number, number, number] = [-0.0005, 0, 0.001];
const pulses = [
  { start: 0, duration: 0.25, d: 0.002 },
  { start: 0.25, duration: 0.25, d: -0.002 },
  { start: 0.5, duration: 0.25, d: 0.002 },
];
const round = (x: number) => Number(x.toPrecision(7));
const st = (s?: { F: number[]; v: number[] }) => (s ? { p: round(s.F[0]), s: round(s.F[1]), v1: round(s.v[0]), v3: round(s.v[1]) } : undefined);

/** Скорости характеристик в состоянии (p, s). */
function speeds(p: number, s: number) {
  const T = (a: number, b: number) => {
    const t = traction(m, [a, 0, b]);
    return [t[0], t[2]];
  };
  const h = 1e-9;
  const a = (T(p + h, s)[0] - T(p - h, s)[0]) / (2 * h);
  const b = (T(p, s + h)[0] - T(p, s - h)[0]) / (2 * h);
  const d = (T(p, s + h)[1] - T(p, s - h)[1]) / (2 * h);
  const D = Math.sqrt(((a - d) / 2) ** 2 + b * b);
  return { cf: Math.sqrt((a + d) / 2 + D), cs: Math.sqrt((a + d) / 2 - D) };
}

/** Поле p − p0 и s − s0 на сетке (t, x) по кадрам. */
function field(r: HybridResult, nx: number, xMax: number) {
  const P: number[][] = [];
  const S: number[][] = [];
  for (const fr of r.frames) {
    const rowP: number[] = [];
    const rowS: number[] = [];
    for (let i = 0; i < nx; i++) {
      const x = (xMax * (i + 0.5)) / nx;
      const q = hybridStateAt(fr, F0, x);
      rowP.push(round(q.F[0] - F0[0]));
      rowS.push(round(q.F[1] - F0[2]));
    }
    P.push(rowP);
    S.push(rowS);
  }
  return { t: r.frames.map((f) => round(f.t)), xMax, nx, P, S };
}

const full = solveHybrid({ material: m, F0, pulses, tEnd: 10, frames: 401, untilQuiet: true, cellsPerL: 400, maxCells: 2400 });
const zoom = solveHybrid({ material: m, F0, pulses, tEnd: 1.6, length: 1.8, frames: 321, cellsPerL: 800 });
const profT = [0.2, 0.45, 0.7, 1.0, 1.5];
const prof = solveHybrid({ material: m, F0, pulses, tEnd: 1.5, length: 1.8, frames: 31, cellsPerL: 800 });
const profiles = profT.map((t) => {
  const fr = prof.frames.reduce((a, b) => (Math.abs(b.t - t) < Math.abs(a.t - t) ? b : a));
  const xs: number[] = [];
  const p: number[] = [];
  const s: number[] = [];
  const v1: number[] = [];
  const v3: number[] = [];
  for (let i = 0; i <= 900; i++) {
    const x = (1.8 * i) / 900;
    const q = hybridStateAt(fr, F0, x);
    xs.push(round(x));
    p.push(round(q.F[0]));
    s.push(round(q.F[1]));
    v1.push(round(q.v[0]));
    v3.push(round(q.v[1]));
  }
  return { t: round(fr.t), x: xs, p, s, v1, v3, fronts: fr.fronts.map((f) => ({ x: round(f.x), fam: f.fam, kind: f.kind, G: round(f.G) })) };
});

const summary = (r: HybridResult) => ({
  tEnd: r.tEnd,
  xView: r.xView,
  N: r.N,
  steps: r.steps,
  events: r.events,
  untracked: r.untracked,
  maxFronts: r.maxFronts,
  kinematicError: r.kinematicError,
  work: r.work,
  energyGain: r.energyGain,
  shockDissipation: r.shockDissipation,
});

const out = {
  material: m,
  F0,
  pulses,
  speedCurve: Array.from({ length: 161 }, (_, i) => {
    const p = -0.004 + (0.008 * i) / 160;
    return { p: round(p), ...speeds(p, F0[2]) };
  }),
  speeds: {
    p0: speeds(F0[0], F0[2]),
    compressed: speeds(-0.00237, 0.00094),
    tension: speeds(0.00158, 0.00125),
  },
  full: summary(full),
  zoomSummary: summary(zoom),
  log: full.log.map((e) => ({
    t: round(e.t),
    x: round(e.x),
    type: e.type,
    d: e.d,
    left: st(e.left),
    right: st(e.right),
    incoming: e.incoming?.map((w) => ({ ...w, G: round(w.G), strength: round(w.strength) })),
    waves: e.waves?.map((w) => ({ fam: w.fam, dir: w.dir, kind: w.kind, G: round(w.G), strength: round(w.strength), tracked: w.tracked, L: st(w.L), R: st(w.R) })),
  })),
  segmentsFull: full.segments.map((s) => ({ ...s, x0: round(s.x0), t0: round(s.t0), x1: round(s.x1), t1: round(s.t1) })),
  segmentsZoom: zoom.segments.map((s) => ({ ...s, x0: round(s.x0), t0: round(s.t0), x1: round(s.x1), t1: round(s.t1) })),
  fieldFull: field(full, 360, full.xView),
  fieldZoom: field(zoom, 360, 1.8),
  profiles,
};
const path = process.argv[2];
if (path) writeFileSync(path, JSON.stringify(out));
else console.log(JSON.stringify({ ...out, fieldFull: undefined, fieldZoom: undefined, profiles: undefined }, null, 1));
