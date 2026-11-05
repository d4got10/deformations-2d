/**
 * Сверка с одноосной задачей статьи Лаптевой, Дудко, Рагозиной (2022, sources/lapteva-dudko-ragozina-2022-cyclic-uniaxial.pdf).
 * Модель Ляховского–Мясникова = наша при α = β = 0; λ = 19.8, μ = 26.9, ν = 15.3 ГПа, ρ = 7600 кг/м³.
 * Единицы: время — мс, длина — L = c1·1 мс (c1 = √((λ+2μ)/ρ) = 3111.9 м/с), скорость — c1.
 * Запуск: node scripts/uniaxial-case.ts
 */
import { normalize } from '../src/core/material.ts';
import { solveBimodular } from '../src/core/bimodular1d.ts';
import { stateAt, trackFronts } from '../src/core/fronts.ts';
import { hybridStateAt, solveHybrid } from '../src/core/hybrid.ts';

const GPa = { lambda: 19.8, mu: 26.9, nu: 15.3, alpha: 0, beta: 0 };
const rho = 7600;
const c1 = Math.sqrt(((GPa.lambda + 2 * GPa.mu) * 1e9) / rho);
const m = normalize(GPa);
const a = Math.sqrt(1 + 2 * m.nu);
const b = Math.sqrt(1 - 2 * m.nu);
console.log('c1 =', c1.toFixed(2), 'м/с; a =', (a * c1).toFixed(2), 'b =', (b * c1).toFixed(2), 'м/с');

/** узлы u0(t): t в мс, u0 в мм (×1e-3 м); последний наклон продолжается */
function pulsesFrom(ts: number[], us: number[]) {
  const times: number[] = [];
  const vel: number[] = [];
  for (let j = 1; j < ts.length; j++) {
    times.push(ts[j - 1]);
    vel.push((us[j] - us[j - 1]) / (ts[j] - ts[j - 1]) / c1); // мм/мс = м/с
  }
  // наш знак: d > 0 — граница вдавливается (u0 растёт вдоль x) — совпадает со знаком u̇0
  const pulses = times.map((t, k) => ({ start: t, duration: (k + 1 < times.length ? times[k + 1] : 1e3) - t, d: vel[k] }));
  return { times, vel, pulses };
}

const cases: Record<string, { ts: number[]; us: number[]; tEnd: number; probe: number[] }> = {
  'I(a)': { ts: [0, 0.2, 0.5, 0.75, 0.9], us: [0, -3, 0, 4, 0], tEnd: 2.5, probe: [0.76, 1.2, 1.9, 2.5] },
  'I(b)': { ts: [0, 0.2, 0.5, 0.75, 0.9], us: [0, -3, -1, 4, 0], tEnd: 3.0, probe: [1.0, 1.4, 1.7, 3.0] },
  II: { ts: [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3], us: [0, -0.3, 0.4, -0.3, 0.4, -0.3, 0.4, 0], tEnd: 2.8, probe: [0.4, 2.8] },
};
for (const [name, cs] of Object.entries(cases)) {
  const { times, vel, pulses } = pulsesFrom(cs.ts, cs.us);
  console.log(`\n=== ${name}: скорости границы, м/с:`, vel.map((v) => (v * c1).toFixed(3)).join(', '));
  const ex = solveBimodular({ a, b, times, velocities: vel, tEnd: cs.tEnd });
  for (const e of ex.events) console.log(`  t=${e.t.toFixed(4)} x=${e.x.toFixed(4)} ${e.type}: ${e.incoming.join(' + ')} → ${e.outgoing.join(' ')}`);
  const xMax = (a * cs.tEnd) * 1.02;
  const fr = trackFronts({ material: m, F0: [0, 0, 0], pulses, tEnd: cs.tEnd, frames: 2, xMax: xMax * 1.1 });
  const hy = solveHybrid({ material: m, F0: [0, 0, 0], pulses, tEnd: cs.tEnd, frames: 2, length: xMax, cellsPerL: 400 / xMax * 1, maxCells: 3200 });
  const err = (f: (x: number) => number) => {
    let e = 0;
    let n = 0;
    const K = 20000;
    for (let i = 0; i < K; i++) {
      const x = (xMax * (i + 0.5)) / K;
      const q = ex.stateAt(x, cs.tEnd).e;
      e += Math.abs(f(x) - q);
      n += Math.abs(q);
    }
    return e / n;
  };
  console.log(`  метод фронтов: событий ${fr.events}, предупреждений ${fr.warnings.length}, отн. L1 по e: ${err((x) => stateAt(fr.frames.at(-1)!, x).F[0]).toExponential(2)}`);
  console.log(`  гибрид: событий ${hy.events}, предупреждений ${hy.warnings.length}, отн. L1 по e: ${err((x) => hybridStateAt(hy.frames.at(-1)!, [0, 0, 0], x).F[0]).toExponential(2)}, N = ${hy.N}`);
  console.log('  гибрид:', hy.log.map((e) => `${e.type}@${e.t.toFixed(4)}`).join(' '));
}

// ---------- данные для рисунка статьи (docs/model/figures/uniaxial.json)
if (process.argv[2]) {
  const { writeFileSync } = await import('node:fs');
  const out: Record<string, unknown> = { a, b, c1 };
  for (const [name, cs] of Object.entries(cases)) {
    const { times, vel, pulses } = pulsesFrom(cs.ts, cs.us);
    const ex = solveBimodular({ a, b, times, velocities: vel, tEnd: cs.tEnd });
    const X = a * cs.tEnd;
    const nt = 400;
    const nx = 300;
    const zone: number[][] = [];
    for (let j = 0; j < nx; j++) {
      const row: number[] = [];
      for (let i = 0; i < nt; i++) {
        const t = (cs.tEnd * (i + 0.5)) / nt;
        const x = (X * (j + 0.5)) / nx;
        {
          const q = ex.stateAt(x, t);
          // e = 0 позади ведущего фронта — жёсткий слой (он может и покоиться), впереди — невозмущённая среда
          const lead = Math.max(0, ...ex.segments.filter((g) => g.t0 <= t && t <= g.t1).map((g) => g.x0 + ((g.x1 - g.x0) * (t - g.t0)) / (g.t1 - g.t0 || 1)));
          row.push(q.e < -1e-12 ? -1 : q.e > 1e-12 ? 1 : x < lead ? 2 : 0);
        }
      }
      zone.push(row);
    }
    // проверка нашего гибрида на этом же варианте — отрезки его выделенных фронтов
    const hy = solveHybrid({ material: m, F0: [0, 0, 0], pulses, tEnd: cs.tEnd, frames: 2, length: X * 1.02, cellsPerL: 400 / X });
    out[name] = { tEnd: cs.tEnd, X, zone, segments: ex.segments, events: ex.events, hybridSegments: hy.segments };
  }
  const tau: { u: number; t1: number; t2: number }[] = [];
  for (let k = 0; k <= 120; k++) {
    const u2 = -2.25 + (6.25 * k) / 120;
    const { times, vel } = pulsesFrom([0, 0.2, 0.5, 0.75, 0.9], [0, -3, u2, 4, 0]);
    // время инцидентов 1 и 2 в первичной картине (аналитически: Σ догоняет γᵇ; ξᵃ догоняет Σ)
    const v1 = vel[0];
    const e1 = -v1 / b;
    const dv = vel[1] - v1;
    const B2 = (a * a - b * b) * e1;
    const G = (-B2 + Math.sqrt(B2 * B2 + 4 * dv * dv * a * a)) / (2 * dv);
    void times;
    tau.push({ u: u2, t1: (0.2 * G) / (G - b), t2: (0.5 * a - 0.2 * G) / (a - G) });
  }
  out.tau = tau;
  writeFileSync(process.argv[2], JSON.stringify(out));
  console.log('записано', process.argv[2]);
}
