import { describe, expect, it } from 'vitest';
import { stateAt, trackFronts } from '../src/core/fronts.ts';
import { hybridStateAt, solveHybrid } from '../src/core/hybrid.ts';
import { profileAt, solveImpact } from '../src/core/impact.ts';
import { normalize } from '../src/core/material.ts';
import { MATERIAL_1, MATERIAL_2 } from '../src/core/presets.ts';

const m1 = normalize(MATERIAL_1);
const m2 = normalize(MATERIAL_2);
const F0_1: [number, number, number] = [-0.00005, 0, 0.001];
const F0_2: [number, number, number] = [-0.005, 0, 0.001];

/** L1-отличие профиля p (и s) гибрида в момент последнего кадра от функции exact(x). */
function l1(res: ReturnType<typeof solveHybrid>, F0: [number, number, number], exact: (x: number) => { F: number[] }, len: number) {
  const fr = res.frames.at(-1)!;
  const K = 20000;
  let ep = 0;
  let es = 0;
  for (let i = 0; i < K; i++) {
    const x = (len * (i + 0.5)) / K;
    const a = hybridStateAt(fr, F0, x);
    const b = exact(x);
    ep += (Math.abs(a.F[0] - b.F[0]) * len) / K;
    es += (Math.abs(a.F[1] - b.F[1]) * len) / K;
  }
  return { ep, es };
}

describe('гибридный метод (ударные волны выделены, вееры на сетке)', () => {
  it('удар сжатием: обе ударные волны совпадают с автомодельным решением §2.3 на любой сетке', () => {
    const ex = solveImpact({ material: m2, F0: F0_2, bc: { type: 'velocity', d: 0.002 } })!;
    const r = solveHybrid({ material: m2, F0: F0_2, pulses: [{ start: 0, duration: 10, d: 0.002 }], tEnd: 0.8, frames: 2, cellsPerL: 100, length: 1 });
    const fr = r.frames.at(-1)!;
    expect(fr.fronts.map((f) => f.kind)).toEqual(['shock', 'shock']);
    expect(fr.fronts[0].G).toBeCloseTo(ex.waves[1].xiFront, 9);
    expect(fr.fronts[1].G).toBeCloseTo(ex.waves[0].xiFront, 9);
    const { ep, es } = l1(r, F0_2, (x) => { const q = profileAt(ex, x / fr.t); return { F: [q.F[0], q.F[2]] }; }, 1);
    expect(ep).toBeLessThan(1e-9);
    expect(es).toBeLessThan(1e-12);
    expect(Math.abs(r.kinematicError)).toBeLessThan(1e-14);
    expect(r.untracked).toBe(0);
  });

  it('широкий веер (растяжение) считается на сетке и сходится к точному решению при измельчении', () => {
    const d = -0.0045;
    const ex = solveImpact({ material: m2, F0: F0_2, bc: { type: 'velocity', d } })!;
    expect(ex.waves.map((w) => w.kind)).toEqual(['fan', 'fan']);
    const err = (cpl: number) => {
      const r = solveHybrid({ material: m2, F0: F0_2, pulses: [{ start: 0, duration: 10, d }], tEnd: 1.5, frames: 2, cellsPerL: cpl, length: 1.8 });
      // на мелкой сетке вееры шире её размазывания — разрывов не выделяется, решение гладкое
      if (cpl >= 400) expect(r.frames.at(-1)!.fronts).toHaveLength(0);
      // баланс точен, кроме момента, когда тонкий веер отдаётся сетке (узкая область ~h/2)
      expect(Math.abs(r.kinematicError)).toBeLessThan(1e-9);
      const t = r.frames.at(-1)!.t;
      return l1(r, F0_2, (x) => { const q = profileAt(ex, x / t); return { F: [q.F[0], q.F[2]] }; }, 1.8).ep;
    };
    const e1 = err(100);
    const e2 = err(400);
    expect(e2).toBeLessThan(0.5 * e1);
    // L1 по p при скачке p ≈ 4e-3 на участке 1.8: ~1 % — погрешность у краёв веера и точки его рождения
    expect(e2).toBeLessThan(1e-4);
  });

  it('тонкий веер (уже ячейки) ведётся одним разрывом с погрешностью порядка своей ширины', () => {
    const d = -0.002;
    const ex = solveImpact({ material: m2, F0: F0_2, bc: { type: 'velocity', d } })!;
    const r = solveHybrid({ material: m2, F0: F0_2, pulses: [{ start: 0, duration: 10, d }], tEnd: 1, frames: 2, cellsPerL: 200, length: 1.2 });
    const fr = r.frames.at(-1)!;
    expect(fr.fronts.map((f) => f.kind)).toEqual(['simple', 'simple']);
    // скорость разрыва — внутри веера
    expect(fr.fronts[1].G).toBeLessThanOrEqual(ex.waves[0].xiFront);
    expect(fr.fronts[1].G).toBeGreaterThanOrEqual(ex.waves[0].xiBack!);
    const { ep } = l1(r, F0_2, (x) => { const q = profileAt(ex, x / fr.t); return { F: [q.F[0], q.F[2]] }; }, 1.2);
    expect(ep).toBeLessThan(0.002 * (ex.waves[0].xiFront - ex.waves[0].xiBack!) * fr.t);
  });

  for (const [name, m, F0, tol] of [
    ['материал 1', m1, F0_1, 1e-5],
    ['материал 2', m2, F0_2, 3e-4],
  ] as const)
    it(`${name}: серия из трёх ударов совпадает с методом фронтов, балансы выполнены`, () => {
      const pulses = [
        { start: 0, duration: 0.08, d: 0.002 },
        { start: 0.25, duration: 0.08, d: 0.002 },
        { start: 0.5, duration: 0.08, d: 0.002 },
      ];
      const ref = trackFronts({ material: m, F0, pulses, tEnd: 0.92, frames: 5, untilQuiet: true, fanStep: 6.25e-5 });
      const h = solveHybrid({ material: m, F0, pulses, tEnd: ref.tEnd, length: ref.xView, frames: 5, cellsPerL: 200 });
      expect(h.warnings).toEqual([]);
      expect(h.untracked).toBe(0);
      let e = 0;
      let mx = 0;
      const K = 4000;
      for (const k of [2, 4])
        for (let i = 0; i < K; i++) {
          const x = (ref.xView * (i + 0.5)) / K;
          const b = stateAt(ref.frames[k], x);
          e += Math.abs(hybridStateAt(h.frames[k], F0, x).F[0] - b.F[0]);
          mx = Math.max(mx, Math.abs(b.F[0] - F0[0]));
        }
      expect(e / (2 * K) / mx).toBeLessThan(tol);
      expect(Math.abs(h.kinematicError)).toBeLessThan(1e-11);
      // работа границы = прирост энергии + диссипация на ударных волнах; численная — ничтожна
      expect(Math.abs(h.work - h.energyGain - h.shockDissipation)).toBeLessThan(1e-6 * h.work);
      // число точных задач — десятки-сотни, а не тысячи встреч лучей
      expect(h.events).toBeLessThan(400);
    });
});
