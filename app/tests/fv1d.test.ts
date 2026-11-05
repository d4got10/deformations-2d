import { describe, expect, it } from 'vitest';
import { type Pulse, boundaryVelocity, cell, simulate } from '../src/core/fv1d.ts';
import { profileAt, solveImpact } from '../src/core/impact.ts';
import { normalize } from '../src/core/material.ts';
import type { Vec3 } from '../src/core/material1d.ts';
import { MATERIAL_1, MATERIAL_2 } from '../src/core/presets.ts';

const CASES = [
  ['материал 1', normalize(MATERIAL_1), [-0.00005, 0, 0.001] as Vec3],
  ['материал 2', normalize(MATERIAL_2), [-0.005, 0, 0.001] as Vec3],
] as const;

/** Средняя по отрезку x ∈ [0, xmax] ошибка p относительно точного автомодельного решения, начавшегося в t0. */
function meanError(res: ReturnType<typeof simulate>, exact: NonNullable<ReturnType<typeof solveImpact>>, t0: number, xmax: number) {
  const fr = res.frames.at(-1)!;
  let e = 0;
  let n = 0;
  for (let i = 0; i < res.N; i++) {
    const x = (i + 0.5) * res.dx;
    if (x > xmax) break;
    e += Math.abs(cell(fr, res.N, i).F[0] - profileAt(exact, x / (fr.t - t0)).F[0]);
    n++;
  }
  return e / n;
}

describe('график скорости границы', () => {
  it('сумма активных импульсов, полуинтервалы [start, start + duration)', () => {
    const p: Pulse[] = [
      { start: 0, duration: 0.2, d: 0.001 },
      { start: 0.1, duration: 0.3, d: 0.002 },
    ];
    expect(boundaryVelocity(p, 0.05)).toBe(0.001);
    expect(boundaryVelocity(p, 0.15)).toBe(0.003);
    expect(boundaryVelocity(p, 0.2)).toBe(0.002);
    expect(boundaryVelocity(p, 0.4)).toBe(0);
  });
});

describe('численный решатель (конечные объёмы) против точного решения главы 2', () => {
  for (const [name, m, F0] of CASES) {
    describe(name, () => {
      const single = solveImpact({ material: m, F0, bc: { type: 'velocity', d: 0.001 } })!;

      it('один удар: сходимость к автомодельному решению при измельчении сетки', () => {
        const run = (N: number) => simulate({ material: m, F0, pulses: [{ start: 0, duration: 10, d: 0.001 }], N, tEnd: 0.5, frames: 2 });
        const e1 = meanError(run(250), single, 0, 1);
        const e2 = meanError(run(500), single, 0, 1);
        expect(e2).toBeLessThan(0.8 * e1);
        expect(e2).toBeLessThan(0.01 * Math.abs(single.zones.at(-1)!.F[0] - F0[0]));
      });

      it('кинематический баланс ∫(p − p0)dx = −∫v*dt выполняется точно', () => {
        const r = simulate({
          material: m,
          F0,
          pulses: [
            { start: 0, duration: 0.15, d: 0.002 },
            { start: 0.3, duration: 0.1, d: 0.003 },
          ],
          N: 300,
          tEnd: 0.6,
          frames: 3,
        });
        expect(Math.abs(r.kinematicError)).toBeLessThan(1e-15);
        // перемещение границы = ∫v*dt
        expect(r.frames.at(-1)!.ub).toBeCloseTo(0.002 * 0.15 + 0.003 * 0.1, 14);
      });

      it('энергия: работа границы = прирост энергии + рассеяние ≥ 0; рассеяние численное — убывает с сеткой', () => {
        const frac = (N: number) => {
          const r = simulate({ material: m, F0, pulses: [{ start: 0, duration: 0.2, d: 0.002 }], N, tEnd: 0.5, frames: 2 });
          expect(r.work).toBeGreaterThanOrEqual(r.energyGain);
          return (r.work - r.energyGain) / r.work;
        };
        const f1 = frac(200);
        const f2 = frac(400);
        expect(f2).toBeLessThan(0.8 * f1);
        expect(f2).toBeLessThan(0.08);
      });

      const zII = single.zones.at(-1)!;
      it('второй удар по движущейся сжатой среде: совпадение с точным решением (v0 ≠ 0) до встречи волн', () => {
        const r = simulate({
          material: m,
          F0,
          pulses: [
            { start: 0, duration: 10, d: 0.001 },
            { start: 0.3, duration: 10, d: 0.002 },
          ],
          N: 800,
          tEnd: 0.42,
          frames: 2,
        });
        const ex = solveImpact({ material: m, F0: zII.F, v0: zII.v, bc: { type: 'velocity', d: 0.003 } })!;
        expect(meanError(r, ex, 0.3, 0.2)).toBeLessThan(0.03 * Math.abs(ex.zones.at(-1)!.F[0] - zII.F[0]));
      });

      it('остановка границы (разгрузка): совпадение с точным решением', () => {
        const r = simulate({ material: m, F0, pulses: [{ start: 0, duration: 0.3, d: 0.001 }], N: 800, tEnd: 0.42, frames: 2 });
        const ex = solveImpact({ material: m, F0: zII.F, v0: zII.v, bc: { type: 'velocity', d: 0 } })!;
        expect(meanError(r, ex, 0.3, 0.2)).toBeLessThan(0.03 * Math.abs(zII.F[0] - F0[0]));
      });
    });
  }

  it('разномодульность: у материала 1 повторное сжатие — веер, разгрузка — составная волна «ударная + веер»; у материала 2 — ударная и веер', () => {
    const solve = (m: (typeof CASES)[number][1], F0: Vec3) => {
      const z = solveImpact({ material: m, F0, bc: { type: 'velocity', d: 0.001 } })!.zones.at(-1)!;
      const up = solveImpact({ material: m, F0: z.F, v0: z.v, bc: { type: 'velocity', d: 0.003 } })!;
      const down = solveImpact({ material: m, F0: z.F, v0: z.v, bc: { type: 'velocity', d: 0 } })!;
      return { up, down };
    };
    const a = solve(CASES[0][1], CASES[0][2]);
    expect(a.up.waves[0].kind).toBe('fan');
    // одна ударная волна разгрузки нарушала бы (2.45) (D < 0): решение — ударная волна, к тылу
    // которой примыкает веер; в профиле волны это скачок (две точки с одним ξ), затем гладкий участок
    const w = a.down.waves[0];
    expect(w.kind).toBe('fan');
    expect(w.profile![1].xi).toBe(w.profile![0].xi);
    expect(Math.abs(w.profile![1].F[0] - w.profile![0].F[0])).toBeGreaterThan(1e-5);
    expect(a.down.warnings.some((s) => s.includes('(2.45)') || s.includes('Лакс'))).toBe(false);
    const b = solve(CASES[1][1], CASES[1][2]);
    expect([b.up.waves[0].kind, b.down.waves[0].kind]).toEqual(['shock', 'fan']);
  });
});
