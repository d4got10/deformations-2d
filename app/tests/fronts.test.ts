import { describe, expect, it } from 'vitest';
import { boundaryProblem, riemannProblem, stateAt, trackFronts } from '../src/core/fronts.ts';
import { cell, simulate } from '../src/core/fv1d.ts';
import { solveImpact } from '../src/core/impact.ts';
import { normalize } from '../src/core/material.ts';
import { traction } from '../src/core/material1d.ts';
import { MATERIAL_1, MATERIAL_2 } from '../src/core/presets.ts';

const CASES = [
  ['материал 1', normalize(MATERIAL_1), -0.00005],
  ['материал 2', normalize(MATERIAL_2), -0.005],
] as const;
const PULSES = [
  { start: 0, duration: 0.08, d: 0.002 },
  { start: 0.25, duration: 0.08, d: 0.002 },
  { start: 0.5, duration: 0.08, d: 0.002 },
];

describe('метод волновых фронтов', () => {
  for (const [name, m, p0] of CASES) {
    describe(name, () => {
      const c = { m };
      const a = { F: [p0, 0.001] as [number, number], v: [0, 0] as [number, number] };

      it('краевая задача (удар) совпадает с автомодельным решением §2.3', () => {
        const bp = boundaryProblem(c, a, 0.001, 2.5e-4)!;
        const ex = solveImpact({ material: m, F0: [p0, 0, 0.001], bc: { type: 'velocity', d: 0.001 } })!;
        expect(bp.fast.pieces).toHaveLength(1);
        expect(bp.slow.pieces).toHaveLength(1);
        expect(bp.fast.pieces[0].G).toBeCloseTo(ex.waves[0].xiFront, 10);
        expect(bp.slow.pieces[0].G).toBeCloseTo(ex.waves[1].xiFront, 9);
        const zb = ex.zones.at(-1)!;
        expect(bp.slow.behind.F[0]).toBeCloseTo(zb.F[0], 13);
        expect(bp.slow.behind.F[1]).toBeCloseTo(zb.F[2], 13);
      });

      it('остановка границы (разгрузка) совпадает с автомодельным решением с v0 ≠ 0', () => {
        const z = solveImpact({ material: m, F0: [p0, 0, 0.001], bc: { type: 'velocity', d: 0.002 } })!.zones.at(-1)!;
        const ex = solveImpact({ material: m, F0: z.F, v0: z.v, bc: { type: 'velocity', d: 0 } })!;
        const bp = boundaryProblem(c, { F: [z.F[0], z.F[2]], v: [z.v[0], z.v[2]] }, 0, 2.5e-4)!;
        const zb = ex.zones.at(-1)!;
        // для простой волны (веер) пучок слабых фронтов даёт конечное состояние с погрешностью O(ε·δ²)
        expect(Math.abs(bp.slow.behind.F[0] - zb.F[0])).toBeLessThan(1e-8);
        expect(Math.abs(bp.slow.behind.v[0])).toBeLessThan(1e-12);
        // тип волны тот же: у материала 1 — ударная, у материала 2 — простая
        expect(bp.fast.pieces.length > 1).toBe(ex.waves[0].kind === 'fan');
      });

      it('каждый фронт точно удовлетворяет условиям на разрывах (2.40)–(2.42)', () => {
        const z = solveImpact({ material: m, F0: [p0, 0, 0.001], bc: { type: 'velocity', d: 0.002 } })!.zones.at(-1)!;
        const bp = boundaryProblem(c, { F: [z.F[0], z.F[2]], v: [z.v[0], z.v[2]] }, 0, 2.5e-4)!;
        let ahead = { F: [z.F[0], z.F[2]], v: [z.v[0], z.v[2]] };
        for (const piece of [...bp.fast.pieces, ...bp.slow.pieces]) {
          const b = piece.behind;
          const Ta = traction(m, [ahead.F[0], 0, ahead.F[1]]);
          const Tb = traction(m, [b.F[0], 0, b.F[1]]);
          const G = piece.G;
          expect(Math.abs(Ta[0] - Tb[0] - G * G * (ahead.F[0] - b.F[0]))).toBeLessThan(1e-14);
          expect(Math.abs(Ta[2] - Tb[2] - G * G * (ahead.F[1] - b.F[1]))).toBeLessThan(1e-14);
          expect(Math.abs(b.v[0] - ahead.v[0] - G * (ahead.F[0] - b.F[0]))).toBeLessThan(1e-14);
          ahead = b;
        }
      });

      it('распад разрыва: одинаковые состояния — волн нет; состояния по обе стороны ударной волны — та же волна', () => {
        const same = riemannProblem(c, a, a, 2.5e-4)!;
        for (const w of [same.LF, same.LS, same.RS, same.RF]) expect(Math.abs(w.behind.F[0] - a.F[0])).toBeLessThan(1e-13);
        const ex = solveImpact({ material: m, F0: [p0, 0, 0.001], bc: { type: 'velocity', d: 0.001 } })!;
        const z0 = ex.zones[0];
        const z1 = ex.zones[1];
        const rp = riemannProblem(c, { F: [z1.F[0], z1.F[2]], v: [z1.v[0], z1.v[2]] }, { F: [z0.F[0], z0.F[2]], v: [z0.v[0], z0.v[2]] }, 2.5e-4)!;
        expect(rp.RF.pieces).toHaveLength(1);
        expect(rp.RF.pieces[0].G).toBeCloseTo(ex.waves[0].xiFront, 9);
        expect(Math.abs(rp.LF.behind.F[0] - z1.F[0])).toBeLessThan(1e-12);
      });

      const res = trackFronts({ material: m, F0: [p0, 0, 0.001], pulses: PULSES, tEnd: 0.92, frames: 47 });
      it('серия из трёх ударов отслеживается до конца без сбоев', () => {
        expect(res.warnings).toEqual([]);
        expect(res.events).toBeGreaterThan(20);
        expect(res.events).toBeLessThan(2000);
      });

      it('балансы: кинематика и энергия (работа = прирост энергии + диссипация на фронтах)', () => {
        expect(Math.abs(res.kinematicError)).toBeLessThan(1e-9);
        expect(Math.abs(res.work - res.energyGain - res.dissipation)).toBeLessThan(0.01 * (res.work - res.energyGain));
        expect(res.dissipation).toBeGreaterThan(0);
      });

      it('конечные объёмы сходятся к решению методом фронтов при измельчении сетки', () => {
        const diff = (N: number) => {
          const fv = simulate({ material: m, F0: [p0, 0, 0.001], pulses: PULSES, N, tEnd: 0.92, frames: 47 });
          let e = 0;
          for (const k of [23, 46]) for (let i = 0; i < N; i++) e += Math.abs(stateAt(res.frames[k], (i + 0.5) / N).F[0] - cell(fv.frames[k], N, i).F[0]) / N;
          return e / 2;
        };
        const d1 = diff(400);
        const d2 = diff(800);
        expect(d2).toBeLessThan(0.8 * d1);
        expect(d2).toBeLessThan(0.05 * 0.002);
      });
    });
  }
});

describe('растяжение границы (рывок наружу)', () => {
  for (const [name, m, p0] of CASES) {
    it(`${name}: растяжение уходит простой волной, ударной волны растяжения нет (§2.1, аналог теоремы Цемплена)`, () => {
      const ex = solveImpact({ material: m, F0: [p0, 0, 0.001], bc: { type: 'velocity', d: -0.002 } })!;
      expect(ex.waves[0].kind).toBe('fan');
      const bp = boundaryProblem({ m }, { F: [p0, 0.001], v: [0, 0] }, -0.002, 2.5e-4)!;
      expect(bp.fast.pieces.length).toBeGreaterThan(1);
      expect(Math.abs(bp.slow.behind.F[0] - ex.zones.at(-1)!.F[0])).toBeLessThan(1e-8);
    });

    it(`${name}: серия «сжатие → растяжение → сжатие» отслеживается фронтами, конечные объёмы к ней сходятся`, () => {
      const pulses = [
        { start: 0, duration: 0.08, d: 0.002 },
        { start: 0.25, duration: 0.08, d: -0.002 },
        { start: 0.5, duration: 0.08, d: 0.002 },
      ];
      const res = trackFronts({ material: m, F0: [p0, 0, 0.001], pulses, tEnd: 0.92, frames: 47 });
      expect(res.warnings).toEqual([]);
      expect(Math.abs(res.kinematicError)).toBeLessThan(1e-8);
      expect(res.mergeResidual).toBeLessThan(1e-4);
      const diff = (N: number) => {
        const fv = simulate({ material: m, F0: [p0, 0, 0.001], pulses, N, tEnd: 0.92, frames: 47 });
        let e = 0;
        for (const k of [23, 46]) for (let i = 0; i < N; i++) e += Math.abs(stateAt(res.frames[k], (i + 0.5) / N).F[0] - cell(fv.frames[k], N, i).F[0]) / N;
        return e / 2;
      };
      expect(diff(800)).toBeLessThan(0.8 * diff(400));
    });
  }

  it('материал 1: при растяжении семейство невыпукло — точный решатель предупреждает о составной волне', () => {
    const ex = solveImpact({ material: CASES[0][1], F0: [-0.00005, 0, 0.001], bc: { type: 'velocity', d: -0.002 } })!;
    expect(ex.warnings.some((w) => w.includes('составная волна'))).toBe(true);
  });
});

describe('разные скорости сжатия и растяжения (материал 2, слабое предсжатие)', () => {
  const m = normalize(MATERIAL_2);
  const F0 = [-0.0005, 0, 0.001] as [number, number, number];
  const pulses = [
    { start: 0, duration: 0.1, d: 0.002 },
    { start: 0.3, duration: 0.1, d: -0.002 },
  ];
  const res = trackFronts({ material: m, F0, pulses, tEnd: 0.92, frames: 47 });

  it('удар — ударная волна сжатия ≈ 1.07 c1, рывок — веер растяжения до ≈ 0.91 c1 (точный решатель)', () => {
    const up = solveImpact({ material: m, F0, bc: { type: 'velocity', d: 0.002 } })!;
    const down = solveImpact({ material: m, F0, bc: { type: 'velocity', d: -0.002 } })!;
    expect(up.waves[0].kind).toBe('shock');
    expect(up.waves[0].xiFront).toBeCloseTo(1.07, 2);
    expect(down.waves[0].kind).toBe('fan');
    expect(down.waves[0].xiBack).toBeCloseTo(0.914, 2);
    expect(down.zones.at(-1)!.F[0]).toBeGreaterThan(0); // среда действительно растянута
  });

  it('метод фронтов отрабатывает сценарий без сбоев (при необходимости — с огрублёнными порогами)', () => {
    expect(res.warnings).toEqual([]);
    const fastSpeeds = res.segments.filter((sg) => sg.fam === 'fast' && sg.x1 > sg.x0 && sg.strength > 1e-5).map((sg) => (sg.x1 - sg.x0) / (sg.t1 - sg.t0));
    expect(Math.max(...fastSpeeds)).toBeGreaterThan(1.06);
    expect(Math.min(...fastSpeeds)).toBeLessThan(0.92);
  });

  it('конечные объёмы согласуются с решением методом фронтов', () => {
    const fv = simulate({ material: m, F0, pulses, N: 1600, tEnd: 0.92, frames: 47 });
    let e = 0;
    for (const k of [23, 46]) for (let i = 0; i < 1600; i++) e += Math.abs(stateAt(res.frames[k], (i + 0.5) / 1600).F[0] - cell(fv.frames[k], 1600, i).F[0]) / 1600;
    expect(e / 2).toBeLessThan(0.02 * 0.002);
  });
});
