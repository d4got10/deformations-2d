import { describe, expect, it } from 'vitest';
import { laxCheck } from '../src/core/acoustic.ts';
import { c2OverC1, normalize, scaleNonlinear } from '../src/core/material.ts';
import {
  type MovingLoadProblem,
  displacement,
  residual10,
  solveMovingLoad,
  sweep,
  zoneIndex,
} from '../src/core/movingLoad.ts';
import { solveLinear } from '../src/core/newton.ts';
import { BASE_CASE, DISSERTATION, LOAD_SWEEP_CASE, MATERIAL_1, MATERIAL_2 } from '../src/core/presets.ts';

const deg = (r: number) => (r * 180) / Math.PI;
const P1: MovingLoadProblem = { material: MATERIAL_1, ...BASE_CASE };
const P2: MovingLoadProblem = { material: MATERIAL_2, ...BASE_CASE };

describe('линейная алгебра', () => {
  it('Гаусс с выбором ведущего элемента', () => {
    const x = solveLinear(
      [
        [0, 2, 1],
        [1, 1, 1],
        [2, 1, 0],
      ],
      [7, 6, 4],
    )!;
    expect(x[0]).toBeCloseTo(1, 14);
    expect(x[1]).toBeCloseTo(2, 14);
    expect(x[2]).toBeCloseTo(3, 14);
  });
});

describe('уровень 2: предельный переход к линейной упругости (ν = α = β = 0)', () => {
  for (const M of [1.2, 1.7, 2.5]) {
    it(`M = ${M}: фронты — конусы Маха, скачки поляризованы как P- и S-волны, диссипация 0`, () => {
      const p = { ...P1, M, material: scaleNonlinear(MATERIAL_1, 0) };
      const s = solveMovingLoad(p)!;
      expect(s).not.toBeNull();
      expect(s.fronts[0].thetaDeg).toBeCloseTo(deg(Math.asin(1 / M)), 9);
      expect(s.fronts[1].thetaDeg).toBeCloseTo(deg(Math.asin(c2OverC1(MATERIAL_1) / M)), 9);
      expect(s.fronts[0].G).toBeCloseTo(1, 10);
      const [z0, z1, z2] = s.zones;
      // P-волна: скачок ∇u ∝ ν⊗ν ⇒ Δk/Δa = tg θ1; S-волна: ∝ τ⊗ν ⇒ Δk/Δa = −ctg θ2
      expect((z1.grad.k - z0.grad.k) / (z1.grad.a - z0.grad.a)).toBeCloseTo(s.fronts[0].xi, 8);
      expect((z2.grad.k - z1.grad.k) / (z2.grad.a - z1.grad.a)).toBeCloseTo(-1 / s.fronts[1].xi, 8);
      for (const f of s.fronts) expect(Math.abs(f.dissipation)).toBeLessThan(1e-18);
    });
  }
});

describe('уровень 3: внутренняя согласованность решения', () => {
  for (const [name, p] of [
    ['материал 1', P1],
    ['материал 2', P2],
  ] as const) {
    const s = solveMovingLoad(p)!;
    it(`${name}: полная система (3.36) из 10 уравнений выполнена`, () => {
      expect(s.residual).toBeLessThan(1e-15);
      const z = s.zones;
      const v = [z[1].grad.a, z[1].grad.b, z[1].grad.k, z[1].grad.n, z[2].grad.a, z[2].grad.b, z[2].grad.k, z[2].grad.n, s.fronts[0].xi, s.fronts[1].xi];
      for (const r of residual10(p, v)) expect(Math.abs(r)).toBeLessThan(1e-15);
    });
    it(`${name}: граничные условия (3.31) в зоне II`, () => {
      expect(s.zones[2].stress.s11).toBeCloseTo(p.s11, 14);
      expect(s.zones[2].stress.s12).toBeCloseTo(p.s21, 14);
    });
    it(`${name}: перемещения непрерывны на фронтах (сетка не рвётся)`, () => {
      for (let i = 0; i < 2; i++) {
        const xi = s.fronts[i].xi;
        for (const eta of [-0.3, -1, -7]) {
          const x1 = -eta * xi;
          const ua = displacement(s.zones[i].grad, x1, eta);
          const ub = displacement(s.zones[i + 1].grad, x1, eta);
          expect(ua[0]).toBeCloseTo(ub[0], 15);
          expect(ua[1]).toBeCloseTo(ub[1], 15);
        }
      }
    });
    it(`${name}: G = M sin θ (3.24), θ1 > θ2, точки разнесены по зонам верно`, () => {
      for (const f of s.fronts) expect(f.G).toBeCloseTo(p.M * Math.sin(Math.atan(f.xi)), 14);
      expect(s.fronts[0].xi).toBeGreaterThan(s.fronts[1].xi);
      expect(zoneIndex(s, 1, 0.5)).toBe(0);
      expect(zoneIndex(s, 0.01, -1)).toBe(2);
      const mid = Math.tan(0.5 * (Math.atan(s.fronts[0].xi) + Math.atan(s.fronts[1].xi)));
      expect(zoneIndex(s, mid, -1)).toBe(1);
    });
    it(`${name}: дополнительное сжатие (3.19) в зонах I, II`, () => {
      for (const z of s.zones.slice(1)) {
        const J1 = z.grad.a + z.grad.n;
        expect(J1).toBeLessThan(0);
        expect(Math.abs(J1)).toBeGreaterThan(Math.abs(p.a0));
      }
    });
  }

  it('M ≤ 1 и a0 ≥ 0 отвергаются даже при старте от соседнего решения', () => {
    const s = solveMovingLoad(P1)!;
    const z = s.zones;
    const x0 = [z[1].grad.a, z[1].grad.k, z[2].grad.a, z[2].grad.k, s.fronts[0].xi, s.fronts[1].xi];
    expect(() => solveMovingLoad({ ...P1, M: 0.9 }, x0)).toThrow();
    expect(() => solveMovingLoad({ ...P1, a0: 0.001 }, x0)).toThrow();
  });

  it('решение не зависит от способа старта (от линейного решения или от соседней точки)', () => {
    const near = solveMovingLoad({ ...P1, M: 1.72 })!;
    const z = near.zones;
    const fromNeighbour = solveMovingLoad(P1, [z[1].grad.a, z[1].grad.k, z[2].grad.a, z[2].grad.k, near.fronts[0].xi, near.fronts[1].xi])!;
    const fresh = solveMovingLoad(P1)!;
    expect(fromNeighbour.fronts[0].xi).toBeCloseTo(fresh.fronts[0].xi, 12);
    expect(fromNeighbour.fronts[1].xi).toBeCloseTo(fresh.fronts[1].xi, 12);
  });
});

describe('уровень 4: сверка с диссертацией (§3.2, стр. 106–111)', () => {
  const s1 = solveMovingLoad(P1)!;
  const s2 = solveMovingLoad(P2)!;
  const d1 = DISSERTATION.material1;
  const d2 = DISSERTATION.material2;

  it('материал 2: θ1, θ2 совпадают с текстом в пределах 0.1°, G — в пределах 0.3 %', () => {
    expect(Math.abs(s2.fronts[0].thetaDeg - d2.theta1)).toBeLessThan(0.1);
    expect(Math.abs(s2.fronts[1].thetaDeg - d2.theta2)).toBeLessThan(0.1);
    expect(Math.abs(s2.fronts[0].G / d2.G1 - 1)).toBeLessThan(0.003);
    expect(Math.abs(s2.fronts[1].G / d2.G2 - 1)).toBeLessThan(0.003);
  });

  it('материал 1: θ1 и G1 совпадают с текстом (0.05°, 0.1 %)', () => {
    expect(Math.abs(s1.fronts[0].thetaDeg - d1.theta1)).toBeLessThan(0.05);
    expect(Math.abs(s1.fronts[0].G / d1.G1 - 1)).toBeLessThan(0.001);
  });

  it('материал 1: θ2 — ИЗВЕСТНОЕ РАСХОЖДЕНИЕ ≈ +0.3° (см. docs/VERIFICATION.md, вопрос 1)', () => {
    const diff = s1.fronts[1].thetaDeg - d1.theta2;
    expect(diff).toBeGreaterThan(0.2);
    expect(diff).toBeLessThan(0.35);
  });

  it('рис. 21: напряжения и сдвиги по зонам совпадают с графиком в пределах точности чтения', () => {
    const f = DISSERTATION.fig21;
    for (let i = 0; i < 3; i++) {
      const z = s1.zones[i];
      expect(Math.abs(z.stress.s11 - f.s11[i])).toBeLessThan(6e-4);
      expect(Math.abs(z.stress.s12 - f.s21[i])).toBeLessThan(3e-4);
      expect(Math.abs(z.grad.b - f.b[i])).toBeLessThan(3e-4);
      expect(Math.abs(z.grad.k - f.k[i])).toBeLessThan(5e-4);
    }
  });

  it('рис. 22б,в (материал 1, M = 1.4): θ1 меняется в пределах 46.05…46.25°, тренды как у автора', () => {
    const base = { material: MATERIAL_1, ...LOAD_SWEEP_CASE };
    const b = sweep(base, [0.006, 0.01, 0.02, 0.04, 0.06], (q, v) => ({ ...q, s11: -v }));
    const c = sweep(base, [0.0005, 0.002, 0.004, 0.006], (q, v) => ({ ...q, s21: v }));
    expect(b.failedAt).toBeNull();
    expect(c.failedAt).toBeNull();
    const t1b = b.points.map((p) => p.solution.fronts[0].thetaDeg);
    const t1c = c.points.map((p) => p.solution.fronts[0].thetaDeg);
    expect(t1b[0]).toBeCloseTo(46.25, 1);
    expect(t1b.at(-1)!).toBeCloseTo(46.02, 1);
    for (let i = 1; i < t1b.length; i++) expect(t1b[i]).toBeLessThan(t1b[i - 1]); // |σ11*| ↑ ⇒ θ1 ↓
    for (let i = 1; i < t1c.length; i++) expect(t1c[i]).toBeGreaterThan(t1c[i - 1]); // σ21* ↑ ⇒ θ1 ↑
  });

  it('рис. 24б,в (материал 2, M = 1.4): тренды противоположны материалу 1, как у автора', () => {
    const base = { material: MATERIAL_2, ...LOAD_SWEEP_CASE };
    const b = sweep(base, [0.006, 0.01, 0.02, 0.04], (q, v) => ({ ...q, s11: -v, s21: 0.3 * v }));
    const c = sweep(base, [0.0005, 0.002, 0.004, 0.006], (q, v) => ({ ...q, s21: v }));
    for (const i of [0, 1] as const) {
      const tb = b.points.map((p) => p.solution.fronts[i].thetaDeg);
      const tc = c.points.map((p) => p.solution.fronts[i].thetaDeg);
      expect(tb.at(-1)!).toBeGreaterThan(tb[0]);
      expect(tc.at(-1)!).toBeLessThan(tc[0]);
    }
  });

  it('θ(M) монотонно убывает и при M = 1.7 совпадает с базовым случаем (рис. 22а — см. вопрос 2)', () => {
    const Ms = [2.5, 2.2, 2.0, 1.8, 1.7, 1.5, 1.3, 1.2];
    const sw = sweep(P1, Ms, (q, M) => ({ ...q, M }));
    expect(sw.failedAt).toBeNull();
    for (let i = 1; i < sw.points.length; i++)
      for (const j of [0, 1]) expect(sw.points[i].solution.fronts[j].thetaDeg).toBeGreaterThan(sw.points[i - 1].solution.fronts[j].thetaDeg);
    expect(sw.points[4].solution.fronts[0].thetaDeg).toBeCloseTo(s1.fronts[0].thetaDeg, 10);
  });
});

describe('допустимость ударных волн: (3.30) и условие Лакса согласованы', () => {
  for (const [name, p] of [
    ['материал 1', P1],
    ['материал 2', P2],
  ] as const) {
    it(`${name}: знак диссипации совпадает с выполнением условия Лакса на каждом фронте`, () => {
      const s = solveMovingLoad(p)!;
      const m = normalize(p.material);
      s.fronts.forEach((f, i) => {
        const lax = laxCheck(m, s.zones[i].grad, s.zones[i + 1].grad, Math.atan(f.xi), f.G, i === 0 ? 1 : 0);
        expect(lax.ok).toBe(f.dissipation > 0);
      });
    });
  }
  it('ИЗВЕСТНОЕ РАСХОЖДЕНИЕ: материал 2, Σ2 термодинамически недопустима (у автора — допустима, вопрос 3)', () => {
    const s = solveMovingLoad(P2)!;
    expect(s.fronts[1].dissipation).toBeLessThan(0);
    expect(s.warnings.some((w) => w.includes('Σ2'))).toBe(true);
  });
});
