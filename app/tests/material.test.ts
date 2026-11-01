import { describe, expect, it } from 'vitest';
import { type Grad, energy, invariants, normalize, stress } from '../src/core/material.ts';
import { MATERIAL_1, MATERIAL_2 } from '../src/core/presets.ts';

let seed = 12345;
const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647 - 0.5;
const randomGrad = (s = 0.02): Grad => ({ a: s * rnd(), b: s * rnd(), k: s * rnd(), n: s * rnd() });

describe('модель материала (1.9)–(1.10), (1.59)–(1.60)', () => {
  for (const [name, m] of [
    ['материал 1', normalize(MATERIAL_1)],
    ['материал 2', normalize(MATERIAL_2)],
  ] as const) {
    it(`${name}: σ = ∂W/∂e (гиперупругая согласованность (1.10) с потенциалом (1.9))`, () => {
      for (let t = 0; t < 200; t++) {
        const g = randomGrad();
        const s = stress(m, g);
        // dW = σ11 de11 + σ22 de22 + 2 σ12 de12, малые приращения
        const h = 1e-8;
        const dW = (dg: Partial<Grad>) =>
          (energy(m, { ...g, ...Object.fromEntries(Object.entries(dg).map(([k, v]) => [k, g[k as keyof Grad] + v])) }) -
            energy(m, { ...g, ...Object.fromEntries(Object.entries(dg).map(([k, v]) => [k, g[k as keyof Grad] - v])) })) /
          (2 * h);
        expect(dW({ a: h })).toBeCloseTo(s.s11, 9);
        expect(dW({ n: h })).toBeCloseTo(s.s22, 9);
        expect(dW({ b: h })).toBeCloseTo(s.s12, 9); // ∂W/∂u1,2 = σ12 (e12 = (b+k)/2, множитель 2·½)
        expect(dW({ k: h })).toBeCloseTo(s.s12, 9);
      }
    });

    it(`${name}: совпадение с явными формулами плоской деформации (1.59)`, () => {
      for (let t = 0; t < 100; t++) {
        const g = randomGrad();
        const { J1, J2, J3 } = invariants(g);
        const r = Math.sqrt(J2);
        const phi = J1 / r;
        const eta = J3 / r ** 3;
        const Mm = 3 * m.alpha * phi - m.nu / phi;
        const N = (3 * m.alpha - m.nu) * phi - m.nu / phi - m.alpha * phi ** 3 - m.beta * eta;
        const K = ((3 * m.beta) / 4) * ((g.b + g.k) / r);
        const F = ((3 * m.beta) / 2 - m.nu / 2) * phi - (m.alpha / 2) * phi ** 3 - (m.beta / 2) * eta;
        const U1 = g.a ** 2 / r;
        const U2 = g.n ** 2 / r;
        const s = stress(m, g);
        const s11 = (m.lambda + 2 * m.mu + N) * g.a + (m.lambda + Mm) * g.n + K * (g.b + g.k) + 3 * m.beta * U1;
        const s22 = (m.lambda + 2 * m.mu + N) * g.n + (m.lambda + Mm) * g.a + K * (g.b + g.k) + 3 * m.beta * U2;
        const s12 = (m.mu + F) * (g.b + g.k);
        const s33 = (m.lambda + Mm) * (g.a + g.n);
        expect(s.s11).toBeCloseTo(s11, 13);
        expect(s.s22).toBeCloseTo(s22, 13);
        expect(s.s12).toBeCloseTo(s12, 13);
        expect(s.s33).toBeCloseTo(s33, 13);
      }
    });

    it(`${name}: начальное напряжённое состояние (3.2)`, () => {
      const a0 = -0.005;
      const k0 = 0.001;
      const k = Math.sign(a0);
      const p = 0.5 * (k0 / a0) ** 2;
      const phi0 = k / Math.sqrt(1 + p);
      const eta0 = (k * (1 + 1.5 * p)) / Math.sqrt((1 + p) ** 3);
      const den = Math.sqrt(a0 ** 2 + 0.5 * k0 ** 2);
      const s = stress(m, { a: a0, b: 0, k: k0, n: 0 });
      expect(s.s11).toBeCloseTo(
        (m.lambda + 2 * m.mu - m.nu / phi0 - m.alpha * phi0 ** 3 - m.beta * eta0) * a0 +
          (3 * m.beta * (a0 ** 2 + 0.25 * k0 ** 2)) / den + 3 * m.alpha * phi0 * a0 - m.nu * phi0 * a0,
        13,
      );
      expect(s.s22).toBeCloseTo((m.lambda - m.nu / phi0 + 3 * m.alpha * phi0) * a0 + ((3 * m.beta) / 4) * (k0 ** 2 / den), 13);
      expect(s.s12).toBeCloseTo((m.mu + ((3 * m.beta) / 2 - m.nu / 2) * phi0 - (m.alpha / 2) * phi0 ** 3 - (m.beta / 2) * eta0) * k0, 13);
    });

    it(`${name}: объективность — поворот деформации поворачивает тензор напряжений`, () => {
      for (let t = 0; t < 50; t++) {
        const g = randomGrad();
        const th = Math.PI * rnd();
        const c = Math.cos(th);
        const s_ = Math.sin(th);
        const e = [
          [g.a, (g.b + g.k) / 2],
          [(g.b + g.k) / 2, g.n],
        ];
        const R = [
          [c, -s_],
          [s_, c],
        ];
        const rot = (T: number[][]) => {
          const out = [
            [0, 0],
            [0, 0],
          ];
          for (let i = 0; i < 2; i++)
            for (let j = 0; j < 2; j++)
              for (let k = 0; k < 2; k++) for (let l = 0; l < 2; l++) out[i][j] += R[i][k] * R[j][l] * T[k][l];
          return out;
        };
        const er = rot(e);
        const sr = stress(m, { a: er[0][0], b: er[0][1], k: er[1][0], n: er[1][1] });
        const s0 = stress(m, g);
        const s0r = rot([
          [s0.s11, s0.s12],
          [s0.s12, s0.s22],
        ]);
        expect(sr.s11).toBeCloseTo(s0r[0][0], 13);
        expect(sr.s12).toBeCloseTo(s0r[0][1], 13);
        expect(sr.s22).toBeCloseTo(s0r[1][1], 13);
        expect(sr.s33).toBeCloseTo(s0.s33, 13);
      }
    });
  }

  it('неаналитичность в нуле: напряжения при J2 = 0 не определены', () => {
    expect(() => stress(normalize(MATERIAL_1), { a: 0, b: 0, k: 0, n: 0 })).toThrow();
  });

  it('при ν = α = β = 0 — закон Гука', () => {
    const m = { lambda: 0.4, mu: 0.3, nu: 0, alpha: 0, beta: 0 };
    const g = { a: 0.01, b: 0.002, k: -0.003, n: -0.004 };
    const s = stress(m, g);
    expect(s.s11).toBeCloseTo(0.4 * (g.a + g.n) + 0.6 * g.a, 15);
    expect(s.s12).toBeCloseTo(0.3 * (g.b + g.k), 15);
  });
});
