import { describe, expect, it } from 'vitest';
import { type ImpactProblem, type ImpactSolution, conservationCheck, displacementTable, profileAt, solveImpact } from '../src/core/impact.ts';
import { normalize, stress } from '../src/core/material.ts';
import { type Vec3, energy1d, shearStiffness, traction, traction229 } from '../src/core/material1d.ts';
import { DISSERTATION_CH2, IMPACT_PRESETS, MATERIAL_1, MATERIAL_2 } from '../src/core/presets.ts';

const m1 = normalize(MATERIAL_1);
const m2 = normalize(MATERIAL_2);
const prob = (name: keyof typeof IMPACT_PRESETS): ImpactProblem => {
  const p = IMPACT_PRESETS[name];
  return { material: normalize(p.material), F0: p.F0, bc: p.bc };
};
const PRESET_NAMES = Object.keys(IMPACT_PRESETS) as (keyof typeof IMPACT_PRESETS)[];

let seed = 777;
const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647 - 0.5;

describe('одномерная модель материала', () => {
  it('T = ∂W/∂F (σ_i1 — производные потенциала (1.9) по u_i,1)', () => {
    for (const m of [m1, m2])
      for (let t = 0; t < 100; t++) {
        const F: Vec3 = [-0.01 * Math.abs(rnd()) - 1e-4, 0.01 * rnd(), 0.01 * rnd()];
        const T = traction(m, F);
        for (let j = 0; j < 3; j++) {
          const h = 1e-8;
          const a = [...F] as Vec3;
          const b = [...F] as Vec3;
          a[j] += h;
          b[j] -= h;
          expect((energy1d(m, a) - energy1d(m, b)) / (2 * h)).toBeCloseTo(T[j], 9);
        }
      }
  });

  it('совпадает с плоской формулой (1.59) при s = 0 (g ↔ k)', () => {
    for (const m of [m1, m2]) {
      const F: Vec3 = [-0.004, 0.0023, 0];
      const pl = stress(m, { a: F[0], b: 0, k: F[1], n: 0 });
      const T = traction(m, F);
      expect(T[0]).toBeCloseTo(pl.s11, 15);
      expect(T[1]).toBeCloseTo(pl.s12, 15);
    }
  });

  it('ОПЕЧАТКА (2.29): со знаком «+3β/2» σ21 расходится с (1.9), со знаком «−3β/2» совпадает', () => {
    for (const m of [m1, m2]) {
      const F: Vec3 = [-0.006, 0.003, 0.002];
      const T = traction(m, F);
      const printed = traction229(m, F, 1);
      const fixed = traction229(m, F, -1);
      expect(fixed[0]).toBeCloseTo(T[0], 15);
      expect(fixed[1]).toBeCloseTo(T[1], 15);
      expect(fixed[2]).toBeCloseTo(T[2], 15);
      expect(printed[0]).toBeCloseTo(T[0], 15); // σ11 напечатан верно
      expect(Math.abs(printed[1] / T[1] - 1)).toBeGreaterThan(0.009); // σ21 — нет
    }
  });

  it('касательное напряжение коллинеарно сдвигу: (σ21, σ31) = Φ (g, s)', () => {
    const F: Vec3 = [-0.005, 0.002, -0.003];
    const T = traction(m1, F);
    const phi = shearStiffness(m1, F);
    expect(T[1]).toBeCloseTo(phi * F[1], 16);
    expect(T[2]).toBeCloseTo(phi * F[2], 16);
  });
});

function checkSolution(sol: ImpactSolution) {
  const { zones, waves } = sol;
  const law = (F: Vec3) => traction(sol.problem.material, F);
  waves.forEach((w, i) => {
    const a = zones[i];
    const b = zones[i + 1];
    if (w.kind !== 'fan') {
      // (2.40)–(2.42): [T] = ξ²[F], [v̂] = −ξ[F]
      const Ta = law(a.F);
      const Tb = law(b.F);
      for (let k = 0; k < 3; k++) {
        expect(Math.abs(Ta[k] - Tb[k] - w.xiFront ** 2 * (a.F[k] - b.F[k]))).toBeLessThan(1e-13);
        expect(Math.abs(a.v[k] - b.v[k] + w.xiFront * (a.F[k] - b.F[k]))).toBeLessThan(1e-13);
      }
    } else {
      const pr = w.profile!;
      // края веера совпадают с соседними зонами, ξ убывает к границе
      for (let k = 0; k < 3; k++) {
        expect(pr[0].F[k]).toBeCloseTo(a.F[k], 14);
        expect(pr.at(-1)!.F[k]).toBeCloseTo(b.F[k], 14);
      }
      for (let j = 1; j < pr.length; j++) expect(pr[j].xi).toBeLessThan(pr[j - 1].xi);
    }
  });
  for (let i = 1; i < waves.length; i++) expect(waves[i].xiFront).toBeLessThanOrEqual(waves[i - 1].xiBack + 1e-12);
}

describe('уровень 2: линейная упругость (ν = α = β = 0)', () => {
  it('косой удар: скорости 1 и c2/c1, диссипация 0, нейтральная волна со скоростью c2/c1', () => {
    const m = { ...m1, nu: 0, alpha: 0, beta: 0 };
    const sol = solveImpact({ material: m, F0: [-0.002, 0, 0.001], bc: { type: 'strain', F: [-0.004, 0.002, 0.003] } })!;
    expect(sol).not.toBeNull();
    const c2 = Math.sqrt(m.mu);
    const fast = sol.waves.find((w) => w.family === 'fast')!;
    const slow = sol.waves.find((w) => w.family === 'slow')!;
    const neu = sol.waves.find((w) => w.family === 'neutral')!;
    expect(fast.xiFront).toBeCloseTo(1, 9);
    expect(slow.xiFront).toBeCloseTo(c2, 9);
    expect(neu.xiFront).toBeCloseTo(c2, 9);
    for (const w of sol.waves) expect(Math.abs(w.dissipation)).toBeLessThan(1e-15);
  });
});

describe('уровень 3: согласованность решений для всех случаев диссертации', () => {
  for (const name of PRESET_NAMES) {
    describe(name, () => {
      const sol = solveImpact(prob(name))!;
      it('решение найдено и допустимо', () => {
        expect(sol).not.toBeNull();
        expect(sol.warnings).toEqual([]);
      });
      it('условия на разрывах и в веерах', () => checkSolution(sol));
      it('граничное условие выполнено', () => {
        const zb = sol.zones.at(-1)!;
        const bc = sol.problem.bc;
        if (bc.type === 'velocity') {
          expect(zb.v[0]).toBeCloseTo(bc.d, 13);
          expect(Math.abs(zb.v[1]) + Math.abs(zb.v[2])).toBeLessThan(1e-13);
        } else for (let k = 0; k < 3; k++) expect(zb.F[k]).toBeCloseTo(bc.F[k], 11);
      });
      it('интегральные законы сохранения: кинематика, импульс, энергия', () => {
        const c = conservationCheck(sol);
        const scale = Math.max(...zonesAbs(sol));
        for (let k = 0; k < 3; k++) {
          expect(Math.abs(c.kinematic.lhs[k] - c.kinematic.rhs[k])).toBeLessThan(1e-6 * scale);
          expect(Math.abs(c.momentum.lhs[k] - c.momentum.rhs[k])).toBeLessThan(1e-6 * scale);
        }
        expect(Math.abs(c.energy.lhs - c.energy.rhs)).toBeLessThan(1e-4 * Math.abs(c.energy.rhs));
      });
      it('перемещение границы: c1·t·w(0) = v(0)·t', () => {
        const w = displacementTable(sol)(0);
        const v = sol.zones.at(-1)!.v;
        for (let k = 0; k < 3; k++) expect(Math.abs(w[k] - v[k])).toBeLessThan(1e-6 * Math.max(...zonesAbs(sol)));
      });
      it('ударные волны термодинамически допустимы (2.45), нейтральная — без диссипации', () => {
        for (const w of sol.waves) {
          if (w.kind === 'shock') {
            expect(w.dissipation).toBeGreaterThan(0);
            expect(w.lax!.ok).toBe(true);
          }
          if (w.kind === 'rotation') expect(Math.abs(w.dissipation)).toBeLessThan(1e-15);
        }
      });
    });
  }
});

function zonesAbs(sol: ImpactSolution): number[] {
  return sol.zones.flatMap((z) => [...z.F, ...z.v].map(Math.abs));
}

describe('уровень 4: сверка с диссертацией (глава 2)', () => {
  const D = DISSERTATION_CH2;
  it('§2.3, материал 1 (рис. 11): ξ1 = 1.0067 (±0.0005), ξ2 = 0.5671 (±0.0002); две ударные волны', () => {
    const s = solveImpact(prob('§2.3 нормальный удар, материал 1 (рис. 11)'))!;
    expect(s.pattern).toBe('SS');
    expect(Math.abs(s.waves[0].xiFront - D.fig11.xi1)).toBeLessThan(0.0005);
    expect(Math.abs(s.waves[1].xiFront - D.fig11.xi2)).toBeLessThan(0.0002);
    // рис. 11: p1 ≈ −1.0e-3, s1 < s0 < s2, h1 > 0
    expect(s.zones[1].F[0]).toBeCloseTo(-1.04e-3, 4);
    expect(s.zones[1].F[2]).toBeLessThan(0.001);
    expect(s.zones[2].F[2]).toBeGreaterThan(0.001);
    expect(s.zones[1].v[2]).toBeGreaterThan(0);
  });

  it('§2.3, материал 2 (рис. 13): ξ1 = 1.08 (±0.002), ξ2 = 0.6799 (±0.001); h1 < 0, s1 > s0 > s2', () => {
    const s = solveImpact(prob('§2.3 нормальный удар, материал 2 (рис. 13)'))!;
    expect(Math.abs(s.waves[0].xiFront - D.fig13.xi1)).toBeLessThan(0.002);
    expect(Math.abs(s.waves[1].xiFront - D.fig13.xi2)).toBeLessThan(0.001);
    expect(s.zones[1].v[2]).toBeCloseTo(-1.0e-6, 7);
    expect(s.zones[1].F[2]).toBeGreaterThan(0.001);
    expect(s.zones[2].F[2]).toBeLessThan(0.001);
  });

  it('рис. 13: при p0 = −0.0005 (как на оси рисунка) ξ1 ≈ 1.061 ≠ 1.08 — подпись оси ошибочна', () => {
    const s = solveImpact({ ...prob('§2.3 нормальный удар, материал 2 (рис. 13)'), F0: [-0.0005, 0, 0.001] })!;
    expect(Math.abs(s.waves[0].xiFront - 1.08)).toBeGreaterThan(0.015);
  });

  it('напечатанная (2.29) для §2.3, материала 1 даёт веер ξ ≈ 0.545 вместо ударной волны 0.5671 — автор считал по верной формуле', () => {
    const s = solveImpact({ ...prob('§2.3 нормальный удар, материал 1 (рис. 11)'), stressLaw: 'printed229' })!;
    expect(s.pattern).toBe('SR');
    expect(Math.abs(s.waves[1].xiFront - D.fig11.xi2)).toBeGreaterThan(0.02);
  });

  it('§2.4.1, материал 1 (рис. 16): три ударные волны Σ1 > Σ2 > Σ3, скорости в пределах 0.0006', () => {
    const s = solveImpact(prob('§2.4.1 косой удар, материал 1 (рис. 16)'))!;
    expect(s.waves.map((w) => w.label)).toEqual(['Σ1', 'Σ2', 'Σ3']);
    expect(s.waves.map((w) => w.kind)).toEqual(['shock', 'shock', 'rotation']);
    const [x1, x2, x3] = s.waves.map((w) => w.xiFront);
    expect(Math.abs(x1 - D.fig16figure.xi1)).toBeLessThan(0.0006);
    expect(Math.abs(x2 - D.fig16figure.xi2)).toBeLessThan(0.0003);
    expect(Math.abs(x3 - D.fig16figure.xi3)).toBeLessThan(0.0003);
    // рис. 16: |p1| > |p0|, s1 < s0; в зоне II g = 0, s = |τ*|
    expect(s.zones[1].F[2]).toBeLessThan(0.001);
    expect(s.zones[2].F[1]).toBeCloseTo(0, 14);
    expect(s.zones[2].F[2]).toBeCloseTo(Math.hypot(0.005, 0.005), 12);
  });

  it('§2.4.2, материал 2 (рис. 18): та же структура — Σ1, затем Σ3, затем простая волна', () => {
    const s = solveImpact(prob('§2.4.2 косой удар, материал 2 (рис. 18)'))!;
    expect(s.waves.map((w) => w.kind)).toEqual(['shock', 'rotation', 'fan']);
  });

  it('§2.4.2 с β = −0.3 (как напечатано): скорость Σ1 = 1.081, а не 1.0424', () => {
    const s = solveImpact(prob('§2.4.2 косой удар, материал 2 (рис. 18)'))!;
    expect(s.waves[0].xiFront).toBeCloseTo(1.0812, 3);
  });

  it('§2.4.2 с β = +0.3: все скорости совпадают с рис. 18 в пределах 0.0012', () => {
    const s = solveImpact(prob('§2.4.2, материал 2 с β = +0.3 (гипотеза)'))!;
    expect(s.waves.map((w) => w.kind)).toEqual(['shock', 'rotation', 'fan']);
    expect(Math.abs(s.waves[0].xiFront - D.fig18.xi1)).toBeLessThan(0.0012);
    expect(Math.abs(s.waves[1].xiFront - D.fig18.xi3)).toBeLessThan(0.0012);
    expect(Math.abs(s.waves[2].xiFront - D.fig18.fanFront)).toBeLessThan(0.0012);
    expect(Math.abs(s.waves[2].xiBack - D.fig18.fanBack)).toBeLessThan(0.0012);
  });

  it('§2.4.2: ширина веера растёт с интенсивностью воздействия (стр. 87)', () => {
    const base = prob('§2.4.2, материал 2 с β = +0.3 (гипотеза)');
    const width = (k: number) => {
      const s = solveImpact({ ...base, bc: { type: 'strain', F: [-0.01, 0.006 * k, 0.005 * k] } })!;
      const f = s.waves.find((w) => w.kind === 'fan')!;
      return f.xiFront - f.xiBack;
    };
    expect(width(1.2)).toBeGreaterThan(width(1));
    expect(width(1)).toBeGreaterThan(width(0.8));
  });

  it('рис. 12д,е / 14д,е: зависимость скоростей волн от d2 ∈ [0.001; 0.01]', () => {
    const run = (name: keyof typeof IMPACT_PRESETS) =>
      [0.001, 0.0025, 0.005, 0.0075, 0.01].map((d) => solveImpact({ ...prob(name), bc: { type: 'velocity', d } })!);
    const a = run('§2.3 нормальный удар, материал 1 (рис. 11)');
    const b = run('§2.3 нормальный удар, материал 2 (рис. 13)');
    const M2a = a.map((s) => s.waves[1].xiFront);
    const M2b = b.map((s) => s.waves[1].xiFront);
    const M1b = b.map((s) => s.waves[0].xiFront);
    for (let i = 1; i < 5; i++) {
      expect(M2a[i]).toBeLessThan(M2a[i - 1]); // материал 1: M2 убывает
      expect(M2b[i]).toBeGreaterThan(M2b[i - 1]); // материал 2: M2 растёт
      expect(M1b[i]).toBeGreaterThanOrEqual(M1b[i - 1]); // материал 2: M1 растёт
    }
    expect(M2a[0]).toBeLessThan(D.fig12.M2[1] + 0.0005);
    expect(M2a.at(-1)!).toBeGreaterThan(D.fig12.M2[0] + 0.0005);
    // материал 1: d1, h1 > 0 и растут; материал 2: h1 < 0, |h1| растёт
    for (let i = 1; i < 5; i++) {
      expect(a[i].zones[1].v[2]).toBeGreaterThan(a[i - 1].zones[1].v[2]);
      expect(b[i].zones[1].v[2]).toBeLessThan(b[i - 1].zones[1].v[2]);
    }
  });
});

describe('профиль и перемещения', () => {
  it('profileAt возвращает зоны и значения внутри веера', () => {
    const s = solveImpact(prob('§2.4.2, материал 2 с β = +0.3 (гипотеза)'))!;
    const fanW = s.waves[2];
    const mid = profileAt(s, 0.5 * (fanW.xiFront + fanW.xiBack)).F;
    const a = s.zones[2].F;
    const b = s.zones[3].F;
    expect(mid[1]).toBeGreaterThan(Math.min(a[1], b[1]));
    expect(mid[1]).toBeLessThan(Math.max(a[1], b[1]));
    expect(profileAt(s, 2).F).toEqual(s.zones[0].F);
    expect(profileAt(s, 0).F).toEqual(s.zones.at(-1)!.F);
  });
});
