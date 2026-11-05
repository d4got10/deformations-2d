/**
 * Сверка с одноосной задачей: Лаптева А.А., Дудко О.В., Рагозина В.Е. «О динамике одномерных
 * деформаций при циклическом одноосном нагружении разномодульного упругого полупространства»
 * (Изв. Сарат. ун-та, 2022; sources/lapteva-dudko-ragozina-2022-cyclic-uniaxial.pdf).
 * Модель Ляховского–Мясникова = наша при α = β = 0; λ = 19.8, μ = 26.9, ν = 15.3 ГПа, ρ = 7600 кг/м³.
 * Единицы: время — мс, длина — c1·1 мс, скорость — c1 = √((λ+2μ)/ρ).
 */
import { describe, expect, it } from 'vitest';
import { solveBimodular } from '../src/core/bimodular1d.ts';
import { stateAt, trackFronts } from '../src/core/fronts.ts';
import { hybridStateAt, solveHybrid } from '../src/core/hybrid.ts';
import { normalize } from '../src/core/material.ts';

const GPa = { lambda: 19.8, mu: 26.9, nu: 15.3, alpha: 0, beta: 0 };
const rho = 7600;
const c1 = Math.sqrt(((GPa.lambda + 2 * GPa.mu) * 1e9) / rho);
const m = normalize(GPa);
const a = Math.sqrt(1 + 2 * m.nu);
const b = Math.sqrt(1 - 2 * m.nu);

/** Узлы u0(t) (t — мс, u0 — мм); последний наклон продолжается. */
function boundary(ts: number[], us: number[]) {
  const times = ts.slice(0, -1);
  const velocities = ts.slice(1).map((t, j) => (us[j + 1] - us[j]) / (t - ts[j]) / c1); // мм/мс = м/с
  const pulses = times.map((t, k) => ({ start: t, duration: (times[k + 1] ?? 1e3) - t, d: velocities[k] }));
  return { times, velocities, pulses };
}

const CASES = {
  'I(a)': { ts: [0, 0.2, 0.5, 0.75, 0.9], us: [0, -3, 0, 4, 0], tEnd: 2.5 },
  'I(b)': { ts: [0, 0.2, 0.5, 0.75, 0.9], us: [0, -3, -1, 4, 0], tEnd: 3.0 },
  II: { ts: [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3], us: [0, -0.3, 0.4, -0.3, 0.4, -0.3, 0.4, 0], tEnd: 2.8 },
} as const;

describe('одноосная задача Лаптевой–Дудко–Рагозиной (2022)', () => {
  it('скорости волн сжатия и растяжения (4): a = 3702.77 м/с, b = 2378.63 м/с', () => {
    expect(a * c1).toBeCloseTo(3702.77, 1);
    expect(b * c1).toBeCloseTo(2378.63, 1);
  });

  it('первое столкновение во вторичной картине (рис. 4 статьи): (a) — инцидент 1 при τ5 ≈ 1.27 мс, (b) — инцидент 2 при τ5 ≈ 1.32 мс', () => {
    const first = (k: 'I(a)' | 'I(b)') => {
      const bc = boundary([...CASES[k].ts], [...CASES[k].us]);
      return solveBimodular({ a, b, ...bc, tEnd: CASES[k].tEnd }).events.find((e) => e.type === 'collision')!;
    };
    const ea = first('I(a)');
    expect(ea.incoming).toEqual(['Σ', 'γᵇ']);
    expect(ea.outgoing).toEqual(['ξ⁻ᵃ', 'γᵃ']);
    expect(ea.t).toBeCloseTo(1.2684, 4);
    const eb = first('I(b)');
    expect(eb.incoming).toEqual(['ξᵃ', 'Σ']);
    expect(eb.outgoing).toEqual(['ξ⁻ᵃ', 'Σ']);
    expect(eb.t).toBeCloseTo(1.3211, 4);
  });

  it('первичная волновая картина: γᵇ (0), Σ (t1), ξᵃ (t2), жёсткий слой [γᵇ, γᵃ] (t3)', () => {
    for (const k of ['I(a)', 'I(b)'] as const) {
      const bc = boundary([...CASES[k].ts], [...CASES[k].us]);
      const ev = solveBimodular({ a, b, ...bc, tEnd: 0.8 }).events;
      expect(ev.map((e) => e.outgoing.join(' '))).toEqual(['γᵇ', 'Σ', 'ξᵃ', 'γᵇ γᵃ']);
    }
  });

  it('исходы встреч — как в табл. 3 статьи (включая обе ветви инцидентов 4 и 5)', () => {
    const outcomes = new Set<string>();
    for (const k of ['I(a)', 'I(b)'] as const) {
      const bc = boundary([...CASES[k].ts], [...CASES[k].us]);
      for (const e of solveBimodular({ a, b, ...bc, tEnd: CASES[k].tEnd }).events) if (e.type !== 'boundary') outcomes.add(`${e.incoming.join('+')} → ${e.outgoing.join(' ')}`);
    }
    for (const row of [
      'Σ+γᵇ → ξ⁻ᵃ γᵃ', // 1
      'ξᵃ+Σ → ξ⁻ᵃ Σ', // 2
      'ξᵃ+ξ⁻ᵃ → ξ⁻ᵃ ξᵃ', // 3
      'γᵃ+ξ⁻ᵃ → γ⁻ᵃ ξᵃ', // 4(a)
      'γᵃ+ξ⁻ᵃ → γ⁻ᵇ γᵇ γᵃ', // 4(b)
      'γᵇ+γ⁻ᵃ → ξ⁻ᵇ γᵇ γᵃ', // 5(a)
      'γᵇ+γ⁻ᵃ → Σ⁻ ξᵃ', // 5(b)
      'γᵇ+γ⁻ᵇ → ξ⁻ᵇ ξᵇ', // 6
      'ξᵇ+Σ⁻ → ξ⁻ᵇ γᵇ γᵃ', // 7
      'ξ⁻ᵇ → ξᵇ', // 8 (отражение)
    ])
      expect(outcomes).toContain(row);
    // в табл. 3 нет встречи ξᵇ и ξ⁻ᵇ, хотя в варианте I(b) она происходит при t ≈ 2.68 мс
    expect(outcomes).toContain('ξᵇ+ξ⁻ᵇ → ξ⁻ᵇ ξᵇ');
  });

  it('критерии выбора исхода (10): для инцидента 4 — знак V (как в статье); для 5 — знак V + b·C_L (в статье ветви (a)/(b) перепутаны)', () => {
    const rows: { k: string; kind: string; V: number; N: number; VbC: number; out: string }[] = [];
    for (const k of ['I(a)', 'I(b)'] as const) {
      const bc = boundary([...CASES[k].ts], [...CASES[k].us]);
      const r = solveBimodular({ a, b, ...bc, tEnd: CASES[k].tEnd });
      for (const e of r.events) {
        const inc = e.incoming.join('+');
        if (e.type !== 'collision' || (inc !== 'γᵃ+ξ⁻ᵃ' && inc !== 'γᵇ+γ⁻ᵃ')) continue;
        const L = r.stateAt(e.x - 1e-6, e.t - 1e-9);
        const R = r.stateAt(e.x + 1e-6, e.t - 1e-9);
        const V = R.v - L.v + a * R.e;
        rows.push({ k, kind: inc, V, N: (V * V - b * b * L.e * L.e) / (V - a * L.e), VbC: V + b * L.e, out: e.outgoing.join(' ') });
      }
    }
    for (const q of rows) {
      if (q.kind === 'γᵃ+ξ⁻ᵃ') expect(q.out).toBe(q.V < 0 ? 'γ⁻ᵃ ξᵃ' : 'γ⁻ᵇ γᵇ γᵃ'); // статья: V < 0 ⇒ 4(a)
      else {
        expect(q.out).toBe(q.VbC > 0 ? 'ξ⁻ᵇ γᵇ γᵃ' : 'Σ⁻ ξᵃ'); // e* = (V + bC_L)/(2b)
        // статья: N < 0 ⇒ 5(a), N > 0 ⇒ 5(b) — в обоих вариантах получается наоборот
        expect(q.out).toBe(q.N < 0 ? 'Σ⁻ ξᵃ' : 'ξ⁻ᵇ γᵇ γᵃ');
      }
    }
    expect(rows.filter((q) => q.kind === 'γᵇ+γ⁻ᵃ').map((q) => q.out).sort()).toEqual(['Σ⁻ ξᵃ', 'ξ⁻ᵇ γᵇ γᵃ']);
  });

  for (const [name, cs] of Object.entries(CASES))
    it(`${name}: метод фронтов и гибрид совпадают с точным решением`, () => {
      const bc = boundary([...cs.ts], [...cs.us]);
      const ex = solveBimodular({ a, b, ...bc, tEnd: cs.tEnd });
      const X = a * cs.tEnd * 1.02;
      const err = (f: (x: number) => number) => {
        let e = 0;
        let n = 0;
        for (let i = 0; i < 8000; i++) {
          const x = (X * (i + 0.5)) / 8000;
          const q = ex.stateAt(x, cs.tEnd).e;
          e += Math.abs(f(x) - q);
          n += Math.abs(q);
        }
        return e / n;
      };
      const fr = trackFronts({ material: m, F0: [0, 0, 0], pulses: bc.pulses, tEnd: cs.tEnd, frames: 2, xMax: X * 1.1 });
      expect(fr.warnings).toEqual([]);
      expect(fr.events).toBe(ex.events.length);
      expect(err((x) => stateAt(fr.frames.at(-1)!, x).F[0])).toBeLessThan(1e-10);
      const hy = solveHybrid({ material: m, F0: [0, 0, 0], pulses: bc.pulses, tEnd: cs.tEnd, frames: 2, length: X, cellsPerL: 400 / X });
      expect(hy.warnings).toEqual([]);
      expect(hy.events).toBe(ex.events.length);
      expect(err((x) => hybridStateAt(hy.frames.at(-1)!, [0, 0, 0], x).F[0])).toBeLessThan(1e-6);
      expect(Math.abs(hy.kinematicError)).toBeLessThan(1e-14);
    }, 30000);
});
