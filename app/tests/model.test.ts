/**
 * Численная проверка формул статьи docs/model/model.tex там, где нужен анализ, а не алгебра
 * (алгебраические тождества проверены в Lean: docs/model/lean).
 */
import { describe, expect, it } from 'vitest';
import { boundaryProblem, riemannProblem } from '../src/core/fronts.ts';
import { linearBoundaryState, solveHybrid, sqrtSym } from '../src/core/hybrid.ts';
import { normalize } from '../src/core/material.ts';
import { energy1d, stress3, traction } from '../src/core/material1d.ts';
import { MATERIAL_1, MATERIAL_2 } from '../src/core/presets.ts';

const m1 = normalize(MATERIAL_1);
const m2 = normalize(MATERIAL_2);
type V2 = [number, number];
const T2 = (m: typeof m1, F: V2): V2 => {
  const t = traction(m, [F[0], 0, F[1]]);
  return [t[0], t[2]];
};
const W = (m: typeof m1, F: V2) => energy1d(m, [F[0], 0, F[1]]);
const acoustic = (m: typeof m1, F: V2): [number, number, number] => {
  const h = 1e-8;
  const a = (T2(m, [F[0] + h, F[1]])[0] - T2(m, [F[0] - h, F[1]])[0]) / (2 * h);
  const b = (T2(m, [F[0], F[1] + h])[0] - T2(m, [F[0], F[1] - h])[0]) / (2 * h);
  const d = (T2(m, [F[0], F[1] + h])[1] - T2(m, [F[0], F[1] - h])[1]) / (2 * h);
  return [a, b, d];
};
const states: V2[] = [
  [-0.0005, 0.001],
  [-0.00237, 0.00094],
  [0.00158, 0.00125],
  [0.0001, 0.001],
  [-0.005, 0.001],
];

describe('формулы статьи: материал', () => {
  for (const [name, m] of [['материал 1', m1], ['материал 2', m2]] as const)
    it(`${name}: T = ∂W/∂F (гиперупругость, (2.2)) и T = (σ11, σ31) из общего закона (1.10)`, () => {
      for (const F of states) {
        const h = 1e-9;
        const T = T2(m, F);
        const dW0 = (W(m, [F[0] + h, F[1]]) - W(m, [F[0] - h, F[1]])) / (2 * h);
        const dW1 = (W(m, [F[0], F[1] + h]) - W(m, [F[0], F[1] - h])) / (2 * h);
        expect(Math.abs(dW0 - T[0])).toBeLessThan(1e-9 * (Math.abs(T[0]) + 1e-6));
        expect(Math.abs(dW1 - T[1])).toBeLessThan(1e-7 * (Math.abs(T[0]) + 1e-6));
        const S = stress3(m, [F[0], 0, F[1]]);
        expect(S.s11).toBeCloseTo(T[0], 15);
        expect(S.s13).toBeCloseTo(T[1], 15);
        // акустическая матрица симметрична (вторые производные потенциала) и положительно определена
        const h2 = 1e-8;
        const a12 = (T2(m, [F[0], F[1] + h2])[0] - T2(m, [F[0], F[1] - h2])[0]) / (2 * h2);
        const a21 = (T2(m, [F[0] + h2, F[1]])[1] - T2(m, [F[0] - h2, F[1]])[1]) / (2 * h2);
        expect(Math.abs(a12 - a21)).toBeLessThan(1e-6);
        const [a, b, d] = acoustic(m, F);
        expect(a).toBeGreaterThan(0);
        expect(a * d - b * b).toBeGreaterThan(0);
      }
    });

  it('материал 2: скорость продольных волн при переходе через p = 0 меняется на ≈17 % (разномодульность)', () => {
    const cf = (F: V2) => {
      const [a, b, d] = acoustic(m2, F);
      return Math.sqrt((a + d) / 2 + Math.sqrt(((a - d) / 2) ** 2 + b * b));
    };
    const comp = cf([-0.01, 0.001]);
    const ten = cf([0.01, 0.001]);
    expect(comp).toBeGreaterThan(1.07);
    expect(ten).toBeLessThan(0.92);
  });
});

describe('формулы статьи: условия на разрывах и диссипация', () => {
  it('ударные волны краевой задачи: [T] = G²[F], [v] = −G[F], D = [T·v] + G[E] = G([W] − ⟨T⟩·[F]) ≥ 0', () => {
    const c = { m: m2 };
    const res = boundaryProblem(c, { F: [-0.0005, 0.001], v: [0, 0] }, 0.002, 2.5e-4)!;
    let ahead = { F: [-0.0005, 0.001] as V2, v: [0, 0] as V2 };
    let n = 0;
    for (const w of [res.fast, res.slow])
      for (const pc of w.pieces) {
        const b = pc.behind;
        const G = pc.G;
        const dF = [ahead.F[0] - b.F[0], ahead.F[1] - b.F[1]]; // [·] = справа − слева: справа — впереди
        const dv = [ahead.v[0] - b.v[0], ahead.v[1] - b.v[1]];
        const Ta = T2(m2, ahead.F);
        const Tb = T2(m2, b.F);
        const dT = [Ta[0] - Tb[0], Ta[1] - Tb[1]];
        for (let k = 0; k < 2; k++) {
          expect(Math.abs(dT[k] - G * G * dF[k])).toBeLessThan(1e-12);
          expect(Math.abs(dv[k] + G * dF[k])).toBeLessThan(1e-12);
        }
        const E = (s: typeof ahead) => W(m2, s.F) + 0.5 * (s.v[0] ** 2 + s.v[1] ** 2);
        const D1 = Ta[0] * ahead.v[0] + Ta[1] * ahead.v[1] - Tb[0] * b.v[0] - Tb[1] * b.v[1] + G * (E(ahead) - E(b));
        const D2 = G * (W(m2, ahead.F) - W(m2, b.F) - (0.5 * (Ta[0] + Tb[0]) * dF[0] + 0.5 * (Ta[1] + Tb[1]) * dF[1]));
        expect(Math.abs(D1 - D2)).toBeLessThan(1e-13 * Math.max(1e-12, Math.abs(D1)) + 1e-18);
        if (pc.kind === 'shock') expect(D1).toBeGreaterThan(0);
        ahead = b;
        n++;
      }
    expect(n).toBe(2);
  });

  it('волны разрежения как разрывы недопустимы: «ударная волна растяжения» дала бы D < 0', () => {
    const c = { m: m2 };
    // рывок наружу: быстрая волна — веер; если бы её крайние состояния соединял разрыв, D < 0
    const res = boundaryProblem(c, { F: [-0.0005, 0.001], v: [0, 0] }, -0.002, 1e-3)!;
    const a = { F: [-0.0005, 0.001] as V2, v: [0, 0] as V2 };
    const b = res.fast.behind;
    expect(res.fast.pieces.every((p) => p.kind === 'simple')).toBe(true);
    const dF = [a.F[0] - b.F[0], a.F[1] - b.F[1]];
    const Ta = T2(m2, a.F);
    const Tb = T2(m2, b.F);
    const G = Math.sqrt(((Ta[0] - Tb[0]) * dF[0] + (Ta[1] - Tb[1]) * dF[1]) / (dF[0] ** 2 + dF[1] ** 2));
    const D = G * (W(m2, a.F) - W(m2, b.F) - (0.5 * (Ta[0] + Tb[0]) * dF[0] + 0.5 * (Ta[1] + Tb[1]) * dF[1]));
    expect(D).toBeLessThan(0);
  });
});

describe('формулы статьи: линейные пределы точных задач', () => {
  it('C = √A: C² = A; граничное состояние F_b = F − C⁻¹(v_b − v) — предел краевой задачи §2.3', () => {
    for (const F of states) {
      const A = acoustic(m2, F);
      const [c11, c12, c22] = sqrtSym(...A);
      expect(c11 * c11 + c12 * c12).toBeCloseTo(A[0], 12);
      expect(c11 * c12 + c12 * c22).toBeCloseTo(A[1], 12);
      expect(c12 * c12 + c22 * c22).toBeCloseTo(A[2], 12);
      const dv = 1e-7;
      const lin = linearBoundaryState(A, [F[0], F[1], 0, 0], dv);
      const ex = boundaryProblem({ m: m2 }, { F, v: [0, 0] }, dv, 1e-3)!.slow.behind;
      // отличие — второго порядка по Δv (нелинейность), а не первого
      expect(Math.abs(lin[0] - ex.F[0])).toBeLessThan(1e-3 * Math.abs(ex.F[0] - F[0]));
      expect(Math.abs(lin[1] - ex.F[1])).toBeLessThan(1e-3 * Math.abs(ex.F[0] - F[0]) + 1e-12);
    }
  });

  it('распад малого разрыва: амплитуды волн — проекции на собственные векторы (явная формула статьи)', () => {
    const F: V2 = [-0.002, 0.001];
    const dF: V2 = [3e-8, -2e-8];
    const dv: V2 = [1e-8, 4e-8];
    const rp = riemannProblem({ m: m2 }, { F, v: [0, 0] }, { F: [F[0] + dF[0], F[1] + dF[1]], v: dv }, 1e-3)!;
    const [a, b, d] = acoustic(m2, F);
    const D = Math.sqrt(((a - d) / 2) ** 2 + b * b);
    const lf = (a + d) / 2 + D;
    const ls = (a + d) / 2 - D;
    const nf = Math.hypot(b, lf - a);
    const rf: V2 = [b / nf, (lf - a) / nf];
    const rs: V2 = [-rf[1], rf[0]];
    const cf = Math.sqrt(lf);
    const cs = Math.sqrt(ls);
    const P = (u: V2, r: V2) => u[0] * r[0] + u[1] * r[1];
    const amp = {
      LF: (P(dF, rf) + P(dv, rf) / cf) / 2,
      RF: (P(dv, rf) / cf - P(dF, rf)) / 2,
      LS: (P(dF, rs) + P(dv, rs) / cs) / 2,
      RS: (P(dv, rs) / cs - P(dF, rs)) / 2,
    };
    const jump = (ahead: { F: V2 }, w: { behind: { F: V2 } }) => Math.hypot(w.behind.F[0] - ahead.F[0], w.behind.F[1] - ahead.F[1]);
    const l = { F };
    const r = { F: [F[0] + dF[0], F[1] + dF[1]] as V2 };
    expect(jump(l, rp.LF)).toBeCloseTo(Math.abs(amp.LF), 12);
    expect(jump(rp.LF.behind, rp.LS)).toBeCloseTo(Math.abs(amp.LS), 12);
    expect(jump(r, rp.RF)).toBeCloseTo(Math.abs(amp.RF), 12);
    expect(jump(rp.RF.behind, rp.RS)).toBeCloseTo(Math.abs(amp.RS), 12);
  });
});

describe('разбор варианта статьи (материал 2, p0 = −0.0005: удар — рывок — удар)', () => {
  it('структура волн и балансы', () => {
    const pulses = [
      { start: 0, duration: 0.25, d: 0.002 },
      { start: 0.25, duration: 0.25, d: -0.002 },
      { start: 0.5, duration: 0.25, d: 0.002 },
    ];
    const r = solveHybrid({ material: m2, F0: [-0.0005, 0, 0.001], pulses, tEnd: 1.4, length: 1.5, frames: 2, cellsPerL: 300 });
    expect(r.warnings).toEqual([]);
    const ev = r.log.filter((e) => e.type === 'boundary' || e.type === 'collision');
    const kinds = (e: (typeof ev)[number]) => e.waves!.map((w) => `${w.fam}:${w.kind}`).join(',');
    // t = 0: удар — две ударные волны сжатия
    expect(ev[0].t).toBe(0);
    expect(kinds(ev[0])).toBe('slow:shock,fast:shock');
    // t = 0.25: рывок в растяжение — быстрый веер (через p = 0) и медленная ударная волна, обе выделены
    expect(ev[1].t).toBeCloseTo(0.25, 12);
    expect(kinds(ev[1])).toBe('slow:shock,fast:simple');
    expect(ev[1].waves!.every((w) => w.tracked)).toBe(true);
    expect(ev[1].waves![1].L.F[0]).toBeGreaterThan(0); // за веером — растяжение
    // t = 0.5: удар из растяжения в сжатие — снова две ударные волны
    expect(kinds(ev[2])).toBe('slow:shock,fast:shock');
    // быстрая ударная волна 3-го удара догоняет медленную ударную волну рывка: медленная волна
    // попадает из растяжения в сжатие и становится веером (смена типа), отражения слабые
    const col = ev.find((e) => e.type === 'collision')!;
    expect(col.t).toBeGreaterThan(0.7);
    expect(col.t).toBeLessThan(0.77);
    expect(col.incoming!.map((w) => `${w.fam}:${w.kind}`)).toEqual(['fast:shock', 'slow:shock']);
    const right = col.waves!.filter((w) => w.dir > 0).map((w) => `${w.fam}:${w.kind}`);
    expect(right).toEqual(['slow:simple', 'fast:shock']);
    for (const w of col.waves!.filter((q) => q.dir < 0)) expect(w.strength).toBeLessThan(1e-4 * col.incoming![0].strength * 1e2);
    // t = 0.75: остановка — два веера разгрузки (на сетке)
    const stop = ev.find((e) => e.type === 'boundary' && Math.abs(e.t - 0.75) < 1e-9)!;
    expect(kinds(stop)).toBe('slow:simple,fast:simple');
    expect(Math.abs(r.kinematicError)).toBeLessThan(1e-12);
  });
});
