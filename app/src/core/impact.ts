/**
 * Глава 2. Одномерные автомодельные задачи об ударе по границе x1 = 0 предварительно
 * деформированного разномодульного полупространства x1 > 0.
 *
 * Всё зависит от ξ = x1/(c1 t). Состояние: F = (p, g, s) = (u1,1, u2,1, u3,1),
 * скорость v̂ = v/c1 = (d, f, h). Напряжения отнесены к λ+2μ, скорости волн — к c1.
 *
 * Волны (по убыванию скорости от фронта к границе):
 *  - квазипродольная (fast): меняет p и величину сдвига |τ|, направление сдвига сохраняет;
 *  - квазипоперечная (slow): то же, но несёт в основном сдвиг;
 *  - нейтральная (волна поворота): поворачивает вектор сдвига (g, s), сохраняя p и |τ|
 *    ([p] = 0, [g² + s²] = 0, стр. 79); скорость ξ² = Φ(p, |τ|).
 * Квазипродольная и квазипоперечная волны бывают ударными (условия (2.40)–(2.42):
 * [T] = ξ²[F], [v̂] = −ξ[F]) или простыми (центрированный веер: ξ² — собственное значение
 * акустического тензора, (2.35)). Тип выбирается по допустимости (условие Лакса), а не задаётся.
 *
 * Граничные условия:
 *  - 'velocity' (§2.3): граница движется со скоростью v* = c1·(d*, 0, 0) — нормальный удар;
 *  - 'strain'   (§2.4): на границе заданы деформации (p*, g*, s*) ⇔ напряжения σ*_i1 — косой удар.
 */
import type { Material } from './material.ts';
import { type Vec3, energy1d, shearStiffness, traction, traction229 } from './material1d.ts';
import { newton } from './newton.ts';
import { boundaryProblem } from './fronts.ts';

export type Family = 'fast' | 'slow';
export type WaveKind = 'shock' | 'fan' | 'rotation';

export interface ImpactProblem {
  /** Модули, нормированные на λ+2μ. */
  material: Material;
  /** Предварительная деформация (покой). */
  F0: Vec3;
  /** Скорость среды впереди (по умолчанию покой). Нужна, например, для повторного удара по уже движущейся среде. */
  v0?: Vec3;
  bc: { type: 'velocity'; d: number } | { type: 'strain'; F: Vec3 };
  /** Какую формулу напряжений использовать: из потенциала (1.9) или напечатанную (2.29). */
  stressLaw?: 'model' | 'printed229';
  /** Внутреннее: не строить решение по волновым кривым fronts.ts (его вызывает сам fronts.ts). */
  noWaveCurves?: boolean;
}

export interface ProfilePoint {
  xi: number;
  F: Vec3;
  v: Vec3;
}

export interface Wave {
  kind: WaveKind;
  family: Family | 'neutral';
  /** Обозначение как в диссертации: Σ1, Σ2, Σ3 или «простая волна». */
  label: string;
  /** Для ударной и нейтральной волны xiFront = xiBack. */
  xiFront: number;
  xiBack: number;
  /** Профиль веера от переднего фронта к заднему. */
  profile?: ProfilePoint[];
  /** Диссипация (2.45), отнесённая к (λ+2μ)c1; для веера — 0. */
  dissipation: number;
  /** Условие Лакса для ударной волны: λ(впереди) < ξ² < λ(позади). */
  lax?: { ahead: number; behind: number; ok: boolean };
}

export interface Zone {
  name: string;
  F: Vec3;
  v: Vec3;
  T: Vec3;
}

export interface ImpactSolution {
  problem: ImpactProblem;
  /** Зоны от невозмущённой (0) к границе. zones[i] — перед waves[i], zones[i+1] — за ней. */
  zones: Zone[];
  waves: Wave[];
  /** Какие типы волн выбраны для квазипродольной и квазипоперечной: 'S' — ударная, 'R' — простая. */
  pattern: string;
  residual: number;
  warnings: string[];
}

// ---------- вспомогательные функции в плоскости поляризации ----------

type Vec2 = [number, number];

interface Plane {
  m: Material;
  e: Vec2; // единичное направление сдвига в плоскости (g, s)
  law: (m: Material, F: Vec3) => Vec3;
}

const to3 = (P: Plane, X: Vec2): Vec3 => [X[0], X[1] * P.e[0], X[1] * P.e[1]];

function T2(P: Plane, X: Vec2): Vec2 {
  const T = P.law(P.m, to3(P, X));
  return [T[0], T[1] * P.e[0] + T[2] * P.e[1]];
}

/** Акустическая матрица ∂T/∂X (2×2) и её собственные пары: [медленная, быстрая]. */
function eig(P: Plane, X: Vec2): { lam: Vec2; vec: [Vec2, Vec2] } {
  const H: number[][] = [
    [0, 0],
    [0, 0],
  ];
  for (let j = 0; j < 2; j++) {
    const h = 1e-7 * Math.max(Math.abs(X[j]), 1e-3);
    const xp = [...X] as Vec2;
    const xm = [...X] as Vec2;
    xp[j] += h;
    xm[j] -= h;
    const tp = T2(P, xp);
    const tm = T2(P, xm);
    for (let i = 0; i < 2; i++) H[i][j] = (tp[i] - tm[i]) / (2 * h);
  }
  const a = H[0][0];
  const d = H[1][1];
  const b = 0.5 * (H[0][1] + H[1][0]);
  const tr = (a + d) / 2;
  const disc = Math.sqrt(((a - d) / 2) ** 2 + b * b);
  const l1 = tr - disc;
  const l2 = tr + disc;
  const vecFor = (l: number): Vec2 => {
    const v: Vec2 = Math.abs(b) > 1e-14 ? [b, l - a] : Math.abs(l - a) < Math.abs(l - d) ? [1, 0] : [0, 1];
    const n = Math.hypot(v[0], v[1]);
    return [v[0] / n, v[1] / n];
  };
  return { lam: [l1, l2], vec: [vecFor(l1), vecFor(l2)] };
}

const famIdx = (f: Family) => (f === 'slow' ? 0 : 1);
export const speed = (lam: number) => Math.sqrt(Math.max(lam, 0));

/**
 * Простая волна: интегральная кривая семейства fam от состояния впереди Xa в направлении
 * убывания собственного значения (к границе ξ уменьшается), длина дуги σ ≥ 0.
 */
function fan(P: Plane, Xa: Vec2, fam: Family, sigma: number, n = 160): { Xb: Vec2; dv: Vec2; prof: { xi: number; X: Vec2; dv: Vec2 }[] } {
  const k = famIdx(fam);
  const e0 = eig(P, Xa);
  let dir = e0.vec[k];
  // направление, в котором λ убывает
  const h = 1e-6;
  const lp = eig(P, [Xa[0] + h * dir[0], Xa[1] + h * dir[1]]).lam[k];
  if (lp > e0.lam[k]) dir = [-dir[0], -dir[1]];
  const field = (X: Vec2, prev: Vec2): { dX: Vec2; xi: number } => {
    const E = eig(P, X);
    let r = E.vec[k];
    if (r[0] * prev[0] + r[1] * prev[1] < 0) r = [-r[0], -r[1]];
    return { dX: r, xi: speed(E.lam[k]) };
  };
  let X: Vec2 = [...Xa];
  let dv: Vec2 = [0, 0];
  const ds = sigma / n;
  const prof = [{ xi: speed(e0.lam[k]), X: [...X] as Vec2, dv: [0, 0] as Vec2 }];
  for (let i = 0; i < n; i++) {
    const k1 = field(X, dir);
    const k2 = field([X[0] + 0.5 * ds * k1.dX[0], X[1] + 0.5 * ds * k1.dX[1]], k1.dX);
    const k3 = field([X[0] + 0.5 * ds * k2.dX[0], X[1] + 0.5 * ds * k2.dX[1]], k2.dX);
    const k4 = field([X[0] + ds * k3.dX[0], X[1] + ds * k3.dX[1]], k3.dX);
    const dX: Vec2 = [
      (ds / 6) * (k1.dX[0] + 2 * k2.dX[0] + 2 * k3.dX[0] + k4.dX[0]),
      (ds / 6) * (k1.dX[1] + 2 * k2.dX[1] + 2 * k3.dX[1] + k4.dX[1]),
    ];
    // dv = −ξ dX (непрерывность скорости в веере)
    const xiAvg = (k1.xi + 2 * k2.xi + 2 * k3.xi + k4.xi) / 6;
    dv = [dv[0] - xiAvg * dX[0], dv[1] - xiAvg * dX[1]];
    X = [X[0] + dX[0], X[1] + dX[1]];
    dir = k4.dX;
    prof.push({ xi: speed(eig(P, X).lam[k]), X: [...X] as Vec2, dv: [...dv] as Vec2 });
  }
  return { Xb: X, dv, prof };
}

// ---------- решение цепочки «квазипродольная + квазипоперечная» ----------

interface ChainResult {
  pattern: string;
  x: number[];
  residual: number;
  /** Состояния и волны в плоскости: X0 → X1 → X2. */
  X: [Vec2, Vec2, Vec2];
  dv: [Vec2, Vec2];
  waves: { fam: Family; kind: 'shock' | 'fan'; xiFront: number; xiBack: number; prof?: { xi: number; X: Vec2; dv: Vec2 }[] }[];
}

function evalChain(P: Plane, X0: Vec2, pattern: string, x: number[]) {
  let i = 0;
  let Xa = X0;
  const res: number[] = [];
  const states: Vec2[] = [X0];
  const dvs: Vec2[] = [];
  const waves: ChainResult['waves'] = [];
  for (const [w, fam] of [
    [pattern[0], 'fast'],
    [pattern[1], 'slow'],
  ] as [string, Family][]) {
    if (w === 'S') {
      const Xb: Vec2 = [x[i], x[i + 1]];
      const xi = x[i + 2];
      i += 3;
      const Ta = T2(P, Xa);
      const Tb = T2(P, Xb);
      res.push(Ta[0] - Tb[0] - xi * xi * (Xa[0] - Xb[0]), Ta[1] - Tb[1] - xi * xi * (Xa[1] - Xb[1]));
      dvs.push([xi * (Xa[0] - Xb[0]), xi * (Xa[1] - Xb[1])]);
      waves.push({ fam, kind: 'shock', xiFront: xi, xiBack: xi });
      Xa = Xb;
    } else {
      const sigma = x[i];
      i += 1;
      if (!(sigma > 0)) throw new Error('fan: σ ≤ 0');
      const f = fan(P, Xa, fam, sigma);
      dvs.push(f.dv);
      waves.push({ fam, kind: 'fan', xiFront: f.prof[0].xi, xiBack: f.prof.at(-1)!.xi, prof: f.prof });
      Xa = f.Xb;
    }
    states.push(Xa);
  }
  return { res, states, dvs, waves };
}

function solveChain(P: Plane, X0: Vec2, end: (X2: Vec2, dvTotal: Vec2) => Vec2, guess: { X1: Vec2; X2: Vec2; xi1: number; xi2: number; s1: number; s2: number }, pattern: string): ChainResult | null {
  const x0: number[] = [];
  const scale: number[] = [];
  if (pattern[0] === 'S') x0.push(guess.X1[0], guess.X1[1], guess.xi1), scale.push(1e-3, 1e-3, 1);
  else x0.push(guess.s1), scale.push(1e-3);
  if (pattern[1] === 'S') x0.push(guess.X2[0], guess.X2[1], guess.xi2), scale.push(1e-3, 1e-3, 1);
  else x0.push(guess.s2), scale.push(1e-3);
  const F = (x: number[]) => {
    const c = evalChain(P, X0, pattern, x);
    const dvT: Vec2 = [c.dvs[0][0] + c.dvs[1][0], c.dvs[0][1] + c.dvs[1][1]];
    return [...c.res, ...end(c.states[2], dvT)];
  };
  const r = newton(F, x0, { scale, tol: 1e-15, maxIter: 80 });
  // Веер интегрируется численно (RK4 + разностный акустический тензор), поэтому невязку
  // ниже ~1e-14 получить нельзя; 1e-11 соответствует относительной точности ~1e-8.
  if (!r.converged && !(r.residual < 1e-11)) return null;
  const c = evalChain(P, X0, pattern, r.x);
  return { pattern, x: r.x, residual: r.residual, X: c.states as [Vec2, Vec2, Vec2], dv: c.dvs as [Vec2, Vec2], waves: c.waves };
}

function admissible(P: Plane, ch: ChainResult): string[] {
  const issues: string[] = [];
  ch.waves.forEach((w, i) => {
    const k = famIdx(w.fam);
    if (w.kind === 'shock') {
      const la = eig(P, ch.X[i]).lam[k];
      const lb = eig(P, ch.X[i + 1]).lam[k];
      const x2 = w.xiFront ** 2;
      if (!(la < x2 + 1e-12 && x2 < lb + 1e-12)) issues.push(`${w.fam}: Лакс нарушен`);
      if (!(w.xiFront > 0)) issues.push(`${w.fam}: ξ ≤ 0`);
    } else if (!(w.xiFront > w.xiBack)) issues.push(`${w.fam}: веер вырожден`);
  });
  if (!(ch.waves[0].xiBack >= ch.waves[1].xiFront - 1e-12)) issues.push('порядок волн нарушен');
  return issues;
}

// ---------- сборка полного решения ----------

const norm2 = (v: Vec2) => Math.hypot(v[0], v[1]);

export function solveImpact(prob: ImpactProblem, forcePattern?: string): ImpactSolution | null {
  const m = prob.material;
  const law = prob.stressLaw === 'printed229' ? (mm: Material, F: Vec3) => traction229(mm, F, 1) : traction;
  const F0 = prob.F0;
  const tau0v: Vec2 = [F0[1], F0[2]];
  const tau0 = norm2(tau0v);
  let eStar: Vec2 | null = null;
  let tauStar = 0;
  if (prob.bc.type === 'strain') {
    const tv: Vec2 = [prob.bc.F[1], prob.bc.F[2]];
    tauStar = norm2(tv);
    if (tauStar > 0) eStar = [tv[0] / tauStar, tv[1] / tauStar];
  }
  const e0: Vec2 = tau0 > 0 ? [tau0v[0] / tau0, tau0v[1] / tau0] : (eStar ?? [0, 1]);
  const P: Plane = { m, e: e0, law };
  const X0: Vec2 = [F0[0], tau0];

  let end: (X2: Vec2, dv: Vec2) => Vec2;
  let target: Vec2 | null = null;
  if (prob.bc.type === 'strain') {
    target = [prob.bc.F[0], tauStar];
    const t = target;
    end = (X2) => [X2[0] - t[0], X2[1] - t[1]];
  } else {
    // у границы v = v0 + Δv = (d, 0, 0); в плоскости поляризации: (v_n, v_τ)
    const d = prob.bc.d;
    const v0 = prob.v0 ?? [0, 0, 0];
    const need: Vec2 = [d - v0[0], -(v0[1] * e0[0] + v0[2] * e0[1])];
    if (Math.abs(v0[1] * e0[1] - v0[2] * e0[0]) > 1e-12)
      throw new Error('Скорость среды впереди должна лежать в плоскости предварительного сдвига');
    end = (_X2, dv) => [dv[0] - need[0], dv[1] - need[1]];
  }

  // Начальное приближение — линеаризация у X0
  const E0 = eig(P, X0);
  const rs = E0.vec[0];
  const rf = E0.vec[1];
  const xf = speed(E0.lam[1]);
  const xs = speed(E0.lam[0]);
  let A: number;
  let B: number;
  if (target) {
    const dX: Vec2 = [target[0] - X0[0], target[1] - X0[1]];
    const det = rf[0] * rs[1] - rs[0] * rf[1];
    A = (dX[0] * rs[1] - rs[0] * dX[1]) / det;
    B = (rf[0] * dX[1] - dX[0] * rf[1]) / det;
  } else {
    // Δv = −ξf A rf − ξs B rs (правило Крамера)
    const d = (prob.bc as { d: number }).d;
    const v0 = prob.v0 ?? [0, 0, 0];
    const need: Vec2 = [d - v0[0], -(v0[1] * e0[0] + v0[2] * e0[1])];
    const a11 = -xf * rf[0];
    const a12 = -xs * rs[0];
    const a21 = -xf * rf[1];
    const a22 = -xs * rs[1];
    const det = a11 * a22 - a12 * a21;
    A = (need[0] * a22 - a12 * need[1]) / det;
    B = (a11 * need[1] - need[0] * a21) / det;
  }
  const X1g: Vec2 = [X0[0] + A * rf[0], X0[1] + A * rf[1]];
  const X2g: Vec2 = target ?? [X1g[0] + B * rs[0], X1g[1] + B * rs[1]];
  const guess = { X1: X1g, X2: X2g, xi1: xf, xi2: xs, s1: Math.abs(A), s2: Math.abs(B) };

  const patterns = forcePattern ? [forcePattern] : ['SS', 'SR', 'RS', 'RR'];
  let best: { ch: ChainResult; issues: string[] } | null = null;
  for (const pat of patterns) {
    let ch: ChainResult | null = null;
    try {
      ch = solveChain(P, X0, end, guess, pat);
    } catch {
      ch = null;
    }
    if (!ch) continue;
    const issues = admissible(P, ch);
    if (!best || (best.issues.length > 0 && issues.length === 0)) best = { ch, issues };
    if (issues.length === 0) break;
  }
  if ((!best || best.issues.length > 0) && prob.bc.type === 'velocity' && !forcePattern && law === traction && !prob.noWaveCurves) {
    // Перебор «ударная волна или веер на семейство» не дал допустимого решения (например, при s0 = 0
    // и переходе p через ноль быстрая волна — два разрыва с полкой p = 0, а медленная вырождена).
    // Тогда решение строится по волновым кривым метода фронтов (fronts.ts), который годится при любом
    // состоянии, включая ненапряжённое.
    const v0 = prob.v0 ?? [0, 0, 0];
    const ch = chainFromWaveCurves(P, X0, [v0[0], v0[1] * e0[0] + v0[2] * e0[1]], prob.bc.d);
    if (ch) {
      const sol = assemble(prob, P, X0, eStar, ch, []);
      sol.warnings.unshift('Решение построено по волновым кривым (метод фронтов): перебор «ударная волна или гладкий веер» не дал допустимого решения.');
      return sol;
    }
  }
  if (!best) return null;
  return assemble(prob, P, X0, eStar, best.ch, best.issues);
}

/**
 * Краевая задача §2.3 по волновым кривым (fronts.ts → boundaryProblem) в виде цепочки
 * «квазипродольная + квазипоперечная». Волна семейства — либо одна ударная, либо «веер» с
 * профилем-лестницей: каждый слабый скачок — две точки с одинаковым ξ. Так же представляются
 * составные волны и вырожденные вееры (разрывы со скоростями c₋, c₊ и полка p = 0 между ними).
 */
function chainFromWaveCurves(P: Plane, X0: Vec2, v0: Vec2, d: number): ChainResult | null {
  const res = boundaryProblem({ m: P.m }, { F: X0, v: v0 }, d, 1e-4);
  if (!res) return null;
  const X: Vec2[] = [X0];
  const dv: Vec2[] = [];
  const waves: ChainResult['waves'] = [];
  let a = { F: X0, v: v0 };
  for (const [w, fam] of [
    [res.fast, 'fast'],
    [res.slow, 'slow'],
  ] as const) {
    const b = w.behind;
    // скачки от фронта к тылу; догоняющие (тыловой быстрее) сливаются в одну ударную волну
    type J = { G: number; A: { F: Vec2; v: Vec2 }; B: { F: Vec2; v: Vec2 } };
    const js: J[] = [];
    let cur0 = a;
    for (const pc of w.pieces) {
      js.push({ G: pc.G, A: cur0, B: pc.behind });
      cur0 = pc.behind;
    }
    for (let k = 0; k + 1 < js.length; ) {
      if (js[k + 1].G > js[k].G + 1e-12) {
        const A = js[k].A;
        const B = js[k + 1].B;
        const Ta = T2(P, A.F);
        const Tb = T2(P, B.F);
        const dF = [A.F[0] - B.F[0], A.F[1] - B.F[1]];
        const G = Math.sqrt(Math.max(((Ta[0] - Tb[0]) * dF[0] + (Ta[1] - Tb[1]) * dF[1]) / (dF[0] ** 2 + dF[1] ** 2 || 1e-300), 0));
        js.splice(k, 2, { G, A, B });
        if (k > 0) k--;
      } else k++;
    }
    if (js.length === 1) {
      waves.push({ fam, kind: 'shock', xiFront: js[0].G, xiBack: js[0].G });
    } else {
      const prof: { xi: number; X: Vec2; dv: Vec2 }[] = [];
      for (const j of js) {
        prof.push({ xi: j.G, X: [...j.A.F] as Vec2, dv: [j.A.v[0] - a.v[0], j.A.v[1] - a.v[1]] });
        prof.push({ xi: j.G, X: [...j.B.F] as Vec2, dv: [j.B.v[0] - a.v[0], j.B.v[1] - a.v[1]] });
      }
      waves.push({ fam, kind: 'fan', xiFront: prof[0].xi, xiBack: prof.at(-1)!.xi, prof });
    }
    dv.push([b.v[0] - a.v[0], b.v[1] - a.v[1]]);
    X.push([...b.F] as Vec2);
    a = b;
  }
  const err = Math.hypot(a.v[0] - d, a.v[1]);
  return { pattern: waves.map((w) => (w.kind === 'shock' ? 'S' : 'R')).join(''), x: [], residual: err, X: X as [Vec2, Vec2, Vec2], dv: dv as [Vec2, Vec2], waves };
}

function assemble(prob: ImpactProblem, P: Plane, X0: Vec2, eStar: Vec2 | null, ch: ChainResult, issues: string[]): ImpactSolution {
  const m = P.m;
  const law = P.law;
  const zones: Zone[] = [];
  const waves: Wave[] = [];
  const warnings = issues.map((s) => `Решение неполностью допустимо: ${s}.`);
  const mk = (name: string, F: Vec3, v: Vec3): Zone => ({ name, F, v, T: law(m, F) });
  const e0 = P.e;
  const to3e = (X: Vec2, e: Vec2): Vec3 => [X[0], X[1] * e[0], X[1] * e[1]];
  const dv3 = (dv: Vec2, e: Vec2): Vec3 => [dv[0], dv[1] * e[0], dv[1] * e[1]];
  const add = (a: Vec3, b: Vec3): Vec3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const romans = ['0', 'I', 'II', 'III', 'IV'];

  // Нужна ли нейтральная волна и где она стоит
  let rot: 'none' | 'after-fast' | 'after-slow' = 'none';
  let angle = 0;
  if (eStar && prob.bc.type === 'strain') {
    angle = Math.atan2(e0[0] * eStar[1] - e0[1] * eStar[0], e0[0] * eStar[0] + e0[1] * eStar[1]);
    if (Math.abs(angle) > 1e-12) {
      const xiA = speed(shearStiffness(m, to3e(ch.X[2], e0))); // у границы
      const xiB = speed(shearStiffness(m, to3e(ch.X[1], e0))); // между Σ1 и квазипоперечной
      const slowBack = ch.waves[1].xiBack;
      const slowFront = ch.waves[1].xiFront;
      if (xiA <= slowBack + 1e-12) rot = 'after-slow';
      else if (xiB >= slowFront - 1e-12 && xiB <= ch.waves[0].xiBack + 1e-12) rot = 'after-fast';
      else {
        rot = xiA - slowBack < slowFront - xiB ? 'after-slow' : 'after-fast';
        warnings.push('Нейтральная волна не укладывается между соседними волнами — порядок волн не согласован.');
      }
    }
  }

  let F: Vec3 = prob.F0;
  let v: Vec3 = prob.v0 ? [...prob.v0] : [0, 0, 0];
  let e = e0;
  zones.push(mk('0', F, v));
  const pushZone = () => zones.push(mk(romans[zones.length] ?? String(zones.length), F, v));
  const lawLam = (Pp: Plane, X: Vec2, k: 0 | 1) => eig(Pp, X).lam[k];

  const addRotation = (X: Vec2) => {
    const Fa = F;
    const va = v;
    const xi = speed(shearStiffness(m, Fa));
    const Fb = to3e(X, eStar!);
    const vb = add(va, [xi * (Fa[0] - Fb[0]), xi * (Fa[1] - Fb[1]), xi * (Fa[2] - Fb[2])]);
    F = Fb;
    v = vb;
    e = eStar!;
    waves.push({ kind: 'rotation', family: 'neutral', label: 'Σ3', xiFront: xi, xiBack: xi, dissipation: dissipation(m, law, mk('a', Fa, va), mk('b', Fb, vb), xi) });
    pushZone();
  };

  ch.waves.forEach((w, i) => {
    const Pe: Plane = { ...P, e };
    const Fa = F;
    const va = v;
    if (w.kind === 'shock') {
      F = to3e(ch.X[i + 1], e);
      v = add(va, dv3(ch.dv[i], e));
      const k = famIdx(w.fam);
      const la = lawLam(Pe, ch.X[i], k as 0 | 1);
      const lb = lawLam(Pe, ch.X[i + 1], k as 0 | 1);
      waves.push({
        kind: 'shock',
        family: w.fam,
        label: w.fam === 'fast' ? 'Σ1' : 'Σ2',
        xiFront: w.xiFront,
        xiBack: w.xiBack,
        dissipation: dissipation(m, law, mk('a', Fa, va), mk('b', F, v), w.xiFront),
        // волна нулевой амплитуды (например, квазипоперечная при s = 0) условию не подлежит
        lax: { ahead: speed(la), behind: speed(lb), ok: (la < w.xiFront ** 2 + 1e-12 && w.xiFront ** 2 < lb + 1e-12) || Math.hypot(ch.X[i + 1][0] - ch.X[i][0], ch.X[i + 1][1] - ch.X[i][1]) < 1e-12 },
      });
    } else {
      const prof = w.prof!.map((q) => ({ xi: q.xi, F: to3e(q.X, e), v: add(va, dv3(q.dv, e)) }));
      F = prof.at(-1)!.F;
      v = prof.at(-1)!.v;
      waves.push({ kind: 'fan', family: w.fam, label: 'простая волна', xiFront: w.xiFront, xiBack: w.xiBack, profile: prof, dissipation: 0 });
      // Невыпуклость: если скорость внутри веера не убывает монотонно, настоящее решение —
      // составная волна (веер + примыкающая ударная), а не веер по всей интегральной кривой.
      if (prof.some((q, j) => j > 0 && q.xi > prof[j - 1].xi + 1e-9))
        warnings.push(
          `Скорость внутри простой волны (${w.fam === 'fast' ? 'квазипродольной' : 'квазипоперечной'}) немонотонна — семейство невыпукло; точное решение здесь — составная волна «веер + ударная» (её строит метод фронтов на вкладке «Серия ударов»).`,
        );
    }
    pushZone();
    if (i === 0 && rot === 'after-fast') addRotation(ch.X[1]);
  });
  if (rot === 'after-slow') addRotation(ch.X[2]);

  for (const w of waves) if (w.kind === 'shock' && w.dissipation < -1e-13) warnings.push(`На ${w.label} нарушено термодинамическое условие (2.45): D = ${w.dissipation.toExponential(2)}.`);
  if (prob.bc.type === 'strain') {
    const Fb = zones.at(-1)!.F;
    const err = Math.hypot(Fb[0] - prob.bc.F[0], Fb[1] - prob.bc.F[1], Fb[2] - prob.bc.F[2]);
    if (err > 1e-10) warnings.push(`Граничное условие выполнено с погрешностью ${err.toExponential(1)}.`);
  }
  return { problem: prob, zones, waves, pattern: ch.pattern, residual: ch.residual, warnings };
}

/** Диссипация на разрыве (2.45): D = [T·v̂] + ξ[W + ½|v̂|²], [·] = впереди − позади. */
function dissipation(m: Material, _law: (m: Material, F: Vec3) => Vec3, a: Zone, b: Zone, xi: number): number {
  const flux = (z: Zone) => z.T[0] * z.v[0] + z.T[1] * z.v[1] + z.T[2] * z.v[2];
  const E = (z: Zone) => energy1d(m, z.F) + 0.5 * (z.v[0] ** 2 + z.v[1] ** 2 + z.v[2] ** 2);
  return flux(a) - flux(b) + xi * (E(a) - E(b));
}

// ---------- профили и интегральные проверки ----------

/** Значения F, v̂ в точке ξ ≥ 0. */
export function profileAt(sol: ImpactSolution, xi: number): { F: Vec3; v: Vec3 } {
  const { zones, waves } = sol;
  for (let i = 0; i < waves.length; i++) {
    const w = waves[i];
    if (xi > w.xiFront) return zones[i];
    if (w.kind === 'fan' && xi >= w.xiBack) {
      const pr = w.profile!;
      for (let j = 1; j < pr.length; j++)
        if (xi >= pr[j].xi) {
          const t = (xi - pr[j].xi) / (pr[j - 1].xi - pr[j].xi || 1);
          const lerp = (a: Vec3, b: Vec3): Vec3 => [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])];
          return { F: lerp(pr[j].F, pr[j - 1].F), v: lerp(pr[j].v, pr[j - 1].v) };
        }
      return { F: pr.at(-1)!.F, v: pr.at(-1)!.v };
    }
  }
  return zones.at(-1)!;
}

/**
 * Безразмерные перемещения относительно предварительно деформированного состояния:
 * Δu(x1, t) = c1 t · w(ξ),  w(ξ) = −∫_ξ^∞ (F − F0) dξ'.  Возвращает функцию w(ξ).
 */
export function displacementTable(sol: ImpactSolution, n = 1200): (xi: number) => Vec3 {
  const xmax = sol.waves[0].xiFront * 1.0001;
  const xs: number[] = [];
  for (let i = 0; i <= n; i++) xs.push((xmax * i) / n);
  // добавляем точные положения разрывов
  for (const w of sol.waves) {
    xs.push(w.xiFront, w.xiBack, w.xiFront + 1e-12, w.xiBack - 1e-12);
    // узлы профиля веера — чтобы интеграл внутри него был точным
    if (w.profile) for (let j = 0; j < w.profile.length - 1; j++) xs.push(w.profile[j].xi, 0.5 * (w.profile[j].xi + w.profile[j + 1].xi));
  }
  xs.sort((a, b) => a - b);
  const uniq = xs.filter((x, i) => x >= 0 && (i === 0 || x - xs[i - 1] > 1e-13));
  const W: Vec3[] = new Array(uniq.length);
  W[uniq.length - 1] = [0, 0, 0];
  const F0 = sol.problem.F0;
  for (let i = uniq.length - 2; i >= 0; i--) {
    const a = uniq[i];
    const b = uniq[i + 1];
    const fa = profileAt(sol, a + 1e-14).F;
    const fb = profileAt(sol, b - 1e-14).F;
    W[i] = [0, 1, 2].map((k) => W[i + 1][k] - 0.5 * (b - a) * (fa[k] - F0[k] + fb[k] - F0[k])) as Vec3;
  }
  return (xi: number) => {
    if (xi >= xmax) return [0, 0, 0];
    let lo = 0;
    let hi = uniq.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (uniq[mid] <= xi) lo = mid;
      else hi = mid;
    }
    const t = (xi - uniq[lo]) / (uniq[hi] - uniq[lo] || 1);
    return [0, 1, 2].map((k) => W[lo][k] + t * (W[hi][k] - W[lo][k])) as Vec3;
  };
}

/**
 * Интегральные законы сохранения для автомодельного решения (независимая проверка сборки,
 * включая вееры):
 *  кинематика: −∫(F − F0)dξ = v̂(0);
 *  импульс:    ∫v̂ dξ = T(F0) − T(F(0));
 *  энергия:    ∫(W + ½|v̂|² − W0)dξ = −T(F(0))·v̂(0) − ΣD.
 */
export function conservationCheck(sol: ImpactSolution, n = 4000) {
  const m = sol.problem.material;
  const law = sol.problem.stressLaw === 'printed229' ? (mm: Material, F: Vec3) => traction229(mm, F, 1) : traction;
  const F0 = sol.problem.F0;
  const W0 = energy1d(m, F0);
  const T0 = law(m, F0);
  const xmax = sol.waves[0].xiFront;
  const breaks = [0, ...sol.waves.flatMap((w) => [w.xiBack, w.xiFront]), xmax].sort((a, b) => a - b);
  const kin: Vec3 = [0, 0, 0];
  const mom: Vec3 = [0, 0, 0];
  let en = 0;
  for (let s = 0; s < breaks.length - 1; s++) {
    const a = breaks[s];
    const b = breaks[s + 1];
    if (b - a < 1e-14) continue;
    const k = Math.max(2, Math.round((n * (b - a)) / xmax));
    for (let i = 0; i < k; i++) {
      // правило средней точки внутри гладкого участка
      const x = a + ((i + 0.5) * (b - a)) / k;
      const { F, v } = profileAt(sol, x);
      const h = (b - a) / k;
      for (let j = 0; j < 3; j++) {
        kin[j] -= h * (F[j] - F0[j]);
        mom[j] += h * v[j];
      }
      en += h * (energy1d(m, F) + 0.5 * (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) - W0);
    }
  }
  const zb = sol.zones.at(-1)!;
  const Tb = zb.T;
  const D = sol.waves.reduce((acc, w) => acc + w.dissipation, 0);
  return {
    kinematic: { lhs: kin, rhs: zb.v },
    momentum: { lhs: mom, rhs: [T0[0] - Tb[0], T0[1] - Tb[1], T0[2] - Tb[2]] as Vec3 },
    energy: { lhs: en, rhs: -(Tb[0] * zb.v[0] + Tb[1] * zb.v[1] + Tb[2] * zb.v[2]) - D },
  };
}
