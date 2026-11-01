/**
 * Метод волновых фронтов для серии нормальных ударов (по построениям главы 2 диссертации).
 *
 * Решение — кусочно-постоянное: между плоскими фронтами состояние (p, s, v1, v3) постоянно.
 * На каждом фронте выполняются условия на разрывах (2.40)–(2.42):
 *     [T] = G²[F],   [v] = −G[F]      ([·] = справа − слева, G — скорость со знаком),
 * где F = (p, s) = (u1,1, u3,1), T = (σ11, σ31)/(λ+2μ), скорости — в единицах c1.
 *
 * События:
 *  1) изменение скорости границы — краевая задача §2.3 для состояния у границы
 *     (рождаются квазипродольная и квазипоперечная волны, бегущие вглубь);
 *  2) встреча двух фронтов — задача о распаде разрыва между состояниями слева и справа
 *     (по две волны уходят влево и вправо);
 *  3) фронт, бегущий к границе, доходит до неё — снова задача §2.3 (отражение).
 * Тип волны выбирается по условию Лакса: сжимающая (λ растёт к «тылу») — ударная волна,
 * иначе — простая волна, которая представляется пучком слабых фронтов, каждый из которых
 * тоже точно удовлетворяет условиям на разрыве (стандартная техника «front tracking»).
 *
 * Ограничение: нормальный удар, предварительный сдвиг в плоскости x1x3 (g0 = 0) — тогда
 * все волны поляризованы в этой плоскости и нейтральная волна не возникает.
 */
import type { Material } from './material.ts';
import { type Vec3, energy1d, traction } from './material1d.ts';
import { solveLinear } from './newton.ts';
import type { Pulse } from './fv1d.ts';
import { boundaryVelocity } from './fv1d.ts';

type V2 = [number, number];

export type Fam = 'fast' | 'slow';

export interface State {
  F: V2; // (p, s)
  v: V2; // (v1, v3) / c1
}

export interface Front {
  x: number;
  /** Скорость со знаком (вправо > 0), в единицах c1. */
  speed: number;
  fam: Fam;
  /** 'shock' — ударная волна; 'simple' — один из фронтов пучка простой волны. */
  kind: 'shock' | 'simple';
  born: number;
  xBorn: number;
}

export interface Segment {
  x0: number;
  t0: number;
  x1: number;
  t1: number;
  fam: Fam;
  kind: 'shock' | 'simple';
  strength: number;
}

export interface FrontsProblem {
  material: Material;
  F0: Vec3; // (p0, 0, s0)
  pulses: Pulse[];
  tEnd: number;
  /** Сколько кадров (моментов времени) сохранить для отображения. */
  frames: number;
  /** Максимальный скачок деформации на одном фронте пучка простой волны. */
  fanStep?: number;
  /** Правая граница отслеживания (фронты дальше удаляются). */
  xMax?: number;
  /** Порог отбрасывания слабых фронтов (скачок деформации). */
  dropTol?: number;
}

export interface FrontsFrame {
  t: number;
  /** Положения фронтов (по возрастанию) и состояния между ними (states.length = xs.length + 1). */
  xs: number[];
  states: State[];
  ub: number;
}

export interface FrontsResult {
  frames: FrontsFrame[];
  segments: Segment[];
  events: number;
  maxFronts: number;
  warnings: string[];
  /** Балансы на конечный момент: кинематика, энергия (работа границы, прирост энергии, диссипация на фронтах). */
  kinematicError: number;
  work: number;
  energyGain: number;
  dissipation: number;
}

// ---------------------------------------------------------------- механика фронтов

interface Ctx {
  m: Material;
}

function T2(c: Ctx, X: V2): V2 {
  const T = traction(c.m, [X[0], 0, X[1]]);
  return [T[0], T[2]];
}

/** Акустическая матрица и собственные пары [медленная, быстрая]. */
function eig(c: Ctx, X: V2): { lam: V2; vec: [V2, V2] } {
  const H = [
    [0, 0],
    [0, 0],
  ];
  for (let j = 0; j < 2; j++) {
    const h = 1e-7 * Math.max(Math.abs(X[j]), 1e-4);
    const xp = [...X] as V2;
    const xm = [...X] as V2;
    xp[j] += h;
    xm[j] -= h;
    const a = T2(c, xp);
    const b = T2(c, xm);
    for (let i = 0; i < 2; i++) H[i][j] = (a[i] - b[i]) / (2 * h);
  }
  const A = H[0][0];
  const D = H[1][1];
  const B = 0.5 * (H[0][1] + H[1][0]);
  const tr = (A + D) / 2;
  const disc = Math.sqrt(((A - D) / 2) ** 2 + B * B);
  const vec = (l: number): V2 => {
    const v: V2 = Math.abs(B) > 1e-15 ? [B, l - A] : Math.abs(l - A) < Math.abs(l - D) ? [1, 0] : [0, 1];
    const n = Math.hypot(v[0], v[1]);
    return [v[0] / n, v[1] / n];
  };
  return { lam: [tr - disc, tr + disc], vec: [vec(tr - disc), vec(tr + disc)] };
}

const fi = (f: Fam) => (f === 'slow' ? 0 : 1);

/** Собственный вектор семейства, ориентированный в сторону роста λ (сжимающее направление). */
function orientedR(c: Ctx, X: V2, f: Fam): { r: V2; lam: number; gnl: boolean } {
  const E = eig(c, X);
  const k = fi(f);
  let r = E.vec[k];
  // шаг должен быть много больше шага внутри eig (иначе знак производной λ тонет в шуме)
  const h = 1e-3 * Math.max(Math.hypot(X[0], X[1]), 1e-5);
  const lp = eig(c, [X[0] + h * r[0], X[1] + h * r[1]]).lam[k];
  const lm = eig(c, [X[0] - h * r[0], X[1] - h * r[1]]).lam[k];
  const g = (lp - lm) / (2 * h);
  if (g < 0) r = [-r[0], -r[1]];
  return { r, lam: E.lam[k], gnl: Math.abs(g) > 1e-9 * Math.max(E.lam[k], 1e-6) };
}

/**
 * Ударная адиабата (условия (2.40)–(2.42)) из состояния впереди Xa для семейства f:
 * Xb = Xa + ε(r + η t), где η находится из коллинеарности [T] и [F]; G² = [T]·[F]/|[F]|².
 * ε > 0 — сжимающий (допустимый по Лаксу) скачок, ε < 0 — слабый скачок пучка простой волны.
 */
function hugoniot(c: Ctx, Xa: V2, f: Fam, eps: number, R?: { r: V2 }): { Xb: V2; G2: number } {
  const r = (R ?? orientedR(c, Xa, f)).r;
  const t: V2 = [-r[1], r[0]];
  const Ta = T2(c, Xa);
  const D = (eta: number): { d: V2; q: V2 } => {
    const q: V2 = [r[0] + eta * t[0], r[1] + eta * t[1]];
    if (Math.abs(eps) < 1e-9) {
      // предел ε → 0: акустическая матрица на q (для столь слабых скачков ошибка O(ε²) пренебрежима,
      // а прямое вычисление (T(Xa + εq) − T(Xa))/ε теряет точность из-за вычитания)
      const h = 1e-8;
      const Tp = T2(c, [Xa[0] + h * q[0], Xa[1] + h * q[1]]);
      const Tm = T2(c, [Xa[0] - h * q[0], Xa[1] - h * q[1]]);
      return { d: [(Tp[0] - Tm[0]) / (2 * h), (Tp[1] - Tm[1]) / (2 * h)], q };
    }
    const Tb = T2(c, [Xa[0] + eps * q[0], Xa[1] + eps * q[1]]);
    return { d: [(Tb[0] - Ta[0]) / eps, (Tb[1] - Ta[1]) / eps], q };
  };
  const cross = (eta: number) => {
    const { d, q } = D(eta);
    return d[0] * q[1] - d[1] * q[0];
  };
  let eta = 0;
  for (let it = 0; it < 30; it++) {
    const f0 = cross(eta);
    const h = 1e-7;
    const df = (cross(eta + h) - cross(eta - h)) / (2 * h);
    if (df === 0) break;
    const step = f0 / df;
    eta -= step;
    if (Math.abs(step) < 1e-15) break;
  }
  const { d, q } = D(eta);
  const G2 = (d[0] * q[0] + d[1] * q[1]) / (q[0] * q[0] + q[1] * q[1]);
  return { Xb: [Xa[0] + eps * q[0], Xa[1] + eps * q[1]], G2 };
}

interface WaveOut {
  /** Фронты волны от переднего (ближе к «впереди») к заднему, с состояниями позади каждого. */
  pieces: { G: number; behind: State; kind: 'shock' | 'simple' }[];
  behind: State;
}

/**
 * Волна семейства f, бегущая в направлении dir (+1 вправо, −1 влево), из состояния впереди a.
 * ε ≥ 0 (или линейно вырожденное семейство) — одна ударная волна; ε < 0 — простая волна,
 * разбитая на n слабых фронтов.
 */
function wave(c: Ctx, a: State, f: Fam, dir: 1 | -1, eps: number, nFan: number): WaveOut {
  const R = orientedR(c, a.F, f);
  // число фронтов пучка задаётся извне и не меняется внутри итераций Ньютона (иначе невязка разрывна)
  const n = eps >= 0 || !R.gnl ? 1 : nFan;
  const pieces: WaveOut['pieces'] = [];
  let cur = a;
  for (let j = 0; j < n; j++) {
    const Rj = j === 0 ? R : orientedR(c, cur.F, f);
    const { Xb, G2 } = hugoniot(c, cur.F, f, eps / n, Rj);
    const G = Math.sqrt(Math.max(G2, 0));
    // [v] = −G_signed [F] с [·] = справа − слева  ⇒  v_b = v_a + dir·G·(F_a − F_b)
    const v: V2 = [cur.v[0] + dir * G * (cur.F[0] - Xb[0]), cur.v[1] + dir * G * (cur.F[1] - Xb[1])];
    const behind = { F: Xb, v };
    pieces.push({ G: dir * G, behind, kind: n === 1 && eps >= 0 ? 'shock' : n === 1 ? 'shock' : 'simple' });
    cur = behind;
  }
  return { pieces, behind: cur };
}

/** Число фронтов пучка для простой волны амплитуды eps. */
const fanCount = (eps: number, fanStep: number) => Math.max(1, Math.ceil(Math.abs(eps) / fanStep));

/**
 * Решить систему с фиксированными числами фронтов пучков; если после схождения амплитуды
 * простых волн сильно отличаются от исходных, числа фронтов пересчитываются и решение повторяется.
 */
function solveWithFans(F: (e: number[], n: number[]) => number[], g: number[], fanStep: number): { e: number[]; n: number[] } | null {
  let n = g.map((x) => fanCount(x, fanStep));
  let x0 = g;
  for (let pass = 0; pass < 3; pass++) {
    const e = newtonSolve((x) => F(x, n), x0, 1e-5);
    if (!e) return null;
    const n2 = e.map((x) => fanCount(x, fanStep));
    if (n2.every((v, i) => v === n[i] || e[i] >= 0)) return { e, n };
    n = n2;
    x0 = e;
  }
  return { e: x0, n };
}

/**
 * Порог приёма решения: невязка в деформациях/скоростях. Ниже ~1e-12 точность ограничена
 * численным дифференцированием внутри ударной адиабаты; 1e-11 — это ~1e-8 от характерной деформации.
 */
const ACCEPT = 1e-11;

function newtonSolve(F: (x: number[]) => number[], x0: number[], scale: number): number[] | null {
  let x = [...x0];
  let fx = F(x);
  let r = Math.hypot(...fx);
  for (let it = 0; it < 40; it++) {
    if (r < 1e-17) return x;
    const J: number[][] = fx.map(() => new Array(x.length).fill(0));
    for (let j = 0; j < x.length; j++) {
      const h = 1e-7 * Math.max(Math.abs(x[j]), scale);
      const xp = [...x];
      const xm = [...x];
      xp[j] += h;
      xm[j] -= h;
      const fp = F(xp);
      const fm = F(xm);
      for (let i = 0; i < fx.length; i++) J[i][j] = (fp[i] - fm[i]) / (2 * h);
    }
    const dx = solveLinear(J, fx.map((v) => -v));
    if (!dx) return null;
    let t = 1;
    let ok = false;
    for (let k = 0; k < 20; k++) {
      const xn = x.map((v, i) => v + t * dx[i]);
      const fn = F(xn);
      const rn = Math.hypot(...fn);
      if (Number.isFinite(rn) && rn < r) {
        x = xn;
        fx = fn;
        r = rn;
        ok = true;
        break;
      }
      t /= 2;
    }
    if (!ok) return r < ACCEPT ? x : null;
  }
  return r < ACCEPT ? x : null;
}

/** Линейное приближение амплитуд по собственным векторам в состоянии X. */
function basis(c: Ctx, X: V2) {
  const s = orientedR(c, X, 'slow');
  const f = orientedR(c, X, 'fast');
  return { rs: s.r, rf: f.r, cs: Math.sqrt(s.lam), cf: Math.sqrt(f.lam) };
}

/**
 * Краевая задача §2.3: у границы задана скорость (d, 0); впереди состояние a (может двигаться).
 * Возвращает волны, бегущие вглубь (быстрая, затем медленная).
 */
export function boundaryProblem(c: Ctx, a: State, d: number, fanStep: number): { fast: WaveOut; slow: WaveOut } | null {
  const run = (e: number[], n: number[]) => {
    const w1 = wave(c, a, 'fast', 1, e[0], n[0]);
    const w2 = wave(c, w1.behind, 'slow', 1, e[1], n[1]);
    return { w1, w2 };
  };
  // линеаризация: Δv = −cf α rf − cs β rs = (d − v1, −v3)
  const B = basis(c, a.F);
  const need: V2 = [d - a.v[0], -a.v[1]];
  const A = [
    [-B.cf * B.rf[0], -B.cs * B.rs[0]],
    [-B.cf * B.rf[1], -B.cs * B.rs[1]],
  ];
  const g = solveLinear(A, need) ?? [0, 0];
  const sol = solveWithFans(
    (e, n) => {
      const { w2 } = run(e, n);
      return [w2.behind.v[0] - d, w2.behind.v[1]];
    },
    g,
    fanStep,
  );
  if (!sol) return null;
  const { w1, w2 } = run(sol.e, sol.n);
  return { fast: w1, slow: w2 };
}

/**
 * Распад разрыва между состояниями l и r (встреча фронтов): влево уходят быстрая и медленная
 * волны, вправо — медленная и быстрая; в середине общее состояние.
 */
export function riemannProblem(c: Ctx, l: State, r: State, fanStep: number) {
  const run = (e: number[], n: number[]) => {
    const LF = wave(c, l, 'fast', -1, e[0], n[0]);
    const LS = wave(c, LF.behind, 'slow', -1, e[1], n[1]);
    const RF = wave(c, r, 'fast', 1, e[3], n[3]);
    const RS = wave(c, RF.behind, 'slow', 1, e[2], n[2]);
    return { LF, LS, RS, RF };
  };
  // линеаризация в среднем состоянии
  const X: V2 = [(l.F[0] + r.F[0]) / 2, (l.F[1] + r.F[1]) / 2];
  const B = basis(c, X);
  // F: l + α1 rf + α2 rs = r + β1 rf + β2 rs;  v: l.v + cf α1 rf + cs α2 rs = r.v − cf β1 rf − cs β2 rs
  const M = [
    [B.rf[0], B.rs[0], -B.rs[0], -B.rf[0]],
    [B.rf[1], B.rs[1], -B.rs[1], -B.rf[1]],
    [B.cf * B.rf[0], B.cs * B.rs[0], B.cs * B.rs[0], B.cf * B.rf[0]],
    [B.cf * B.rf[1], B.cs * B.rs[1], B.cs * B.rs[1], B.cf * B.rf[1]],
  ];
  const rhs = [r.F[0] - l.F[0], r.F[1] - l.F[1], r.v[0] - l.v[0], r.v[1] - l.v[1]];
  const g = solveLinear(M, rhs) ?? [0, 0, 0, 0];
  const sol = solveWithFans(
    (e, n) => {
      const w = run(e, n);
      const a = w.LS.behind;
      const b = w.RS.behind;
      return [a.F[0] - b.F[0], a.F[1] - b.F[1], a.v[0] - b.v[0], a.v[1] - b.v[1]];
    },
    g,
    fanStep,
  );
  if (!sol) return null;
  return run(sol.e, sol.n);
}

// ---------------------------------------------------------------- отслеживание фронтов

interface Tracked extends Front {
  strength: number;
}

const strengthOf = (a: State, b: State) => Math.hypot(b.F[0] - a.F[0], b.F[1] - a.F[1]);

export function trackFronts(prob: FrontsProblem): FrontsResult {
  const c: Ctx = { m: prob.material };
  const fanStep = prob.fanStep ?? 2.5e-4;
  const xMax = prob.xMax ?? 1.05;
  // Фронты со скачком деформации меньше dropTol не создаются (поглощаются соседним состоянием).
  // Сдвиговые волны в этих задачах ~1e-6, продольные ~1e-3; отражения при встречах — 1e-10…1e-12.
  const dropTol = prob.dropTol ?? 1e-9;
  const F0: V2 = [prob.F0[0], prob.F0[2]];
  let fronts: Tracked[] = [];
  let states: State[] = [{ F: F0, v: [0, 0] }];
  const segments: Segment[] = [];
  const warnings: string[] = [];
  let t = 0;
  let ub = 0;
  let events = 0;
  let maxFronts = 0;
  let work = 0;
  let dissipation = 0;
  const W0 = energy1d(c.m, prob.F0);
  const frames: FrontsFrame[] = [];
  const nF = Math.max(2, prob.frames);

  const closeSeg = (f: Tracked, tEnd: number) => {
    segments.push({ x0: f.xBorn, t0: f.born, x1: f.xBorn + f.speed * (tEnd - f.born), t1: tEnd, fam: f.fam, kind: f.kind, strength: f.strength });
  };
  const advance = (dt: number) => {
    for (const f of fronts) f.x += f.speed * dt;
    const d = boundaryVelocity(prob.pulses, t + 1e-13);
    // работа границы: мощность −σ11(0)·v1 (σ — у границы, состояние states[0])
    const T0 = T2(c, states[0].F);
    work += -T0[0] * d * dt;
    ub += d * dt;
    // диссипация на ударных фронтах (2.45): D = [T·v] + G[W + ½v²] (впереди − позади, по направлению движения)
    for (let i = 0; i < fronts.length; i++) {
      const f = fronts[i];
      const L = states[i];
      const R = states[i + 1];
      const [ahead, behind] = f.speed > 0 ? [R, L] : [L, R];
      const Ta = T2(c, ahead.F);
      const Tb = T2(c, behind.F);
      const E = (s: State) => energy1d(c.m, [s.F[0], 0, s.F[1]]) + 0.5 * (s.v[0] ** 2 + s.v[1] ** 2);
      const flux = (Tt: V2, s: State) => Tt[0] * s.v[0] + Tt[1] * s.v[1];
      const sgn = f.speed > 0 ? 1 : -1;
      dissipation += dt * (sgn * (flux(Ta, ahead) - flux(Tb, behind)) + Math.abs(f.speed) * (E(ahead) - E(behind)));
    }
    t += dt;
  };
  /** Вставить волны (в порядке слева направо) в точку x между фронтами i−1 и i. */
  const insertWaves = (pos: number, x: number, left: State, waves: { G: number; behind: State; kind: 'shock' | 'simple'; fam: Fam }[], right: State) => {
    // waves упорядочены слева направо; состояния между ними
    const newF: Tracked[] = [];
    const newS: State[] = [left];
    let prev = left;
    for (const w of waves) {
      const next = w.behind; // состояние справа от фронта (см. сборку ниже)
      const s = strengthOf(prev, next);
      if (s < dropTol) {
        newS[newS.length - 1] = next; // слишком слабый фронт — поглощается
        prev = next;
        continue;
      }
      newF.push({ x, speed: w.G, fam: w.fam, kind: w.kind, born: t, xBorn: x, strength: s });
      newS.push(next);
      prev = next;
    }
    newS[newS.length - 1] = right; // последнее состояние совпадает с исходным справа (с точностью решения)
    fronts.splice(pos, 0, ...newF);
    states.splice(pos, 1, ...newS);
  };

  /** Волны краевой задачи (бегут вправо), упорядоченные слева направо, с состоянием справа от каждого фронта. */
  const boundaryWaves = (a: State, d: number) => {
    const res = boundaryProblem(c, a, d, fanStep);
    if (!res) return null;
    // от границы: состояние у границы = res.slow.behind; затем фронты медленной волны (от заднего к переднему),
    // затем быстрой. Для фронта k: справа от него — состояние «впереди» этого фронта.
    const seq: { G: number; behind: State; kind: 'shock' | 'simple'; fam: Fam }[] = [];
    const chain = (w: WaveOut, ahead0: State, fam: Fam) => {
      const ahead = [ahead0, ...w.pieces.map((p) => p.behind)];
      const out: typeof seq = [];
      for (let j = w.pieces.length - 1; j >= 0; j--) out.push({ G: w.pieces[j].G, behind: ahead[j], kind: w.pieces[j].kind, fam });
      return out;
    };
    seq.push(...chain(res.slow, res.fast.behind, 'slow'), ...chain(res.fast, a, 'fast'));
    return { left: res.slow.behind, seq };
  };

  const snapshot = () => frames.push({ t, xs: fronts.map((f) => f.x), states: states.map((s) => ({ F: [...s.F] as V2, v: [...s.v] as V2 })), ub });

  // начальная краевая задача
  const changes = [...new Set(prob.pulses.flatMap((p) => [p.start, p.start + p.duration]))].filter((x) => x >= 0 && x < prob.tEnd).sort((a, b) => a - b);
  let nextFrame = 0;
  const appliedBc = new Set<number>();
  const applyBoundary = () => {
    const d = boundaryVelocity(prob.pulses, t + 1e-13);
    const bw = boundaryWaves(states[0], d);
    if (!bw) {
      warnings.push(`t = ${t.toFixed(4)}: краевая задача не решена`);
      return;
    }
    insertWaves(0, 0, bw.left, bw.seq, states[0]);
    events++;
  };

  while (t < prob.tEnd - 1e-13 && events < 20000) {
    // ближайшее событие
    let tNext = prob.tEnd;
    let kind: 'end' | 'bc' | 'collide' | 'hit' | 'frame' = 'end';
    let which = -1;
    const tFrame = (nextFrame * prob.tEnd) / (nF - 1);
    if (tFrame < tNext - 1e-13 && tFrame >= t - 1e-13) {
      tNext = tFrame;
      kind = 'frame';
    }
    for (const b of changes)
      if (b >= t - 1e-13 && b < tNext - 1e-13 && !(b <= t + 1e-13 && appliedBc.has(b))) {
        tNext = b;
        kind = 'bc';
      }
    for (let i = 0; i + 1 < fronts.length; i++) {
      const a = fronts[i];
      const b = fronts[i + 1];
      if (a.speed > b.speed + 1e-12) {
        const tc = t + Math.max(0, b.x - a.x) / (a.speed - b.speed);
        if (tc < tNext - 1e-14) {
          tNext = tc;
          kind = 'collide';
          which = i;
        }
      }
    }
    if (fronts.length && fronts[0].speed < 0) {
      const th = t + fronts[0].x / -fronts[0].speed;
      if (th < tNext - 1e-14) {
        tNext = th;
        kind = 'hit';
      }
    }
    advance(Math.max(0, tNext - t));
    if (kind === 'frame' || Math.abs(t - tFrame) < 1e-13) {
      if (nextFrame < nF) {
        snapshot();
        nextFrame++;
      }
      if (kind === 'frame') continue;
    }
    if (kind === 'bc') {
      appliedBc.add(tNext);
      applyBoundary();
    } else if (kind === 'collide') {
      const a = fronts[which];
      const b = fronts[which + 1];
      const x = 0.5 * (a.x + b.x);
      closeSeg(a, t);
      closeSeg(b, t);
      const l = states[which];
      const r = states[which + 2];
      fronts.splice(which, 2);
      states.splice(which, 3, l); // между l и r будет вставлено решение задачи о распаде разрыва
      const rp = riemannProblem(c, l, r, fanStep);
      if (!rp) {
        warnings.push(`t = ${t.toFixed(4)}: распад разрыва не решён — расчёт остановлен`);
        states.splice(which + 1, 0, r);
        fronts.splice(which, 0, { x, speed: 0, fam: 'fast', kind: 'shock', born: t, xBorn: x, strength: strengthOf(l, r) });
        break;
      }
      // слева направо: быстрая влево (от переднего к заднему фронту), медленная влево, медленная вправо (от заднего к переднему), быстрая вправо
      const seq: { G: number; behind: State; kind: 'shock' | 'simple'; fam: Fam }[] = [];
      const leftChain = (w: WaveOut, fam: Fam) => w.pieces.map((p) => ({ G: p.G, behind: p.behind, kind: p.kind, fam }));
      seq.push(...leftChain(rp.LF, 'fast'), ...leftChain(rp.LS, 'slow'));
      const rightChain = (w: WaveOut, ahead0: State, fam: Fam) => {
        const ahead = [ahead0, ...w.pieces.map((p) => p.behind)];
        const out: typeof seq = [];
        for (let j = w.pieces.length - 1; j >= 0; j--) out.push({ G: w.pieces[j].G, behind: ahead[j], kind: w.pieces[j].kind, fam });
        return out;
      };
      seq.push(...rightChain(rp.RS, rp.RF.behind, 'slow'), ...rightChain(rp.RF, r, 'fast'));
      insertWaves(which, x, l, seq, r);
      events++;
    } else if (kind === 'hit') {
      const f = fronts[0];
      closeSeg(f, t);
      fronts.shift();
      states.shift();
      applyBoundary();
    }
    // удалить ушедшие за пределы фронты
    while (fronts.length && fronts.at(-1)!.x > xMax) {
      closeSeg(fronts.at(-1)!, t);
      fronts.pop();
      states.pop();
    }
    maxFronts = Math.max(maxFronts, fronts.length);
  }
  if (events >= 20000) warnings.push('Достигнут предел числа событий (20000).');
  for (const f of fronts) closeSeg(f, t);
  while (frames.length < nF) snapshot();

  // кинематический баланс: ∫(p − p0)dx + ∫v*dt (по последнему кадру, интеграл до xMax)
  const last = frames.at(-1)!;
  let kin = 0;
  for (let i = 0; i < last.states.length; i++) {
    const xa = i === 0 ? 0 : last.xs[i - 1];
    const xb = i < last.xs.length ? last.xs[i] : xMax;
    kin += (last.states[i].F[0] - F0[0]) * (xb - xa);
  }
  let en = 0;
  for (let i = 0; i < last.states.length; i++) {
    const xa = i === 0 ? 0 : last.xs[i - 1];
    const xb = i < last.xs.length ? last.xs[i] : xMax;
    const s = last.states[i];
    en += (energy1d(c.m, [s.F[0], 0, s.F[1]]) - W0 + 0.5 * (s.v[0] ** 2 + s.v[1] ** 2)) * (xb - xa);
  }
  return { frames, segments, events, maxFronts, warnings, kinematicError: kin + last.ub, work, energyGain: en, dissipation };
}

/** Состояние в точке x кадра. */
export function stateAt(fr: FrontsFrame, x: number): State {
  let lo = 0;
  let hi = fr.xs.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (fr.xs[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return fr.states[lo];
}
