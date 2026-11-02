/**
 * Гибридный метод для серии нормальных ударов: ударные волны выделяются точно, простые волны
 * (вееры) считаются на сетке как гладкое решение.
 *
 * Зачем. В методе волновых фронтов (fronts.ts) веер — пучок слабых разрывов; каждая встреча двух
 * лучей решается как задача о распаде разрыва и порождает новые слабые фронты. Решение получается
 * ступенчатым, число событий растёт как произведение числа лучей, а точность — первого порядка по
 * шагу веера. Здесь ступенек нет: гладкие области — это сетка второго порядка.
 *
 * Как устроено (лагранжева координата x̂ = x1/L, система w_t + f(w)_x = 0, w = (p, s, v1, v3),
 * f = −(v1, v3, T1(F), T3(F))):
 *  • Ударные волны — подвижные внутренние границы. На каждом шаге состояние по обе стороны
 *    находится из условий на разрыве (2.40)–(2.42) [f] = G[w] и соотношений на характеристиках,
 *    приходящих к фронту из соседних областей (по Лаксу их 5, неизвестных 9: w₋, w₊, G). Свободны
 *    лишь амплитуды трёх уходящих от фронта характеристических волн — так на фронте рождаются
 *    отражённые и прошедшие гладкие волны (например, при пересечении ударной волны с веером).
 *  • Между фронтами — конечные объёмы на равномерной сетке; объёмы, примыкающие к фронту, имеют
 *    подвижную границу (длина от h/2 до 3h/2), поток через неё f(w±) − G w± (по условиям на разрыве
 *    он одинаков с обеих сторон, поэтому схема консервативна). Реконструкция MUSCL (minmod),
 *    предиктор Хэнкока, поток по характеристикам |∂f/∂w| = diag(√A, √A) (A — акустическая матрица):
 *    каждое семейство сносится со своей скоростью, медленные волны не размываются быстрыми.
 *  • Встреча двух ударных волн, удар волны о границу и изменение скорости границы решаются точно
 *    (задача о распаде разрыва и краевая задача §2.3 из fronts.ts); ударные волны решения становятся
 *    новыми выделенными фронтами, простые волны уходят на сетку. Узкая область между только что
 *    родившимися (или сходящимися) фронтами, пока она уже h/2, хранится двумя крайними состояниями.
 *  • Ударная волна, ослабевшая ниже порога или переставшая быть допустимой (её «съел» веер),
 *    перестаёт выделяться — дальше её ведёт сетка.
 */
import type { Material } from './material.ts';
import { type Vec3, energy1d } from './material1d.ts';
import { type Pulse, boundaryVelocity, tractionInto } from './fv1d.ts';
import { type Fam, type Segment, type State, boundaryProblem, riemannProblem, trackFronts } from './fronts.ts';
import { solveLinear } from './newton.ts';

export interface HybridProblem {
  material: Material;
  F0: Vec3; // (p0, 0, s0)
  pulses: Pulse[];
  tEnd: number;
  frames: number;
  /** Ячеек на единицу длины L (по умолчанию 400). */
  cellsPerL?: number;
  /** Длина расчётного отрезка (по умолчанию 1). */
  length?: number;
  /** Предел общего числа ячеек (на длинном отрезке сетка огрубляется; время растёт как квадрат числа ячеек). */
  maxCells?: number;
  cfl?: number;
  /** Шаг разбиения вееров внутри точных задач о распаде разрыва (только чтобы распознать составные волны). */
  fanStep?: number;
  /** Ударные волны слабее (скачок деформации) не выделяются. */
  trackTol?: number;
  /** Считать «до затухания»: длина отрезка и конечное время — как у метода фронтов в этом режиме. */
  untilQuiet?: boolean;
  maxTime?: number;
  maxSteps?: number;
}

/** Кадр: кусочно-линейное решение по узлам; на выделенном фронте два узла с одинаковым x. */
export interface HybridFrame {
  t: number;
  ub: number;
  xs: Float64Array;
  /** (p − p0, s − s0, v1, v3) в узлах: [4·k + j]. */
  W: Float32Array;
  /** Выделенные фронты в этот момент. */
  fronts: { x: number; G: number; fam: Fam; kind: 'shock' | 'simple' }[];
}

export interface HybridResult {
  frames: HybridFrame[];
  segments: Segment[];
  steps: number;
  /** Точные задачи (удары о границу, встречи ударных волн, изменения скорости границы). */
  events: number;
  maxFronts: number;
  /** Сколько ударных волн перестали выделяться (ослабли или стали недопустимыми). */
  untracked: number;
  warnings: string[];
  kinematicError: number;
  work: number;
  energyGain: number;
  /** Диссипация на выделенных ударных волнах (2.45); остаток work − energyGain − она — численная. */
  shockDissipation: number;
  tEnd: number;
  xView: number;
  cellsPerL: number;
  N: number;
  tQuiet?: number;
  horizonReached?: boolean;
}

type V4 = [number, number, number, number];

interface HFront {
  x: number;
  G: number;
  fam: Fam;
  /** Номер характеристического семейства по возрастанию скорости: −c_f, −c_s, +c_s, +c_f. */
  m: number;
  wm: V4;
  wp: V4;
  /** 'shock' — ударная волна; 'simple' — тонкий веер, который до конца расчёта уже ячейки сетки. */
  kind: 'shock' | 'simple';
  /** Момент рождения (от него считается, насколько разошёлся бы тонкий веер). */
  born: number;
  segX: number;
  segT: number;
  segG: number;
  px: number;
  pt: number;
}

interface Region {
  passive: boolean;
  /** Узкая область: состояния у левого и правого концов (для соседних фронтов) и содержимое ∫w dx. */
  wl: V4;
  wr: V4;
  Q: V4;
  /** Активная область: края объёмов (xs[0], xs[n] — стенки) и средние (4 на объём). */
  n: number;
  xs: Float64Array;
  U: Float64Array;
}

const mOf = (fam: Fam, dir: number) => (dir < 0 ? (fam === 'fast' ? 0 : 1) : fam === 'slow' ? 2 : 3);
const toState = (w: V4): State => ({ F: [w[0], w[1]], v: [w[2], w[3]] });
const fromState = (s: State): V4 => [s.F[0], s.F[1], s.v[0], s.v[1]];
const minmod = (a: number, b: number) => (a * b <= 0 ? 0 : Math.abs(a) < Math.abs(b) ? a : b);

export function solveHybrid(prob: HybridProblem, onProgress?: (fraction: number) => void): HybridResult {
  if (!prob.untilQuiet) return run(prob, prob.length ?? 1, prob.tEnd, onProgress);
  // те же отрезок и время, что у метода фронтов «до затухания» (грубый шаг веера — нужны лишь размеры)
  const fr = trackFronts({ material: prob.material, F0: prob.F0, pulses: prob.pulses, tEnd: prob.tEnd, frames: 2, untilQuiet: true, maxTime: prob.maxTime, fanStep: Math.max(prob.fanStep ?? 0, 1e-3) });
  onProgress?.(0.05);
  const res = run(prob, fr.xView, fr.tEnd, (f) => onProgress?.(0.05 + 0.95 * f));
  res.tQuiet = fr.tQuiet;
  res.horizonReached = fr.horizonReached;
  return res;
}

function run(prob: HybridProblem, L: number, tEnd: number, onProgress?: (fraction: number) => void): HybridResult {
  const mat = prob.material;
  const c = { m: mat };
  const N = Math.max(8, Math.min(Math.round((prob.cellsPerL ?? 400) * L), prob.maxCells ?? Infinity));
  const h = L / N;
  const cfl = prob.cfl ?? 0.4;
  const fanStep = prob.fanStep ?? 2.5e-4;
  const trackTol = prob.trackTol ?? 1e-9;
  const p0 = prob.F0[0];
  const s0 = prob.F0[2];
  const U0: V4 = [p0, s0, 0, 0];
  const warnings: string[] = [];

  // ------------------------------------------------------------ материал
  const tb = new Float64Array(3);
  /** T = (σ11, σ31) в out[o], out[o+1]. */
  const Tq = (p: number, s: number, out: Float64Array | number[], o: number) => {
    tractionInto(mat, p, 0, s, tb, 0);
    out[o] = tb[0];
    out[o + 1] = tb[2];
  };
  const Ta = new Float64Array(4);
  /** Акустическая матрица ∂T/∂F (симметричная): [a11, a12, a22]. */
  const acoustic = (p: number, s: number, out: Float64Array | number[], o: number) => {
    const hp = 1e-7 * Math.max(Math.abs(p), 1e-4);
    const hs = 1e-7 * Math.max(Math.abs(s), Math.abs(p), 1e-4);
    Tq(p + hp, s, Ta, 0);
    Tq(p - hp, s, Ta, 2);
    const a11 = (Ta[0] - Ta[2]) / (2 * hp);
    const a21 = (Ta[1] - Ta[3]) / (2 * hp);
    Tq(p, s + hs, Ta, 0);
    Tq(p, s - hs, Ta, 2);
    const a12 = (Ta[0] - Ta[2]) / (2 * hs);
    const a22 = (Ta[1] - Ta[3]) / (2 * hs);
    out[o] = a11;
    out[o + 1] = 0.5 * (a12 + a21);
    out[o + 2] = a22;
  };
  const Am = [0, 0, 0];
  /** Собственные пары акустической матрицы: скорости c_s ≤ c_f и векторы. */
  const eig = (w: V4 | Float64Array, o = 0) => {
    acoustic(w[o], w[o + 1], Am, 0);
    const [a, b, d] = Am;
    const tr = (a + d) / 2;
    const disc = Math.sqrt(((a - d) / 2) ** 2 + b * b);
    const vec = (l: number): [number, number] => {
      const v: [number, number] = Math.abs(b) > 1e-300 ? [b, l - a] : Math.abs(l - a) < Math.abs(l - d) ? [1, 0] : [0, 1];
      const n = Math.hypot(v[0], v[1]);
      return [v[0] / n, v[1] / n];
    };
    const ls = Math.max(tr - disc, 0);
    const lf = Math.max(tr + disc, 0);
    return { cs: Math.sqrt(ls), cf: Math.sqrt(lf), rs: vec(tr - disc), rf: vec(tr + disc), A: [a, b, d] as [number, number, number] };
  };
  /** Характеристическая скорость семейства m. */
  const sigma = (e: ReturnType<typeof eig>, m: number) => [-e.cf, -e.cs, e.cs, e.cf][m];
  /** Направление в пространстве w простой волны семейства j: (r, −σ r). */
  const waveDir = (e: ReturnType<typeof eig>, j: number): V4 => {
    const r = j === 0 || j === 3 ? e.rf : e.rs;
    const sg = sigma(e, j);
    return [r[0], r[1], -sg * r[0], -sg * r[1]];
  };

  // численный поток по характеристикам между состояниями L и R
  const fT = new Float64Array(4);
  /** a, b, d — акустическая матрица у грани (среднее по соседним объёмам). */
  const numFlux = (L: Float64Array | number[], lo: number, R: Float64Array | number[], ro: number, out: Float64Array, o: number, a: number, b: number, d: number) => {
    const det = Math.max(a * d - b * b, 0);
    const sq = Math.sqrt(det);
    const tt = Math.sqrt(Math.max(a + d + 2 * sq, 1e-300));
    const c11 = (a + sq) / tt;
    const c12 = b / tt;
    const c22 = (d + sq) / tt;
    Tq(L[lo], L[lo + 1], fT, 0);
    Tq(R[ro], R[ro + 1], fT, 2);
    const d0 = R[ro] - L[lo];
    const d1 = R[ro + 1] - L[lo + 1];
    const d2 = R[ro + 2] - L[lo + 2];
    const d3 = R[ro + 3] - L[lo + 3];
    out[o] = -0.5 * (L[lo + 2] + R[ro + 2]) - 0.5 * (c11 * d0 + c12 * d1);
    out[o + 1] = -0.5 * (L[lo + 3] + R[ro + 3]) - 0.5 * (c12 * d0 + c22 * d1);
    out[o + 2] = -0.5 * (fT[0] + fT[2]) - 0.5 * (c11 * d2 + c12 * d3);
    out[o + 3] = -0.5 * (fT[1] + fT[3]) - 0.5 * (c12 * d2 + c22 * d3);
  };
  const exactFlux = (w: V4): V4 => {
    Tq(w[0], w[1], fT, 0);
    return [-w[2], -w[3], -fT[0], -fT[1]];
  };
  /** Состояние у границы x = 0 со скоростью (d, 0): по приходящим к границе характеристикам. */
  const boundaryState = (w: V4, d: number): V4 => {
    const e = eig(w);
    const [a, b, dd] = e.A;
    const det = Math.max(a * dd - b * b, 1e-300);
    const sq = Math.sqrt(det);
    const tt = Math.sqrt(a + dd + 2 * sq);
    const c11 = (a + sq) / tt;
    const c12 = b / tt;
    const c22 = (dd + sq) / tt;
    const cdet = c11 * c22 - c12 * c12;
    const dv0 = d - w[2];
    const dv1 = -w[3];
    // F_b = F + C⁻¹ (v_b − v)
    return [w[0] + (c22 * dv0 - c12 * dv1) / cdet, w[1] + (-c12 * dv0 + c11 * dv1) / cdet, d, 0];
  };

  /** Расхождение характеристик семейства m по обе стороны разрыва (> 0 — это веер). */
  const fanSpread = (wm: V4, wp: V4, m: number) => sigma(eig(wp), m) - sigma(eig(wm), m);
  /**
   * Веер, родившийся в момент born, ведётся одним разрывом, пока его ширина spread·τ (τ — время
   * после рождения) меньше того, на сколько сетка размазала бы его сама: max(h/2, ½·h^⅔·τ^⅓) — так
   * схема второго порядка размазывает почти линейную волну (скорость характеристик почти одна).
   * Сразу после рождения так ведётся любой веер: центрированный веер из точки сетка всё равно не
   * разрешит, а соотношения на соседнем фронте не порождают веер своего же семейства (составная
   * волна при невыпуклости). Если и к концу расчёта веер останется таким узким — он ведётся до конца.
   * При измельчении сетки вееры отдаются сетке раньше и решение сходится к точному.
   */
  const numWidth = (tau: number) => Math.max(0.5 * h, 0.5 * Math.cbrt(h * h * Math.max(tau, 0)));
  const thin = (spread: number, born: number) => spread <= 0 || spread * (t - born) <= numWidth(t - born) || spread * (tEnd - born) <= numWidth(tEnd - born);

  // ------------------------------------------------------------ фронты и области
  let fronts: HFront[] = [];
  const uniform = (a: number, b: number, w: V4): Region => {
    const r: Region = { passive: false, wl: [...w], wr: [...w], Q: [0, 0, 0, 0], n: 1, xs: Float64Array.of(a, b), U: Float64Array.from(w) };
    return r;
  };
  let regions: Region[] = [];
  const leftWall = (i: number) => (i === 0 ? 0 : fronts[i - 1].x);
  const rightWall = (i: number) => (i === fronts.length ? L : fronts[i].x);

  /** Наклоны (minmod по соседям внутри области; у стенок — односторонние). */
  const slopesOf = (r: Region, S: Float64Array) => {
    const n = r.n;
    for (let i = 0; i < n; i++) {
      const ci = 0.5 * (r.xs[i] + r.xs[i + 1]);
      for (let k = 0; k < 4; k++) {
        const u = r.U[4 * i + k];
        let sl = 0;
        if (n > 1) {
          const dl = i > 0 ? (u - r.U[4 * (i - 1) + k]) / (ci - 0.5 * (r.xs[i - 1] + r.xs[i])) : NaN;
          const dr = i < n - 1 ? (r.U[4 * (i + 1) + k] - u) / (0.5 * (r.xs[i + 1] + r.xs[i + 2]) - ci) : NaN;
          sl = Number.isNaN(dl) ? dr : Number.isNaN(dr) ? dl : minmod(dl, dr);
        }
        S[4 * i + k] = sl;
      }
    }
  };
  let scratchS = new Float64Array(4 * (N + 16));
  const ensure = (n: number) => {
    if (scratchS.length < 4 * (n + 4)) scratchS = new Float64Array(4 * (n + 4) * 2);
  };
  /** Значение реконструкции области в точке x (по ближайшему объёму). */
  const endValue = (r: Region, right: boolean): V4 => {
    if (r.passive) return right ? [...r.wr] : [...r.wl];
    ensure(r.n);
    slopesOf(r, scratchS);
    const i = right ? r.n - 1 : 0;
    const x = right ? r.xs[r.n] : r.xs[0];
    const ci = 0.5 * (r.xs[i] + r.xs[i + 1]);
    return [0, 1, 2, 3].map((k) => r.U[4 * i + k] + scratchS[4 * i + k] * (x - ci)) as V4;
  };

  /** Привести область к сетке: узкая — пассивная; иначе края объёмов — узлы сетки не ближе h/2 к стенкам. */
  const normalize = (i: number) => {
    const r = regions[i];
    const a = leftWall(i);
    const b = rightWall(i);
    const len = b - a;
    if (r.passive) {
      if (len >= 0.6 * h) {
        const w: V4 = [0, 1, 2, 3].map((k) => r.Q[k] / len) as V4;
        regions[i] = uniform(a, b, w);
        normalize(i);
      }
      return;
    }
    r.xs[0] = a;
    r.xs[r.n] = b;
    if (len < 0.5 * h) {
      const wl = endValue(r, false);
      const wr = endValue(r, true);
      // среднее сохраняется: (wl + wr)/2 = среднее по области для линейной реконструкции одного объёма
      const Q: V4 = [0, 0, 0, 0];
      for (let j = 0; j < r.n; j++) for (let k = 0; k < 4; k++) Q[k] += r.U[4 * j + k] * (r.xs[j + 1] - r.xs[j]);
      const tot = len > 0 ? (Q.map((v) => v / len) as V4) : Q;
      if (r.n > 1 && len > 0) {
        // линейный профиль с тем же средним
        const mid = wl.map((v, k) => 0.5 * (v + wr[k]));
        for (let k = 0; k < 4; k++) {
          wl[k] += tot[k] - mid[k];
          wr[k] += tot[k] - mid[k];
        }
      }
      regions[i] = { passive: true, wl, wr, Q, n: 0, xs: new Float64Array(0), U: new Float64Array(0) };
      return;
    }
    const j0 = Math.floor((a + 0.5 * h) / h) + 1;
    const j1 = Math.ceil((b - 0.5 * h) / h) - 1;
    const want = Math.max(0, j1 - j0 + 1);
    const n = want + 1;
    let same = n === r.n;
    if (same) for (let j = 1; j < n && same; j++) same = Math.abs(r.xs[j] - (j0 + j - 1) * h) < 1e-12 * h;
    if (same) return;
    // консервативная перекладка по линейной реконструкции
    ensure(r.n);
    slopesOf(r, scratchS);
    const xs = new Float64Array(n + 1);
    xs[0] = a;
    xs[n] = b;
    for (let j = 1; j < n; j++) xs[j] = (j0 + j - 1) * h;
    const U = new Float64Array(4 * n);
    let o = 0;
    for (let j = 0; j < n; j++) {
      const xa = xs[j];
      const xb = xs[j + 1];
      while (o < r.n - 1 && r.xs[o + 1] <= xa) o++;
      let q = o;
      while (q < r.n && r.xs[q] < xb) {
        const lo = Math.max(xa, r.xs[q]);
        const hi = Math.min(xb, r.xs[q + 1]);
        if (hi > lo) {
          const cq = 0.5 * (r.xs[q] + r.xs[q + 1]);
          const mid = 0.5 * (lo + hi);
          for (let k = 0; k < 4; k++) U[4 * j + k] += (hi - lo) * (r.U[4 * q + k] + scratchS[4 * q + k] * (mid - cq));
        }
        q++;
      }
      for (let k = 0; k < 4; k++) U[4 * j + k] /= xb - xa;
    }
    regions[i] = { passive: false, wl: r.wl, wr: r.wr, Q: [0, 0, 0, 0], n, xs, U };
  };

  /** Слить области i и i+1 (фронт i перестал выделяться или исчез). */
  const mergeAt = (i: number) => {
    const A = regions[i];
    const B = regions[i + 1];
    const x = fronts[i].x;
    const a = leftWall(i);
    const b = rightWall(i + 1);
    const vols: { xa: number; xb: number; w: V4 }[] = [];
    const extra: V4 = [0, 0, 0, 0]; // содержимое узких областей нулевой длины
    const add = (r: Region, ra: number, rb: number) => {
      if (r.passive) {
        if (rb > ra) vols.push({ xa: ra, xb: rb, w: r.Q.map((q) => q / (rb - ra)) as V4 });
        else for (let k = 0; k < 4; k++) extra[k] += r.Q[k];
      } else
        for (let j = 0; j < r.n; j++) {
          const xa = j === 0 ? ra : r.xs[j];
          const xb = j === r.n - 1 ? rb : r.xs[j + 1];
          if (xb > xa) vols.push({ xa, xb, w: [r.U[4 * j], r.U[4 * j + 1], r.U[4 * j + 2], r.U[4 * j + 3]] });
        }
    };
    add(A, a, x);
    add(B, x, b);
    closeSeg(fronts[i]);
    fronts.splice(i, 1);
    let R: Region;
    if (!vols.length) R = { passive: true, wl: A.passive ? A.wl : endValue(A, false), wr: B.passive ? B.wr : endValue(B, true), Q: extra, n: 0, xs: new Float64Array(0), U: new Float64Array(0) };
    else {
      const n = vols.length;
      const xs = new Float64Array(n + 1);
      const U = new Float64Array(4 * n);
      vols.forEach((v, j) => {
        xs[j] = v.xa;
        for (let k = 0; k < 4; k++) U[4 * j + k] = v.w[k];
      });
      xs[n] = vols[n - 1].xb;
      xs[0] = a;
      xs[n] = b;
      R = { passive: false, wl: [...U0], wr: [...U0], Q: [0, 0, 0, 0], n, xs, U };
      const l0 = xs[1] - xs[0];
      for (let k = 0; k < 4; k++) U[k] += extra[k] / l0;
    }
    regions.splice(i, 2, R);
    normalize(i);
  };

  // ------------------------------------------------------------ отрезки траекторий (для диаграммы)
  const segments: Segment[] = [];
  let t = 0;
  const segTol = 2e-4 * L;
  const strengthF = (f: HFront) => Math.hypot(f.wp[0] - f.wm[0], f.wp[1] - f.wm[1]);
  const closeSeg = (f: HFront) => {
    if (t > f.segT) segments.push({ x0: f.segX, t0: f.segT, x1: f.x, t1: t, fam: f.fam, kind: f.kind, strength: strengthF(f) });
  };
  const trackSeg = (f: HFront) => {
    if (Math.abs(f.x - (f.segX + f.segG * (t - f.segT))) > segTol && f.pt > f.segT) {
      segments.push({ x0: f.segX, t0: f.segT, x1: f.px, t1: f.pt, fam: f.fam, kind: f.kind, strength: strengthF(f) });
      f.segG = (f.x - f.px) / Math.max(t - f.pt, 1e-300);
      f.segX = f.px;
      f.segT = f.pt;
    }
    f.px = f.x;
    f.pt = t;
  };

  // ------------------------------------------------------------ точные задачи: новые фронты
  type Piece = { G: number; behind: State; kind: 'shock' | 'simple'; fam: Fam };
  /** Скорость разрыва L | R по условиям на разрыве (МНК): G² = [T]·[F]/|[F]|². */
  const rhSpeed = (Lw: V4, Rw: V4, sign: number) => {
    const TL = [0, 0];
    const TR = [0, 0];
    Tq(Lw[0], Lw[1], TL, 0);
    Tq(Rw[0], Rw[1], TR, 0);
    const dF = [Rw[0] - Lw[0], Rw[1] - Lw[1]];
    const G2 = ((TR[0] - TL[0]) * dF[0] + (TR[1] - TL[1]) * dF[1]) / (dF[0] ** 2 + dF[1] ** 2 || 1e-300);
    return sign * Math.sqrt(Math.max(G2, 0));
  };
  let merges = 0;
  /**
   * Последовательность волн (слева направо, у каждой — состояние справа) → выделяемые ударные волны
   * и пары крайних состояний промежутков между ними (gaps.length = новых фронтов + 1).
   */
  const buildFronts = (x: number, left: State, seq: Piece[], parents?: { fam: Fam; dir: number; s: number }[]) => {
    const ws: { G: number; fam: Fam; kind: 'shock' | 'simple'; L: V4; R: V4 }[] = [];
    let prev = fromState(left);
    for (const p of seq) {
      const next = fromState(p.behind);
      if (Math.hypot(next[0] - prev[0], next[1] - prev[1]) < 1e-14) {
        prev = next;
        continue;
      }
      ws.push({ G: p.G, fam: p.fam, kind: p.kind, L: prev, R: next });
      prev = next;
    }
    // составная волна (невыпуклость): соседние разрывы одного семейства в обратном порядке скоростей
    // сливаются в один ударный (как в fronts.ts)
    for (let k = 0; k + 1 < ws.length; ) {
      const A = ws[k];
      const B = ws[k + 1];
      if (A.fam === B.fam && Math.sign(A.G) === Math.sign(B.G) && A.G > B.G + 1e-13) {
        ws.splice(k, 2, { G: rhSpeed(A.L, B.R, Math.sign(A.G)), fam: A.fam, kind: 'shock', L: A.L, R: B.R });
        merges++;
        if (k > 0) k--;
      } else k++;
    }
    // лучи одного веера (подряд идущие простые волны одного семейства и направления) — одна волна
    const waves: typeof ws = [];
    for (let k = 0; k < ws.length; ) {
      const A = ws[k];
      let q = k + 1;
      if (A.kind === 'simple') while (q < ws.length && ws[q].kind === 'simple' && ws[q].fam === A.fam && Math.sign(ws[q].G) === Math.sign(A.G)) q++;
      const B = ws[q - 1];
      waves.push(q === k + 1 ? A : { G: rhSpeed(A.L, B.R, Math.sign(A.G)), fam: A.fam, kind: 'simple', L: A.L, R: B.R });
      k = q;
    }
    // какие волны выделять: ударные и тонкие вееры (сильнее порога)
    const str = (w: (typeof waves)[number]) => Math.hypot(w.R[0] - w.L[0], w.R[1] - w.L[1]);
    const maxParent = Math.max(0, ...(parents ?? []).map((q) => q.s));
    const track = waves.map((w, k) => {
      if (str(w) < trackTol) return false;
      // При встрече: продолжения встретившихся волн выделяются всегда, а новые (отражения и
      // переходы в другое семейство) — только если они не слабее 1 % встретившихся. Слабые волны
      // отлично считает сетка, а их выделение порождало бы лавину всё более слабых встреч.
      if (parents && !parents.some((q) => q.fam === w.fam && q.dir === Math.sign(w.G)) && str(w) < 0.01 * maxParent) return false;
      if (w.kind === 'shock') return true;
      // Веер ведётся одним разрывом, если до конца расчёта он останется уже размазывания сетки, или
      // если к нему примыкает ударная волна того же семейства (составная волна при невыпуклости:
      // соотношения на фронте не порождают веер своего же семейства) — тогда, пока он узкий.
      const adj = [waves[k - 1], waves[k + 1]].some((q) => q && q.kind === 'shock' && q.fam === w.fam && Math.sign(q.G) === Math.sign(w.G));
      const spread = fanSpread(w.L, w.R, mOf(w.fam, Math.sign(w.G)));
      return spread <= 0 || spread * (tEnd - t) <= numWidth(tEnd - t) || adj;
    });
    const nf: HFront[] = [];
    const gaps: [V4, V4][] = [];
    let gl = fromState(left);
    waves.forEach((w, k) => {
      if (!track[k]) return;
      const m = mOf(w.fam, Math.sign(w.G));
      gaps.push([gl, w.L]);
      nf.push({ x, G: w.G, fam: w.fam, m, wm: w.L, wp: w.R, kind: w.kind, born: t, segX: x, segT: t, segG: w.G, px: x, pt: t });
      gl = w.R;
    });
    gaps.push([gl, prev]);
    return { nf, gaps };
  };
  const passive = (wl: V4, wr: V4): Region => ({ passive: true, wl: [...wl], wr: [...wr], Q: [0, 0, 0, 0], n: 0, xs: new Float64Array(0), U: new Float64Array(0) });
  /** Содержимое исчезающей узкой области — соседней (чтобы балансы сохранялись точно). */
  const dumpQ = (Q: V4, i: number, atRight: boolean) => {
    const r = regions[i];
    if (!r) return;
    if (r.passive) for (let k = 0; k < 4; k++) r.Q[k] += Q[k];
    else {
      const j = atRight ? r.n - 1 : 0;
      const len = r.xs[j + 1] - r.xs[j];
      if (len > 0) for (let k = 0; k < 4; k++) r.U[4 * j + k] += Q[k] / len;
    }
  };

  let events = 0;
  let maxFronts = 0;
  let untracked = 0;
  /** Новые волны у границы x = 0 из состояния a (впереди) при скорости границы d. */
  const boundaryWaves = (a: V4, d: number) => {
    const res = boundaryProblem(c, toState(a), d, fanStep);
    if (!res) {
      warnings.push(`t = ${t.toFixed(4)}: краевая задача не решена`);
      return null;
    }
    const seq: Piece[] = [];
    const chain = (w: (typeof res)['fast'], ahead0: State, fam: Fam) => {
      const ahead = [ahead0, ...w.pieces.map((p) => p.behind)];
      for (let j = w.pieces.length - 1; j >= 0; j--) seq.push({ G: w.pieces[j].G, behind: ahead[j], kind: w.pieces[j].kind, fam });
    };
    chain(res.slow, res.fast.behind, 'slow');
    chain(res.fast, toState(a), 'fast');
    return buildFronts(0, res.slow.behind, seq);
  };
  /** Вставить рождённые у границы фронты (слева от всех). */
  const insertAtBoundary = (b: ReturnType<typeof buildFronts>) => {
    const k = b.nf.length;
    if (!k) {
      // только простые волны: у пассивной области обновляем левое состояние, активную ведёт граничное условие
      if (regions[0].passive) regions[0].wl = b.gaps[0][0];
      return;
    }
    fronts.splice(0, 0, ...b.nf);
    regions.splice(0, 0, ...b.gaps.slice(0, k).map(([wl, wr]) => passive(wl, wr)));
    if (regions[k].passive) regions[k].wl = b.gaps[k][0];
  };
  const applyBoundary = (dNew: number, dOld: number) => {
    const a = regions[0].passive ? boundaryState(regions[0].wl, dOld) : boundaryState(endValue(regions[0], false), dOld);
    const b = boundaryWaves(a, dNew);
    if (!b) return;
    insertAtBoundary(b);
    events++;
  };
  const collide = (i: number) => {
    const A = fronts[i];
    const B = fronts[i + 1];
    const x = 0.5 * (A.x + B.x);
    A.x = B.x = x;
    const l = A.wm;
    const r = B.wp;
    const rp = riemannProblem(c, toState(l), toState(r), fanStep);
    if (!rp) {
      // не решилась: обе волны перестают выделяться (их ведёт сетка)
      warnings.push(`t = ${t.toFixed(4)}: встреча ударных волн не решена точно — дальше их ведёт сетка`);
      mergeAt(i + 1);
      mergeAt(i);
      untracked += 2;
      return;
    }
    const seq: Piece[] = [];
    for (const p of rp.LF.pieces) seq.push({ ...p, fam: 'fast' });
    for (const p of rp.LS.pieces) seq.push({ ...p, fam: 'slow' });
    const rightChain = (w: (typeof rp)['RF'], ahead0: State, fam: Fam) => {
      const ahead = [ahead0, ...w.pieces.map((p) => p.behind)];
      for (let j = w.pieces.length - 1; j >= 0; j--) seq.push({ G: w.pieces[j].G, behind: ahead[j], kind: w.pieces[j].kind, fam });
    };
    rightChain(rp.RS, rp.RF.behind, 'slow');
    rightChain(rp.RF, toState(r), 'fast');
    const b = buildFronts(x, toState(l), seq, [A, B].map((f) => ({ fam: f.fam, dir: f.m >= 2 ? 1 : -1, s: strengthF(f) })));
    const k = b.nf.length;
    events++;
    if (!k) {
      // все волны — простые: области сливаются, разрыв l | r дальше ведёт сетка
      regions[i + 1].wl = l;
      regions[i + 1].wr = r;
      mergeAt(i);
      mergeAt(i);
      return;
    }
    closeSeg(A);
    closeSeg(B);
    dumpQ(regions[i + 1].Q, i, true);
    fronts.splice(i, 2, ...b.nf);
    regions.splice(i + 1, 1, ...b.gaps.slice(1, k).map(([wl, wr]) => passive(wl, wr)));
    if (regions[i].passive) regions[i].wr = b.gaps[0][1];
    if (regions[i + k].passive) regions[i + k].wl = b.gaps[k][0];
  };
  const hit = () => {
    const f = fronts[0];
    closeSeg(f);
    const a = f.wp;
    fronts.shift();
    const gone = regions.shift()!;
    dumpQ(gone.Q, 0, false);
    const b = boundaryWaves(a, dCur);
    events++;
    if (!b) return;
    insertAtBoundary(b);
  };

  // ------------------------------------------------------------ решение на выделенном фронте
  /**
   * По данным wL (слева) и wR (справа) найти G и состояния w₋, w₊: [f] = G[w], а слева/справа
   * меняются только амплитуды уходящих от фронта характеристических волн.
   */
  const solveFront = (f: HFront, wL: V4, wR: V4): boolean => {
    const eL = eig(wL);
    const eR = eig(wR);
    const dirs: { side: -1 | 1; e: V4 }[] = [];
    for (let j = 0; j < 4; j++) if (j < f.m) dirs.push({ side: -1, e: waveDir(eL, j) });
    for (let j = 0; j < 4; j++) if (j > f.m) dirs.push({ side: 1, e: waveDir(eR, j) });
    const z = [0, 0, 0, f.G];
    const states = () => {
      const wm = [...wL] as V4;
      const wp = [...wR] as V4;
      dirs.forEach((d, j) => {
        const w = d.side < 0 ? wm : wp;
        for (let k = 0; k < 4; k++) w[k] += z[j] * d.e[k];
      });
      return { wm, wp };
    };
    let ok = false;
    for (let it = 0; it < 30; it++) {
      const { wm, wp } = states();
      const fm = exactFlux(wm);
      const fp = exactFlux(wp);
      const G = z[3];
      const R = [0, 1, 2, 3].map((k) => fp[k] - fm[k] - G * (wp[k] - wm[k]));
      const scale = Math.hypot(fp[0] - fm[0], fp[1] - fm[1], fp[2] - fm[2], fp[3] - fm[3]);
      // предел точности — округление в напряжениях (~1e-16 от их величины)
      const floor = 1e-14 * (Math.abs(fm[2]) + Math.abs(fm[3]) + Math.abs(fp[2]) + Math.abs(fp[3]));
      if (Math.hypot(...R) <= 1e-12 * scale + floor) {
        ok = true;
        break;
      }
      const Am_ = eig(wm).A;
      const Ap = eig(wp).A;
      const J: number[][] = [[], [], [], []];
      dirs.forEach((d, j) => {
        const A = d.side < 0 ? Am_ : Ap;
        // Df(w)·e = −(e_v, A e_F)
        const Df = [-d.e[2], -d.e[3], -(A[0] * d.e[0] + A[1] * d.e[1]), -(A[1] * d.e[0] + A[2] * d.e[1])];
        for (let k = 0; k < 4; k++) J[k][j] = d.side < 0 ? -Df[k] + G * d.e[k] : Df[k] - G * d.e[k];
      });
      for (let k = 0; k < 4; k++) J[k][3] = -(wp[k] - wm[k]);
      const dz = solveLinear(J, R.map((v) => -v));
      if (!dz) return false;
      for (let k = 0; k < 4; k++) z[k] += dz[k];
    }
    if (!ok || !Number.isFinite(z[3])) return false;
    const { wm, wp } = states();
    // Ударная волна рождается допустимой (из точной задачи о распаде разрыва) и меняется непрерывно;
    // при невыпуклости (материал 1 около p = 0) простое условие Лакса σ(w₋) ≥ G ≥ σ(w₊) к ней
    // неприменимо, поэтому её перестаём вести, только если она ослабла до порога или решение
    // «перескочило» на другую ветвь (скорость резко изменилась). Тонкий веер ведём, пока он к концу
    // расчёта остаётся уже половины ячейки; иначе его ведёт сетка.
    const G = z[3];
    if (Math.sign(G) !== (f.m >= 2 ? 1 : -1) || Math.abs(G - f.G) > 0.02 * Math.abs(f.G)) return false;
    if (f.kind === 'simple' && !thin(fanSpread(wm, wp, f.m), f.born)) return false;
    if (Math.hypot(wp[0] - wm[0], wp[1] - wm[1]) < trackTol) return false;
    f.G = G;
    f.wm = wm;
    f.wp = wp;
    return true;
  };

  // ------------------------------------------------------------ шаг
  const nF = Math.max(2, prob.frames);
  const frames: HybridFrame[] = [];
  let ub = 0;
  let work = 0;
  let shockDiss = 0;
  let outKin = 0;
  let outEnergy = 0;
  const E = (w: V4) => energy1d(mat, [w[0], 0, w[1]]) + 0.5 * (w[2] ** 2 + w[3] ** 2);
  const W0 = energy1d(mat, prob.F0);

  /** Узлы кадра: центры объёмов и концы областей (реконструкция), у фронта — два узла. */
  const snapshot = () => {
    const xs: number[] = [];
    const W: number[] = [];
    const push = (x: number, w: ArrayLike<number>, o = 0) => {
      xs.push(x);
      W.push(w[o] - p0, w[o + 1] - s0, w[o + 2], w[o + 3]);
    };
    regions.forEach((r, i) => {
      const a = leftWall(i);
      const b = rightWall(i);
      if (r.passive) {
        push(a, r.wl);
        push(b, r.wr);
        return;
      }
      ensure(r.n);
      slopesOf(r, scratchS);
      const c0 = 0.5 * (r.xs[0] + r.xs[1]);
      push(a, [0, 1, 2, 3].map((k) => r.U[k] + scratchS[k] * (a - c0)));
      for (let j = 0; j < r.n; j++) push(0.5 * (r.xs[j] + r.xs[j + 1]), r.U, 4 * j);
      const cn = 0.5 * (r.xs[r.n - 1] + r.xs[r.n]);
      push(b, [0, 1, 2, 3].map((k) => r.U[4 * (r.n - 1) + k] + scratchS[4 * (r.n - 1) + k] * (b - cn)));
    });
    frames.push({ t, ub, xs: Float64Array.from(xs), W: Float32Array.from(W), fronts: fronts.map((f) => ({ x: f.x, G: f.G, fam: f.fam, kind: f.kind })) });
  };

  regions = [uniform(0, L, U0)];
  normalize(0);
  const changes = [...new Set(prob.pulses.flatMap((p) => [p.start, p.start + p.duration]))].filter((x) => x >= 0 && x < tEnd).sort((a, b) => a - b);
  const applied = new Set<number>();
  let nextFrame = 0;
  let steps = 0;
  const maxSteps = prob.maxSteps ?? 400000;
  let dCur = 0;

  // рабочие массивы шага
  const pred: Float64Array[] = []; // предсказанные на полшага средние (по областям)
  const slope: Float64Array[] = [];
  const nAct: number[] = [];
  const acm: Float64Array[] = []; // акустические матрицы объёмов [a11, a12, a22]

  const handleEventsAt = () => {
    // встречи (повторяем, пока есть: после встречи возможны новые совпадения)
    for (let guard = 0; guard < 1000; guard++) {
      let done = true;
      for (let i = 0; i + 1 < fronts.length; i++) {
        if (fronts[i + 1].x - fronts[i].x <= 1e-12 * Math.max(1, L) && fronts[i].G > fronts[i + 1].G + 1e-12) {
          collide(i);
          done = false;
          break;
        }
      }
      if (fronts.length && fronts[0].G < 0 && fronts[0].x <= 1e-12 * Math.max(1, L)) {
        hit();
        done = false;
      }
      while (fronts.length && fronts.at(-1)!.x >= L - 1e-12 * L) {
        const i = fronts.length - 1;
        closeSeg(fronts[i]);
        fronts.pop();
        dumpQ(regions.pop()!.Q, regions.length - 1, true);
        done = false;
      }
      if (done) break;
    }
    for (const b of changes)
      if (Math.abs(b - t) < 1e-12 && !applied.has(b)) {
        applied.add(b);
        const dNew = boundaryVelocity(prob.pulses, t + 1e-13);
        applyBoundary(dNew, dCur);
        dCur = dNew;
      }
    for (let i = 0; i < regions.length; i++) normalize(i);
    maxFronts = Math.max(maxFronts, fronts.length);
  };

  // начальный кадр и воздействия в t = 0
  const frameTime = (k: number) => (k * tEnd) / (nF - 1);
  snapshot();
  nextFrame = 1;
  handleEventsAt();

  while (t < tEnd - 1e-13 && steps < maxSteps) {
    // наклоны, скорость звука, активная часть каждой области
    let cmax = 1e-6;
    for (let i = 0; i < regions.length; i++) {
      const r = regions[i];
      if (r.passive) {
        nAct[i] = 0;
        continue;
      }
      if (!slope[i] || slope[i].length < 4 * r.n) slope[i] = new Float64Array(4 * r.n + 64);
      if (!pred[i] || pred[i].length < 4 * r.n) pred[i] = new Float64Array(4 * r.n + 64);
      if (!acm[i] || acm[i].length < 3 * r.n) acm[i] = new Float64Array(3 * r.n + 48);
      // правее последнего возмущённого объёма всё равно исходному состоянию — не считаем
      let last = -1;
      for (let j = r.n - 1; j >= 0; j--) {
        const o = 4 * j;
        if (r.U[o] !== U0[0] || r.U[o + 1] !== U0[1] || r.U[o + 2] !== 0 || r.U[o + 3] !== 0) {
          last = j;
          break;
        }
      }
      nAct[i] = i === regions.length - 1 ? Math.min(r.n, last + 4) : r.n;
      slopesOf(r, slope[i]);
      const Ac = acm[i];
      for (let j = 0; j < Math.min(r.n, nAct[i] + 1); j++) {
        acoustic(r.U[4 * j], r.U[4 * j + 1], Ac, 3 * j);
        const a = Ac[3 * j];
        const b = Ac[3 * j + 1];
        const d = Ac[3 * j + 2];
        cmax = Math.max(cmax, Math.sqrt(Math.max(0, 0.5 * (a + d) + Math.sqrt(0.25 * (a - d) ** 2 + b * b))));
      }
    }
    for (const f of fronts) cmax = Math.max(cmax, Math.abs(f.G));
    let dt = (cfl * h) / cmax;
    const tFrame = frameTime(nextFrame);
    for (const b of [...changes.filter((x) => !applied.has(x)), tFrame, tEnd]) if (b > t + 1e-13 && t + dt > b) dt = b - t;

    // предиктор Хэнкока и решение на фронтах; если встреча раньше — шаг укорачивается
    let tries = 0;
    let failed = -1;
    for (;;) {
      for (let i = 0; i < regions.length; i++) {
        const r = regions[i];
        if (r.passive) continue;
        const S = slope[i];
        const P = pred[i];
        const Ac = acm[i];
        for (let j = 0; j < Math.min(r.n, nAct[i] + 1); j++) {
          const o = 4 * j;
          const q = 3 * j;
          P[o] = r.U[o] + 0.5 * dt * S[o + 2];
          P[o + 1] = r.U[o + 1] + 0.5 * dt * S[o + 3];
          P[o + 2] = r.U[o + 2] + 0.5 * dt * (Ac[q] * S[o] + Ac[q + 1] * S[o + 1]);
          P[o + 3] = r.U[o + 3] + 0.5 * dt * (Ac[q + 1] * S[o] + Ac[q + 2] * S[o + 1]);
        }
      }
      const faceVal = (i: number, right: boolean, x: number): V4 => {
        const r = regions[i];
        if (r.passive) return right ? [...r.wr] : [...r.wl];
        const j = right ? r.n - 1 : 0;
        const cj = 0.5 * (r.xs[j] + r.xs[j + 1]);
        return [0, 1, 2, 3].map((k) => pred[i][4 * j + k] + slope[i][4 * j + k] * (x - cj)) as V4;
      };
      failed = -1;
      for (let i = 0; i < fronts.length; i++) {
        const f = fronts[i];
        const xh = f.x + 0.5 * dt * f.G;
        if (!solveFront(f, faceVal(i, true, xh), faceVal(i + 1, false, xh))) {
          failed = i;
          break;
        }
      }
      if (failed >= 0) break;
      // ближайшая встреча/удар/выход по новым скоростям
      let tc = Infinity;
      for (let i = 0; i + 1 < fronts.length; i++) {
        const a = fronts[i];
        const b = fronts[i + 1];
        if (a.G > b.G + 1e-12) tc = Math.min(tc, Math.max(0, b.x - a.x) / (a.G - b.G));
      }
      if (fronts.length && fronts[0].G < 0) tc = Math.min(tc, fronts[0].x / -fronts[0].G);
      if (fronts.length && fronts.at(-1)!.G > 0) tc = Math.min(tc, (L - fronts.at(-1)!.x) / fronts.at(-1)!.G);
      if (tc < dt - 1e-15 && tries < 2) {
        dt = Math.max(tc, 0);
        tries++;
        continue;
      }
      if (tc < dt) dt = Math.max(tc, 0);
      break;
    }
    if (failed >= 0) {
      // ударная волна ослабла (или её «съел» веер) — дальше её ведёт сетка
      mergeAt(failed);
      untracked++;
      continue;
    }

    // потоки и обновление
    const d = dCur;
    for (let i = 0; i < regions.length; i++) {
      const r = regions[i];
      if (r.passive) continue;
      const n = r.n;
      const na = Math.max(nAct[i], 1);
      const P = pred[i];
      const S = slope[i];
      const flux = new Float64Array(4 * (n + 1));
      const fl = new Float64Array(4);
      const fr = new Float64Array(4);
      // внутренние грани
      for (let j = 0; j < Math.min(na, n - 1); j++) {
        const cj = 0.5 * (r.xs[j] + r.xs[j + 1]);
        const cj1 = 0.5 * (r.xs[j + 1] + r.xs[j + 2]);
        const xf = r.xs[j + 1];
        for (let k = 0; k < 4; k++) {
          fl[k] = P[4 * j + k] + S[4 * j + k] * (xf - cj);
          fr[k] = P[4 * (j + 1) + k] + S[4 * (j + 1) + k] * (xf - cj1);
        }
        const Ac = acm[i];
        numFlux(fl, 0, fr, 0, flux, 4 * (j + 1), 0.5 * (Ac[3 * j] + Ac[3 * j + 3]), 0.5 * (Ac[3 * j + 1] + Ac[3 * j + 4]), 0.5 * (Ac[3 * j + 2] + Ac[3 * j + 5]));
      }
      // левая стенка
      if (i === 0) {
        const c0 = 0.5 * (r.xs[0] + r.xs[1]);
        for (let k = 0; k < 4; k++) fr[k] = P[k] + S[k] * (0 - c0);
        const g: V4 = [fr[0], fr[1], 2 * d - fr[2], -fr[3]];
        numFlux(g, 0, fr, 0, flux, 0, acm[i][0], acm[i][1], acm[i][2]);
        work += dt * d * flux[2];
      } else {
        const f = fronts[i - 1];
        const ff = exactFlux(f.wp);
        for (let k = 0; k < 4; k++) flux[k] = ff[k] - f.G * f.wp[k];
      }
      // правая стенка
      if (i === regions.length - 1) {
        if (na >= n) {
          // свободный выход: поток по среднему в последнем объёме (экстраполяция наклона здесь неустойчива)
          for (let k = 0; k < 4; k++) fl[k] = P[4 * (n - 1) + k];
          const ff = exactFlux(fl as unknown as V4);
          for (let k = 0; k < 4; k++) flux[4 * n + k] = ff[k];
          // ушедшее через правый край (для балансов): ∫(p − p0) и энергия (мощность −T·v)
          outKin += dt * ff[0];
          outEnergy += dt * (ff[2] * fl[2] + ff[3] * fl[3]);
        }
      } else {
        const f = fronts[i];
        const ff = exactFlux(f.wm);
        for (let k = 0; k < 4; k++) flux[4 * n + k] = ff[k] - f.G * f.wm[k];
      }
      // новые длины: стенки-фронты сдвигаются
      const a1 = i === 0 ? 0 : fronts[i - 1].x + dt * fronts[i - 1].G;
      const b1 = i === regions.length - 1 ? L : fronts[i].x + dt * fronts[i].G;
      const upto = i === regions.length - 1 && na < n ? na : n;
      for (let j = 0; j < upto; j++) {
        const lo = r.xs[j + 1] - r.xs[j];
        const xa = j === 0 ? a1 : r.xs[j];
        const xb = j === n - 1 ? b1 : r.xs[j + 1];
        const ln = xb - xa;
        for (let k = 0; k < 4; k++) r.U[4 * j + k] = (lo * r.U[4 * j + k] - dt * (flux[4 * (j + 1) + k] - flux[4 * j + k])) / ln;
      }
      r.xs[0] = a1;
      r.xs[n] = b1;
    }
    // узкие области: содержимое меняется потоками через стенки (как у объёма)
    for (let i = 0; i < regions.length; i++) {
      const r = regions[i];
      if (!r.passive) continue;
      let fl: V4;
      if (i === 0) {
        r.wl = boundaryState(r.wl, d);
        fl = exactFlux(r.wl);
        work += dt * d * fl[2];
      } else {
        const f = fronts[i - 1];
        const ff = exactFlux(f.wp);
        fl = ff.map((v, k) => v - f.G * f.wp[k]) as V4;
      }
      let fr: V4;
      if (i === regions.length - 1) fr = exactFlux(r.wr);
      else {
        const f = fronts[i];
        const ff = exactFlux(f.wm);
        fr = ff.map((v, k) => v - f.G * f.wm[k]) as V4;
      }
      for (let k = 0; k < 4; k++) r.Q[k] += dt * (fl[k] - fr[k]);
    }
    // диссипация на выделенных ударных волнах (2.45)
    for (const f of fronts) {
      const [ahead, behind] = f.G > 0 ? [f.wp, f.wm] : [f.wm, f.wp];
      const tA = [0, 0];
      const tB = [0, 0];
      Tq(ahead[0], ahead[1], tA, 0);
      Tq(behind[0], behind[1], tB, 0);
      const pw = (T: number[], w: V4) => T[0] * w[2] + T[1] * w[3];
      const sg = f.G > 0 ? 1 : -1;
      shockDiss += dt * (sg * (pw(tA, ahead) - pw(tB, behind)) + Math.abs(f.G) * (E(ahead) - E(behind)));
    }
    for (const f of fronts) f.x += dt * f.G;
    // пассивные области между фронтами получают состояния с фронтов
    for (let i = 0; i < regions.length; i++) {
      const r = regions[i];
      if (!r.passive) continue;
      if (i > 0) r.wl = [...fronts[i - 1].wp];
      if (i < fronts.length) r.wr = [...fronts[i].wm];
    }
    // фронты не должны перескочить друг друга из-за округления
    for (let i = 0; i + 1 < fronts.length; i++) if (fronts[i + 1].x < fronts[i].x) fronts[i + 1].x = fronts[i].x;
    if (fronts.length && fronts[0].x < 0) fronts[0].x = 0;
    ub += d * dt;
    t += dt;
    steps++;
    for (const f of fronts) trackSeg(f);
    handleEventsAt();
    if (t >= tFrame - 1e-12 && nextFrame < nF) {
      snapshot();
      nextFrame++;
      if (nextFrame % 10 === 0) onProgress?.(t / tEnd);
    }
  }
  if (steps >= maxSteps) warnings.push(`Достигнут предел числа шагов (${maxSteps}).`);
  for (const f of fronts) closeSeg(f);
  while (frames.length < nF) snapshot();

  // балансы по последнему кадру
  let kin = 0;
  let en = 0;
  regions.forEach((r, i) => {
    const a = leftWall(i);
    const b = rightWall(i);
    if (r.passive) {
      kin += r.Q[0] - p0 * (b - a);
      if (b > a) {
        const w = r.Q.map((q) => q / (b - a));
        en += (energy1d(mat, [w[0], 0, w[1]]) - W0 + 0.5 * (w[2] ** 2 + w[3] ** 2)) * (b - a);
      }
      return;
    }
    for (let j = 0; j < r.n; j++) {
      const len = r.xs[j + 1] - r.xs[j];
      const o = 4 * j;
      kin += (r.U[o] - p0) * len;
      en += (energy1d(mat, [r.U[o], 0, r.U[o + 1]]) - W0 + 0.5 * (r.U[o + 2] ** 2 + r.U[o + 3] ** 2)) * len;
    }
  });
  void merges;
  return {
    frames,
    segments,
    steps,
    events,
    maxFronts,
    untracked,
    warnings,
    kinematicError: kin + ub + outKin,
    work,
    energyGain: en + outEnergy,
    shockDissipation: shockDiss,
    tEnd: t,
    xView: L,
    cellsPerL: N / L,
    N,
  };
}

/** Состояние в точке x кадра (линейная интерполяция по узлам; на фронте — справа, если right). */
export function hybridStateAt(fr: HybridFrame, F0: Vec3, x: number, right = false): State {
  const xs = fr.xs;
  const n = xs.length;
  let lo = 0;
  let hi = n - 1;
  if (x <= xs[0]) return node(fr, F0, 0);
  if (x >= xs[n - 1]) return node(fr, F0, n - 1);
  // первый узел с xs > x (или ≥ x для правой стороны)
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (right ? xs[mid] <= x : xs[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  const j = Math.max(1, lo);
  const xa = xs[j - 1];
  const xb = xs[j];
  const u = xb > xa ? (x - xa) / (xb - xa) : right ? 1 : 0;
  const A = node(fr, F0, j - 1);
  const B = node(fr, F0, j);
  return {
    F: [A.F[0] + u * (B.F[0] - A.F[0]), A.F[1] + u * (B.F[1] - A.F[1])],
    v: [A.v[0] + u * (B.v[0] - A.v[0]), A.v[1] + u * (B.v[1] - A.v[1])],
  };
}

function node(fr: HybridFrame, F0: Vec3, k: number): State {
  const W = fr.W;
  return { F: [F0[0] + W[4 * k], F0[2] + W[4 * k + 1]], v: [W[4 * k + 2], W[4 * k + 3]] };
}
