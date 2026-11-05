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
import { solveImpact } from './impact.ts';
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
  /** Если на встрече один фронт слабее другого в 1/passRatio раз — сильный проходит насквозь. */
  passRatio?: number;
  /** Предел числа событий. */
  maxEvents?: number;
  /**
   * Считать «до затухания»: пока возможна хоть одна встреча фронтов или отражение от границы
   * (tEnd тогда определяется автоматически, обрезки по x нет).
   */
  untilQuiet?: boolean;
  /** Внутреннее: остановиться, как только взаимодействия закончились. */
  stopWhenQuiet?: boolean;
  /** Горизонт режима «до затухания» (в единицах L/c1), по умолчанию 10. */
  maxTime?: number;
}

export interface FrontsFrame {
  t: number;
  /** Положения фронтов (по возрастанию) и состояния между ними (states.length = xs.length + 1). */
  xs: number[];
  /** Скорости фронтов в этот момент (до следующего события положения меняются линейно). */
  speeds: number[];
  states: State[];
  ub: number;
  /**
   * Принадлежность фронтов пучкам простых волн: соседние фронты с одинаковым номером ≥ 0 — лучи
   * одного центрированного пучка (одно семейство, общее место и время рождения); −1 — отдельный фронт.
   * В точном решении пучок — непрерывная волна; набор слабых фронтов — его кусочно-постоянная аппроксимация.
   */
  fan: number[];
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
  /** Сколько раз фронты составной волны (невыпуклость) сливались при рождении, и наибольшая относительная невязка [T] − G²[F] у слитых фронтов. */
  merges: number;
  mergeResidual: number;
  /** Сколько крошечных фронтов поглощено при встрече с сильными (в 1e5 раз и более). */
  absorbed: number;
  /** Сколько встреч не удалось решить и они обработаны как прохождение сильной волны сквозь слабую (≥100×). */
  fallbacks: number;
  /** Расчёт пришлось повторить с грубыми порогами (см. trackFronts). */
  coarsened?: boolean;
  /** Конечное время расчёта и длина участка, на котором находятся все фронты к этому моменту. */
  tEnd: number;
  xView: number;
  /** Момент, после которого фронты больше не встречаются (если считали «до затухания»). */
  tQuiet?: number;
  /** Встречи, не решённые даже запасными путями и обработанные линейным наложением. */
  linearPasses: number;
  /** Когда случилась бы ближайшая ещё возможная встреча после окончания расчёта (Infinity — таких нет). */
  nextInteraction: number;
  /** «До затухания» упёрлось в горизонт: остались далёкие встречи. */
  horizonReached?: boolean;
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
  const sc = Math.max(Math.hypot(X[0], X[1]), 1e-5);
  const dl = (h: number) => eig(c, [X[0] + h * r[0], X[1] + h * r[1]]).lam[k] - eig(c, [X[0] - h * r[0], X[1] - h * r[1]]).lam[k];
  // Порог — выше шума численного дифференцирования в eig (~1e-9·λ): λ однородна нулевой степени,
  // поэтому у истинно нелинейного семейства на шаге 1e-3·|X| она меняется на ~1e-3·λ.
  const tol = 1e-6 * Math.max(E.lam[k], 1e-6);
  const g = dl(1e-3 * sc);
  const gnl = Math.abs(g) > tol;
  // Если локально скорость постоянна (при s = 0 материал кусочно-линеен), направление сжатия
  // определяется на большем шаге — по тому, в какую сторону лежит более жёсткая область.
  let go = g;
  for (let h = 1e-2 * sc; !gnl && Math.abs(go) <= tol && h <= 3 * sc; h *= 3) go = dl(h);
  if (go < 0) r = [-r[0], -r[1]];
  return { r, lam: E.lam[k], gnl };
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
function wave(c: Ctx, a: State, f: Fam, dir: 1 | -1, eps: number, nFan: number, split = false): WaveOut {
  const R = orientedR(c, a.F, f);
  // число фронтов пучка задаётся извне и не меняется внутри итераций Ньютона (иначе невязка разрывна)
  // Разрежение разбивается на слабые скачки и у линейно вырожденного в точке a семейства: при s = 0
  // материал кусочно-линеен (разные модули при p < 0 и p > 0), и веер, проходящий через p = 0, —
  // это два разрыва со скоростями c₋ и c₊ и полка p = 0 между ними, а не один разрыв (он был бы
  // недопустимой «ударной волной разрежения»). Скачки одной скорости потом снова сливаются.
  // split: сжимающая волна тоже разбивается на слабые скачки — нужно при невыпуклости, когда одна
  // ударная волна недопустима (D < 0): догоняющие скачки потом сливаются в ударную волну, остальные
  // образуют веер позади неё (составная волна «ударная + веер»).
  if (split && eps > 0) return splitWave(c, a, f, dir, eps, nFan, R);
  const n = eps >= 0 ? 1 : nFan;
  const pieces: WaveOut['pieces'] = [];
  let cur = a;
  let prevDir: V2 | null = null;
  for (let j = 0; j < n; j++) {
    const Rj = j === 0 ? R : orientedR(c, cur.F, f);
    // Пучок идёт по волновой кривой непрерывно (без разворотов). Если семейство невыпукло
    // (скорость вдоль кривой проходит через экстремум, как у материала 1 при растяжении),
    // дальше экстремума шаги становятся сжимающими скачками: они догоняют друг друга и
    // сливаются в ударную волну — получается составная волна «веер + ударная».
    let e = eps / n;
    if (prevDir && (-Rj.r[0] * prevDir[0] - Rj.r[1] * prevDir[1]) * Math.sign(-e) < 0) e = -e;
    let { Xb, G2 } = hugoniot(c, cur.F, f, e, Rj);
    // Одноосный случай (s = 0): материал кусочно-линеен, излом σ(p) при p = 0. Шаг, пересекающий
    // p = 0, делится ровно в нуле — тогда полка «жёсткого слоя» [γᵇ, γᵃ] точно e = 0, а скорости
    // по обе стороны — точно c₋ и c₊ (без промежуточного шага с секущей скоростью).
    let rest: number | null = null;
    // (только для веера: сжимающий скачок растяжение → сжатие — это одна ударная волна Σ)
    if (eps < 0 && cur.F[1] === 0 && Math.abs(Xb[1]) <= 1e-18 && cur.F[0] * Xb[0] < 0) {
      const frac = cur.F[0] / (cur.F[0] - Xb[0]);
      rest = e * (1 - frac);
      ({ Xb, G2 } = hugoniot(c, cur.F, f, e * frac, Rj));
      Xb = [0, 0];
      const T = T2(c, cur.F);
      G2 = T[0] / cur.F[0];
    }
    const pushStep = (Xb: V2, G2: number, e: number) => {
      const G = Math.sqrt(Math.max(G2, 0));
      // [v] = −G_signed [F] с [·] = справа − слева  ⇒  v_b = v_a + dir·G·(F_a − F_b)
      const v: V2 = [cur.v[0] + dir * G * (cur.F[0] - Xb[0]), cur.v[1] + dir * G * (cur.F[1] - Xb[1])];
      const behind = { F: Xb, v };
      // ступенька разрежения (e < 0 у истинно нелинейного семейства) — фронт простой волны, даже если
      // он один (так выходят лучи пучка из задачи о распаде разрыва после пересечения с другой волной)
      pieces.push({ G: dir * G, behind, kind: (R.gnl || n > 1) && e < 0 ? 'simple' : 'shock' });
      const dl = Math.hypot(Xb[0] - cur.F[0], Xb[1] - cur.F[1]);
      if (dl > 0) prevDir = [(Xb[0] - cur.F[0]) / dl, (Xb[1] - cur.F[1]) / dl];
      cur = behind;
    };
    pushStep(Xb, G2, e);
    if (rest !== null) {
      const R0 = orientedR(c, cur.F, f);
      let e2 = rest;
      if (prevDir && (-R0.r[0] * prevDir[0] - R0.r[1] * prevDir[1]) * Math.sign(-e2) < 0) e2 = -e2;
      const h2 = hugoniot(c, cur.F, f, e2, R0);
      // из нуля: однородность — G² = T(Xb)·Xb/|Xb|² точно
      const Tb = T2(c, h2.Xb);
      const d2 = h2.Xb[0] ** 2 + h2.Xb[1] ** 2;
      pushStep(h2.Xb, d2 > 0 ? (Tb[0] * h2.Xb[0] + Tb[1] * h2.Xb[1]) / d2 : h2.G2, e2);
    }
  }
  // Соседние шаги одной скорости (линейная ветвь материала, например сжатие или растяжение при s = 0)
  // — один контактный разрыв: он не расходится, поэтому ведётся как разрыв, а не как луч веера.
  if (pieces.length > 1) {
    const merged: WaveOut['pieces'] = [];
    let cnt = 0;
    for (const pc of pieces) {
      const last = merged.at(-1);
      if (last && Math.abs(last.G - pc.G) <= 1e-9 * Math.abs(pc.G)) {
        last.behind = pc.behind;
        cnt++;
        last.kind = 'shock';
      } else {
        merged.push({ ...pc });
        cnt = 1;
      }
    }
    void cnt;
    pieces.length = 0;
    pieces.push(...merged);
  }
  // Одноосный случай (s ≡ 0): материал кусочно-линеен, веера нет — разрежение внутри ветви —
  // контактный разрыв, через p = 0 — два контактных разрыва (жёсткий слой). Все они не расходятся
  // и ведутся как разрывы (у ненапряжённого состояния нет определённой скорости характеристик,
  // поэтому оценка «ширины веера» по ней была бы неверна).
  if (a.F[1] === 0 && pieces.every((pc) => pc.behind.F[1] === 0)) for (const pc of pieces) pc.kind = 'shock';
  return { pieces, behind: cur };
}

/**
 * Составная волна «ударная + веер» при невыпуклости (правило Лю). Ударная волна растёт по ударной
 * адиабате из a, пока её скорость не сравняется со скоростью характеристики своего семейства позади
 * (звуковая точка ε*): дальше одна ударная волна стала бы недопустимой (D < 0). Если амплитуда больше
 * ε*, остаток уходит веером от звукового состояния. ε* зависит только от a, поэтому волна непрерывна
 * по амплитуде — это важно для Ньютона. Каждый разрыв точно удовлетворяет условиям на разрыве.
 */
function splitWave(c: Ctx, a: State, f: Fam, dir: 1 | -1, eps: number, nFan: number, R: { r: V2; gnl: boolean }): WaveOut {
  const k = fi(f);
  const h = (e: number) => {
    const { Xb, G2 } = hugoniot(c, a.F, f, e, R);
    return Math.sqrt(Math.max(G2, 0)) - Math.sqrt(Math.max(eig(c, Xb).lam[k], 0));
  };
  const M = 24;
  let eStar = -1;
  let prev = 0;
  for (let j = 1; j <= M; j++) {
    const e = (eps * j) / M;
    if (h(e) >= 0) {
      let lo = prev;
      let hi = e;
      for (let it = 0; it < 40; it++) {
        const mid = 0.5 * (lo + hi);
        if (h(mid) >= 0) hi = mid;
        else lo = mid;
      }
      eStar = hi;
      break;
    }
    prev = e;
  }
  const shockTo = (e: number): State => {
    const { Xb, G2 } = hugoniot(c, a.F, f, e, R);
    const G = Math.sqrt(Math.max(G2, 0));
    return { F: Xb, v: [a.v[0] + dir * G * (a.F[0] - Xb[0]), a.v[1] + dir * G * (a.F[1] - Xb[1])] };
  };
  const shockG = (e: number) => Math.sqrt(Math.max(hugoniot(c, a.F, f, e, R).G2, 0));
  if (eStar < 0) {
    const b = shockTo(eps);
    return { pieces: [{ G: dir * shockG(eps), behind: b, kind: 'shock' }], behind: b };
  }
  const b = shockTo(eStar);
  const pieces: WaveOut['pieces'] = [{ G: dir * shockG(eStar), behind: b, kind: 'shock' }];
  const fan = wave(c, b, f, dir, -(eps - eStar), Math.max(1, nFan));
  pieces.push(...fan.pieces);
  return { pieces, behind: fan.behind };
}

/**
 * Допустимость волны: на каждом сжимающем скачке диссипация (2.45) D = G([W] − ⟨T⟩·[F]) ≥ 0
 * ([·] = справа − слева, G со знаком). Около p = 0 материал сильно нелинеен, и Ньютон может сойтись
 * к «ударной волне разрежения» (D < 0) — такой корень отбрасывается.
 */
function admissibleWave(c: Ctx, a: State, w: WaveOut, dir: 1 | -1, rel = 1e-7): boolean {
  // Проверяется только волна, которая целиком одна ударная волна. У веера (в том числе составного
  // «веер + ударная») отдельные шаги — дискретизация пути, а не разрывы решения: шаг, перешагнувший
  // экстремум скорости, формально даёт D < 0 порядка шага, хотя после слияния всё допустимо.
  if (w.pieces.length !== 1) return true;
  let cur = a;
  for (const pc of w.pieces) {
    if (pc.kind === 'shock') {
      const [L, R] = dir > 0 ? [pc.behind.F, cur.F] : [cur.F, pc.behind.F];
      const TL = T2(c, L);
      const TR = T2(c, R);
      const dF = [R[0] - L[0], R[1] - L[1]];
      const work = 0.5 * ((TL[0] + TR[0]) * dF[0] + (TL[1] + TR[1]) * dF[1]);
      const D = pc.G * (energy1d(c.m, [R[0], 0, R[1]]) - energy1d(c.m, [L[0], 0, L[1]]) - work);
      if (D < -rel * Math.abs(pc.G * work) - 1e-20) return false;
    }
    cur = pc.behind;
  }
  return true;
}

/** Число фронтов пучка для простой волны амплитуды eps. */
const fanCount = (eps: number, fanStep: number) => Math.max(1, Math.ceil(Math.abs(eps) / fanStep));

/**
 * Решить систему с фиксированными числами фронтов пучков; если после схождения амплитуды
 * простых волн сильно отличаются от исходных, числа фронтов пересчитываются и решение повторяется.
 */
function solveWithFans(F: (e: number[], n: number[]) => number[], g: number[], fanStep: number, split = false): { e: number[]; n: number[] } | null {
  let n = g.map((x) => fanCount(x, fanStep));
  let x0 = g;
  for (let pass = 0; pass < 3; pass++) {
    const e = newtonSolve((x) => F(x, n), x0, 1e-5);
    if (!e) return null;
    const n2 = e.map((x) => fanCount(x, fanStep));
    if (n2.every((v, i) => v === n[i] || (e[i] >= 0 && !split))) return { e, n };
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
  let split = false;
  const run = (e: number[], n: number[]) => {
    const w1 = wave(c, a, 'fast', 1, e[0], n[0], split);
    const w2 = wave(c, w1.behind, 'slow', 1, e[1], n[1], split);
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
  const solveFor = (target: V2, guess: number[]) =>
    solveWithFans(
      (e, n) => {
        const { w2 } = run(e, n);
        return [w2.behind.v[0] - target[0], w2.behind.v[1] - target[1]];
      },
      guess,
      fanStep,
      split,
    );
  let rejected = false;
  const okSol = (q: { e: number[]; n: number[] } | null) => {
    if (!q) return null;
    const { w1, w2 } = run(q.e, q.n);
    // при разбиении скачки — дискретизация составной волны: допуск соразмерен их величине
    const rel = split ? 1e-2 : 1e-7;
    if (admissibleWave(c, a, w1, 1, rel) && admissibleWave(c, w1.behind, w2, 1, rel)) return q;
    rejected = true;
    return null;
  };
  const attempt = () => {
    let sol = okSol(solveFor([d, 0], g));
    // сошлось, но к недопустимой ударной волне — нужна составная волна, сразу к разбиению
    if (!sol && rejected && !split) return null;
    if (!sol) {
      // Продолжение по параметру: целевая скорость меняется от скорости среды (волн нет) до заданной.
      // Нужно, когда волна проходит через сильно нелинейную область около p = 0.
      let e = [0, 0];
      let ok = true;
      for (let k = 1; k <= 16 && ok; k++) {
        const lam = k / 16;
        const target: V2 = [(1 - lam) * a.v[0] + lam * d, (1 - lam) * a.v[1]];
        const sk = solveFor(target, e);
        if (!sk) ok = false;
        else e = sk.e;
      }
      if (ok) sol = okSol(solveFor([d, 0], e));
    }
    if (!sol) {
      // Начальное приближение из автомодельного решения главы 2 (impact.ts): амплитуды — проекции
      // скачков между его зонами на собственные векторы семейств.
      try {
        const ex = solveImpact({ material: c.m, F0: [a.F[0], 0, a.F[1]], v0: [a.v[0], 0, a.v[1]], bc: { type: 'velocity', d }, noWaveCurves: true });
        if (ex) {
          const zs = ex.zones.map((z) => [z.F[0], z.F[2]] as V2);
          const iFast = ex.waves.findIndex((w) => w.family === 'fast');
          const iSlow = ex.waves.findIndex((w) => w.family === 'slow');
          const proj = (i: number, fam: Fam) => {
            if (i < 0) return 0;
            const A0 = zs[i];
            const B0 = zs[i + 1];
            const r = orientedR(c, A0, fam).r;
            const sgn = ex.waves[i].kind === 'fan' ? -1 : 1;
            return sgn * Math.abs((B0[0] - A0[0]) * r[0] + (B0[1] - A0[1]) * r[1]);
          };
          sol = okSol(solveFor([d, 0], [proj(iFast, 'fast'), proj(iSlow, 'slow')]));
        }
      } catch {
        /* остаётся null */
      }
    }
    return sol;
  };
  let sol = attempt();
  if (!sol) {
    // невыпуклость: одна ударная волна недопустима — составная волна «ударная + веер»
    split = true;
    sol = attempt();
  }
  if (!sol) return null;
  const { w1, w2 } = run(sol.e, sol.n);
  return { fast: w1, slow: w2 };
}

/**
 * Распад разрыва между состояниями l и r (встреча фронтов): влево уходят быстрая и медленная
 * волны, вправо — медленная и быстрая; в середине общее состояние.
 */
export function riemannProblem(c: Ctx, l: State, r: State, fanStep: number) {
  let split = false;
  const run = (e: number[], n: number[]) => {
    const LF = wave(c, l, 'fast', -1, e[0], n[0], split);
    const LS = wave(c, LF.behind, 'slow', -1, e[1], n[1], split);
    const RF = wave(c, r, 'fast', 1, e[3], n[3], split);
    const RS = wave(c, RF.behind, 'slow', 1, e[2], n[2], split);
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
  const solveFor = (rr: State, guess: number[]) =>
    solveWithFans(
      (e, n) => {
        const LF = wave(c, l, 'fast', -1, e[0], n[0], split);
        const LS = wave(c, LF.behind, 'slow', -1, e[1], n[1], split);
        const RF = wave(c, rr, 'fast', 1, e[3], n[3], split);
        const RS = wave(c, RF.behind, 'slow', 1, e[2], n[2], split);
        const a = LS.behind;
        const b = RS.behind;
        return [a.F[0] - b.F[0], a.F[1] - b.F[1], a.v[0] - b.v[0], a.v[1] - b.v[1]];
      },
      guess,
      fanStep,
      split,
    );
  let rejected = false;
  const admissible = (q: { e: number[]; n: number[] } | null, rr: State) => {
    if (!q) return null;
    const LF = wave(c, l, 'fast', -1, q.e[0], q.n[0], split);
    const LS = wave(c, LF.behind, 'slow', -1, q.e[1], q.n[1], split);
    const RF = wave(c, rr, 'fast', 1, q.e[3], q.n[3], split);
    const RS = wave(c, RF.behind, 'slow', 1, q.e[2], q.n[2], split);
    const rel = split ? 1e-2 : 1e-7;
    if (admissibleWave(c, l, LF, -1, rel) && admissibleWave(c, LF.behind, LS, -1, rel) && admissibleWave(c, rr, RF, 1, rel) && admissibleWave(c, RF.behind, RS, 1, rel)) return q;
    rejected = true;
    return null;
  };
  const attempt = () => {
    let sol = admissible(solveFor(r, g), r);
    if (!sol && rejected && !split) return null;
    if (!sol) {
      // Продолжение по параметру: правое состояние меняется от левого (волн нет) до заданного —
      // решение остаётся на физической ветви.
      let e = [0, 0, 0, 0];
      let ok = true;
      for (let k = 1; k <= 16 && ok; k++) {
        const lam = k / 16;
        const rr: State = { F: [l.F[0] + lam * (r.F[0] - l.F[0]), l.F[1] + lam * (r.F[1] - l.F[1])], v: [l.v[0] + lam * (r.v[0] - l.v[0]), l.v[1] + lam * (r.v[1] - l.v[1])] };
        const sk = solveFor(rr, k === 1 ? g.map((x) => x / 16) : e);
        if (!sk) ok = false;
        else e = sk.e;
      }
      if (ok) sol = admissible(solveFor(r, e), r);
    }
    return sol;
  };
  let sol = attempt();
  if (!sol) {
    // невыпуклость: одна ударная волна недопустима — составная волна «ударная + веер»
    split = true;
    sol = attempt();
  }
  if (!sol) return null;
  return run(sol.e, sol.n);
}

// ---------------------------------------------------------------- отслеживание фронтов

interface Tracked extends Front {
  strength: number;
}

const strengthOf = (a: State, b: State) => Math.hypot(b.F[0] - a.F[0], b.F[1] - a.F[1]);

/**
 * Отслеживание фронтов. Сначала — с точными порогами (слабые фронты от 1e-9, прохождение насквозь
 * при разнице в 1e5 раз). Если событий слишком много (сильная нелинейность около p = 0 у материала 2
 * порождает лавину сопоставимых отражений) или встреча не решается, расчёт повторяется с грубыми
 * порогами (1e-7 и 1e2) — это отмечается в результате (coarsened).
 */
export function trackFronts(prob: FrontsProblem): FrontsResult {
  const fine = runMode({ ...prob, maxEvents: prob.maxEvents ?? 1500 });
  if (fine.events < (prob.maxEvents ?? 1500) && fine.warnings.length === 0) return fine;
  const coarse = runMode({ ...prob, dropTol: Math.max(prob.dropTol ?? 0, 1e-7), passRatio: Math.max(prob.passRatio ?? 0, 1e-2), maxEvents: 20000 });
  coarse.coarsened = true;
  return coarse;
}

/** «До затухания»: пробный проход находит момент последней встречи, затем — расчёт до него с кадрами. */
function runMode(prob: FrontsProblem): FrontsResult {
  if (!prob.untilQuiet) return trackOnce(prob);
  const lastChange = Math.max(0, ...prob.pulses.map((p) => p.start + p.duration));
  const horizon = Math.max(prob.maxTime ?? 10, lastChange * 1.2);
  const probe = trackOnce({ ...prob, tEnd: horizon, frames: 2, xMax: Infinity, stopWhenQuiet: true });
  if (probe.warnings.length) return probe;
  const reached = probe.tQuiet === undefined;
  const tq = Math.max(probe.tQuiet ?? horizon, lastChange);
  const tEnd = reached ? horizon : Math.max(0.3, tq * 1.1, lastChange * 1.15);
  const res = trackOnce({ ...prob, tEnd, xMax: Infinity, untilQuiet: false });
  res.tQuiet = reached ? undefined : tq;
  res.horizonReached = reached;
  return res;
}

function trackOnce(prob: FrontsProblem): FrontsResult {
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
  let merges = 0;
  let absorbed = 0;
  let fallbacks = 0;
  let linearPasses = 0;
  const passRatio = prob.passRatio ?? 1e-5;
  let mergeResidual = 0;
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
    // Составная волна при невыпуклости: соседние фронты одного семейства и направления, которые
    // рождаются «в обратном порядке» (задний быстрее переднего), сразу сливаются в один ударный фронт.
    // Его скорость — из условий на разрыве по МНК: G² = [T]·[F]/|[F]|² (невязка ~ шаг пучка³).
    for (let k = 0; k + 1 < newF.length; ) {
      const A = newF[k];
      const B = newF[k + 1];
      if (A.fam === B.fam && Math.sign(A.speed) === Math.sign(B.speed) && A.speed > B.speed + 1e-13) {
        const L = newS[k];
        const R = newS[k + 2];
        const TL = T2(c, L.F);
        const TR = T2(c, R.F);
        const dF: V2 = [R.F[0] - L.F[0], R.F[1] - L.F[1]];
        const G2 = ((TR[0] - TL[0]) * dF[0] + (TR[1] - TL[1]) * dF[1]) / (dF[0] ** 2 + dF[1] ** 2);
        const G = Math.sign(A.speed) * Math.sqrt(Math.max(G2, 0));
        const resid = Math.hypot(TR[0] - TL[0] - G2 * dF[0], TR[1] - TL[1] - G2 * dF[1]) / Math.max(Math.hypot(TR[0] - TL[0], TR[1] - TL[1]), 1e-30);
        mergeResidual = Math.max(mergeResidual, resid);
        newF.splice(k, 2, { ...A, speed: G, kind: 'shock', strength: strengthOf(L, R) });
        newS.splice(k + 1, 1);
        if (k > 0) k--; // слитый фронт может оказаться быстрее предыдущего
        merges++;
      } else k++;
    }
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

  // лучи одного пучка: соседние фронты простой волны одного семейства, идущие в одну сторону,
  // со скачками одного знака (после пересечений с другими волнами лучи «перерождаются», поэтому
  // время и место рождения не годятся как признак)
  const fanIds = () => {
    const ids: number[] = [];
    let id = -1;
    fronts.forEach((f, i) => {
      const p = fronts[i - 1];
      if (f.kind !== 'simple') return ids.push(-1);
      let same = !!p && p.kind === 'simple' && p.fam === f.fam && Math.sign(p.speed) === Math.sign(f.speed);
      if (same) {
        const a = states[i - 1];
        const b = states[i];
        const c2 = states[i + 1];
        const d1 = [b.F[0] - a.F[0], b.F[1] - a.F[1]];
        const d2 = [c2.F[0] - b.F[0], c2.F[1] - b.F[1]];
        same = d1[0] * d2[0] + d1[1] * d2[1] > 0;
      }
      if (!same) id++;
      ids.push(id);
    });
    return ids;
  };
  const snapshot = () =>
    frames.push({ t, xs: fronts.map((f) => f.x), speeds: fronts.map((f) => f.speed), states: states.map((s) => ({ F: [...s.F] as V2, v: [...s.v] as V2 })), ub, fan: fanIds() });

  // начальная краевая задача
  const changes = [...new Set(prob.pulses.flatMap((p) => [p.start, p.start + p.duration]))].filter((x) => x >= 0 && x < prob.tEnd).sort((a, b) => a - b);
  let nextFrame = 0;
  const appliedBc = new Set<number>();
  const applyBoundary = () => {
    const d = boundaryVelocity(prob.pulses, t + 1e-13);
    const bw = boundaryWaves(states[0], d);
    if (!bw) {
      warnings.push(`t = ${t.toFixed(4)}: краевая задача не решена`);
      if ((globalThis as any).__dbg) console.log('BPFAIL', JSON.stringify({ a: states[0], d }));
      return;
    }
    insertWaves(0, 0, bw.left, bw.seq, states[0]);
    events++;
  };

  const maxEvents = prob.maxEvents ?? 20000;
  let tQuiet: number | undefined;
  while (t < prob.tEnd - 1e-13 && events < maxEvents) {
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
    if (prob.stopWhenQuiet && (kind === 'frame' || kind === 'end')) {
      // ни встреч, ни отражений, ни новых воздействий впереди — взаимодействия закончились
      const pendingBc = changes.some((b) => b > t + 1e-13 || (b >= t - 1e-13 && !appliedBc.has(b)));
      if (!pendingBc) {
        tQuiet = t;
        break;
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
      const m0 = states[which + 1];
      const r = states[which + 2];
      // Встреча сильного фронта с крошечным (в 1e5 раз слабее): сильный проходит насквозь, слабый
      // поглощается (ошибка в балансах ~ скачку слабого). Задача о распаде разрыва для таких пар плохо
      // обусловлена и лишь плодит ещё более слабые отражения.
      if (Math.min(a.strength, b.strength) < passRatio * Math.max(a.strength, b.strength)) {
        const keep = a.strength >= b.strength ? a : b;
        fronts.splice(which, 2, { ...keep, x, born: t, xBorn: x, strength: strengthOf(l, r) });
        states.splice(which, 3, l, r);
        absorbed++;
        events++;
        while (fronts.length && fronts.at(-1)!.x > xMax) {
          closeSeg(fronts.at(-1)!, t);
          fronts.pop();
          states.pop();
        }
        maxFronts = Math.max(maxFronts, fronts.length);
        continue;
      }
      fronts.splice(which, 2);
      states.splice(which, 3, l); // между l и r будет вставлено решение задачи о распаде разрыва
      const rp = riemannProblem(c, l, r, fanStep);
      if (!rp && Math.min(a.strength, b.strength) < 1e-2 * Math.max(a.strength, b.strength)) {
        // запасной путь: задача не решилась, но одна волна намного слабее — сильная проходит насквозь
        const keep = a.strength >= b.strength ? a : b;
        states.splice(which + 1, 0, r); // сейчас states[which] = l; восстанавливаем l | r
        fronts.splice(which, 0, { ...keep, x, born: t, xBorn: x, strength: strengthOf(l, r) });
        fallbacks++;
        events++;
        continue;
      }
      if (!rp) {
        // последний запасной путь: волны проходят друг сквозь друга по линейному наложению
        // (новое среднее состояние = l + (r − m)); точно в линейной теории, ошибка ~ нелинейность × амплитуды
        const mid = { F: [l.F[0] + r.F[0] - m0.F[0], l.F[1] + r.F[1] - m0.F[1]] as V2, v: [l.v[0] + r.v[0] - m0.v[0], l.v[1] + r.v[1] - m0.v[1]] as V2 };
        states.splice(which + 1, 0, mid, r);
        fronts.splice(which, 0, { ...b, x, born: t, xBorn: x, strength: strengthOf(l, mid) }, { ...a, x, born: t, xBorn: x, strength: strengthOf(mid, r) });
        linearPasses++;
        events++;
        continue;
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
  if (events >= maxEvents) warnings.push(`Достигнут предел числа событий (${maxEvents}).`);
  for (const f of fronts) closeSeg(f, t);
  while (frames.length < nF) snapshot();

  // кинематический баланс: ∫(p − p0)dx + ∫v*dt (по последнему кадру; за последним фронтом — исходное состояние)
  const last = frames.at(-1)!;
  const xEndInt = Number.isFinite(xMax) ? xMax : (last.xs.at(-1) ?? 0) + 1;
  let kin = 0;
  for (let i = 0; i < last.states.length; i++) {
    const xa = i === 0 ? 0 : last.xs[i - 1];
    const xb = i < last.xs.length ? last.xs[i] : xEndInt;
    kin += (last.states[i].F[0] - F0[0]) * (xb - xa);
  }
  let en = 0;
  for (let i = 0; i < last.states.length; i++) {
    const xa = i === 0 ? 0 : last.xs[i - 1];
    const xb = i < last.xs.length ? last.xs[i] : xEndInt;
    const s = last.states[i];
    en += (energy1d(c.m, [s.F[0], 0, s.F[1]]) - W0 + 0.5 * (s.v[0] ** 2 + s.v[1] ** 2)) * (xb - xa);
  }
  // ближайшая ещё возможная встреча/отражение после окончания расчёта
  let nextInteraction = Infinity;
  for (let i = 0; i + 1 < fronts.length; i++)
    if (fronts[i].speed > fronts[i + 1].speed + 1e-12) nextInteraction = Math.min(nextInteraction, t + Math.max(0, fronts[i + 1].x - fronts[i].x) / (fronts[i].speed - fronts[i + 1].speed));
  if (fronts.length && fronts[0].speed < 0) nextInteraction = Math.min(nextInteraction, t + fronts[0].x / -fronts[0].speed);
  const xFar = Math.max(0, ...segments.map((sg) => Math.max(sg.x0, sg.x1)));
  const xView = Number.isFinite(xMax) ? Math.min(xMax, 1) : Math.max(1e-3, xFar * 1.04);
  return { frames, segments, events, maxFronts, warnings, kinematicError: kin + last.ub, work, energyGain: en, dissipation, merges, mergeResidual, absorbed, fallbacks, linearPasses, nextInteraction, tEnd: t, xView, tQuiet };
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
