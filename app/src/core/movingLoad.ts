/**
 * §3.2. Движение постоянной нагрузки со сверхсейсмической скоростью по границе
 * предварительно деформированного разномодульного упругого полупространства.
 *
 * Геометрия (рис. 19): x1 — вглубь полупространства, x2 — вдоль границы, нагрузка
 * σ11*, σ21* действует на участке x1 = 0, x2 ≤ st и движется со скоростью s = M·c1.
 * Автомодельная переменная ξ = x1/(st − x2) (3.4). От точки L0 расходятся две ударные
 * волны: квазипродольная Σ1 (ξ = ξ1 = tg θ1) и квазипоперечная Σ2 (ξ = ξ2 = tg θ2),
 * которые делят полупространство на зоны 0 (ξ > ξ1), I (ξ2 < ξ < ξ1), II (0 < ξ < ξ2).
 * В каждой зоне градиент перемещений постоянен: (a, b, k, n) (3.10)–(3.11).
 *
 * Неизвестные (3.36): a1, b1, k1, n1, a2, b2, k2, n2, ξ1, ξ2 — 10 уравнений.
 * Условия непрерывности перемещений (3.28) линейны относительно b и n, поэтому
 * b1, n1, b2, n2 выражены явно, и Ньютону остаются 6 неизвестных: a1, k1, a2, k2, ξ1, ξ2.
 *
 * Условия на разрыве выведены заново из закона сохранения импульса
 * [σ_ij ν_j] = −ρ G [v_i],  [f] = f(впереди) − f(позади),  ν = (cos θ, sin θ),  G = s sin θ,
 * что после деления на cos θ и (λ+2μ) даёт ровно (3.26)/(3.36):
 *   [σ̃_i1] + [σ̃_i2] ξ = M² ξ [b или n],  M² = ρs²/(λ+2μ).
 */
import { type Grad, type Material, type Stress, c2OverC1, energy, normalize, scaleNonlinear, stress } from './material.ts';
import { newton } from './newton.ts';

export interface MovingLoadProblem {
  material: Material;
  /** Предварительная деформация (3.15): a0 = u1,1 < 0, k0 = u2,1 > 0. b0 = n0 = 0 (покой, (3.16)). */
  a0: number;
  k0: number;
  /** Число Маха M = s/c1, c1 = √((λ+2μ)/ρ). */
  M: number;
  /** Нагрузка на границе, отнесённая к λ+2μ: σ̃11* < 0, σ̃21*. */
  s11: number;
  s21: number;
}

export interface Zone {
  name: '0' | 'I' | 'II';
  grad: Grad;
  /** Напряжения, отнесённые к λ+2μ. */
  stress: Stress;
  /** Скорость частиц, отнесённая к s: v/s = (−b, −n) (3.12). */
  v: [number, number];
  /** Упругая энергия, отнесённая к λ+2μ. */
  W: number;
}

export interface Front {
  name: 'Σ1' | 'Σ2';
  xi: number;
  /** Угол фронта с границей x1 = 0, градусы. */
  thetaDeg: number;
  /** Скорость фронта G/c1 = M sin θ (3.24). */
  G: number;
  /**
   * Диссипация на фронте (на единицу площади фронта, отнесённая к (λ+2μ)·s·cosθ):
   *   D = [σ̃_i1 v̂_i + ξ σ̃_i2 v̂_i] + ξ [W̃ + ½ M² |v̂|²] ≥ 0 — термодинамическое условие,
   * эквивалент (3.30) из закона сохранения энергии.
   */
  dissipation: number;
}

export interface MovingLoadSolution {
  problem: MovingLoadProblem;
  zones: [Zone, Zone, Zone];
  fronts: [Front, Front];
  /** Норма невязки системы (3.36) в полной (10-мерной) форме. */
  residual: number;
  iterations: number;
  /** Нарушенные физические ограничения: сжатие (3.18)/(3.19), термодинамика, порядок волн. */
  warnings: string[];
}

/** Восстановить b, n в зонах I и II из условий непрерывности перемещений (3.28). */
function grads(p: MovingLoadProblem, x: number[]): [Grad, Grad, Grad] {
  const [a1, k1, a2, k2, xi1, xi2] = x;
  const g0: Grad = { a: p.a0, b: 0, k: p.k0, n: 0 };
  // (a_ahead − a_behind) ξ − (b_ahead − b_behind) = 0
  const b1 = g0.b - (g0.a - a1) * xi1;
  const n1 = g0.n - (g0.k - k1) * xi1;
  const b2 = b1 - (a1 - a2) * xi2;
  const n2 = n1 - (k1 - k2) * xi2;
  return [g0, { a: a1, b: b1, k: k1, n: n1 }, { a: a2, b: b2, k: k2, n: n2 }];
}

/** Невязки уравнений импульса на фронте: ahead → behind. */
function jumpResidual(m: Material, M2: number, ga: Grad, gb: Grad, xi: number): [number, number] {
  const sa = stress(m, ga);
  const sb = stress(m, gb);
  return [
    sa.s11 - sb.s11 + (sa.s12 - sb.s12) * xi - M2 * xi * (ga.b - gb.b),
    sa.s12 - sb.s12 + (sa.s22 - sb.s22) * xi - M2 * xi * (ga.n - gb.n),
  ];
}

/** Невязка 6 уравнений для неизвестных [a1, k1, a2, k2, ξ1, ξ2]; m — нормированный материал. */
export function residual6(p: MovingLoadProblem, m: Material, x: number[]): number[] {
  const [g0, g1, g2] = grads(p, x);
  const M2 = p.M * p.M;
  const s2 = stress(m, g2);
  return [
    ...jumpResidual(m, M2, g0, g1, x[4]),
    ...jumpResidual(m, M2, g1, g2, x[5]),
    s2.s11 - p.s11,
    s2.s12 - p.s21,
  ];
}

/**
 * Полная система (3.36) из 10 уравнений, как в диссертации, — для независимой проверки решения.
 * Порядок неизвестных: a1, b1, k1, n1, a2, b2, k2, n2, ξ1, ξ2.
 */
export function residual10(p: MovingLoadProblem, z: number[]): number[] {
  const m = normalize(p.material);
  const [a1, b1, k1, n1, a2, b2, k2, n2, xi1, xi2] = z;
  const g0: Grad = { a: p.a0, b: 0, k: p.k0, n: 0 };
  const g1: Grad = { a: a1, b: b1, k: k1, n: n1 };
  const g2: Grad = { a: a2, b: b2, k: k2, n: n2 };
  const [S0, S1, S2] = [g0, g1, g2].map((g) => stress(m, g));
  const M2 = p.M * p.M;
  return [
    S0.s11 - S1.s11 + (S0.s12 - S1.s12) * xi1 + M2 * xi1 * b1,
    S0.s12 - S1.s12 + (S0.s22 - S1.s22) * xi1 + M2 * xi1 * n1,
    S1.s11 - S2.s11 + (S1.s12 - S2.s12) * xi2 - M2 * xi2 * (b1 - b2),
    S1.s12 - S2.s12 + (S1.s22 - S2.s22) * xi2 - M2 * xi2 * (n1 - n2),
    (p.a0 - a1) * xi1 + b1,
    (p.k0 - k1) * xi1 + n1,
    (a1 - a2) * xi2 - (b1 - b2),
    (k1 - k2) * xi2 - (n1 - n2),
    S2.s11 - p.s11,
    S2.s12 - p.s21,
  ];
}

/**
 * Начальное приближение — решение линейной теории упругости (ν = α = β = 0):
 * фронты — конусы Маха sin θ1 = 1/M, sin θ2 = (c2/c1)/M; скачок градиента на Σ1
 * поляризован по нормали ν (продольная волна), на Σ2 — по касательной τ (поперечная).
 * Амплитуды p, q подбираются из граничных условий.
 */
export function linearGuess(p: MovingLoadProblem): number[] {
  const m = normalize(scaleNonlinear(p.material, 0));
  const c2 = c2OverC1(p.material);
  if (p.M <= 1) throw new Error('Нагрузка должна быть сверхсейсмической: M > 1');
  const th1 = Math.asin(1 / p.M);
  const th2 = Math.asin(c2 / p.M);
  const [c1s, s1s] = [Math.cos(th1), Math.sin(th1)];
  const [c2s, s2s] = [Math.cos(th2), Math.sin(th2)];
  const build = (P: number, Q: number) => {
    const a1 = p.a0 + P * c1s * c1s;
    const k1 = p.k0 + P * s1s * c1s;
    const a2 = a1 - Q * s2s * c2s;
    const k2 = k1 + Q * c2s * c2s;
    return [a1, k1, a2, k2, Math.tan(th1), Math.tan(th2)];
  };
  // Линейная задача линейна по (P, Q): σ(II) аффинно зависит от них.
  const f = (P: number, Q: number) => {
    const g2 = grads(p, build(P, Q))[2];
    const s = stress(m, g2);
    return [s.s11 - p.s11, s.s12 - p.s21];
  };
  const f0 = f(0, 0);
  const fP = f(1, 0).map((v, i) => v - f0[i]);
  const fQ = f(0, 1).map((v, i) => v - f0[i]);
  const det = fP[0] * fQ[1] - fQ[0] * fP[1];
  const P = (-f0[0] * fQ[1] + fQ[0] * f0[1]) / det;
  const Q = (-fP[0] * f0[1] + f0[0] * fP[1]) / det;
  return build(P, Q);
}

const admissible = (x: number[]) => x[4] > 0 && x[5] > 0 && x.every(Number.isFinite);

/**
 * Решить задачу. Если задано начальное приближение x0 (например, из соседней точки
 * при изменении параметра), Ньютон стартует с него; иначе — продолжение по параметру
 * от линейной упругости (ε: 0 → 1 масштабирует ν, α, β).
 */
export function solveMovingLoad(p: MovingLoadProblem, x0?: number[]): MovingLoadSolution | null {
  if (!(p.M > 1)) throw new Error('Нагрузка должна быть сверхсейсмической: M > 1 (§3.2).');
  if (!(p.a0 < 0)) throw new Error('Предварительная деформация должна быть сжатием: a0 < 0 (3.15).');
  const mFull = normalize(p.material);
  const scale = [1e-3, 1e-3, 1e-3, 1e-3, 1, 1];
  const run = (m: Material, start: number[]) =>
    newton((x) => residual6(p, m, x), start, { scale, admissible, tol: 1e-14, maxIter: 80 });

  let x: number[] | null = null;
  let iterations = 0;
  if (x0) {
    const r = run(mFull, x0);
    iterations += r.iterations;
    if (r.converged) x = r.x;
  }
  if (!x) {
    let cur = linearGuess(p);
    let eps = 0;
    let step = 0.25;
    while (eps < 1) {
      const next = Math.min(1, eps + step);
      const r = run(normalize(scaleNonlinear(p.material, next)), cur);
      iterations += r.iterations;
      if (r.converged) {
        cur = r.x;
        eps = next;
        step = Math.min(step * 2, 0.5);
      } else {
        step /= 2;
        if (step < 1e-4) return null;
      }
    }
    x = cur;
  }
  return assemble(p, mFull, x, iterations);
}

function assemble(p: MovingLoadProblem, m: Material, x: number[], iterations: number): MovingLoadSolution {
  const gs = grads(p, x);
  const names = ['0', 'I', 'II'] as const;
  const zones = gs.map((g, i) => ({
    name: names[i],
    grad: g,
    stress: stress(m, g),
    v: [-g.b, -g.n] as [number, number],
    W: energy(m, g),
  })) as [Zone, Zone, Zone];
  const M2 = p.M * p.M;
  const front = (name: 'Σ1' | 'Σ2', xi: number, za: Zone, zb: Zone): Front => {
    const theta = Math.atan(xi);
    const flux = (z: Zone) =>
      z.stress.s11 * z.v[0] + z.stress.s12 * z.v[1] + xi * (z.stress.s12 * z.v[0] + z.stress.s22 * z.v[1]);
    const energyDensity = (z: Zone) => z.W + 0.5 * M2 * (z.v[0] ** 2 + z.v[1] ** 2);
    return {
      name,
      xi,
      thetaDeg: (theta * 180) / Math.PI,
      G: p.M * Math.sin(theta),
      dissipation: flux(za) - flux(zb) + xi * (energyDensity(za) - energyDensity(zb)),
    };
  };
  const fronts: [Front, Front] = [front('Σ1', x[4], zones[0], zones[1]), front('Σ2', x[5], zones[1], zones[2])];
  const z10 = [x[0], gs[1].b, x[1], gs[1].n, x[2], gs[2].b, x[3], gs[2].n, x[4], x[5]];
  const res = Math.hypot(...residual10(p, z10));

  const warnings: string[] = [];
  if (!(x[4] > x[5])) warnings.push('Порядок волн нарушен: ожидается θ1 > θ2.');
  const S0 = zones[0].stress;
  if (!(p.s11 < 0 && Math.abs(p.s11) > Math.abs(S0.s11)))
    warnings.push('Нарушено (3.18): нагрузка должна дополнительно сжимать среду, σ11* < 0, |σ11*| > |σ11(0)|.');
  for (const z of zones.slice(1)) {
    const J1 = z.grad.a + z.grad.n;
    if (!(J1 < 0 && Math.abs(J1) > Math.abs(p.a0))) warnings.push(`Нарушено (3.19) в зоне ${z.name}: ожидается дополнительное сжатие.`);
  }
  for (const f of fronts)
    if (f.dissipation < -1e-12) warnings.push(`На ${f.name} нарушено термодинамическое условие (3.30): D = ${f.dissipation.toExponential(2)} < 0.`);
  return { problem: p, zones, fronts, residual: res, iterations, warnings };
}

/** Перемещения точки (x1, η = x2 − st) в зоне с градиентом g: u1 = a x1 + b η, u2 = k x1 + n η. */
export function displacement(g: Grad, x1: number, eta: number): [number, number] {
  return [g.a * x1 + g.b * eta, g.k * x1 + g.n * eta];
}

/** Номер зоны (0, 1, 2) для точки (x1 ≥ 0, η = x2 − st) в подвижной системе координат. */
export function zoneIndex(sol: MovingLoadSolution, x1: number, eta: number): 0 | 1 | 2 {
  if (eta >= 0) return 0;
  const xi = x1 / -eta;
  if (xi > sol.fronts[0].xi) return 0;
  if (xi > sol.fronts[1].xi) return 1;
  return 2;
}

export interface SweepPoint {
  value: number;
  solution: MovingLoadSolution;
}

/**
 * Продолжение по параметру: решения для последовательности значений, каждое — от предыдущего.
 * Останавливается, когда решение перестаёт существовать (Ньютон не сходится) —
 * так находится, например, M_min из рис. 22а.
 */
export function sweep(
  base: MovingLoadProblem,
  values: number[],
  apply: (p: MovingLoadProblem, v: number) => MovingLoadProblem,
): { points: SweepPoint[]; failedAt: number | null } {
  const points: SweepPoint[] = [];
  let prev: number[] | undefined;
  for (const v of values) {
    const p = apply(base, v);
    let sol: MovingLoadSolution | null = null;
    try {
      sol = solveMovingLoad(p, prev);
    } catch {
      sol = null;
    }
    if (!sol) return { points, failedAt: v };
    points.push({ value: v, solution: sol });
    const z = sol.zones;
    prev = [z[1].grad.a, z[1].grad.k, z[2].grad.a, z[2].grad.k, sol.fronts[0].xi, sol.fronts[1].xi];
  }
  return { points, failedAt: null };
}
