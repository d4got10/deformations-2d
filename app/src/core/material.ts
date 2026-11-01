/**
 * Разномодульная изотропная упругая среда, модель Мясникова–Олейникова.
 * Дудко О.В., дисс. 1998, формулы (1.9), (1.10), (1.59), (1.60).
 *
 * Плоская деформация: u3 = 0, градиент перемещений задаётся четырьмя числами
 *   a = u1,1   b = u1,2   k = u2,1   n = u2,2
 * (обозначения §3.1, (3.11)).
 */

/** Упругие модули. Единицы любые, но одинаковые (в диссертации — 10^3 МПа). */
export interface Material {
  lambda: number;
  mu: number;
  nu: number;
  alpha: number;
  beta: number;
}

/** Градиент перемещений в плоскости (a, b, k, n). */
export interface Grad {
  a: number;
  b: number;
  k: number;
  n: number;
}

/** Ненулевые компоненты тензора напряжений Коши при плоской деформации. */
export interface Stress {
  s11: number;
  s12: number; // = s21
  s22: number;
  s33: number;
}

export interface Invariants {
  J1: number;
  J2: number;
  J3: number;
}

/** Инварианты тензора малых деформаций (1.60). */
export function invariants(g: Grad): Invariants {
  const { a, n } = g;
  const e12 = 0.5 * (g.b + g.k);
  return {
    J1: a + n,
    J2: a * a + n * n + 2 * e12 * e12,
    J3: a * a * a + n * n * n + 3 * e12 * e12 * (a + n),
  };
}

/**
 * Напряжения по (1.10):
 *   σ = (λ − ν/φ + 3αφ)·J1·δ + (2μ − νφ − αφ³ − βη)·e + 3β·e²/√J2,
 *   φ = J1/√J2,  η = J3/√J2³.
 * Слагаемое (ν/φ)·J1 записано как ν·√J2, чтобы не делить на J1 (при J1 = 0 оно конечно).
 * При J2 = 0 потенциал неаналитичен — бросаем исключение.
 */
export function stress(m: Material, g: Grad): Stress {
  const { J1, J2, J3 } = invariants(g);
  if (!(J2 > 0)) throw new Error('stress: J2 = 0, потенциал (1.9) неаналитичен в ненапряжённом состоянии');
  const r = Math.sqrt(J2);
  const phi = J1 / r;
  const eta = J3 / (r * r * r);
  const A = m.lambda * J1 - m.nu * r + 3 * m.alpha * phi * J1;
  const B = 2 * m.mu - m.nu * phi - m.alpha * phi ** 3 - m.beta * eta;
  const C = (3 * m.beta) / r;
  const { a, n } = g;
  const e12 = 0.5 * (g.b + g.k);
  return {
    s11: A + B * a + C * (a * a + e12 * e12),
    s22: A + B * n + C * (n * n + e12 * e12),
    s12: B * e12 + C * e12 * (a + n),
    s33: A,
  };
}

/** Упругий потенциал (1.9) (энергия на единицу объёма в отсчётной конфигурации). */
export function energy(m: Material, g: Grad): number {
  const { J1, J2, J3 } = invariants(g);
  if (!(J2 > 0)) return 0;
  const r = Math.sqrt(J2);
  return (
    0.5 * m.lambda * J1 * J1 +
    m.mu * J2 -
    m.nu * J1 * r +
    (m.alpha * J1 ** 3) / r +
    (m.beta * J3) / r
  );
}

/** Модули, отнесённые к λ+2μ (так обезразмерены напряжения в (3.36)). */
export function normalize(m: Material): Material {
  const s = m.lambda + 2 * m.mu;
  return { lambda: m.lambda / s, mu: m.mu / s, nu: m.nu / s, alpha: m.alpha / s, beta: m.beta / s };
}

/** Отношение скоростей поперечных и продольных волн линейной теории c2/c1 = √(μ/(λ+2μ)). */
export function c2OverC1(m: Material): number {
  return Math.sqrt(m.mu / (m.lambda + 2 * m.mu));
}

/** Масштабирование «нелинейной» части модели: ε = 0 — линейная упругость, ε = 1 — полная модель. */
export function scaleNonlinear(m: Material, eps: number): Material {
  return { ...m, nu: m.nu * eps, alpha: m.alpha * eps, beta: m.beta * eps };
}
