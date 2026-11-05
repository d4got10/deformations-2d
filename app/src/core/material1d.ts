/**
 * Модель Мясникова–Олейникова для одномерной деформации (глава 2):
 * перемещения зависят только от x1 и t, но все три компоненты ненулевые.
 * Состояние — вектор F = (p, g, s) = (u1,1, u2,1, u3,1).
 *
 * Тензор малых деформаций: e11 = p, e12 = g/2, e13 = s/2, остальные 0.
 * Напряжения — общая формула (1.10): σ = A·δ + B·e + C·e²,
 *   A = λJ1 − ν√J2 + 3αφJ1,  B = 2μ − νφ − αφ³ − βη,  C = 3β/√J2.
 */
import type { Material } from './material.ts';

export type Vec3 = [number, number, number];

export interface Stress3 {
  s11: number;
  s22: number;
  s33: number;
  s12: number;
  s13: number;
  s23: number;
}

/**
 * Коэффициенты закона (1.10). В ненапряжённом состоянии (J2 = 0) потенциал (1.9) не гладок
 * (однороден второй степени, «конус»), но непрерывен вместе с напряжениями: W(0) = 0, T(0) = 0 —
 * тогда возвращается null, и вызывающие функции дают нули. Жёсткость в нуле зависит от направления
 * (сжатие и растяжение — разные модули), её дают секущие (см. acousticAt в fronts.ts/hybrid.ts).
 */
function coeffs(m: Material, F: Vec3) {
  const [p, g, s] = F;
  const q = g * g + s * s;
  const J1 = p;
  const J2 = p * p + q / 2;
  const J3 = p * p * p + 0.75 * p * q;
  if (!(J2 > 0)) return null;
  const r = Math.sqrt(J2);
  const phi = J1 / r;
  const eta = J3 / (r * r * r);
  return {
    A: m.lambda * J1 - m.nu * r + 3 * m.alpha * phi * J1,
    B: 2 * m.mu - m.nu * phi - m.alpha * phi ** 3 - m.beta * eta,
    C: (3 * m.beta) / r,
    J1,
    J2,
    J3,
    r,
  };
}

/** Все компоненты тензора напряжений. */
export function stress3(m: Material, F: Vec3): Stress3 {
  const [p, g, s] = F;
  const c = coeffs(m, F);
  if (!c) return { s11: 0, s12: 0, s13: 0, s22: 0, s33: 0, s23: 0 };
  const { A, B, C } = c;
  return {
    s11: A + B * p + C * (p * p + (g * g + s * s) / 4),
    s12: (B * g) / 2 + (C * p * g) / 2,
    s13: (B * s) / 2 + (C * p * s) / 2,
    s22: A + (C * g * g) / 4,
    s33: A + (C * s * s) / 4,
    s23: (C * g * s) / 4,
  };
}

/** Вектор напряжений на площадке x1 = const: T = (σ11, σ21, σ31) = ∂W/∂F. */
export function traction(m: Material, F: Vec3): Vec3 {
  const [p, g, s] = F;
  const c = coeffs(m, F);
  if (!c) return [0, 0, 0];
  const { A, B, C } = c;
  const k = (B + C * p) / 2;
  return [A + B * p + C * (p * p + (g * g + s * s) / 4), k * g, k * s];
}

/**
 * Сдвиговая жёсткость Φ: (σ21, σ31) = Φ·(g, s). При одномерной деформации изотропной среды
 * касательное напряжение всегда коллинеарно сдвигу. Φ — квадрат скорости нейтральной волны
 * (волны поворота) в единицах c1², если модули нормированы на λ+2μ.
 */
export function shearStiffness(m: Material, F: Vec3): number {
  const c = coeffs(m, F);
  // в нуле сдвиговая жёсткость зависит от направления; берём предел вдоль чистого сдвига (p = 0)
  if (!c) return m.mu;
  return (c.B + c.C * F[0]) / 2;
}

/** Упругий потенциал (1.9). */
export function energy1d(m: Material, F: Vec3): number {
  const c = coeffs(m, F);
  if (!c) return 0;
  const { J1, J2, J3, r } = c;
  return 0.5 * m.lambda * J1 * J1 + m.mu * J2 - m.nu * J1 * r + (m.alpha * J1 ** 3) / r + (m.beta * J3) / r;
}

/**
 * Явные формулы (2.29) диссертации для σ11, σ21, σ31 при p < 0 (k = −1), θ = (g² + s²)/(2p²).
 * signFix = +1 — как напечатано (… + (ν + 3β/2)θ в σ21), −1 — как следует из (1.9) и (1.59).
 */
export function traction229(m: Material, F: Vec3, signFix: 1 | -1 = 1): Vec3 {
  const [p, g, s] = F;
  const th = (g * g + s * s) / (2 * p * p);
  const d = Math.sqrt((1 + th) ** 3);
  const s11 = (m.lambda + 2 * m.mu + ((m.nu - m.alpha - m.beta) * (2 + 3 * th) + (m.nu - 1.5 * m.beta) * th * th) / d) * p;
  const k = m.mu + (0.5 * (m.nu + m.alpha - 2 * m.beta + (m.nu + signFix * 1.5 * m.beta) * th)) / d;
  return [s11, k * g, k * s];
}
