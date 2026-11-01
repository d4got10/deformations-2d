/**
 * Материалы и расчётные случаи из §3.2 диссертации (стр. 106–111).
 * Модули — в 10^3 МПа.
 */
import type { Material } from './material.ts';
import type { MovingLoadProblem } from './movingLoad.ts';

/** Материал 1 (стр. 106): ν < 0, α + β < ν. */
export const MATERIAL_1: Material = { lambda: 38.7, mu: 34.8, nu: -4.8, alpha: -2.3, beta: -3.3 };

/** Материал 2 (стр. 109): ν > 0, α + β < ν. */
export const MATERIAL_2: Material = { lambda: 4.8, mu: 4.7, nu: 2.1, alpha: 1.2, beta: -0.3 };

/** Базовый случай рис. 21/23: a0 = −0.005, k0 = 0.001, σ̃11* = −0.01, σ̃21* = 0.003, M = 1.7. */
export const BASE_CASE: Omit<MovingLoadProblem, 'material'> = {
  a0: -0.005,
  k0: 0.001,
  M: 1.7,
  s11: -0.01,
  s21: 0.003,
};

/** Случай рис. 22б,в / 24б,в: M = 1.4, a0 = −0.003, k0 = 0.001. */
export const LOAD_SWEEP_CASE: Omit<MovingLoadProblem, 'material'> = {
  a0: -0.003,
  k0: 0.001,
  M: 1.4,
  s11: -0.01,
  s21: 0.003,
};

export const PRESETS = {
  'Материал 1 (ν < 0), рис. 21': { material: MATERIAL_1, ...BASE_CASE },
  'Материал 2 (ν > 0), рис. 23': { material: MATERIAL_2, ...BASE_CASE },
} satisfies Record<string, MovingLoadProblem>;

/**
 * Значения из текста диссертации для сверки.
 * Табличные числа — как напечатаны; «graph» — грубо считаны с рисунков.
 */
export const DISSERTATION = {
  material1: { G1: 1.0087, G2: 0.5593, theta1: 36.4, theta2: 19.2, Mmin: 1.11 },
  material2: { G1: 1.0575, G2: 0.5978, theta1: 38.5, theta2: 20.6, Mmin: 1.005 },
  /** Рис. 21 (материал 1), считано с графика. Индексы: зона 0, I, II. */
  fig21: {
    a: [-0.005, -0.0095, -0.0118],
    b: [0, -0.0014, -0.0028],
    k: [0.001, -0.0002, 0.012],
    n: [0, -0.0012, 0.0028],
    s11: [-0.005, -0.0074, -0.01],
    s21: [0.0003, -0.0005, 0.003],
    s22: [-0.0026, -0.004, -0.0011],
  },
  /** Рис. 22а (материал 1): θ(M), считано с графика и с пометок на полях. */
  fig22a: [
    { M: 1.11, theta1: 63.7, theta2: 29.5 },
    { M: 2.5, theta1: 29.0, theta2: 20.0 },
  ],
} as const;
