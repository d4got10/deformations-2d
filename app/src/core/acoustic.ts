/**
 * Скорости малых возмущений (акустический тензор) в однородно деформированной зоне
 * и проверка условия эволюционности Лакса для ударной волны.
 * Это независимая от (3.30) проверка допустимости разрыва.
 */
import type { Grad, Material } from './material.ts';
import { stress } from './material.ts';

type H = [number, number, number, number]; // H11=a, H12=b, H21=k, H22=n

function sigma(m: Material, h: H): [[number, number], [number, number]] {
  const s = stress(m, { a: h[0], b: h[1], k: h[2], n: h[3] });
  return [
    [s.s11, s.s12],
    [s.s12, s.s22],
  ];
}

/**
 * Квадраты скоростей плоских волн с нормалью ν в зоне с градиентом g, в единицах c1²
 * (модули m должны быть нормированы на λ+2μ). Возвращает [медленная, быстрая].
 * A_ik = Σ_jl ∂σ_ij/∂u_{k,l} ν_j ν_l.
 */
export function waveSpeedsSquared(m: Material, g: Grad, nu: [number, number]): [number, number] {
  const h0: H = [g.a, g.b, g.k, g.n];
  const A = [
    [0, 0],
    [0, 0],
  ];
  for (let k = 0; k < 2; k++)
    for (let l = 0; l < 2; l++) {
      const idx = 2 * k + l;
      const d = 1e-7 * Math.max(1e-3, Math.abs(h0[idx]));
      const hp = [...h0] as H;
      const hm = [...h0] as H;
      hp[idx] += d;
      hm[idx] -= d;
      const sp = sigma(m, hp);
      const sm = sigma(m, hm);
      for (let i = 0; i < 2; i++)
        for (let j = 0; j < 2; j++) A[i][k] += ((sp[i][j] - sm[i][j]) / (2 * d)) * nu[j] * nu[l];
    }
  const tr = A[0][0] + A[1][1];
  const det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
  const disc = Math.sqrt(Math.max(0, (tr * tr) / 4 - det));
  return [tr / 2 - disc, tr / 2 + disc];
}

export interface LaxCheck {
  /** Скорость характеристик того же типа впереди и позади фронта, в единицах c1. */
  ahead: number;
  behind: number;
  G: number;
  /** Условие Лакса: behind > G > ahead. */
  ok: boolean;
}

/** Условие Лакса для фронта типа family (0 — квазипоперечный, 1 — квазипродольный). */
export function laxCheck(m: Material, ahead: Grad, behind: Grad, thetaRad: number, G: number, family: 0 | 1): LaxCheck {
  const nu: [number, number] = [Math.cos(thetaRad), Math.sin(thetaRad)];
  const ca = Math.sqrt(waveSpeedsSquared(m, ahead, nu)[family]);
  const cb = Math.sqrt(waveSpeedsSquared(m, behind, nu)[family]);
  return { ahead: ca, behind: cb, G, ok: cb > G && G > ca };
}
