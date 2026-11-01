/**
 * Нестационарный численный решатель одномерной задачи (глава 2 без предположения автомодельности):
 * серия нормальных ударов по границе x1 = 0 предварительно деформированного полупространства.
 *
 * Безразмерные переменные: x̂ = x1/L, t̂ = c1 t/L, F = (p, g, s) = (u1,1, u2,1, u3,1), v̂ = v/c1,
 * T̂ = σ_i1/(λ+2μ). Система законов сохранения (ρ c1² = λ+2μ):
 *     ∂F/∂t̂ = ∂v̂/∂x̂,     ∂v̂/∂t̂ = ∂T̂(F)/∂x̂.
 * Схема: конечные объёмы, MUSCL (minmod) + поток Русанова, SSP-RK2. Граница x̂ = 0 — заданная скорость
 * v̂*(t̂) = (d(t̂), 0, 0) через зеркальные фиктивные ячейки (v_g = 2v* − v, F_g = F); справа —
 * свободный выход. Ударные волны схема выделяет сама (численная вязкость обеспечивает допустимость).
 */
import type { Material } from './material.ts';
import { type Vec3, energy1d } from './material1d.ts';

/** Импульс: граница движется вглубь со скоростью d·c1 при t̂ ∈ [start, start + duration). */
export interface Pulse {
  start: number;
  duration: number;
  d: number;
}

export interface SeriesProblem {
  /** Модули, нормированные на λ+2μ. */
  material: Material;
  F0: Vec3;
  pulses: Pulse[];
  /** Число ячеек на отрезке x̂ ∈ [0, 1]. */
  N: number;
  /** Конечное время t̂. */
  tEnd: number;
  /** Сколько кадров сохранить. */
  frames: number;
  cfl?: number;
}

export interface SeriesFrame {
  t: number;
  /** p, g, s, v1, v2, v3 по ячейкам: [k*N + i]. */
  U: Float32Array;
  /** Перемещение границы u1(0, t) (в единицах L). */
  ub: number;
}

export interface SeriesResult {
  N: number;
  dx: number;
  frames: SeriesFrame[];
  /** Баланс: кинематика — ∫(F − F0)dx + ∫v* dt (должно быть ≈ 0, точно для схемы). */
  kinematicError: number;
  /** Энергия: работа границы, прирост энергии среды, численно рассеянная энергия. */
  work: number;
  energyGain: number;
  steps: number;
}

export function boundaryVelocity(pulses: Pulse[], t: number): number {
  let d = 0;
  for (const p of pulses) if (t >= p.start && t < p.start + p.duration) d += p.d;
  return d;
}

/** Напряжения T = (σ11, σ21, σ31) без выделения памяти (модель (1.9)). */
function tractionInto(m: Material, p: number, g: number, s: number, out: Float64Array, o: number) {
  const q = g * g + s * s;
  const J2 = Math.max(p * p + q / 2, 1e-30);
  const r = Math.sqrt(J2);
  const phi = p / r;
  const eta = (p * p * p + 0.75 * p * q) / (r * r * r);
  const A = m.lambda * p - m.nu * r + 3 * m.alpha * phi * p;
  const B = 2 * m.mu - m.nu * phi - m.alpha * phi * phi * phi - m.beta * eta;
  const C = (3 * m.beta) / r;
  const k = (B + C * p) / 2;
  out[o] = A + B * p + C * (p * p + q / 4);
  out[o + 1] = k * g;
  out[o + 2] = k * s;
}

const minmod = (a: number, b: number) => (a * b <= 0 ? 0 : Math.abs(a) < Math.abs(b) ? a : b);

export function simulate(prob: SeriesProblem, onProgress?: (fraction: number) => void): SeriesResult {
  const { material: m, F0, N, tEnd } = prob;
  const cfl = prob.cfl ?? 0.4;
  const G = 2; // фиктивные ячейки с каждой стороны
  const M = N + 2 * G;
  const dx = 1 / N;
  const U = Array.from({ length: 6 }, () => new Float64Array(M));
  const U1 = Array.from({ length: 6 }, () => new Float64Array(M));
  const R = Array.from({ length: 6 }, () => new Float64Array(M));
  for (let i = 0; i < M; i++) for (let k = 0; k < 3; k++) U[k][i] = F0[k];
  const T = new Float64Array(6);
  const slope = Array.from({ length: 6 }, () => new Float64Array(M));
  const cf = new Float64Array(M);
  const flux = Array.from({ length: 6 }, () => new Float64Array(M + 1));
  const TL = new Float64Array(3);
  const TR = new Float64Array(3);
  const L = new Float64Array(6);
  const Rr = new Float64Array(6);
  let boundaryT1 = 0; // напряжение σ11 на границе (из потока) — для работы

  const fillGhosts = (W: Float64Array[], vstar: number) => {
    for (let j = 0; j < G; j++) {
      const src = G + j;
      const dst = G - 1 - j;
      for (let k = 0; k < 3; k++) W[k][dst] = W[k][src];
      W[3][dst] = 2 * vstar - W[3][src];
      W[4][dst] = -W[4][src];
      W[5][dst] = -W[5][src];
      const srcR = G + N - 1 - j;
      const dstR = G + N + j;
      for (let k = 0; k < 6; k++) W[k][dstR] = W[k][srcR];
    }
  };

  /** Правая часть dU/dt; возвращает максимальную скорость. */
  const rhs = (W: Float64Array[], vstar: number): number => {
    fillGhosts(W, vstar);
    let amax = 0;
    for (let i = 0; i < M; i++) {
      // оценка скорости быстрых возмущений: √(∂σ11/∂p) с запасом
      const h = 1e-7 * Math.max(Math.abs(W[0][i]), 1e-4);
      tractionInto(m, W[0][i] + h, W[1][i], W[2][i], T, 0);
      tractionInto(m, W[0][i] - h, W[1][i], W[2][i], T, 3);
      cf[i] = 1.1 * Math.sqrt(Math.max((T[0] - T[3]) / (2 * h), 1e-6));
      if (i >= G && i < G + N) amax = Math.max(amax, cf[i]);
    }
    for (let k = 0; k < 6; k++)
      for (let i = 1; i < M - 1; i++) slope[k][i] = minmod(W[k][i] - W[k][i - 1], W[k][i + 1] - W[k][i]);
    // интерфейсы i+1/2 между ячейками i и i+1, i = G−1 … G+N−1
    for (let i = G - 1; i < G + N; i++) {
      for (let k = 0; k < 6; k++) {
        L[k] = W[k][i] + 0.5 * slope[k][i];
        Rr[k] = W[k][i + 1] - 0.5 * slope[k][i + 1];
      }
      tractionInto(m, L[0], L[1], L[2], TL, 0);
      tractionInto(m, Rr[0], Rr[1], Rr[2], TR, 0);
      const a = Math.max(cf[i], cf[i + 1]);
      for (let k = 0; k < 3; k++) {
        flux[k][i] = 0.5 * (-L[3 + k] - Rr[3 + k]) - 0.5 * a * (Rr[k] - L[k]);
        flux[3 + k][i] = 0.5 * (-TL[k] - TR[k]) - 0.5 * a * (Rr[3 + k] - L[3 + k]);
      }
      if (i === G - 1) boundaryT1 = -flux[3][i];
    }
    for (let k = 0; k < 6; k++) for (let i = G; i < G + N; i++) R[k][i] = -(flux[k][i] - flux[k][i - 1]) / dx;
    return amax;
  };

  const frames: SeriesFrame[] = [];
  const nF = Math.max(2, prob.frames);
  let nextFrame = 0;
  let t = 0;
  let ub = 0;
  let kinInt = 0; // ∫ v* dt
  let work = 0;
  const W0 = energy1d(m, F0);
  const energy = (W: Float64Array[]) => {
    let e = 0;
    for (let i = G; i < G + N; i++) e += (energy1d(m, [W[0][i], W[1][i], W[2][i]]) - W0 + 0.5 * (W[3][i] ** 2 + W[4][i] ** 2 + W[5][i] ** 2)) * dx;
    return e;
  };
  const snap = () => {
    const A = new Float32Array(6 * N);
    for (let k = 0; k < 6; k++) for (let i = 0; i < N; i++) A[k * N + i] = U[k][G + i];
    frames.push({ t, U: A, ub });
  };
  const breakTimes = prob.pulses.flatMap((p) => [p.start, p.start + p.duration]).filter((x) => x > 0 && x < tEnd);
  let steps = 0;
  snap();
  nextFrame = 1;
  while (t < tEnd - 1e-12) {
    const v0 = boundaryVelocity(prob.pulses, t + 1e-12);
    const amax = rhs(U, v0);
    let dt = (cfl * dx) / amax;
    const tFrame = (nextFrame * tEnd) / (nF - 1);
    // не перешагивать моменты включения/выключения импульсов и кадры
    for (const b of [...breakTimes, tFrame, tEnd]) if (b > t + 1e-12 && t + dt > b) dt = b - t;
    // SSP-RK2 (граничная скорость постоянна на шаге)
    const tb1 = boundaryT1;
    for (let k = 0; k < 6; k++) for (let i = G; i < G + N; i++) U1[k][i] = U[k][i] + dt * R[k][i];
    rhs(U1, v0);
    const tb2 = boundaryT1;
    for (let k = 0; k < 6; k++) for (let i = G; i < G + N; i++) U[k][i] = 0.5 * (U[k][i] + U1[k][i] + dt * R[k][i]);
    work += dt * v0 * -0.5 * (tb1 + tb2); // мощность: −σ11(0)·v*
    kinInt += dt * v0;
    ub += dt * v0;
    t += dt;
    steps++;
    if (t >= tFrame - 1e-12 && nextFrame < nF) {
      snap();
      nextFrame++;
      onProgress?.(t / tEnd);
    }
  }
  let kin = 0;
  for (let i = G; i < G + N; i++) kin += (U[0][i] - F0[0]) * dx;
  return {
    N,
    dx,
    frames,
    kinematicError: kin + kinInt,
    work,
    energyGain: energy(U),
    steps,
  };
}

/** Значения в ячейке i кадра. */
export function cell(frame: SeriesFrame, N: number, i: number): { F: Vec3; v: Vec3 } {
  const U = frame.U;
  return { F: [U[i], U[N + i], U[2 * N + i]], v: [U[3 * N + i], U[4 * N + i], U[5 * N + i]] };
}
