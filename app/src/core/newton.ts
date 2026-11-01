/** Плотная линейная алгебра и демпфированный метод Ньютона для небольших систем. */

/** Решение A x = b методом Гаусса с выбором ведущего элемента. Возвращает null, если матрица вырождена. */
export function solveLinear(A: number[][], b: number[]): number[] | null {
  const n = b.length;
  const M = A.map((row, i) => [...row, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    if (Math.abs(M[p][c]) < 1e-300) return null;
    [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < n; r++) {
      const f = M[r][c] / M[c][c];
      if (f === 0) continue;
      for (let j = c; j <= n; j++) M[r][j] -= f * M[c][j];
    }
  }
  const x = new Array<number>(n).fill(0);
  for (let r = n - 1; r >= 0; r--) {
    let s = M[r][n];
    for (let j = r + 1; j < n; j++) s -= M[r][j] * x[j];
    x[r] = s / M[r][r];
  }
  return x.every(Number.isFinite) ? x : null;
}

export function norm(v: number[]): number {
  return Math.sqrt(v.reduce((s, x) => s + x * x, 0));
}

/** Якобиан центральными разностями. */
export function jacobian(F: (x: number[]) => number[], x: number[], scale: number[]): number[][] {
  const n = x.length;
  const J: number[][] = [];
  const cols: number[][] = [];
  for (let j = 0; j < n; j++) {
    const h = 1e-6 * scale[j];
    const xp = [...x];
    const xm = [...x];
    xp[j] += h;
    xm[j] -= h;
    const fp = F(xp);
    const fm = F(xm);
    cols.push(fp.map((v, i) => (v - fm[i]) / (2 * h)));
  }
  for (let i = 0; i < cols[0].length; i++) J.push(cols.map((c) => c[i]));
  return J;
}

export interface NewtonResult {
  x: number[];
  residual: number;
  iterations: number;
  converged: boolean;
}

export interface NewtonOptions {
  tol?: number;
  maxIter?: number;
  /** Характерные масштабы неизвестных (для шага разностной производной). */
  scale?: number[];
  /** Допустима ли точка (например, ξ > 0). Недопустимые шаги дробятся. */
  admissible?: (x: number[]) => boolean;
}

/** Метод Ньютона с дроблением шага по норме невязки. */
export function newton(F: (x: number[]) => number[], x0: number[], opt: NewtonOptions = {}): NewtonResult {
  const tol = opt.tol ?? 1e-13;
  const maxIter = opt.maxIter ?? 60;
  const scale = opt.scale ?? x0.map((v) => Math.max(Math.abs(v), 1e-3));
  const ok = opt.admissible ?? (() => true);
  const safeF = (x: number[]): number[] | null => {
    if (!ok(x)) return null;
    try {
      const f = F(x);
      return f.every(Number.isFinite) ? f : null;
    } catch {
      return null;
    }
  };
  let x = [...x0];
  let f = safeF(x);
  if (!f) return { x, residual: Infinity, iterations: 0, converged: false };
  let r = norm(f);
  for (let it = 1; it <= maxIter; it++) {
    if (r < tol) return { x, residual: r, iterations: it - 1, converged: true };
    let J: number[][];
    try {
      J = jacobian((y) => safeF(y) ?? y.map(() => NaN), x, scale);
    } catch {
      return { x, residual: r, iterations: it, converged: false };
    }
    if (!J.every((row) => row.every(Number.isFinite))) return { x, residual: r, iterations: it, converged: false };
    const dx = solveLinear(J, f.map((v) => -v));
    if (!dx) return { x, residual: r, iterations: it, converged: false };
    let t = 1;
    let accepted = false;
    for (let k = 0; k < 30; k++) {
      const xn = x.map((v, i) => v + t * dx[i]);
      const fn = safeF(xn);
      if (fn) {
        const rn = norm(fn);
        if (rn < r * (1 - 1e-4 * t) || rn < tol) {
          x = xn;
          f = fn;
          r = rn;
          accepted = true;
          break;
        }
      }
      t *= 0.5;
    }
    if (!accepted) return { x, residual: r, iterations: it, converged: r < tol };
  }
  return { x, residual: r, iterations: maxIter, converged: r < tol };
}
