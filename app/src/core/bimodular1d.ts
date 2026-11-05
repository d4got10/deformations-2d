/**
 * Точное решение одноосной задачи для разномодульной среды Ляховского–Мясникова (α = β = 0, s ≡ 0):
 * σ(e) = a²e при e < 0, σ(e) = b²e при e > 0 (Лаптева, Дудко, Рагозина, 2022, формулы (2)–(4)).
 * Независимая от fronts.ts/hybrid.ts реализация для сверки.
 *
 * Материал кусочно-линеен, поэтому все волны — разрывы:
 *   ξ (простой разрыв, деформация одного знака по обе стороны) — скорость a (сжатие) или b (растяжение);
 *   γ (полусигнотон, с одной стороны e = 0) — скорость a или b;
 *   Σ (ударная волна, растяжение впереди → сжатие позади) — G² = (b²e₊ − a²e₋)/(e₊ − e₋), b < G < a;
 *   переход сжатие → растяжение — «жёсткий слой» [γᵇ, γᵃ]: два полусигнотона и e = 0 между ними.
 * Единицы безразмерные (как в остальном коде): скорости в c1 = √((λ+2μ)/ρ), v = u̇/c1, e = u,x.
 * Условия на разрыве (6): σ₊ − σ₋ = G²(e₊ − e₋), v₋ = v₊ + G(e₊ − e₋) для волны вправо
 * («+» — впереди, «−» — позади), v₋ = v₊ − G(e₊ − e₋) — для волны влево.
 */

export interface BState {
  e: number;
  v: number;
}

export type BWaveKind = 'xi' | 'gamma' | 'shock';

export interface BFront {
  x0: number;
  t0: number;
  /** Скорость со знаком. */
  G: number;
  kind: BWaveKind;
  /** 'a' — быстрая (сжатие), 'b' — медленная (растяжение), '' — ударная. */
  fam: 'a' | 'b' | '';
  left: BState;
  right: BState;
}

export interface BEvent {
  t: number;
  x: number;
  type: 'boundary' | 'collision' | 'reflection';
  incoming: string[];
  outgoing: string[];
}

const sgnTol = 0;

/** Волна, бегущая от состояния впереди A к состоянию позади с деформацией e (dir = ±1). */
function waveTo(a: number, b: number, A: BState, e: number, dir: 1 | -1): { v: number; pieces: { G: number; kind: BWaveKind; fam: 'a' | 'b' | ''; behind: BState }[] } {
  const pieces: { G: number; kind: BWaveKind; fam: 'a' | 'b' | ''; behind: BState }[] = [];
  const jump = (S: BState, eb: number, G: number, kind: BWaveKind, fam: 'a' | 'b' | ''): BState => {
    const st = { e: eb, v: S.v + dir * G * (S.e - eb) };
    pieces.push({ G: dir * G, kind, fam, behind: st });
    return st;
  };
  if (e === A.e) return { v: A.v, pieces };
  const kindFor = (e1: number, e2: number): BWaveKind => (e1 === 0 || e2 === 0 ? 'gamma' : 'xi');
  if (A.e >= sgnTol && e >= sgnTol) return { v: jump(A, e, b, kindFor(A.e, e), 'b').v, pieces };
  if (A.e <= sgnTol && e <= sgnTol) return { v: jump(A, e, a, kindFor(A.e, e), 'a').v, pieces };
  if (A.e > 0 && e < 0) {
    const G = Math.sqrt((b * b * A.e - a * a * e) / (A.e - e));
    return { v: jump(A, e, G, 'shock', '').v, pieces };
  }
  // сжатие впереди, растяжение позади: жёсткий слой — γᵃ (до e = 0), затем γᵇ
  const z = jump(A, 0, a, 'gamma', 'a');
  return { v: jump(z, e, b, 'gamma', 'b').v, pieces };
}

/** Деформация, при которой волна от A даёт скорость v (σ монотонна ⇒ v(e) монотонна). */
function solveE(a: number, b: number, A: BState, v: number, dir: 1 | -1): number {
  const f = (e: number) => waveTo(a, b, A, e, dir).v - v;
  let lo = -1;
  let hi = 1;
  for (let it = 0; it < 200; it++) {
    const mid = 0.5 * (lo + hi);
    // v(e) убывает по e для волны вправо и возрастает для волны влево
    if ((f(mid) > 0) === (dir > 0)) lo = mid;
    else hi = mid;
  }
  const e = 0.5 * (lo + hi);
  return Math.abs(e) < 1e-15 ? 0 : e;
}

/** Распад разрыва L | R: волны влево (от L) и вправо (от R) с общим средним состоянием. */
function riemann(a: number, b: number, L: BState, R: BState) {
  const g = (e: number) => waveTo(a, b, L, e, -1).v - waveTo(a, b, R, e, 1).v; // возрастает по e
  let lo = -1;
  let hi = 1;
  for (let it = 0; it < 200; it++) {
    const mid = 0.5 * (lo + hi);
    if (g(mid) > 0) hi = mid;
    else lo = mid;
  }
  let e = 0.5 * (lo + hi);
  if (Math.abs(e) < 1e-15) e = 0;
  return { left: waveTo(a, b, L, e, -1), right: waveTo(a, b, R, e, 1), e };
}

export interface BimodularProblem {
  /** Безразмерные скорости характеристик: a — сжатие, b — растяжение (a > b). */
  a: number;
  b: number;
  /** Скорость границы (u̇0/c1) — кусочно-постоянная: значения после моментов times[k]. */
  times: number[];
  velocities: number[];
  tEnd: number;
}

export interface BimodularResult {
  fronts: BFront[];
  /** Все отрезки траекторий (для диаграммы). */
  segments: { x0: number; t0: number; x1: number; t1: number; kind: BWaveKind; fam: string }[];
  events: BEvent[];
  /** Профиль в момент t. */
  stateAt: (x: number, t: number) => BState;
}

const label = (f: { G: number; kind: BWaveKind; fam: string }) => (f.kind === 'shock' ? 'Σ' : f.kind === 'gamma' ? 'γ' : 'ξ') + (f.G < 0 ? '⁻' : '') + (f.fam === 'a' ? 'ᵃ' : f.fam === 'b' ? 'ᵇ' : '');

export function solveBimodular(p: BimodularProblem): BimodularResult {
  const { a, b } = p;
  let fronts: BFront[] = [];
  const segments: BimodularResult['segments'] = [];
  const events: BEvent[] = [];
  const history: { t: number; fronts: BFront[]; boundary: number }[] = [];
  let vb = 0;
  let t = 0;
  const pos = (f: BFront, tt: number) => f.x0 + f.G * (tt - f.t0);
  const close = (f: BFront, tt: number) => segments.push({ x0: f.x0, t0: f.t0, x1: pos(f, tt), t1: tt, kind: f.kind, fam: f.fam });
  const stateNearBoundary = (): BState => (fronts.length ? fronts[0].left : { e: 0, v: 0 });
  const mk = (x: number, tt: number, left: BState, pcs: { G: number; kind: BWaveKind; fam: 'a' | 'b' | ''; behind: BState }[], dir: 1 | -1): BFront[] => {
    // pcs идут от «впереди» к «позади»; возвращаем фронты слева направо
    const out: BFront[] = [];
    if (dir > 0) {
      // впереди — справа: первый кусок самый правый
      let right = left; // здесь left — состояние впереди (справа)
      const arr: BFront[] = [];
      for (const pc of pcs) {
        arr.push({ x0: x, t0: tt, G: pc.G, kind: pc.kind, fam: pc.fam, left: pc.behind, right });
        right = pc.behind;
      }
      out.push(...arr.reverse());
    } else {
      let lft = left; // состояние впереди (слева)
      for (const pc of pcs) {
        out.push({ x0: x, t0: tt, G: pc.G, kind: pc.kind, fam: pc.fam, left: lft, right: pc.behind });
        lft = pc.behind;
      }
    }
    return out;
  };
  const boundary = (vNew: number, type: 'boundary' | 'reflection', incoming: string[]) => {
    const A = stateNearBoundary();
    const e = solveE(a, b, A, vNew, 1);
    const w = waveTo(a, b, A, e, 1);
    const nf = mk(0, t, A, w.pieces, 1);
    fronts = [...nf, ...fronts];
    vb = vNew;
    events.push({ t, x: 0, type, incoming, outgoing: nf.map(label) });
  };
  const changes = p.times.map((tt, k) => ({ t: tt, v: p.velocities[k] }));
  let ci = 0;
  history.push({ t: 0, fronts: [], boundary: 0 });
  for (let guard = 0; guard < 100000; guard++) {
    // ближайшее событие
    let tn = p.tEnd;
    let kind: 'change' | 'collide' | 'hit' | 'end' = 'end';
    let idx = -1;
    if (ci < changes.length && changes[ci].t <= tn) {
      tn = changes[ci].t;
      kind = 'change';
    }
    for (let i = 0; i + 1 < fronts.length; i++) {
      const A = fronts[i];
      const B = fronts[i + 1];
      if (A.G > B.G + 1e-14) {
        const tc = t + Math.max(0, pos(B, t) - pos(A, t)) / (A.G - B.G);
        if (tc < tn - 1e-14) {
          tn = tc;
          kind = 'collide';
          idx = i;
        }
      }
    }
    if (fronts.length && fronts[0].G < 0) {
      const th = t + pos(fronts[0], t) / -fronts[0].G;
      if (th < tn - 1e-14) {
        tn = th;
        kind = 'hit';
      }
    }
    t = tn;
    if (kind === 'end') break;
    if (kind === 'change') {
      boundary(changes[ci].v, 'boundary', []);
      ci++;
    } else if (kind === 'hit') {
      const f = fronts.shift()!;
      close(f, t);
      const A = f.right;
      const e = solveE(a, b, A, vb, 1);
      const w = waveTo(a, b, A, e, 1);
      const nf = mk(0, t, A, w.pieces, 1);
      fronts = [...nf, ...fronts];
      events.push({ t, x: 0, type: 'reflection', incoming: [label(f)], outgoing: nf.map(label) });
    } else {
      const A = fronts[idx];
      const B = fronts[idx + 1];
      const x = pos(A, t);
      close(A, t);
      close(B, t);
      const rp = riemann(a, b, A.left, B.right);
      const lf = mk(x, t, A.left, rp.left.pieces, -1);
      const rf = mk(x, t, B.right, rp.right.pieces, 1);
      fronts.splice(idx, 2, ...lf, ...rf);
      events.push({ t, x, type: 'collision', incoming: [label(A), label(B)], outgoing: [...lf, ...rf].map(label) });
    }
    history.push({ t, fronts: fronts.map((f) => ({ ...f })), boundary: vb });
  }
  for (const f of fronts) close(f, p.tEnd);
  const stateAt = (x: number, tt: number): BState => {
    let h = history[0];
    for (const q of history) if (q.t <= tt) h = q;
    const fs = h.fronts;
    if (!fs.length) return { e: 0, v: 0 };
    for (const f of fs) if (x < pos(f, tt)) return f.left;
    return fs[fs.length - 1].right;
  };
  return { fronts, segments, events, stateAt };
}
