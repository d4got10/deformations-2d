/**
 * Масштабируемый график профиля величины вдоль x1 (canvas).
 *
 * Данные передаются не точками, а «как есть»: кусочно-постоянный профиль метода фронтов — списком
 * положений фронтов и значений между ними, точное решение — функцией. Поэтому при любом увеличении
 * рисуется точная картинка: ступени — по истинным положениям фронтов, функции — с шагом в пиксель.
 *
 * Колесо мыши — масштаб по x1 вокруг курсора (общий для всех графиков группы), Shift + колесо — по
 * вертикали; перетаскивание — сдвиг; двойной щелчок или «⤢» — весь вид. Под курсором — значения.
 */

import { fmt } from './charts.ts';

export interface StepData {
  /** Положения разрывов по возрастанию. */
  edges: number[];
  /** Значения: vals[0] — левее edges[0], vals[k] — между edges[k−1] и edges[k]; длина edges.length + 1. */
  vals: number[];
}

export interface PlotSeries {
  label: string;
  color: string;
  dash?: number[];
  steps?: StepData;
  /** Если задана — рисуется она (с шагом в пиксель), а steps используются только для автомасштаба. */
  fn?: (x: number) => number;
}

/** Общее окно по x1 для группы графиков. */
export interface SharedX {
  win: [number, number] | null;
  full: () => [number, number];
  redraw: () => void;
}

const css = (name: string) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

/** Подпись деления с точностью, достаточной для шага между делениями. */
function tickLabel(v: number, step: number): string {
  if (Math.abs(v) < step * 1e-6) return '0';
  const dec = Math.max(0, Math.ceil(-Math.log10(step) + 1e-9) + (Math.abs(step / 10 ** Math.floor(Math.log10(step)) - 2.5) < 1e-6 ? 1 : 0));
  if (dec <= 6 && Math.abs(v) < 1e5) return v.toFixed(dec);
  const sig = Math.max(1, Math.ceil(Math.log10(Math.abs(v) / step) + 1e-9) + 1);
  return v.toExponential(Math.min(8, sig));
}

function ticks(a: number, b: number): number[] {
  const raw = (b - a) / 5;
  if (!(raw > 0)) return [a];
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((k) => k * mag).find((v) => v >= raw) ?? raw;
  const out: number[] = [];
  for (let v = Math.ceil(a / step) * step; v <= b + 1e-9 * (b - a); v += step) out.push(+v.toPrecision(10));
  return out;
}

/** Значение ступенчатой функции в точке. */
function stepAt(s: StepData, x: number): number {
  let lo = 0;
  let hi = s.edges.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (s.edges[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return s.vals[lo];
}

export function createProfilePlot(parent: HTMLElement, title: string, xLabel: string, shared: SharedX) {
  const wrap = document.createElement('div');
  wrap.className = 'prof-plot';
  const cv = document.createElement('canvas');
  cv.setAttribute('aria-label', `Профиль: ${title}`);
  wrap.appendChild(cv);
  parent.appendChild(wrap);
  let series: PlotSeries[] = [];
  let yWin: [number, number] | null = null;
  let hover: number | null = null; // x курсора в пикселях
  const m = { l: 78, r: 14, t: 34, b: 38 };
  let geom = { w: 1, h: 1, x0: 0, x1: 1, y0: 0, y1: 1 };

  const xwin = () => shared.win ?? shared.full();

  /** Автомасштаб по вертикали по видимой части. */
  function autoY(x0: number, x1: number, pw: number): [number, number] {
    let lo = Infinity;
    let hi = -Infinity;
    for (const s of series) {
      if (s.steps) {
        const { edges, vals } = s.steps;
        for (let k = 0; k < vals.length; k++) {
          const a = k === 0 ? -Infinity : edges[k - 1];
          const b = k === edges.length ? Infinity : edges[k];
          if (b < x0 || a > x1) continue;
          lo = Math.min(lo, vals[k]);
          hi = Math.max(hi, vals[k]);
        }
      }
      if (s.fn)
        for (let i = 0; i <= pw; i += 2) {
          const v = s.fn(x0 + ((x1 - x0) * i) / pw);
          if (Number.isFinite(v)) {
            lo = Math.min(lo, v);
            hi = Math.max(hi, v);
          }
        }
    }
    if (!Number.isFinite(lo)) return [-1, 1];
    if (hi - lo < 1e-14 * Math.max(1, Math.abs(hi))) {
      const d = Math.abs(hi) * 0.1 || 1e-6;
      return [lo - d, hi + d];
    }
    const pad = 0.08 * (hi - lo);
    return [lo - pad, hi + pad];
  }

  function draw(next?: PlotSeries[]) {
    if (next) series = next;
    const dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth;
    const h = cv.clientHeight;
    if (!w || !h) return;
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(h * dpr);
    }
    const ctx = cv.getContext('2d')!;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = css('--panel');
    ctx.fillRect(0, 0, w, h);
    const [x0, x1] = xwin();
    const pw = w - m.l - m.r;
    const ph = h - m.t - m.b;
    const [y0, y1] = yWin ?? autoY(x0, x1, pw);
    geom = { w, h, x0, x1, y0, y1 };
    const X = (x: number) => m.l + ((x - x0) / (x1 - x0)) * pw;
    const Y = (y: number) => m.t + ((y1 - y) / (y1 - y0)) * ph;
    // сетка и оси
    ctx.font = '11px system-ui, sans-serif';
    ctx.strokeStyle = css('--border');
    ctx.fillStyle = css('--muted');
    ctx.lineWidth = 1;
    const ty = ticks(y0, y1);
    const tx = ticks(x0, x1);
    const sy = ty.length > 1 ? ty[1] - ty[0] : Math.abs(y1 - y0);
    const sx = tx.length > 1 ? tx[1] - tx[0] : Math.abs(x1 - x0);
    for (const v of ty) {
      ctx.beginPath();
      ctx.moveTo(m.l, Y(v));
      ctx.lineTo(w - m.r, Y(v));
      ctx.stroke();
      ctx.textAlign = 'right';
      ctx.fillText(tickLabel(v, sy), m.l - 6, Y(v) + 4);
    }
    ctx.textAlign = 'center';
    for (const v of tx) {
      ctx.beginPath();
      ctx.moveTo(X(v), m.t);
      ctx.lineTo(X(v), h - m.b);
      ctx.stroke();
      ctx.fillText(tickLabel(v, sx), X(v), h - m.b + 15);
    }
    ctx.textAlign = 'right';
    ctx.fillText(xLabel, w - m.r, h - 6);
    ctx.textAlign = 'left';
    ctx.strokeStyle = css('--border');
    ctx.strokeRect(m.l, m.t, pw, ph);
    ctx.font = '600 12px system-ui, sans-serif';
    ctx.fillStyle = css('--text');
    ctx.fillText(title, 10, 20);
    const titleW = ctx.measureText(title).width;
    const zoom = (shared.full()[1] - shared.full()[0]) / (x1 - x0);
    if (zoom > 1.001 || yWin) {
      ctx.font = '11px system-ui, sans-serif';
      ctx.fillStyle = css('--muted');
      ctx.fillText(`×${fmt(zoom, 3)}${yWin ? ' · вертикаль вручную' : ''}`, 18 + titleW, 20);
    }
    // данные
    ctx.save();
    ctx.beginPath();
    ctx.rect(m.l, m.t, pw, ph);
    ctx.clip();
    for (const s of series) {
      ctx.strokeStyle = s.color;
      ctx.lineWidth = s.dash ? 1.4 : 1.8;
      ctx.setLineDash(s.dash ?? []);
      ctx.beginPath();
      if (s.steps && !s.fn) {
        const { edges, vals } = s.steps;
        // только видимые ступени
        let k = 0;
        while (k < edges.length && edges[k] < x0) k++;
        ctx.moveTo(X(x0), Y(vals[k]));
        for (; k < edges.length && edges[k] <= x1; k++) {
          ctx.lineTo(X(edges[k]), Y(vals[k]));
          ctx.lineTo(X(edges[k]), Y(vals[k + 1]));
        }
        ctx.lineTo(X(x1), Y(vals[k]));
      } else if (s.fn) {
        for (let i = 0; i <= pw; i++) {
          const x = x0 + ((x1 - x0) * i) / pw;
          const v = s.fn(x);
          if (i === 0) ctx.moveTo(X(x), Y(v));
          else ctx.lineTo(X(x), Y(v));
        }
      }
      ctx.stroke();
    }
    ctx.setLineDash([]);
    ctx.restore();
    // легенда
    if (series.length > 1) {
      ctx.font = '11px system-ui, sans-serif';
      let lx = m.l + 8;
      for (const s of series) {
        ctx.strokeStyle = s.color;
        ctx.setLineDash(s.dash ?? []);
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(lx, m.t + 12);
        ctx.lineTo(lx + 20, m.t + 12);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = css('--muted');
        ctx.fillText(s.label, lx + 24, m.t + 16);
        lx += 34 + ctx.measureText(s.label).width;
      }
    }
    // значения под курсором
    if (hover !== null && hover >= m.l && hover <= w - m.r) {
      const x = x0 + ((hover - m.l) / pw) * (x1 - x0);
      ctx.strokeStyle = css('--muted');
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(hover, m.t);
      ctx.lineTo(hover, h - m.b);
      ctx.stroke();
      ctx.setLineDash([]);
      const lines = [`x1/L = ${fmt(x, 5)}`];
      for (const s of series) {
        const v = s.fn ? s.fn(x) : stepAt(s.steps!, x);
        lines.push(`${series.length > 1 ? s.label + ': ' : ''}${v === 0 ? '0' : Math.abs(v) < 1e-3 || Math.abs(v) >= 1e5 ? v.toExponential(5) : v.toPrecision(6)}`);
        ctx.fillStyle = s.color;
        ctx.beginPath();
        ctx.arc(hover, Y(v), 3, 0, 2 * Math.PI);
        ctx.fill();
      }
      ctx.font = '11px system-ui, sans-serif';
      const bw = Math.max(...lines.map((l) => ctx.measureText(l).width)) + 12;
      const bh = 15 * lines.length + 6;
      let bx = hover + 10;
      if (bx + bw > w - m.r) bx = hover - bw - 10;
      const by = h - m.b - bh - 6;
      ctx.fillStyle = css('--panel');
      ctx.globalAlpha = 0.92;
      ctx.fillRect(bx, by, bw, bh);
      ctx.globalAlpha = 1;
      ctx.strokeStyle = css('--border');
      ctx.strokeRect(bx, by, bw, bh);
      ctx.fillStyle = css('--text');
      lines.forEach((l, i) => ctx.fillText(l, bx + 6, by + 15 + 15 * i));
    }
  }

  // ---- управление
  const full = shared.full;
  function setX(a: number, b: number) {
    const [f0, f1] = full();
    const span = Math.max((f1 - f0) / 1e5, Math.min(b - a, f1 - f0));
    a = Math.max(f0, Math.min(a, f1 - span));
    shared.win = span >= f1 - f0 - 1e-12 ? null : [a, a + span];
  }
  function zoomAt(factor: number, px: number | null, py: number | null, vertical: boolean) {
    const g = geom;
    const pw = g.w - m.l - m.r;
    const ph = g.h - m.t - m.b;
    if (vertical) {
      const cy = py === null ? (g.y0 + g.y1) / 2 : g.y1 - ((py - m.t) / ph) * (g.y1 - g.y0);
      yWin = [cy - (cy - g.y0) / factor, cy + (g.y1 - cy) / factor];
      draw();
      return;
    }
    const cx = px === null ? (g.x0 + g.x1) / 2 : g.x0 + ((px - m.l) / pw) * (g.x1 - g.x0);
    setX(cx - (cx - g.x0) / factor, cx + (g.x1 - cx) / factor);
    shared.redraw();
  }
  cv.style.touchAction = 'none';
  cv.style.cursor = 'crosshair';
  cv.addEventListener(
    'wheel',
    (e) => {
      e.preventDefault();
      const r = cv.getBoundingClientRect();
      const d = e.deltaY || e.deltaX;
      zoomAt(Math.exp(-d * 0.0015), e.clientX - r.left, e.clientY - r.top, e.shiftKey);
    },
    { passive: false },
  );
  const pts = new Map<number, { x: number; y: number }>();
  let pinch0 = 0;
  cv.addEventListener('pointerdown', (e) => {
    cv.setPointerCapture(e.pointerId);
    const r = cv.getBoundingClientRect();
    pts.set(e.pointerId, { x: e.clientX - r.left, y: e.clientY - r.top });
    if (pts.size === 2) {
      const [a, b] = [...pts.values()];
      pinch0 = Math.abs(a.x - b.x);
    }
  });
  cv.addEventListener('pointermove', (e) => {
    const r = cv.getBoundingClientRect();
    const cur = { x: e.clientX - r.left, y: e.clientY - r.top };
    const prev = pts.get(e.pointerId);
    if (!prev) {
      hover = cur.x;
      draw();
      return;
    }
    pts.set(e.pointerId, cur);
    if (pts.size === 2) {
      const [a, b] = [...pts.values()];
      const d = Math.abs(a.x - b.x);
      if (pinch0 > 10 && d > 10) zoomAt(d / pinch0, (a.x + b.x) / 2, null, false);
      pinch0 = d;
      return;
    }
    const g = geom;
    const dx = ((cur.x - prev.x) / (g.w - m.l - m.r)) * (g.x1 - g.x0);
    const dy = ((cur.y - prev.y) / (g.h - m.t - m.b)) * (g.y1 - g.y0);
    hover = cur.x;
    if (yWin) yWin = [g.y0 + dy, g.y1 + dy]; // по вертикали сдвигается, только если масштаб задан вручную
    if (shared.win) {
      setX(g.x0 - dx, g.x1 - dx);
      shared.redraw();
    } else draw();
  });
  const end = (e: PointerEvent) => {
    pts.delete(e.pointerId);
    if (pts.size < 2) pinch0 = 0;
  };
  cv.addEventListener('pointerup', end);
  cv.addEventListener('pointercancel', end);
  cv.addEventListener('pointerleave', () => {
    hover = null;
    draw();
  });
  const reset = () => {
    shared.win = null;
    yWin = null;
    shared.redraw();
  };
  cv.addEventListener('dblclick', reset);
  const box = document.createElement('div');
  box.className = 'zoom-ctrl';
  for (const [label, t, fn] of [
    ['+', 'приблизить по x1', () => zoomAt(1.6, null, null, false)],
    ['−', 'отдалить по x1', () => zoomAt(1 / 1.6, null, null, false)],
    ['↕', 'растянуть по вертикали (Shift + колесо)', () => zoomAt(1.6, null, null, true)],
    ['⤢', 'весь вид', reset],
  ] as [string, string, () => void][]) {
    const b = document.createElement('button');
    b.textContent = label;
    b.title = t;
    b.setAttribute('aria-label', t);
    b.addEventListener('click', fn);
    box.appendChild(b);
  }
  wrap.appendChild(box);
  new ResizeObserver(() => draw()).observe(cv);
  return { el: wrap, draw, resetY: () => (yWin = null) };
}
