/** Простые SVG-графики без внешних библиотек: ступенчатые диаграммы по ξ и кривые θ(параметр). */

const NS = 'http://www.w3.org/2000/svg';

function el<K extends keyof SVGElementTagNameMap>(tag: K, attrs: Record<string, string | number>, parent?: Element) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  parent?.appendChild(e);
  return e;
}

/** Подпись деления с числом знаков, достаточным для шага сетки. */
function tickLabel(v: number, ticks: number[]): string {
  const step = ticks.length > 1 ? Math.abs(ticks[1] - ticks[0]) : Math.abs(v) || 1;
  const a = Math.max(...ticks.map(Math.abs));
  if (a !== 0 && (a < 1e-3 || a >= 1e5)) return v === 0 ? '0' : v.toExponential(Math.max(0, Math.ceil(-Math.log10(step / a))));
  const dec = Math.max(0, Math.ceil(-Math.log10(step) - 1e-9));
  return v.toFixed(dec);
}

function niceTicks(min: number, max: number, count = 4): number[] {
  if (min === max) {
    const d = Math.abs(min) || 1;
    min -= d * 0.5;
    max += d * 0.5;
  }
  const span = max - min;
  const step0 = span / count;
  const mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= step0) ?? 10 * mag;
  const out: number[] = [];
  for (let v = Math.ceil(min / step) * step; v <= max + 1e-12 * span; v += step) out.push(+v.toPrecision(12));
  return out;
}

export function fmt(v: number, digits = 4): string {
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e4) return v.toFixed(Math.max(0, digits - Math.floor(Math.log10(a)) - 1));
  return v.toExponential(2);
}

/**
 * Ступенчатая диаграмма величины по зонам II | I | 0 вдоль ξ, как на рис. 21.
 * values — [зона 0, зона I, зона II].
 */
export function stepChart(title: string, values: [number, number, number], xi1: number, xi2: number, color: string): HTMLElement {
  const W = 260;
  const H = 120;
  const m = { l: 58, r: 10, t: 20, b: 18 };
  const box = document.createElement('div');
  box.className = 'chart';
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': title }, box);
  el('text', { x: 6, y: 13, class: 'title' }, svg).textContent = title;
  const xmax = xi1 * 1.45;
  let ymin = Math.min(...values, 0);
  let ymax = Math.max(...values, 0);
  const pad = (ymax - ymin || Math.abs(ymax) || 1) * 0.15;
  ymin -= pad;
  ymax += pad;
  const X = (x: number) => m.l + (x / xmax) * (W - m.l - m.r);
  const Y = (y: number) => m.t + ((ymax - y) / (ymax - ymin)) * (H - m.t - m.b);
  const ticksY = niceTicks(ymin, ymax, 3);
  for (const t of ticksY) {
    el('line', { x1: m.l - 3, x2: m.l, y1: Y(t), y2: Y(t), class: 'axis' }, svg);
    el('text', { x: m.l - 5, y: Y(t) + 3.5, 'text-anchor': 'end' }, svg).textContent = tickLabel(t, ticksY);
  }
  el('line', { x1: m.l, x2: m.l, y1: m.t, y2: H - m.b, class: 'axis' }, svg);
  el('line', { x1: m.l, x2: W - m.r, y1: Y(0), y2: Y(0), class: 'zero' }, svg);
  for (const [x, lab] of [
    [xi2, 'Σ2'],
    [xi1, 'Σ1'],
  ] as const) {
    el('line', { x1: X(x), x2: X(x), y1: m.t, y2: H - m.b, class: 'zero' }, svg);
    el('text', { x: X(x), y: H - 5, 'text-anchor': 'middle' }, svg).textContent = lab;
  }
  const segs: [number, number, number][] = [
    [0, xi2, values[2]],
    [xi2, xi1, values[1]],
    [xi1, xmax, values[0]],
  ];
  let d = '';
  segs.forEach(([a, b, v], i) => {
    d += `${i === 0 ? 'M' : 'L'}${X(a)},${Y(v)} L${X(b)},${Y(v)} `;
  });
  el('path', { d, fill: 'none', stroke: color, 'stroke-width': 2 }, svg);
  const zoneMid: [number, string, number][] = [
    [xi2 / 2, 'II', values[2]],
    [(xi1 + xi2) / 2, 'I', values[1]],
    [(xi1 + xmax) / 2, '0', values[0]],
  ];
  for (const [x, , v] of zoneMid) {
    const above = Y(v) > (m.t + H - m.b) / 2;
    el('text', { x: X(x), y: Y(v) + (above ? -4 : 12), 'text-anchor': 'middle' }, svg).textContent = fmt(v, 3);
  }
  return box;
}

export interface Series {
  label: string;
  x: number[];
  y: number[];
  color: string;
  dash?: string;
  /** Для каждой точки: допустима ли (недопустимые участки рисуются пунктиром). */
  ok?: boolean[];
}

export interface Marker {
  x: number;
  y: number;
  color: string;
  label: string;
  shape: 'dot' | 'ring';
}

export function lineChart(opts: {
  xLabel: string;
  yLabel: string;
  series: Series[];
  markers?: Marker[];
  vline?: number;
  xRange?: [number, number];
  /** Разрешить отрицательные значения по оси Y (по умолчанию ось начинается не ниже 0). */
  allowNegative?: boolean;
  height?: number;
  width?: number;
}): HTMLElement {
  const W = opts.width ?? 720;
  const H = opts.height ?? 300;
  const m = { l: 66, r: 14, t: 12, b: 38 };
  const box = document.createElement('div');
  box.className = 'chart';
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': opts.yLabel }, box);
  const xs = opts.series.flatMap((s) => s.x);
  const ys = [...opts.series.flatMap((s) => s.y), ...(opts.markers ?? []).map((p) => p.y)];
  if (!xs.length) {
    el('text', { x: W / 2, y: H / 2, 'text-anchor': 'middle' }, svg).textContent = 'Нет решений в этом диапазоне';
    return box;
  }
  const [x0, x1] = opts.xRange ?? [Math.min(...xs), Math.max(...xs)];
  let y0 = Math.min(...ys);
  let y1 = Math.max(...ys);
  if (y1 - y0 < 1e-15 * Math.max(1, Math.abs(y1))) {
    // все значения одинаковы: показать от нуля до значения с запасом
    y0 = Math.min(0, y0);
    y1 = y1 > 0 ? y1 * 1.25 : y1 < 0 ? 0 : 1;
  }
  const pad = (y1 - y0) * 0.08;
  y0 = opts.allowNegative ? y0 - pad : Math.max(0, y0 - pad);
  y1 += pad;
  const X = (x: number) => m.l + ((x - x0) / (x1 - x0 || 1)) * (W - m.l - m.r);
  const Y = (y: number) => m.t + ((y1 - y) / (y1 - y0 || 1)) * (H - m.t - m.b);
  const ticksY = niceTicks(y0, y1, 6);
  for (const t of ticksY) {
    el('line', { x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), class: 'axis', 'stroke-opacity': 0.5 }, svg);
    el('text', { x: m.l - 6, y: Y(t) + 4, 'text-anchor': 'end' }, svg).textContent = tickLabel(t, ticksY);
  }
  const ticksX = niceTicks(x0, x1, 8);
  for (const t of ticksX) {
    el('line', { x1: X(t), x2: X(t), y1: H - m.b, y2: H - m.b + 4, class: 'axis' }, svg);
    el('text', { x: X(t), y: H - m.b + 16, 'text-anchor': 'middle' }, svg).textContent = tickLabel(t, ticksX);
  }
  el('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: 'axis' }, svg);
  el('text', { x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle' }, svg).textContent = opts.xLabel;
  el('text', { x: 12, y: (m.t + H - m.b) / 2, 'text-anchor': 'middle', transform: `rotate(-90 12 ${(m.t + H - m.b) / 2})` }, svg).textContent =
    opts.yLabel;
  if (opts.vline !== undefined)
    el('line', { x1: X(opts.vline), x2: X(opts.vline), y1: m.t, y2: H - m.b, class: 'zero' }, svg);

  for (const s of opts.series) {
    // разбиваем на участки по допустимости
    let start = 0;
    for (let i = 1; i <= s.x.length; i++) {
      const okPrev = s.ok?.[i - 1] ?? true;
      const okCur = i < s.x.length ? (s.ok?.[i] ?? true) : !okPrev;
      if (okCur !== okPrev || i === s.x.length) {
        const pts = [];
        for (let j = start; j < Math.min(i + 1, s.x.length); j++) pts.push(`${X(s.x[j])},${Y(s.y[j])}`);
        if (pts.length > 1)
          el(
            'polyline',
            {
              points: pts.join(' '),
              fill: 'none',
              stroke: s.color,
              'stroke-width': okPrev ? 2 : 1.6,
              'stroke-dasharray': okPrev ? (s.dash ?? 'none') : '1.5 3',
              'stroke-linecap': 'round',
            },
            svg,
          );
        start = i;
      }
    }
  }
  for (const p of opts.markers ?? []) {
    const c = el('circle', { cx: X(p.x), cy: Y(p.y), r: p.shape === 'dot' ? 4 : 5 }, svg);
    if (p.shape === 'dot') c.setAttribute('fill', p.color);
    else {
      c.setAttribute('fill', 'none');
      c.setAttribute('stroke', p.color);
      c.setAttribute('stroke-width', '2');
    }
    el('title', {}, c).textContent = p.label;
  }
  const leg = document.createElement('div');
  leg.className = 'legend-inline';
  for (const s of opts.series) {
    const sp = document.createElement('span');
    sp.style.color = s.color;
    if (s.dash) sp.className = 'dash';
    sp.textContent = s.label;
    leg.appendChild(sp);
  }
  const seen = new Set<string>();
  for (const p of opts.markers ?? []) {
    if (seen.has(p.label)) continue;
    seen.add(p.label);
    const sp = document.createElement('span');
    sp.className = 'dot';
    sp.style.color = p.color;
    sp.textContent = p.label;
    leg.appendChild(sp);
  }
  box.appendChild(leg);
  return box;
}

/**
 * Компактный профиль величины по ξ (как диаграммы рис. 11, 16, 18): кусочно-гладкая функция
 * со скачками на ударных волнах; вертикальные метки — положения волн.
 */
export function profileChart(
  title: string,
  xs: number[],
  ys: number[],
  marks: { x: number; label: string; fan?: [number, number] }[],
  color: string,
): HTMLElement {
  const W = 260;
  const H = 130;
  const m = { l: 60, r: 8, t: 20, b: 22 };
  const box = document.createElement('div');
  box.className = 'chart';
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': title }, box);
  el('text', { x: 6, y: 13, class: 'title' }, svg).textContent = title;
  const xmax = Math.max(...xs);
  let ymin = Math.min(...ys);
  let ymax = Math.max(...ys);
  const span = ymax - ymin || Math.abs(ymax) || 1e-6;
  ymin -= span * 0.15;
  ymax += span * 0.15;
  const X = (x: number) => m.l + (x / xmax) * (W - m.l - m.r);
  const Y = (y: number) => m.t + ((ymax - y) / (ymax - ymin)) * (H - m.t - m.b);
  for (const mk of marks) {
    if (mk.fan) el('rect', { x: X(mk.fan[1]), y: m.t, width: Math.max(1, X(mk.fan[0]) - X(mk.fan[1])), height: H - m.t - m.b, fill: color, 'fill-opacity': 0.12 }, svg);
    el('line', { x1: X(mk.x), x2: X(mk.x), y1: m.t, y2: H - m.b, class: 'zero' }, svg);
  }
  const labelled: number[] = [];
  for (const mk of marks) {
    if (labelled.some((x) => Math.abs(X(x) - X(mk.x)) < 22)) continue;
    labelled.push(mk.x);
    el('text', { x: X(mk.x), y: H - 7, 'text-anchor': 'middle' }, svg).textContent = mk.label;
  }
  const ticksY = niceTicks(ymin, ymax, 3);
  for (const t of ticksY) {
    el('line', { x1: m.l - 3, x2: m.l, y1: Y(t), y2: Y(t), class: 'axis' }, svg);
    el('text', { x: m.l - 5, y: Y(t) + 3.5, 'text-anchor': 'end' }, svg).textContent = tickLabel(t, ticksY);
  }
  el('line', { x1: m.l, x2: m.l, y1: m.t, y2: H - m.b, class: 'axis' }, svg);
  if (ymin < 0 && ymax > 0) el('line', { x1: m.l, x2: W - m.r, y1: Y(0), y2: Y(0), class: 'zero' }, svg);
  el('polyline', { points: xs.map((x, i) => `${X(x)},${Y(ys[i])}`).join(' '), fill: 'none', stroke: color, 'stroke-width': 2, 'stroke-linejoin': 'round' }, svg);
  return box;
}
