import { laxCheck } from '../core/acoustic.ts';
import { type Material, c2OverC1, normalize } from '../core/material.ts';
import { type MovingLoadProblem, type MovingLoadSolution, solveMovingLoad, sweep } from '../core/movingLoad.ts';
import { DISSERTATION, MATERIAL_1, MATERIAL_2, PRESETS } from '../core/presets.ts';
import { type Marker, type Series, fmt, lineChart, stepChart } from './charts.ts';
import { FIELD_LABELS, type FieldKey, type SceneOptions, drawScene, fieldValue } from './scene.ts';

type PresetName = keyof typeof PRESETS;

const state = {
  problem: structuredClone(PRESETS['Материал 1 (ν < 0), рис. 21']) as MovingLoadProblem,
  scene: { magnify: 25, showPrestrain: false, field: 'none', delta: false, showGrid: true, phase: 0 } as SceneOptions,
  playing: true,
  speed: 0.25,
  sweepParam: 'M' as 'M' | 's11' | 's21',
  solution: null as MovingLoadSolution | null,
  error: '',
};

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const canvas = $<HTMLCanvasElement>('scene');

// ---------- Панель управления ----------

interface NumField {
  key: string;
  label: string;
  get: () => number;
  set: (v: number) => void;
  min: number;
  max: number;
  step: number;
  slider?: boolean;
  hint?: string;
}

const controls = $('controls');
const inputs = new Map<string, { range?: HTMLInputElement; num: HTMLInputElement }>();

function numRow(f: NumField, parent: HTMLElement) {
  const row = document.createElement('div');
  row.className = f.slider ? 'row' : 'row';
  const lab = document.createElement('label');
  lab.textContent = f.label;
  lab.htmlFor = `in-${f.key}`;
  if (f.hint) lab.title = f.hint;
  row.appendChild(lab);
  let range: HTMLInputElement | undefined;
  if (f.slider) {
    range = document.createElement('input');
    range.type = 'range';
    Object.assign(range, { min: f.min, max: f.max, step: f.step });
    range.value = String(f.get());
    range.setAttribute('aria-label', f.label);
    row.appendChild(range);
  } else row.appendChild(document.createElement('span'));
  const num = document.createElement('input');
  num.type = 'number';
  num.id = `in-${f.key}`;
  Object.assign(num, { step: f.step });
  num.value = String(f.get());
  row.appendChild(num);
  const apply = (v: number) => {
    if (!Number.isFinite(v)) return;
    f.set(v);
    num.value = String(+v.toPrecision(6));
    if (range) range.value = String(v);
    recompute();
  };
  range?.addEventListener('input', () => apply(+range!.value));
  num.addEventListener('change', () => apply(+num.value));
  parent.appendChild(row);
  inputs.set(f.key, { range, num });
}

function fieldset(title: string) {
  const fs = document.createElement('fieldset');
  const lg = document.createElement('legend');
  lg.textContent = title;
  fs.appendChild(lg);
  controls.appendChild(fs);
  return fs;
}

function buildControls() {
  const fsPreset = fieldset('Расчётный случай');
  const sel = document.createElement('select');
  sel.style.width = '100%';
  for (const name of Object.keys(PRESETS)) sel.add(new Option(name, name));
  sel.add(new Option('Свои параметры', 'custom'));
  sel.id = 'preset';
  sel.addEventListener('change', () => {
    if (sel.value === 'custom') return;
    state.problem = structuredClone(PRESETS[sel.value as PresetName]);
    syncInputs();
    recompute();
  });
  fsPreset.appendChild(sel);

  const p = () => state.problem;
  const mat = (k: keyof Material) => ({
    get: () => p().material[k],
    set: (v: number) => {
      p().material[k] = v;
      sel.value = 'custom';
    },
  });
  const fsMat = fieldset('Упругие модули, 10³ МПа');
  for (const [k, label] of [
    ['lambda', 'λ'],
    ['mu', 'μ'],
    ['nu', 'ν'],
    ['alpha', 'α'],
    ['beta', 'β'],
  ] as const)
    numRow({ key: k, label, ...mat(k), min: -100, max: 100, step: 0.1 }, fsMat);

  const fsPre = fieldset('Предварительная деформация (3.15)');
  numRow({ key: 'a0', label: 'a0 = u1,1', get: () => p().a0, set: (v) => ((p().a0 = v), (sel.value = 'custom')), min: -0.02, max: -0.0005, step: 0.0005, slider: true, hint: 'a0 < 0: среда сжата' }, fsPre);
  numRow({ key: 'k0', label: 'k0 = u2,1', get: () => p().k0, set: (v) => ((p().k0 = v), (sel.value = 'custom')), min: 0, max: 0.01, step: 0.0005, slider: true, hint: 'предварительный сдвиг' }, fsPre);

  const fsLoad = fieldset('Нагрузка');
  numRow({ key: 'M', label: 'M = s/c1', get: () => p().M, set: (v) => ((p().M = v), (sel.value = 'custom')), min: 1.01, max: 3, step: 0.01, slider: true, hint: 'число Маха, M > 1' }, fsLoad);
  numRow({ key: 's11', label: 'σ̃11*', get: () => p().s11, set: (v) => ((p().s11 = v), (sel.value = 'custom')), min: -0.07, max: -0.001, step: 0.0005, slider: true, hint: 'нормальная нагрузка /(λ+2μ), < 0' }, fsLoad);
  numRow({ key: 's21', label: 'σ̃21*', get: () => p().s21, set: (v) => ((p().s21 = v), (sel.value = 'custom')), min: -0.01, max: 0.01, step: 0.0005, slider: true, hint: 'касательная нагрузка /(λ+2μ)' }, fsLoad);

  const fsView = fieldset('Отображение');
  const fsel = document.createElement('select');
  fsel.style.width = '100%';
  fsel.id = 'field';
  for (const [k, v] of Object.entries(FIELD_LABELS)) fsel.add(new Option(v, k));
  fsel.addEventListener('change', () => {
    state.scene.field = fsel.value as FieldKey;
    draw();
  });
  const fl = document.createElement('div');
  fl.className = 'row wide';
  fl.innerHTML = '<label for="field">Цвет</label>';
  fl.appendChild(fsel);
  fsView.appendChild(fl);
  numRow({ key: 'mag', label: 'Увеличение', get: () => state.scene.magnify, set: (v) => (state.scene.magnify = v), min: 1, max: 100, step: 1, slider: true, hint: 'перемещения умножаются на этот множитель' }, fsView);
  numRow({ key: 'speed', label: 'Скорость', get: () => state.speed, set: (v) => (state.speed = v), min: 0, max: 1, step: 0.05, slider: true, hint: 'скорость анимации' }, fsView);
  const check = (label: string, get: () => boolean, set: (v: boolean) => void, id: string) => {
    const l = document.createElement('label');
    l.className = 'check';
    const c = document.createElement('input');
    c.type = 'checkbox';
    c.checked = get();
    c.id = id;
    c.addEventListener('change', () => {
      set(c.checked);
      draw();
    });
    l.append(c, label);
    fsView.appendChild(l);
  };
  check('приращение относительно зоны 0', () => state.scene.delta, (v) => (state.scene.delta = v), 'chk-delta');
  check('показывать предварительную деформацию', () => state.scene.showPrestrain, (v) => (state.scene.showPrestrain = v), 'chk-pre');
  check('материальная сетка', () => state.scene.showGrid, (v) => (state.scene.showGrid = v), 'chk-grid');
  const btn = document.createElement('button');
  btn.id = 'play';
  btn.className = 'primary';
  btn.textContent = 'Пауза';
  btn.addEventListener('click', () => {
    state.playing = !state.playing;
    btn.textContent = state.playing ? 'Пауза' : 'Пуск';
    if (state.playing) requestAnimationFrame(loop);
  });
  fsView.appendChild(btn);
  const hint = document.createElement('p');
  hint.className = 'hint';
  hint.textContent =
    'Система координат движется вместе с передней точкой нагрузки L0, поэтому картина стационарна, а материальные линии протекают через фронты. Масштаб длины произволен: решение автомодельно.';
  fsView.appendChild(hint);
}

function syncInputs() {
  const p = state.problem;
  const vals: Record<string, number> = { ...p.material, a0: p.a0, k0: p.k0, M: p.M, s11: p.s11, s21: p.s21 };
  for (const [k, v] of Object.entries(vals)) {
    const i = inputs.get(k);
    if (!i) continue;
    i.num.value = String(v);
    if (i.range) i.range.value = String(v);
  }
}

// ---------- Расчёт ----------

let lastGuess: number[] | undefined;

function recompute() {
  state.error = '';
  let sol: MovingLoadSolution | null = null;
  try {
    sol = solveMovingLoad(state.problem, lastGuess);
  } catch (e) {
    state.error = (e as Error).message;
  }
  if (sol) {
    const z = sol.zones;
    lastGuess = [z[1].grad.a, z[1].grad.k, z[2].grad.a, z[2].grad.k, sol.fronts[0].xi, sol.fronts[1].xi];
  } else if (!state.error) state.error = 'При этих параметрах автомодельное решение с двумя ударными волнами не найдено.';
  state.solution = sol;
  renderStatus();
  renderResults();
  renderSteps();
  scheduleSweep();
  draw();
}

// ---------- Вывод ----------

function renderStatus() {
  const st = $('status');
  st.innerHTML = '';
  const add = (cls: string, text: string) => {
    const d = document.createElement('div');
    d.className = `msg ${cls}`;
    d.textContent = text;
    st.appendChild(d);
  };
  if (state.error) add('err', state.error);
  for (const w of state.solution?.warnings ?? []) add('warn', w);
}

function table(head: string[], rows: (string | { v: string; bad?: boolean })[][]) {
  const wrap = document.createElement('div');
  wrap.className = 'table-wrap';
  const t = document.createElement('table');
  const tr = t.createTHead().insertRow();
  for (const h of head) {
    const th = document.createElement('th');
    th.textContent = h;
    tr.appendChild(th);
  }
  const tb = t.createTBody();
  for (const r of rows) {
    const row = tb.insertRow();
    for (const c of r) {
      const td = row.insertCell();
      if (typeof c === 'string') td.textContent = c;
      else {
        td.textContent = c.v;
        if (c.bad) td.className = 'bad';
      }
    }
  }
  wrap.appendChild(t);
  return wrap;
}

function renderResults() {
  const box = $('results');
  box.innerHTML = '';
  const s = state.solution;
  if (!s) return;
  box.className = 'tables';
  const m = normalize(s.problem.material);
  box.appendChild(
    table(
      ['Фронт', 'θ, °', 'G / c1', 'ξ = tg θ', 'Диссипация D', 'Лакс: c впереди < G < c позади'],
      s.fronts.map((f, i) => {
        const lax = laxCheck(m, s.zones[i].grad, s.zones[i + 1].grad, Math.atan(f.xi), f.G, i === 0 ? 1 : 0);
        return [
          `${f.name} (${i === 0 ? 'квазипродольная' : 'квазипоперечная'})`,
          f.thetaDeg.toFixed(3),
          f.G.toFixed(4),
          f.xi.toFixed(4),
          { v: f.dissipation.toExponential(2), bad: f.dissipation < 0 },
          { v: `${lax.ahead.toFixed(4)} · ${f.G.toFixed(4)} · ${lax.behind.toFixed(4)}`, bad: !lax.ok },
        ];
      }),
    ),
  );
  box.appendChild(
    table(
      ['Зона', 'a = u1,1', 'b = u1,2', 'k = u2,1', 'n = u2,2', 'σ̃11', 'σ̃12', 'σ̃22', 'σ̃33', 'v1 / s', 'v2 / s'],
      s.zones.map((z) => [
        z.name,
        ...[z.grad.a, z.grad.b, z.grad.k, z.grad.n, z.stress.s11, z.stress.s12, z.stress.s22, z.stress.s33, z.v[0], z.v[1]].map((v) => fmt(v, 4)),
      ]),
    ),
  );
}

function renderSteps() {
  const box = $('steps');
  box.innerHTML = '';
  const s = state.solution;
  if (!s) return;
  const color = getComputedStyle(document.documentElement).getPropertyValue('--accent').trim();
  const z = s.zones;
  const items: [string, (i: 0 | 1 | 2) => number][] = [
    ['a = u1,1', (i) => z[i].grad.a],
    ['b = u1,2', (i) => z[i].grad.b],
    ['k = u2,1', (i) => z[i].grad.k],
    ['n = u2,2', (i) => z[i].grad.n],
    ['σ11 / (λ+2μ)', (i) => z[i].stress.s11],
    ['σ21 / (λ+2μ)', (i) => z[i].stress.s12],
    ['σ22 / (λ+2μ)', (i) => z[i].stress.s22],
  ];
  for (const [title, f] of items) box.appendChild(stepChart(title, [f(0), f(1), f(2)], s.fronts[0].xi, s.fronts[1].xi, color));
}

// ---------- θ(параметр) ----------

let sweepTimer = 0;
function scheduleSweep() {
  clearTimeout(sweepTimer);
  sweepTimer = window.setTimeout(renderSweep, 120);
}

function grid(from: number, to: number, step: number) {
  const out: number[] = [];
  const n = Math.round((to - from) / step);
  for (let i = 0; i <= n; i++) out.push(+(from + i * step).toFixed(6));
  return out;
}

/** Продолжение по параметру в обе стороны от текущего значения. */
function sweepBoth(values: number[], current: number, apply: (p: MovingLoadProblem, v: number) => MovingLoadProblem) {
  const lower = values.filter((v) => v <= current).sort((a, b) => b - a);
  const upper = values.filter((v) => v > current).sort((a, b) => a - b);
  const base = state.problem;
  const a = sweep(base, [current, ...lower], apply).points.slice(1);
  const b = sweep(base, [current, ...upper], apply).points.slice(1);
  return [...a.reverse(), ...b];
}

function buildSweepControls() {
  const box = $('sweep-controls');
  box.innerHTML = '<label for="sweep-param">Параметр по оси X:</label>';
  const sel = document.createElement('select');
  sel.id = 'sweep-param';
  sel.add(new Option('число Маха M', 'M'));
  sel.add(new Option('нормальная нагрузка |σ̃11*|', 's11'));
  sel.add(new Option('касательная нагрузка σ̃21*', 's21'));
  sel.addEventListener('change', () => {
    state.sweepParam = sel.value as typeof state.sweepParam;
    renderSweep();
  });
  box.appendChild(sel);
  const note = document.createElement('span');
  note.className = 'hint';
  note.textContent = 'Пунктир из точек — участок, где ударная волна термодинамически недопустима (D < 0).';
  box.appendChild(note);
}

function renderSweep() {
  const box = $('sweep');
  box.innerHTML = '';
  const p = state.problem;
  if (!state.solution) return;
  const cs = getComputedStyle(document.documentElement);
  const c1 = cs.getPropertyValue('--front1').trim();
  const c2 = cs.getPropertyValue('--front2').trim();
  const muted = cs.getPropertyValue('--muted').trim();
  let values: number[];
  let current: number;
  let apply: (q: MovingLoadProblem, v: number) => MovingLoadProblem;
  let xLabel: string;
  if (state.sweepParam === 'M') {
    values = grid(1.005, 3, 0.005);
    current = p.M;
    apply = (q, v) => ({ ...q, M: v });
    xLabel = 'M = s / c1';
  } else if (state.sweepParam === 's11') {
    values = grid(0.001, 0.07, 0.001);
    current = -p.s11;
    apply = (q, v) => ({ ...q, s11: -v });
    xLabel = '|σ̃11*|';
  } else {
    values = grid(-0.01, 0.01, 0.00025);
    current = p.s21;
    apply = (q, v) => ({ ...q, s21: v });
    xLabel = 'σ̃21*';
  }
  const pts = sweepBoth(values, current, apply);
  const xs = pts.map((q) => q.value);
  const series: Series[] = [
    { label: 'θ1 (Σ1)', x: xs, y: pts.map((q) => q.solution.fronts[0].thetaDeg), color: c1, ok: pts.map((q) => q.solution.fronts[0].dissipation >= 0) },
    { label: 'θ2 (Σ2)', x: xs, y: pts.map((q) => q.solution.fronts[1].thetaDeg), color: c2, ok: pts.map((q) => q.solution.fronts[1].dissipation >= 0) },
  ];
  if (state.sweepParam === 'M') {
    const r = c2OverC1(p.material);
    const lx = xs.length ? xs : values;
    series.push(
      { label: 'θ1 линейной теории', x: lx, y: lx.map((M) => (Math.asin(1 / M) * 180) / Math.PI), color: muted, dash: '5 4' },
      { label: 'θ2 линейной теории', x: lx, y: lx.map((M) => (Math.asin(Math.min(1, r / M)) * 180) / Math.PI), color: muted, dash: '2 3' },
    );
  }
  const markers: Marker[] = [
    { x: current, y: state.solution.fronts[0].thetaDeg, color: c1, label: 'текущее решение', shape: 'ring' },
    { x: current, y: state.solution.fronts[1].thetaDeg, color: c2, label: 'текущее решение', shape: 'ring' },
  ];
  // Точки из текста диссертации (только для её параметров)
  const same = (m: Material, ref: Material) => (Object.keys(ref) as (keyof Material)[]).every((k) => m[k] === ref[k]);
  const base = p.a0 === -0.005 && p.k0 === 0.001 && p.s11 === -0.01 && p.s21 === 0.003;
  const ref = same(p.material, MATERIAL_1) ? DISSERTATION.material1 : same(p.material, MATERIAL_2) ? DISSERTATION.material2 : null;
  if (ref && base && state.sweepParam === 'M')
    markers.push(
      { x: 1.7, y: ref.theta1, color: c1, label: 'диссертация (текст)', shape: 'dot' },
      { x: 1.7, y: ref.theta2, color: c2, label: 'диссертация (текст)', shape: 'dot' },
    );
  box.appendChild(lineChart({ xLabel, yLabel: 'θ, °', series, markers, vline: current }));
}

// ---------- Анимация ----------

let lastT = 0;
function draw() {
  const r = drawScene(canvas, state.solution, state.scene);
  const leg = $('legend');
  leg.innerHTML = '';
  const s = state.solution;
  if (!s || !r || state.scene.field === 'none') return;
  const title = document.createElement('div');
  title.textContent = (state.scene.delta ? 'Δ ' : '') + FIELD_LABELS[state.scene.field];
  title.style.marginBottom = '4px';
  leg.appendChild(title);
  const vmax = Math.max(Math.abs(r.min), Math.abs(r.max)) || 1;
  for (const [i, name] of [
    [0, '0'],
    [1, 'I'],
    [2, 'II'],
  ] as const) {
    const v = fieldValue(s, i, state.scene.field, state.scene.delta);
    const t = v / vmax;
    const a = 0.12 + 0.68 * Math.abs(t);
    const row = document.createElement('div');
    row.innerHTML = `<span class="sw" style="background:${t >= 0 ? `rgba(214,72,52,${a})` : `rgba(48,104,214,${a})`}"></span>зона ${name}: ${fmt(v, 4)}`;
    leg.appendChild(row);
  }
}

function loop(t: number) {
  const dt = Math.min(0.05, (t - lastT) / 1000 || 0);
  lastT = t;
  if (state.playing) {
    state.scene.phase += dt * state.speed;
    draw();
    requestAnimationFrame(loop);
  }
}

buildControls();
buildSweepControls();
recompute();
new ResizeObserver(() => draw()).observe(canvas);
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => recompute());
requestAnimationFrame(loop);

// Для автотестов интерфейса
(window as unknown as { __app: typeof state }).__app = state;
