/**
 * Страница «Серия ударов»: несколько нормальных ударов (импульсов скорости границы) с интервалами.
 * Решение нестационарное (волны разных ударов догоняют друг друга) — считается численно (fv1d.ts),
 * а до встречи волн сравнивается с точным автомодельным решением главы 2.
 */
import { type FrontsFrame, type FrontsResult, type Segment, stateAt } from '../core/fronts.ts';
import { type Pulse, type SeriesFrame, type SeriesProblem, type SeriesResult, boundaryVelocity, cell } from '../core/fv1d.ts';
import { type HybridFrame, type HybridResult, hybridStateAt } from '../core/hybrid.ts';
import { profileAt, solveImpact } from '../core/impact.ts';
import { type Material, normalize } from '../core/material.ts';
import { type Vec3, traction } from '../core/material1d.ts';
import { MATERIAL_1, MATERIAL_2 } from '../core/presets.ts';
import { fmt, lineChart } from './charts.ts';
import { type PlotSeries, type SharedX, createProfilePlot } from './profilePlot.ts';

type FieldKey = 'p' | 's' | 'v1' | 's11';
const FIELD_LABELS: Record<FieldKey, string> = {
  p: 'p − p0 (продольная деформация)',
  s: 's − s0 (сдвиг)',
  v1: 'v1 / c1 (скорость частиц)',
  s11: 'σ11 / (λ+2μ) − σ11⁰',
};

const MATERIALS: Record<string, { material: Material; F0: Vec3 }> = {
  'Материал 1 (ν < 0), как в §2.3': { material: MATERIAL_1, F0: [-0.00005, 0, 0.001] },
  'Материал 2 (ν > 0), как в §2.3': { material: MATERIAL_2, F0: [-0.005, 0, 0.001] },
  'Материал 2, слабое предсжатие p0 = −0.0005 (растяжение достижимо)': { material: MATERIAL_2, F0: [-0.0005, 0, 0.001] },
  'Материал 2, предварительное растяжение p0 = +0.001': { material: MATERIAL_2, F0: [0.001, 0, 0.001] },
  'Материал 2, ненапряжённая среда p0 = s0 = 0': { material: MATERIAL_2, F0: [0, 0, 0] },
};

const PULSE_PRESETS: Record<string, Pulse[]> = {
  'Три одинаковых удара': [
    { start: 0, duration: 0.08, d: 0.002 },
    { start: 0.25, duration: 0.08, d: 0.002 },
    { start: 0.5, duration: 0.08, d: 0.002 },
  ],
  'Нарастающие удары': [
    { start: 0, duration: 0.08, d: 0.001 },
    { start: 0.2, duration: 0.08, d: 0.002 },
    { start: 0.4, duration: 0.08, d: 0.003 },
  ],
  'Ступеньки скорости (без остановки)': [
    { start: 0, duration: 1, d: 0.001 },
    { start: 0.2, duration: 1, d: 0.001 },
    { start: 0.4, duration: 1, d: 0.001 },
  ],
  'Сжатие → растяжение → сжатие': [
    { start: 0, duration: 0.08, d: 0.002 },
    { start: 0.25, duration: 0.08, d: -0.002 },
    { start: 0.5, duration: 0.08, d: 0.002 },
  ],
  'Рывок наружу (растяжение)': [{ start: 0, duration: 0.15, d: -0.002 }],
  'Удар, затем рывок (сравнить скорости сжатия и растяжения)': [
    { start: 0, duration: 0.1, d: 0.002 },
    { start: 0.3, duration: 0.1, d: -0.002 },
  ],
  'Удар и сразу рывок (знакопеременный импульс)': [
    { start: 0, duration: 0.08, d: 0.002 },
    { start: 0.08, duration: 0.08, d: -0.002 },
  ],
  'Один удар (проверка по §2.3)': [{ start: 0, duration: 1, d: 0.001 }],
};

const matNames = Object.keys(MATERIALS);
const st = {
  matName: matNames[0],
  material: { ...MATERIALS[matNames[0]].material } as Material,
  F0: [...MATERIALS[matNames[0]].F0] as Vec3,
  pulses: PULSE_PRESETS['Три одинаковых удара'].map((p) => ({ ...p })),
  N: 800,
  /** Время показа (для метода фронтов «до затухания» определяется расчётом). */
  tEnd: 0.92,
  /** Время расчёта по умолчанию (конечные объёмы и метод фронтов без «до затухания»). */
  tBase: 0.92,
  /** Длина показанного участка в единицах L. */
  xView: 1,
  /** Метод фронтов: считать, пока фронты ещё могут встретиться. */
  untilQuiet: true,
  /** Для КО в режиме «до затухания»: до какого момента и на какой длине считали. */
  fvMeta: null as FvMeta | null,
  /** Шаг веера: скачок деформации на один луч простой волны (метод фронтов). */
  fanStep: 2.5e-4,
  /** Окно диаграммы x–t при масштабировании (null — весь расчёт). */
  xtWin: null as null | { x0: number; x1: number; t0: number; t1: number },
  /** Точка, выбранная щелчком на диаграмме x–t (с подписью). */
  xtPick: null as null | { x: number; t: number; lines: string[] },
  field: 'p' as FieldKey,
  magnify: 100,
  /** Увеличение поперечных (сдвиговых) перемещений — при нормальном ударе они в сотни раз меньше. */
  magnifyT: 1000,
  /** Усиливать поперечные (сдвиговые) перемещения отдельно; по умолчанию масштаб общий — видно чистое сжатие. */
  boostShear: false,
  /** Пучки простых волн на профилях — непрерывно (линейно между лучами), а не ступеньками. */
  smoothFans: true,
  view: '2d' as '2d' | '3d',
  playing: true,
  speed: 0.12,
  t: 0.6,
  res: null as SeriesResult | null,
  /**
   * 'fronts' — метод волновых фронтов (как в диссертации), 'hybrid' — ударные волны выделяются точно,
   * вееры считаются на сетке (гладкие), 'fv' — конечные объёмы (для сравнения).
   */
  method: 'fronts' as 'fronts' | 'hybrid' | 'fv',
  /** Гибрид: ячеек на L. */
  hybCells: 400,
  hyb: null as HybridResult | null,
  showFronts: true,
  /** Рисовать и слабые фронты (отражения при встречах), с яркостью по логарифму скачка. */
  weakFronts: false,
  fronts: null as FrontsResult | null,
  ms: 0,
  computing: 0,
  error: '',
};

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const css = (name: string) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

// ---------------------------------------------------------------- управление

function fieldset(parent: HTMLElement, title: string) {
  const fs = document.createElement('fieldset');
  const lg = document.createElement('legend');
  lg.textContent = title;
  fs.appendChild(lg);
  parent.appendChild(fs);
  return fs;
}

function numInput(value: number, step: string, onChange: (v: number) => void, label: string) {
  const i = document.createElement('input');
  i.type = 'number';
  i.step = step;
  i.value = String(value);
  i.setAttribute('aria-label', label);
  i.addEventListener('change', () => Number.isFinite(+i.value) && onChange(+i.value));
  return i;
}

function row(parent: HTMLElement, label: string, input: HTMLElement, id?: string) {
  const r = document.createElement('div');
  r.className = 'row wide';
  const l = document.createElement('label');
  l.textContent = label;
  if (id) {
    input.id = id;
    l.htmlFor = id;
  }
  r.append(l, input);
  parent.appendChild(r);
}

let pulsesBox: HTMLElement;
let matInputs: Record<string, HTMLInputElement> = {};
let preInputs: HTMLInputElement[] = [];

function buildControls() {
  const c = $('s-controls');
  const fsM = fieldset(c, 'Материал');
  const sel = document.createElement('select');
  sel.id = 's-material';
  sel.style.width = '100%';
  for (const n of matNames) sel.add(new Option(n, n));
  sel.add(new Option('Свои параметры', 'custom'));
  sel.addEventListener('change', () => {
    if (sel.value === 'custom') return;
    st.matName = sel.value;
    st.material = { ...MATERIALS[sel.value].material };
    st.F0 = [...MATERIALS[sel.value].F0] as Vec3;
    for (const [k, i] of Object.entries(matInputs)) i.value = String(st.material[k as keyof Material]);
    preInputs[0].value = String(st.F0[0]);
    preInputs[1].value = String(st.F0[2]);
    schedule();
  });
  fsM.appendChild(sel);
  for (const [k, l] of [
    ['lambda', 'λ'],
    ['mu', 'μ'],
    ['nu', 'ν'],
    ['alpha', 'α'],
    ['beta', 'β'],
  ] as [keyof Material, string][]) {
    const i = numInput(st.material[k], 'any', (v) => ((st.material[k] = v), (sel.value = 'custom'), schedule()), l);
    matInputs[k] = i;
    row(fsM, `${l}, 10³ МПа`, i, `s-in-${k}`);
  }
  const fsP = fieldset(c, 'Предварительная деформация');
  preInputs = [
    numInput(st.F0[0], '0.00005', (v) => ((st.F0[0] = v), schedule()), 'p0'),
    numInput(st.F0[2], '0.0005', (v) => ((st.F0[2] = v), schedule()), 's0'),
  ];
  row(fsP, 'p0 = u1,1', preInputs[0], 's-in-p0');
  row(fsP, 's0 = u3,1', preInputs[1], 's-in-s0');

  const fsU = fieldset(c, 'Воздействия на границу');
  const pre = document.createElement('select');
  pre.id = 's-pulse-preset';
  pre.style.width = '100%';
  pre.add(new Option('Шаблон…', ''));
  for (const n of Object.keys(PULSE_PRESETS)) pre.add(new Option(n, n));
  pre.addEventListener('change', () => {
    if (!pre.value) return;
    st.pulses = PULSE_PRESETS[pre.value].map((p) => ({ ...p }));
    pre.value = '';
    renderPulses();
    schedule();
  });
  fsU.appendChild(pre);
  const hint = document.createElement('p');
  hint.className = 'hint';
  hint.textContent =
    'С момента t граница движется со скоростью d·c1 в течение Δt, затем стоит. d > 0 — удар (граница вдавливается, сжатие), d < 0 — рывок наружу (растяжение). Время — в единицах L/c1 (L — длина показанного участка).';
  fsU.appendChild(hint);
  pulsesBox = document.createElement('div');
  pulsesBox.className = 'pulses';
  fsU.appendChild(pulsesBox);
  const add = document.createElement('button');
  add.id = 's-add-pulse';
  add.textContent = '+ добавить удар';
  add.addEventListener('click', () => {
    const last = st.pulses.at(-1);
    st.pulses.push({ start: last ? +(last.start + 0.2).toFixed(3) : 0, duration: 0.08, d: 0.002 });
    renderPulses();
    schedule();
  });
  fsU.appendChild(add);
  renderPulses();

  const fsN = fieldset(c, 'Расчёт и отображение');
  const msel = document.createElement('select');
  msel.add(new Option('волновые фронты (как в диссертации)', 'fronts'));
  msel.add(new Option('ударные волны + сетка (гладкие вееры)', 'hybrid'));
  msel.add(new Option('конечные объёмы (для сравнения)', 'fv'));
  msel.value = st.method;
  msel.addEventListener('change', () => {
    st.method = msel.value as typeof st.method;
    schedule();
  });
  row(fsN, 'Метод', msel, 's-method');
  const ql = document.createElement('label');
  ql.className = 'check';
  const qcb = document.createElement('input');
  qcb.type = 'checkbox';
  qcb.id = 's-until-quiet';
  qcb.checked = st.untilQuiet;
  qcb.addEventListener('change', () => {
    st.untilQuiet = qcb.checked;
    schedule();
  });
  ql.append(qcb, 'считать, пока фронты могут встретиться (конечные объёмы — до того же момента)');
  fsN.appendChild(ql);
  const fanSel = document.createElement('select');
  for (const [v, l] of [
    [5e-4, 'грубо — 5·10⁻⁴ на луч'],
    [2.5e-4, 'обычно — 2.5·10⁻⁴ (по умолчанию)'],
    [1e-4, 'мелко — 1·10⁻⁴ (в 2–3 раза дольше)'],
    [5e-5, 'очень мелко — 5·10⁻⁵ (в 5–10 раз дольше)'],
  ] as [number, string][])
    fanSel.add(new Option(l, String(v)));
  fanSel.value = String(st.fanStep);
  fanSel.title = 'Простая волна (веер) в методе фронтов представлена набором слабых фронтов-лучей; чем мельче шаг, тем ближе к непрерывному вееру, но тем больше встреч фронтов';
  fanSel.addEventListener('change', () => {
    st.fanStep = +fanSel.value;
    schedule();
  });
  row(fsN, 'Шаг веера', fanSel, 's-fan-step');
  const fl = document.createElement('label');
  fl.className = 'check';
  const fcb = document.createElement('input');
  fcb.type = 'checkbox';
  fcb.id = 's-show-fronts';
  fcb.checked = st.showFronts;
  fcb.addEventListener('change', () => {
    st.showFronts = fcb.checked;
    drawAll();
  });
  fl.append(fcb, 'фронты на диаграмме x–t');
  fsN.appendChild(fl);
  const wl = document.createElement('label');
  wl.className = 'check';
  const wcb = document.createElement('input');
  wcb.type = 'checkbox';
  wcb.id = 's-weak-fronts';
  wcb.checked = st.weakFronts;
  wcb.addEventListener('change', () => {
    st.weakFronts = wcb.checked;
    drawAll();
  });
  wl.append(wcb, 'показать слабые фронты (отражения при встречах; яркость — по логарифму скачка)');
  fsN.appendChild(wl);
  const sl = document.createElement('label');
  sl.className = 'check';
  const scb = document.createElement('input');
  scb.type = 'checkbox';
  scb.id = 's-smooth-fans';
  scb.checked = st.smoothFans;
  scb.addEventListener('change', () => {
    st.smoothFans = scb.checked;
    drawAll();
  });
  sl.append(scb, 'простые волны (пучки) на профилях — плавно, а не набором ступенек');
  fsN.appendChild(sl);
  const nsel = document.createElement('select');
  for (const n of [400, 800, 1600]) nsel.add(new Option(`${n} ячеек на L${n === 800 ? ' (по умолчанию)' : ''}`, String(n)));
  nsel.value = String(st.N);
  nsel.addEventListener('change', () => {
    st.N = +nsel.value;
    schedule();
  });
  row(fsN, 'Сетка (КО)', nsel, 's-grid');
  const hsel = document.createElement('select');
  for (const n of [200, 400, 800]) hsel.add(new Option(`${n} ячеек на L${n === 400 ? ' (по умолчанию)' : n === 800 ? ' (в 4 раза дольше)' : ''}`, String(n)));
  hsel.value = String(st.hybCells);
  hsel.title = 'Гибридный метод: сетка, на которой считаются вееры между выделенными ударными волнами';
  hsel.addEventListener('change', () => {
    st.hybCells = +hsel.value;
    schedule();
  });
  row(fsN, 'Сетка (гибрид)', hsel, 's-hyb-grid');
  const fsel = document.createElement('select');
  for (const [k, v] of Object.entries(FIELD_LABELS)) fsel.add(new Option(v, k));
  fsel.addEventListener('change', () => {
    st.field = fsel.value as FieldKey;
    if (st.res) computeFieldMax();
    heatCache = null;
    drawAll();
  });
  row(fsN, 'Цвет', fsel, 's-field');
  const mag = numInput(st.magnify, '1', (v) => ((st.magnify = v), drawAll()), 'увеличение');
  row(fsN, 'Увеличение', mag, 's-mag');
  const magT = numInput(st.magnifyT, '1', (v) => ((st.magnifyT = v), drawAll()), 'увеличение сдвига');
  const bl = document.createElement('label');
  bl.className = 'check';
  const bcb = document.createElement('input');
  bcb.type = 'checkbox';
  bcb.id = 's-boost-shear';
  bcb.checked = st.boostShear;
  bl.append(bcb, 'усилить сдвиг u3 отдельно (сдвиг в ~1000 раз слабее сжатия)');
  fsN.appendChild(bl);
  const magTRow = document.createElement('div');
  fsN.appendChild(magTRow);
  row(magTRow, '× сдвиг', magT, 's-magT');
  magTRow.hidden = !st.boostShear;
  bcb.addEventListener('change', () => {
    st.boostShear = bcb.checked;
    magTRow.hidden = !st.boostShear;
    autoMagnify();
    drawAll();
  });
  const tIn = document.createElement('input');
  tIn.type = 'range';
  tIn.min = '0';
  tIn.max = '1';
  tIn.step = '0.001';
  tIn.addEventListener('input', () => {
    st.playing = false;
    st.t = +tIn.value * st.tEnd;
    syncPlay();
    drawAll();
  });
  row(fsN, 'Время', tIn, 's-time');
  const btn = document.createElement('button');
  btn.id = 's-play';
  btn.className = 'primary';
  btn.addEventListener('click', () => {
    st.playing = !st.playing;
    if (st.playing && st.t >= st.tEnd) st.t = 0;
    syncPlay();
    if (st.playing) requestAnimationFrame(loop);
  });
  fsN.appendChild(btn);
  syncPlay();
}

function renderPulses() {
  pulsesBox.innerHTML = '';
  const head = document.createElement('div');
  head.className = 'pulse-row head';
  head.innerHTML = '<span>№</span><span>t начала</span><span>Δt</span><span>d (± )</span><span></span>';
  pulsesBox.appendChild(head);
  st.pulses.forEach((p, i) => {
    const r = document.createElement('div');
    r.className = 'pulse-row';
    const n = document.createElement('span');
    n.textContent = String(i + 1);
    r.append(
      n,
      numInput(p.start, '0.01', (v) => ((p.start = Math.max(0, v)), schedule()), `начало удара ${i + 1}`),
      numInput(p.duration, '0.01', (v) => ((p.duration = Math.max(0.001, v)), schedule()), `длительность удара ${i + 1}`),
      numInput(p.d, '0.0005', (v) => ((p.d = v), schedule()), `скорость удара ${i + 1}`),
    );
    const del = document.createElement('button');
    del.textContent = '×';
    del.title = 'удалить';
    del.setAttribute('aria-label', `удалить удар ${i + 1}`);
    del.addEventListener('click', () => {
      st.pulses.splice(i, 1);
      renderPulses();
      schedule();
    });
    r.appendChild(del);
    pulsesBox.appendChild(r);
  });
}

function syncPlay() {
  const b = document.getElementById('s-play');
  if (b) b.textContent = st.playing ? 'Пауза' : 'Пуск';
  const t = document.getElementById('s-time') as HTMLInputElement | null;
  if (t) t.value = String(st.t / st.tEnd);
}

// ---------------------------------------------------------------- расчёт

let worker: Worker | null = null;
let reqId = 0;
let timer = 0;

function schedule() {
  clearTimeout(timer);
  if (!st.computing) busySince = performance.now();
  st.computing = 0.001; // расчёт запланирован
  updateBusy();
  timer = window.setTimeout(compute, 250);
}

/**
 * Индикация пересчёта: закреплённая плашка над графиками (метод, процент, если решатель его
 * сообщает, и сколько секунд идёт расчёт), а прежние графики приглушаются — видно, что они устарели.
 */
let busySince = 0;
let busyTimer = 0;
const METHOD_NAMES = { fronts: 'метод фронтов', hybrid: 'гибрид (ударные волны + сетка)', fv: 'конечные объёмы' } as const;
function updateBusy() {
  const box = document.getElementById('s-busy');
  const sec = document.querySelector('#page-series section');
  if (!box || !sec) return;
  const on = st.computing > 0;
  box.hidden = !on;
  sec.classList.toggle('busy', on && !!st.res);
  if (!on) {
    clearInterval(busyTimer);
    busyTimer = 0;
    return;
  }
  const el = ((performance.now() - busySince) / 1000).toFixed(1);
  const pct = st.computing > 0.002 ? ` ${Math.round(st.computing * 100)} %` : '';
  $('s-busy-text').textContent = `Пересчёт: ${METHOD_NAMES[st.method]}…${pct} · ${el} с`;
  if (!busyTimer) busyTimer = window.setInterval(updateBusy, 200);
}

function problem(): SeriesProblem {
  return { material: normalize(st.material), F0: [...st.F0] as Vec3, pulses: st.pulses.map((p) => ({ ...p })), N: st.N, tEnd: st.tBase, frames: 240 };
}

async function compute() {
  st.error = '';
  const id = ++reqId;
  if (!st.computing) busySince = performance.now();
  st.computing = 0.001;
  renderStatus();
  const prob = problem();
  const done = (res: SeriesResult, ms: number, fronts: FrontsResult | null = null, meta: FvMeta | null = null, hyb: HybridResult | null = null) => {
    if (id !== reqId) return;
    st.res = res;
    st.fronts = fronts;
    st.hyb = hyb;
    st.fvMeta = meta;
    st.tEnd = fronts ? fronts.tEnd : hyb ? hyb.tEnd : meta ? meta.tEnd : st.tBase;
    st.xView = fronts ? fronts.xView : hyb ? hyb.xView : meta ? meta.length : 1;
    st.xtWin = null;
    profShared.win = null;
    st.xtPick = null;
    st.t = Math.min(st.t, st.tEnd);
    st.ms = ms;
    st.computing = 0;
    autoMagnify();
    renderStatus();
    renderStatic();
    drawAll();
  };
  try {
    worker?.terminate();
    worker = new Worker(new URL('./seriesWorker.ts', import.meta.url), { type: 'module' });
    worker.onmessage = (e) => {
      const msg = e.data;
      if (msg.id !== reqId) return;
      if (msg.type === 'progress') {
        st.computing = msg.f;
        renderStatus();
      } else if (msg.type === 'done') {
        if (msg.method === 'fronts') done(framesFromFronts(msg.res), msg.ms, msg.res);
        else if (msg.method === 'hybrid') done(framesFromHybrid(msg.res), msg.ms, null, null, msg.res);
        else done(msg.res, msg.ms, null, msg.meta);
      }
      else {
        st.error = msg.message;
        st.computing = 0;
        renderStatus();
      }
    };
    if (st.method === 'fronts') worker.postMessage({ id, method: 'fronts', prob: { material: prob.material, F0: prob.F0, pulses: prob.pulses, tEnd: prob.tEnd, frames: FRONT_FRAMES, untilQuiet: st.untilQuiet, fanStep: st.fanStep } });
    else if (st.method === 'hybrid') worker.postMessage({ id, method: 'hybrid', prob: hybridProblem(prob) });
    else worker.postMessage({ id, method: 'fv', prob, match: fvMatch() });
  } catch {
    // запасной путь — в основном потоке
    const t0 = performance.now();
    if (st.method === 'fronts') {
      const { trackFronts } = await import('../core/fronts.ts');
      const fr = trackFronts({ material: prob.material, F0: prob.F0, pulses: prob.pulses, tEnd: prob.tEnd, frames: FRONT_FRAMES, untilQuiet: st.untilQuiet, fanStep: st.fanStep });
      done(framesFromFronts(fr), performance.now() - t0, fr);
    } else if (st.method === 'hybrid') {
      const { solveHybrid } = await import('../core/hybrid.ts');
      const hr = solveHybrid(hybridProblem(prob));
      done(framesFromHybrid(hr), performance.now() - t0, null, null, hr);
    } else {
      const match = fvMatch();
      if (match) {
        const { simulateMatched } = await import('../core/fvMatched.ts');
        const m = simulateMatched(prob, match);
        done(m.res, performance.now() - t0, null, { tEnd: m.tEnd, length: m.length, cellsPerL: m.cellsPerL, capped: m.capped });
      } else {
        const { simulate } = await import('../core/fv1d.ts');
        done(simulate(prob), performance.now() - t0);
      }
    }
  }
}

interface FvMeta {
  tEnd: number;
  length: number;
  cellsPerL: number;
  capped: boolean;
}

/** Предел общего числа ячеек КО на длинном отрезке (время растёт как квадрат числа ячеек). */
const FV_MAX_CELLS = 3200;

/** КО в режиме «до затухания» считаются до того же момента и на той же длине, что и метод фронтов. */
function fvMatch() {
  return st.untilQuiet ? { fanStep: st.fanStep, frames: FRONT_FRAMES, maxCells: FV_MAX_CELLS } : null;
}

/** Сколько моментов сохраняет метод фронтов (по ним строится точный фон диаграммы x–t). */
const FRONT_FRAMES = 1201;
/** Гибрид хранит кусочно-линейные профили (~ число ячеек узлов на кадр) — моментов меньше. */
const HYB_FRAMES = 601;
/** Предел общего числа ячеек гибрида на длинном отрезке (расчёт «до затухания»). */
const HYB_MAX_CELLS = 2400;

function hybridProblem(prob: SeriesProblem) {
  return { material: prob.material, F0: prob.F0, pulses: prob.pulses, tEnd: prob.tEnd, frames: HYB_FRAMES, cellsPerL: st.hybCells, maxCells: HYB_MAX_CELLS, untilQuiet: st.untilQuiet, fanStep: st.fanStep };
}

/** Кадр гибрида, ближайший к моменту t. */
function hybFrameAt(t: number): HybridFrame {
  const fr = st.hyb!.frames;
  const k = Math.max(0, Math.min(fr.length - 1, Math.round((t / st.hyb!.tEnd) * (fr.length - 1))));
  return fr[k];
}

/** Значения (p − p0, s − s0, v1, v3) кусочно-линейного профиля кадра в точке x → out. */
function sampleW(f: HybridFrame, x: number, out: Float64Array, o: number) {
  const xs = f.xs;
  const W = f.W;
  const n = xs.length;
  let lo = 1;
  let hi = n - 1;
  if (x <= xs[0]) {
    for (let k = 0; k < 4; k++) out[o + k] = W[k];
    return;
  }
  if (x >= xs[n - 1]) {
    for (let k = 0; k < 4; k++) out[o + k] = W[4 * (n - 1) + k];
    return;
  }
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  const xa = xs[lo - 1];
  const xb = xs[lo];
  const u = xb > xa ? (x - xa) / (xb - xa) : 1;
  const a = 4 * (lo - 1);
  for (let k = 0; k < 4; k++) out[o + k] = W[a + k] + u * (W[a + 4 + k] - W[a + k]);
}

/** Скорости переноса для интерполяции по времени: покой и характеристики ±c_f, ±c_s в исходном состоянии. */
let hybInterpInfo: { id: number; speeds: number[]; scale: number[] } | null = null;
function hybInfo() {
  if (hybInterpInfo?.id === reqId) return hybInterpInfo;
  const m = normalize(st.material);
  const [p, , s] = st.F0;
  const T = (pp: number, ss: number) => {
    const t = traction(m, [pp, 0, ss]);
    return [t[0], t[2]];
  };
  const h = 1e-7 * Math.max(Math.abs(p), 1e-4);
  const a = (T(p + h, s)[0] - T(p - h, s)[0]) / (2 * h);
  const b = (T(p, s + h)[0] - T(p, s - h)[0]) / (2 * h);
  const d = (T(p, s + h)[1] - T(p, s - h)[1]) / (2 * h);
  const tr = (a + d) / 2;
  const disc = Math.sqrt(((a - d) / 2) ** 2 + b * b);
  const cf = Math.sqrt(Math.max(tr + disc, 0));
  const cs = Math.sqrt(Math.max(tr - disc, 0));
  const scale = [1e-30, 1e-30, 1e-30, 1e-30];
  for (const f of st.hyb!.frames) for (let i = 0; i < f.W.length; i++) scale[i & 3] = Math.max(scale[i & 3], Math.abs(f.W[i]));
  hybInterpInfo = { id: reqId, speeds: [0, cf, -cf, cs, -cs], scale };
  return hybInterpInfo;
}

const iwA = new Float64Array(4);
const iwB = new Float64Array(4);
/**
 * Решение гибрида в момент t между сохранёнными моментами t0 < t < t1. Простое смешивание кадров
 * давало бы двойное изображение, а ближайший кадр — ступеньки на диаграмме t–x. Волны переносятся
 * со скоростями характеристик, поэтому значение в точке x берётся из кадра t0 в точке x − c(t − t0)
 * и из кадра t1 в точке x + c(t1 − t), где c — та из скоростей {0, ±c_f, ±c_s}, при которой эти
 * два значения ближе всего (волна, проходящая через точку), и смешивается линейно по времени.
 */
function hybInterp(t: number) {
  const fr = st.hyb!.frames;
  const n = fr.length;
  const pos = Math.max(0, Math.min(n - 1, (t / st.hyb!.tEnd) * (n - 1)));
  const k = Math.min(n - 2, Math.floor(pos));
  const f0 = fr[k];
  const f1 = fr[k + 1];
  const tt = Math.max(f0.t, Math.min(f1.t, t));
  const w = f1.t > f0.t ? (tt - f0.t) / (f1.t - f0.t) : 0;
  const tau0 = tt - f0.t;
  const tau1 = f1.t - tt;
  const { speeds, scale } = hybInfo();
  return {
    f0,
    f1,
    w,
    at(x: number, out: Float64Array) {
      if (w < 1e-9 || w > 1 - 1e-9) return sampleW(w < 0.5 ? f0 : f1, x, out, 0);
      let best = Infinity;
      for (const c of speeds) {
        sampleW(f0, x - c * tau0, iwA, 0);
        sampleW(f1, x + c * tau1, iwB, 0);
        let mis = 0;
        for (let q = 0; q < 4; q++) mis += Math.abs(iwA[q] - iwB[q]) / scale[q];
        if (mis < best - 1e-12) {
          best = mis;
          for (let q = 0; q < 4; q++) out[q] = iwA[q] + w * (iwB[q] - iwA[q]);
        }
        // точка, где ничего не движется (или почти): дальше не ищем
        if (mis < 1e-3) break;
      }
    },
  };
}

/** Состояние гибрида в момент t как функция x. */
function hybSampler(t: number) {
  const I = hybInterp(t);
  const out = new Float64Array(4);
  return (x: number) => {
    I.at(x, out);
    return { F: [st.F0[0] + out[0], 0, st.F0[2] + out[1]] as Vec3, v: [out[2], 0, out[3]] as Vec3 };
  };
}

/** Профили гибрида в кадры на равномерной сетке (для общих графиков: автомасштаб, граница и т. п.). */
function framesFromHybrid(hr: HybridResult): SeriesResult {
  const N = 1000;
  const dx = hr.xView / N;
  const frames = hr.frames
    .filter((_, k) => k % 3 === 0 || k === hr.frames.length - 1)
    .map((f) => {
      const U = new Float32Array(6 * N);
      fillHybrid(f, U, N, dx);
      return { t: f.t, U, ub: f.ub };
    });
  return { N, dx, frames, kinematicError: hr.kinematicError, work: hr.work, energyGain: hr.energyGain, steps: hr.steps };
}

/** Значения кусочно-линейного профиля в центрах ячеек (один проход по узлам). */
function fillHybrid(f: HybridFrame, U: Float32Array, N: number, dx: number) {
  const xs = f.xs;
  const W = f.W;
  const n = xs.length;
  let j = 1;
  for (let i = 0; i < N; i++) {
    const x = (i + 0.5) * dx;
    while (j < n - 1 && xs[j] < x) j++;
    const xa = xs[j - 1];
    const xb = xs[j];
    const u = xb > xa ? Math.max(0, Math.min(1, (x - xa) / (xb - xa))) : 1;
    const o = 4 * (j - 1);
    const val = (k: number) => W[o + k] + u * (W[o + 4 + k] - W[o + k]);
    U[i] = st.F0[0] + val(0);
    U[N + i] = 0;
    U[2 * N + i] = st.F0[2] + val(1);
    U[3 * N + i] = val(2);
    U[4 * N + i] = 0;
    U[5 * N + i] = val(3);
  }
}

/** Точные перемещения в узлах x_i = i·dx для кусочно-линейного профиля гибрида (интеграл трапециями по узлам). */
function hybridDisplacements(f: HybridFrame, N: number, dx: number) {
  const xs = f.xs;
  const W = f.W;
  const n = xs.length;
  const P1 = new Float64Array(N + 1);
  const P3 = new Float64Array(N + 1);
  // интеграл от узла 0 до x
  let j = 0;
  let a1 = 0;
  let a3 = 0;
  const lin = (k: number, x: number, c: number) => {
    const xa = xs[k];
    const xb = xs[k + 1];
    const u = xb > xa ? (x - xa) / (xb - xa) : 0;
    return W[4 * k + c] + u * (W[4 * k + 4 + c] - W[4 * k + c]);
  };
  for (let i = 0; i <= N; i++) {
    const x = Math.min(i * dx, xs[n - 1]);
    // целые участки до x
    while (j < n - 1 && xs[j + 1] <= x) {
      const len = xs[j + 1] - xs[j];
      a1 += 0.5 * len * (W[4 * j] + W[4 * j + 4]);
      a3 += 0.5 * len * (W[4 * j + 1] + W[4 * j + 5]);
      j++;
    }
    let p1 = a1;
    let p3 = a3;
    if (j < n - 1 && x > xs[j]) {
      const len = x - xs[j];
      p1 += 0.5 * len * (W[4 * j] + lin(j, x, 0));
      p3 += 0.5 * len * (W[4 * j + 1] + lin(j, x, 1));
    }
    P1[i] = p1;
    P3[i] = p3;
  }
  const du1 = new Float64Array(N + 1);
  const du3 = new Float64Array(N + 1);
  for (let i = 0; i <= N; i++) {
    du1[i] = P1[i] - P1[N];
    du3[i] = P3[i] - P3[N];
  }
  return { du1, du3 };
}

/** Перевести кусочно-постоянное решение метода фронтов в кадры на сетке (для общих графиков). */
function framesFromFronts(fr: FrontsResult): SeriesResult {
  const N = 1000;
  const dx = fr.xView / N;
  // для сечения, профилей и 3D хватает каждого 5-го момента; точная диаграмма x–t строится по всем
  const frames = fr.frames.filter((_, k) => k % 5 === 0 || k === fr.frames.length - 1).map((f) => {
    const U = new Float32Array(6 * N);
    for (let i = 0; i < N; i++) {
      const s = stateAt(f, (i + 0.5) * dx);
      U[i] = s.F[0];
      U[N + i] = 0;
      U[2 * N + i] = s.F[1];
      U[3 * N + i] = s.v[0];
      U[4 * N + i] = 0;
      U[5 * N + i] = s.v[1];
    }
    return { t: f.t, U, ub: f.ub };
  });
  return { N, dx, frames, kinematicError: fr.kinematicError, work: fr.work, energyGain: fr.energyGain, steps: fr.events };
}

function autoMagnify() {
  if (!st.res) return;
  const r = st.res;
  let mx = 1e-9;
  for (const fr of r.frames) for (let i = 0; i < r.N; i++) mx = Math.max(mx, Math.abs(fr.U[i] - st.F0[0]));
  st.magnify = Math.max(1, Math.round(0.3 / mx));
  const m = document.getElementById('s-mag') as HTMLInputElement | null;
  if (m) m.value = String(st.magnify);
  // сдвиг: максимальное поперечное перемещение ×Kt ≈ 0.05 L
  let us = 1e-12;
  for (const fr of r.frames) {
    let acc = 0;
    for (let i = r.N - 1; i >= 0; i--) {
      acc -= (fr.U[2 * r.N + i] - st.F0[2]) * r.dx;
      us = Math.max(us, Math.abs(acc));
    }
  }
  st.magnifyT = st.boostShear ? Math.max(1, Math.min(1e6, Math.round(0.05 / us))) : st.magnify;
  const mt = document.getElementById('s-magT') as HTMLInputElement | null;
  if (mt) mt.value = String(st.magnifyT);
}

// ---------------------------------------------------------------- вывод

function renderStatus() {
  const box = $('s-status');
  box.innerHTML = '';
  const add = (cls: string, text: string) => {
    const d = document.createElement('div');
    d.className = `msg ${cls}`;
    d.textContent = text;
    box.appendChild(d);
  };
  if (st.error) add('err', st.error);
  if (st.computing) add('info', `Расчёт… ${Math.round(st.computing * 100)}%`);
  updateBusy();
  const r = st.res;
  const fr = st.fronts;
  const hy = st.hyb;
  if (r && hy && !st.computing) {
    const numDiss = hy.work - hy.energyGain - hy.shockDissipation;
    add(
      'info',
      `Гибридный метод: ударные волны выделены точно (условия на разрывах (2.40)–(2.42) на каждом шаге), вееры посчитаны на сетке ${Math.round(hy.cellsPerL)} ячеек на L (всего ${hy.N}) как гладкое решение — без лучей. ` +
        `${hy.steps} шагов, ${hy.events} точных задач (удары, встречи и отражения ударных волн), до ${hy.maxFronts} выделенных разрывов одновременно, ${st.ms.toFixed(0)} мс. ` +
        `Кинематический баланс ${hy.kinematicError.toExponential(1)}. Энергия: работа границы ${fmt(hy.work, 4)} = прирост энергии ${fmt(hy.energyGain, 4)} + диссипация на ударных волнах (2.45) ${fmt(hy.shockDissipation, 3)} + численная ${numDiss.toExponential(1)}.`,
    );
    add(
      'info',
      'Узкий веер, который сетка всё равно не разрешила бы (уже численного размазывания, max(h/2, ½·h^⅔·t^⅓)), ведётся одним разрывом — на диаграмме тонкая линия; как только он становится шире, его ведёт сетка.' +
        (hy.untracked ? ` Отдано сетке разрывов: ${hy.untracked}.` : ''),
    );
    if (hy.cellsPerL < st.hybCells - 0.5)
      add('warn', `Отрезок длинный (до x1/L = ${fmt(hy.xView, 3)}), поэтому, чтобы расчёт занял разумное время, сетка огрублена до ${Math.round(hy.cellsPerL)} ячеек на L (вместо ${st.hybCells}). Ударные волны от этого не размываются — только вееры.`);
    if (st.untilQuiet && hy.tQuiet !== undefined) add('info', `Как и у метода фронтов, расчёт идёт до момента, после которого волны больше не встречаются: t c1/L = ${fmt(hy.tEnd, 3)}, на участке до x1/L = ${fmt(hy.xView, 3)}.`);
    for (const w of hy.warnings) add('err', w);
  } else if (r && fr && !st.computing) {
    add(
      'info',
      `Метод волновых фронтов: ${fr.events} событий (удары, остановки, встречи и отражения фронтов), до ${fr.maxFronts} фронтов одновременно, шаг веера ${st.fanStep.toExponential(1)}, ${st.ms.toFixed(0)} мс. ` +
        `На каждом фронте точно выполнены условия на разрывах (2.40)–(2.42); простые волны представлены пучками слабых фронтов. ` +
        `Кинематический баланс ${fr.kinematicError.toExponential(1)}. Энергия: работа границы ${fmt(fr.work, 4)} = прирост энергии ${fmt(fr.energyGain, 4)} + диссипация на ударных фронтах (2.45) ${fmt(fr.dissipation, 3)} (расхождение ${(fr.work - fr.energyGain - fr.dissipation).toExponential(1)}).`,
    );
    if (st.untilQuiet) {
      const far = Number.isFinite(fr.nextInteraction) && fr.nextInteraction < 1e3;
      add(
        'info',
        fr.horizonReached
          ? `Расчёт до горизонта t = ${fmt(fr.tEnd, 3)}: взаимодействия ещё продолжаются.`
          : `Последняя встреча фронтов — при t c1/L ≈ ${fmt(fr.tQuiet ?? fr.tEnd, 3)}; дальше все фронты разбегаются (каждый следующий не медленнее предыдущего), к границе ничего не бежит. Показан участок до x1/L ≈ ${fmt(fr.xView, 3)}.` +
            (far ? ` За горизонтом (t = 10) осталась далёкая встреча почти параллельных фронтов при t ≈ ${fmt(fr.nextInteraction, 3)}.` : ''),
      );
    }
    if (fr.linearPasses) add('warn', `${fr.linearPasses} встреч(а) не решены точно и обработаны линейным наложением волн (ошибка ~ нелинейность × амплитуды).`);
    if (fr.coarsened)
      add(
        'warn',
        `Сильная нелинейность (переход через p = 0) порождает лавину сопоставимых отражений, поэтому расчёт выполнен с огрублёнными порогами: фронты слабее 1e-7 отбрасываются, а при встрече волн, различающихся в 100 раз и более, сильная проходит насквозь (поглощено ${fr.absorbed}). Погрешность баланса — порядка этих порогов.`,
      );
    for (const w of fr.warnings) add('err', w);
  } else if (r && !st.computing) {
    const diss = r.work > 0 ? (r.work - r.energyGain) / r.work : 0;
    add(
      'info',
      `Численное решение: ${r.N} ячеек, ${r.steps} шагов, ${st.ms.toFixed(0)} мс. ` +
        `Кинематический баланс ∫(p − p0)dx + ∫v*dt = ${r.kinematicError.toExponential(1)}. ` +
        `Работа границы ${fmt(r.work, 4)}, прирост энергии среды ${fmt(r.energyGain, 4)}, рассеяно ${(diss * 100).toFixed(2)} % (ударные волны + численная вязкость; убывает с измельчением сетки).`,
    );    const mt = st.fvMeta;
    if (mt)
      add(
        mt.capped ? 'warn' : 'info',
        `Как и метод фронтов, расчёт идёт до момента, после которого фронты больше не встречаются: t c1/L = ${fmt(mt.tEnd, 3)}, на участке до x1/L = ${fmt(mt.length, 3)}.` +
          (mt.capped
            ? ` Чтобы расчёт занял разумное время, сетка огрублена до ${Math.round(mt.cellsPerL)} ячеек на L (вместо ${st.N}; всего ${r.N}) — фронты размазаны сильнее, чем на коротком отрезке.`
            : ` ${Math.round(mt.cellsPerL)} ячеек на L.`),
      );
  }
}

function fieldOf(F: Vec3, v: Vec3): number {
  switch (st.field) {
    case 'p':
      return F[0] - st.F0[0];
    case 's':
      return F[2] - st.F0[2];
    case 'v1':
      return v[0];
    case 's11':
      return traction(normalize(st.material), F)[0] - traction(normalize(st.material), st.F0)[0];
  }
}

/** Ближайший сохранённый кадр (для КО — единственный источник; для фронтов — для редких перерисовок). */
function storedFrameAt(t: number) {
  const r = st.res!;
  const k = Math.max(0, Math.min(r.frames.length - 1, Math.round((t / st.tEnd) * (r.frames.length - 1))));
  return r.frames[k];
}

/**
 * Фронты в момент t: берётся ближайший сохранённый момент и фронты сдвигаются со своими скоростями.
 * Сдвигать можно только до первой встречи фронтов (или удара о границу): дальше рождаются новые
 * волны, а фронты, «прошедшие сквозь» друг друга, поставили бы состояния не на свои места —
 * при большом увеличении смещений это выглядит как прыжок сетки. Поэтому сначала пробуем шаг
 * вперёд от предыдущего момента, затем назад от следующего, а если встреча есть в обе стороны —
 * встретившиеся фронты стоят в точке встречи до следующего момента, остальные движутся точно.
 */
function frontsAt(t: number): { f: FrontsFrame; xs: number[] } {
  const fr = st.fronts!;
  const nF = fr.frames.length;
  const k = Math.max(0, Math.min(nF - 1, Math.floor((t / fr.tEnd) * (nF - 1) + 1e-9)));
  const move = (f: FrontsFrame, dt: number) => f.xs.map((x, i) => x + f.speeds[i] * dt);
  const valid = (xs: number[]) => {
    for (let i = 0; i < xs.length; i++) if (xs[i] < -1e-12 || (i > 0 && xs[i] < xs[i - 1] - 1e-12)) return false;
    return true;
  };
  const f = fr.frames[k];
  const fwd = move(f, t - f.t);
  if (valid(fwd)) return { f, xs: fwd };
  if (k + 1 < nF) {
    const g = fr.frames[k + 1];
    const back = move(g, t - g.t);
    if (valid(back)) return { f: g, xs: back };
  }
  // Встреча есть в обе стороны: сдвигаем вперёд и исправляем порядок только у встретившихся фронтов
  // (они стоят в точке встречи, пока не наступит следующий сохранённый момент). Остальные фронты
  // остаются точными — иначе останавливалась бы вся картина и фон диаграммы шёл бы «лесенкой».
  return { f, xs: monotone(fwd) };
}

/** Ближайшая неубывающая последовательность (объединение соседних нарушителей), x ≥ 0. */
function monotone(xs: number[]): number[] {
  const val: number[] = [];
  const cnt: number[] = [];
  for (const x0 of xs) {
    let v = Math.max(0, x0);
    let c = 1;
    while (val.length && val[val.length - 1] > v) {
      const pv = val.pop()!;
      const pc = cnt.pop()!;
      v = (pv * pc + v * c) / (pc + c);
      c += pc;
    }
    val.push(v);
    cnt.push(c);
  }
  const out: number[] = [];
  val.forEach((v, i) => {
    for (let k = 0; k < cnt[i]; k++) out.push(v);
  });
  return out;
}

let exactFrame: { t: number; id: number; f: SeriesFrame; du1: Float64Array; du3: Float64Array } | null = null;

/**
 * Точные перемещения в узлах x_i = i·dx для метода фронтов: Δu(x) = −∫_x^{xView} (F − F0) dx',
 * интеграл кусочно-постоянной функции берётся по истинным положениям фронтов, а не по ячейкам.
 * (Суммирование по ячейкам ошибается на ΔF·dx/2 на каждом фронте — после умножения на
 * коэффициент увеличения смещений картинка дрожит при движении фронтов через ячейки.)
 */
function exactDisplacements(xs: number[], states: { F: number[] }[], N: number, dx: number) {
  const P1 = new Float64Array(N + 1); // ∫_0^{x_i} (p − p0)
  const P3 = new Float64Array(N + 1); // ∫_0^{x_i} (s − s0)
  let j = 0;
  let a1 = 0;
  let a3 = 0;
  let last = 0;
  for (let i = 0; i <= N; i++) {
    const x = i * dx;
    while (j < xs.length && xs[j] < x) {
      const xj = Math.max(last, xs[j]);
      a1 += (xj - last) * (states[j].F[0] - st.F0[0]);
      a3 += (xj - last) * (states[j].F[1] - st.F0[2]);
      last = xj;
      j++;
    }
    a1 += (x - last) * (states[j].F[0] - st.F0[0]);
    a3 += (x - last) * (states[j].F[1] - st.F0[2]);
    last = x;
    P1[i] = a1;
    P3[i] = a3;
  }
  const du1 = new Float64Array(N + 1);
  const du3 = new Float64Array(N + 1);
  for (let i = 0; i <= N; i++) {
    du1[i] = P1[i] - P1[N];
    du3[i] = P3[i] - P3[N];
  }
  return { du1, du3 };
}

/**
 * Кадр в момент t. Для метода фронтов строится точно в момент t (фронты сдвигаются со своими
 * скоростями от ближайшего предыдущего события), поэтому анимация плавная при любой длине расчёта.
 */
function frameAt(t: number): SeriesFrame {
  const fr = st.fronts;
  if (st.hyb) {
    if (exactFrame && exactFrame.t === t && exactFrame.id === reqId) return exactFrame.f;
    const r = st.res!;
    const I = hybInterp(t);
    const U = exactFrame?.f.U.length === 6 * r.N ? exactFrame.f.U : new Float32Array(6 * r.N);
    const N = r.N;
    const out = new Float64Array(4);
    for (let i = 0; i < N; i++) {
      I.at((i + 0.5) * r.dx, out);
      U[i] = st.F0[0] + out[0];
      U[N + i] = 0;
      U[2 * N + i] = st.F0[2] + out[1];
      U[3 * N + i] = out[2];
      U[4 * N + i] = 0;
      U[5 * N + i] = out[3];
    }
    // перемещения непрерывны и меняются плавно — линейно между точными перемещениями соседних кадров
    const A = hybridDisplacements(I.f0, N, r.dx);
    const B = hybridDisplacements(I.f1, N, r.dx);
    const du1 = A.du1.map((v, i) => v + I.w * (B.du1[i] - v));
    const du3 = A.du3.map((v, i) => v + I.w * (B.du3[i] - v));
    const ub = I.f0.ub + I.w * (I.f1.ub - I.f0.ub);
    exactFrame = { t, id: reqId, f: { t, U, ub }, du1, du3 };
    return exactFrame.f;
  }
  if (!fr) return storedFrameAt(t);
  if (exactFrame && exactFrame.t === t && exactFrame.id === reqId) return exactFrame.f;
  const r = st.res!;
  const N = r.N;
  const { f, xs } = frontsAt(t);
  const dt = t - f.t;
  const U = exactFrame?.f.U.length === 6 * N ? exactFrame.f.U : new Float32Array(6 * N);
  let j = 0;
  for (let i = 0; i < N; i++) {
    const x = (i + 0.5) * r.dx;
    while (j < xs.length && xs[j] < x) j++;
    const sv = f.states[j];
    U[i] = sv.F[0];
    U[N + i] = 0;
    U[2 * N + i] = sv.F[1];
    U[3 * N + i] = sv.v[0];
    U[4 * N + i] = 0;
    U[5 * N + i] = sv.v[1];
  }
  // перемещение границы: скорость границы постоянна между событиями
  const ub = f.ub + f.states[0].v[0] * dt;
  const { du1, du3 } = exactDisplacements(xs, f.states, N, r.dx);
  exactFrame = { t, id: reqId, f: { t, U, ub }, du1, du3 };
  return exactFrame.f;
}

let exactCache: { key: string; img: HTMLCanvasElement } | null = null;

/**
 * Точный фон диаграммы x–t для метода фронтов: решение кусочно-постоянное, поэтому состояние в
 * точке (x, t) берётся из ближайшего предыдущего сохранённого момента, в котором фронты сдвинуты
 * со своими скоростями до t. Растр строится под текущее окно с разрешением экрана
 * (время — по горизонтали, x — по вертикали).
 */
/** Перестройка фона t–x: когда окно последний раз менялось и сколько заняла последняя чёткая перестройка. */
const xtPerf = { winKey: '', changed: 0, fullMs: 0, timer: 0 };

function drawExactXt(ctx: CanvasRenderingContext2D, m: { l: number; r: number; t: number; b: number }, w: number, h: number, win: { x0: number; x1: number; t0: number; t1: number }) {
  const PW = Math.max(1, Math.round(w - m.l - m.r));
  const PH = Math.max(1, Math.round(h - m.t - m.b));
  // Пока окно меняется (колесо, перетаскивание), фон строится в пониженном разрешении, если чёткий
  // строится заметно долго; через ~0.25 с после остановки — перестраивается чётко.
  const now = performance.now();
  const winKey = [win.x0, win.x1, win.t0, win.t1].join('|');
  if (winKey !== xtPerf.winKey) {
    xtPerf.winKey = winKey;
    xtPerf.changed = now;
  }
  const interacting = now - xtPerf.changed < 250 && xtPerf.fullMs > 30;
  const sc = interacting ? Math.min(4, Math.max(2, Math.ceil(Math.sqrt(xtPerf.fullMs / 25)))) : 1;
  const pw = Math.ceil(PW / sc);
  const ph = Math.ceil(PH / sc);
  const key = [winKey, pw, ph, st.field, reqId, css('--panel'), fieldMax].join('|');
  if (!exactCache || exactCache.key !== key) {
    const tStart = performance.now();
    const img = document.createElement('canvas');
    img.width = pw;
    img.height = ph;
    const ic = img.getContext('2d')!;
    const data = ic.createImageData(pw, ph);
    const base = parseColor(css('--panel'));
    const hv = new Float64Array(4);
    const a = new Float64Array(4);
    const b = new Float64Array(4);
    const scale = st.hyb ? hybInfo().scale : null;
    /** Значение кадра в точке x проходом по узлам (x растёт от строки к строке); возвращает новый указатель. */
    const walk = (f: HybridFrame, x: number, j: number, out: Float64Array) => {
      const xs = f.xs;
      const W = f.W;
      const n = xs.length;
      while (j < n - 1 && xs[j] < x) j++;
      const xa = xs[j - 1];
      const xb = xs[j];
      const u = x <= xa ? 0 : x >= xb ? 1 : xb > xa ? (x - xa) / (xb - xa) : 1;
      const o = 4 * (j - 1);
      for (let k = 0; k < 4; k++) out[k] = W[o + k] + u * (W[o + 4 + k] - W[o + k]);
      return j;
    };
    for (let col = 0; col < pw && st.hyb; col++) {
      // гибрид: между сохранёнными моментами — перенос вдоль характеристик (без ступенек по времени);
      // там, где ничего не движется (соседние кадры совпадают), — быстрый проход по узлам
      const t = win.t0 + ((col + 0.5) / pw) * (win.t1 - win.t0);
      const I = hybInterp(t);
      let j0 = 1;
      let j1 = 1;
      for (let row = ph - 1; row >= 0; row--) {
        const x = win.x0 + ((ph - 1 - row + 0.5) / ph) * (win.x1 - win.x0);
        j0 = walk(I.f0, x, j0, a);
        j1 = walk(I.f1, x, j1, b);
        let mis = 0;
        for (let k = 0; k < 4; k++) mis += Math.abs(a[k] - b[k]) / scale![k];
        if (mis < 1e-3) for (let k = 0; k < 4; k++) hv[k] = a[k] + I.w * (b[k] - a[k]);
        else I.at(x, hv);
        const F: Vec3 = [st.F0[0] + hv[0], 0, st.F0[2] + hv[1]];
        const [R, G, B] = diverging(fieldOf(F, [hv[2], 0, hv[3]]) / fieldMax, base);
        const q = 4 * (row * pw + col);
        data.data[q] = R;
        data.data[q + 1] = G;
        data.data[q + 2] = B;
        data.data[q + 3] = 255;
      }
    }
    for (let col = 0; col < pw && !st.hyb; col++) {
      const t = win.t0 + ((col + 0.5) / pw) * (win.t1 - win.t0);
      const { f, xs } = frontsAt(t);
      // цвет каждого состояния этого момента
      const cols = f.states.map((sv) => diverging(fieldOf([sv.F[0], 0, sv.F[1]], [sv.v[0], 0, sv.v[1]]) / fieldMax, base));
      let j = 0;
      for (let row = ph - 1; row >= 0; row--) {
        const x = win.x0 + ((ph - 1 - row + 0.5) / ph) * (win.x1 - win.x0);
        while (j < xs.length && xs[j] < x) j++;
        const [R, G, B] = cols[j];
        const o = 4 * (row * pw + col);
        data.data[o] = R;
        data.data[o + 1] = G;
        data.data[o + 2] = B;
        data.data[o + 3] = 255;
      }
    }
    ic.putImageData(data, 0, 0);
    exactCache = { key, img };
    const ms = performance.now() - tStart;
    if (sc === 1) xtPerf.fullMs = ms;
    (window as unknown as { __xtMs: number }).__xtMs = ms;
  }
  if (sc > 1) {
    clearTimeout(xtPerf.timer);
    xtPerf.timer = window.setTimeout(() => drawXt(), 270);
  }
  ctx.imageSmoothingEnabled = sc > 1;
  ctx.drawImage(exactCache.img, m.l, m.t, w - m.l - m.r, h - m.t - m.b);
}

/** Текущее окно диаграммы x–t. */
function xtWindow() {
  return st.xtWin ?? { x0: 0, x1: st.xView, t0: 0, t1: st.tEnd };
}
let xtGeom: { m: { l: number; r: number; t: number; b: number }; w: number; h: number } | null = null;

/** (x, t) → пиксели: время по горизонтали, x по вертикали (граница x = 0 внизу). */
function xtToPx(x: number, t: number): [number, number] {
  const g = xtGeom!;
  const win = xtWindow();
  return [
    g.m.l + ((t - win.t0) / (win.t1 - win.t0)) * (g.w - g.m.l - g.m.r),
    g.h - g.m.b - ((x - win.x0) / (win.x1 - win.x0)) * (g.h - g.m.t - g.m.b),
  ];
}

/** Деления оси на отрезке [a, b]. */
function rangeTicks(a: number, b: number): number[] {
  const raw = (b - a) / 5;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((k) => k * mag).find((v) => v >= raw) ?? raw;
  const out: number[] = [];
  for (let v = Math.ceil(a / step) * step; v <= b + 1e-12 * (b - a); v += step) out.push(+v.toPrecision(10));
  return out;
}

/** Ограничить окно пределами расчёта и минимальным размером. */
function clampWin(wn: { x0: number; x1: number; t0: number; t1: number }) {
  const minX = st.xView / 5000;
  const minT = st.tEnd / 5000;
  let { x0, x1, t0, t1 } = wn;
  if (x1 - x0 < minX) [x0, x1] = [(x0 + x1) / 2 - minX / 2, (x0 + x1) / 2 + minX / 2];
  if (t1 - t0 < minT) [t0, t1] = [(t0 + t1) / 2 - minT / 2, (t0 + t1) / 2 + minT / 2];
  const sx = Math.min(x1 - x0, st.xView);
  const sT = Math.min(t1 - t0, st.tEnd);
  x0 = Math.max(0, Math.min(x0, st.xView - sx));
  t0 = Math.max(0, Math.min(t0, st.tEnd - sT));
  const full = sx >= st.xView - 1e-12 && sT >= st.tEnd - 1e-12;
  st.xtWin = full ? null : { x0, x1: x0 + sx, t0, t1: t0 + sT };
}

/** Перевести координаты экрана в (x, t) диаграммы. */
function xtAt(px: number, py: number) {
  const g = xtGeom!;
  const win = xtWindow();
  const ft = (px - g.m.l) / (g.w - g.m.l - g.m.r);
  const fx = (g.h - g.m.b - py) / (g.h - g.m.t - g.m.b);
  return { x: win.x0 + fx * (win.x1 - win.x0), t: win.t0 + ft * (win.t1 - win.t0), fx, ft };
}

function zoomXt(factor: number, px?: number, py?: number) {
  const g = xtGeom;
  if (!g) return;
  const win = xtWindow();
  const c = px === undefined ? { x: (win.x0 + win.x1) / 2, t: (win.t0 + win.t1) / 2 } : xtAt(px, py!);
  clampWin({ x0: c.x - (c.x - win.x0) / factor, x1: c.x + (win.x1 - c.x) / factor, t0: c.t - (c.t - win.t0) / factor, t1: c.t + (win.t1 - c.t) / factor });
  drawXt();
}

function setupXtZoom() {
  const cv = $<HTMLCanvasElement>('s-xt');
  cv.style.touchAction = 'none';
  cv.style.cursor = 'crosshair';
  cv.addEventListener(
    'wheel',
    (e) => {
      if (!st.res) return;
      e.preventDefault();
      const r = cv.getBoundingClientRect();
      zoomXt(Math.exp(-e.deltaY * 0.0015), e.clientX - r.left, e.clientY - r.top);
    },
    { passive: false },
  );
  // перетаскивание — сдвиг, щипок — масштаб, щелчок без движения — выбрать момент времени
  const pts = new Map<number, { x: number; y: number }>();
  let moved = false;
  let pinch0 = 0;
  cv.addEventListener('pointerdown', (e) => {
    cv.setPointerCapture(e.pointerId);
    const r = cv.getBoundingClientRect();
    pts.set(e.pointerId, { x: e.clientX - r.left, y: e.clientY - r.top });
    moved = false;
    if (pts.size === 2) {
      const [a, b] = [...pts.values()];
      pinch0 = Math.hypot(a.x - b.x, a.y - b.y);
    }
  });
  cv.addEventListener('pointermove', (e) => {
    const prev = pts.get(e.pointerId);
    if (!prev || !xtGeom || !st.res) return;
    const r = cv.getBoundingClientRect();
    const cur = { x: e.clientX - r.left, y: e.clientY - r.top };
    if (Math.hypot(cur.x - prev.x, cur.y - prev.y) > 2) moved = true;
    if (pts.size === 2) {
      pts.set(e.pointerId, cur);
      const [a, b] = [...pts.values()];
      const d = Math.hypot(a.x - b.x, a.y - b.y);
      if (pinch0 > 0) zoomXt(d / pinch0, (a.x + b.x) / 2, (a.y + b.y) / 2);
      pinch0 = d;
      return;
    }
    const g = xtGeom;
    const win = xtWindow();
    const dt = ((cur.x - prev.x) / (g.w - g.m.l - g.m.r)) * (win.t1 - win.t0);
    const dx = ((cur.y - prev.y) / (g.h - g.m.t - g.m.b)) * (win.x1 - win.x0);
    pts.set(e.pointerId, cur);
    if (st.xtWin) {
      clampWin({ x0: win.x0 + dx, x1: win.x1 + dx, t0: win.t0 - dt, t1: win.t1 - dt });
      cv.style.cursor = 'grabbing';
      drawXt();
    }
  });
  const end = (e: PointerEvent) => {
    const p = pts.get(e.pointerId);
    pts.delete(e.pointerId);
    cv.style.cursor = 'crosshair';
    if (pts.size < 2) pinch0 = 0;
    if (!moved && p && xtGeom && st.res && pts.size === 0) {
      const q = xtAt(p.x, p.y);
      if (q.fx >= 0 && q.fx <= 1 && q.ft >= 0 && q.ft <= 1) pickXt(p.x, p.y, q.x, q.t);
    }
  };
  cv.addEventListener('pointerup', end);
  cv.addEventListener('pointercancel', end);
  cv.addEventListener('dblclick', () => {
    st.xtWin = null;
    drawXt();
  });
  // кнопки
  const box = document.createElement('div');
  box.className = 'zoom-ctrl';
  for (const [label, title, fn] of [
    ['+', 'приблизить', () => zoomXt(1.6)],
    ['−', 'отдалить', () => zoomXt(1 / 1.6)],
    ['⤢', 'весь вид', () => ((st.xtWin = null), drawXt())],
  ] as [string, string, () => void][]) {
    const b = document.createElement('button');
    b.textContent = label;
    b.title = title;
    b.setAttribute('aria-label', title);
    b.addEventListener('click', fn);
    box.appendChild(b);
  }
  cv.parentElement!.appendChild(box);
}

/** Состояние среды в точке (x, t): точное для метода фронтов, из ближайшего кадра для КО. */
function stateAtXt(x: number, t: number): { F: Vec3; v: Vec3 } {
  if (st.hyb) return hybSampler(t)(x);
  if (st.fronts) {
    const { f, xs } = frontsAt(t);
    let j = 0;
    while (j < xs.length && xs[j] < x) j++;
    const sv = f.states[j];
    return { F: [sv.F[0], 0, sv.F[1]], v: [sv.v[0], 0, sv.v[1]] };
  }
  const r = st.res!;
  const i = Math.max(0, Math.min(r.N - 1, Math.floor(x / r.dx)));
  return cell(frameAt(t), r.N, i);
}

const num = (v: number) => (Math.abs(v) < 1e-10 ? '0' : Math.abs(v) < 1e-3 || Math.abs(v) >= 1e4 ? v.toExponential(3) : v.toFixed(5));

/** Строки подписи со значениями деформаций, скоростей и напряжений в состоянии. */
function stateLines(S: { F: Vec3; v: Vec3 }, title?: string): string[] {
  const m = normalize(st.material);
  const T = traction(m, S.F);
  const T0 = traction(m, st.F0);
  const out = [
    `p = u1,1 = ${num(S.F[0])}   (p − p0 = ${num(S.F[0] - st.F0[0])})`,
    `s = u3,1 = ${num(S.F[2])}   (s − s0 = ${num(S.F[2] - st.F0[2])})`,
    `v1/c1 = ${num(S.v[0])},  v3/c1 = ${num(S.v[2])}`,
    `σ11/(λ+2μ) = ${num(T[0])}   (Δ = ${num(T[0] - T0[0])})`,
    `σ31/(λ+2μ) = ${num(T[2])}`,
  ];
  return title ? [title, ...out.map((l) => '  ' + l)] : out;
}

const FAM_NAME = { fast: 'квазипродольный', slow: 'квазипоперечный' } as const;

/** Траектории фронтов (метод фронтов или выделенные разрывы гибрида). */
function segsOf(): Segment[] {
  return st.fronts?.segments ?? st.hyb?.segments ?? [];
}

/**
 * Щелчок по диаграмме: если рядом точка события (встреча фронтов, рождение на границе, отражение) —
 * берётся она; если рядом фронт — показываются его характеристики; иначе — просто (x, t).
 * Время анимации переводится на выбранный момент.
 */
function pickXt(px: number, py: number, x: number, t: number) {
  const toPx = xtToPx;
  let pick = { x, t, lines: [`t c1/L = ${fmt(t, 4)},  x1/L = ${fmt(x, 4)}`, ...stateLines(stateAtXt(x, t))] };
  const segs = segsOf();
  if (segs.length && st.showFronts) {
    // 1) события — начала отрезков фронтов (рождение на границе, встреча, отражение)
    const maxStr = Math.max(...segs.map((sg) => sg.strength), 1e-30);
    const visible = (sg: Segment) => st.weakFronts || sg.strength / maxStr >= 1e-4;
    let best: Segment | null = null;
    let bestD = 9; // пикселей
    for (const sg of segs) {
      if (!visible(sg) || sg.t0 <= 1e-12) continue;
      const [ex, ey] = toPx(sg.x0, sg.t0);
      const d = Math.hypot(ex - px, ey - py);
      if (d < bestD) {
        bestD = d;
        best = sg;
      }
    }
    if (best) {
      const at = segs.filter((sg) => Math.abs(sg.t0 - best!.t0) < 1e-9 && Math.abs(sg.x0 - best!.x0) < 1e-9 && visible(sg));
      const where = best.x0 < 1e-9 ? 'у границы x1 = 0' : 'встреча фронтов';
      const born = at.map((sg) => `${FAM_NAME[sg.fam]} ${sg.x1 >= sg.x0 ? '→' : '←'}`);
      pick = {
        x: best.x0,
        t: best.t0,
        lines: [`${where}: t c1/L = ${fmt(best.t0, 4)}`, `x1/L = ${fmt(best.x0, 4)}`, `выходят: ${[...new Set(born)].join(', ')}`],
      };
    } else {
      // 2) ближайший фронт (отрезок)
      let fb: Segment | null = null;
      let fd = 6;
      for (const sg of segs) {
        if (!visible(sg)) continue;
        const [ax, ay] = toPx(sg.x0, sg.t0);
        const [bx, by] = toPx(sg.x1, sg.t1);
        const L2 = (bx - ax) ** 2 + (by - ay) ** 2 || 1;
        const u = Math.max(0, Math.min(1, ((px - ax) * (bx - ax) + (py - ay) * (by - ay)) / L2));
        const d = Math.hypot(ax + u * (bx - ax) - px, ay + u * (by - ay) - py);
        if (d < fd) {
          fd = d;
          fb = sg;
        }
      }
      if (fb) {
        const G = (fb.x1 - fb.x0) / (fb.t1 - fb.t0);
        const tt = Math.max(fb.t0, Math.min(fb.t1, t));
        pick = {
          x: fb.x0 + G * (tt - fb.t0),
          t: tt,
          lines: [
            `${FAM_NAME[fb.fam]} фронт${fb.kind === 'simple' ? (st.hyb ? ' (тонкий веер — уже ячейки, ведётся одним разрывом)' : ' (пучок простой волны)') : ' (ударный)'}, ${G >= 0 ? 'вглубь' : 'к границе'}`,
            `скорость ${fmt(Math.abs(G), 4)} c1, скачок деформации ${fb.strength.toExponential(1)}`,
            `t c1/L = ${fmt(tt, 4)}, x1/L = ${fmt(fb.x0 + G * (tt - fb.t0), 4)}`,
            `рождён при t = ${fmt(fb.t0, 4)}`,
          ],
        };
        // состояния по обе стороны фронта (немного отступив от него)
        const eps = 1e-6 * Math.max(1, st.xView);
        const ahead = stateAtXt(pick.x + eps, tt);
        const behind = stateAtXt(pick.x - eps, tt);
        pick.lines.push(...stateLines(ahead, 'перед фронтом (x больше):'), ...stateLines(behind, 'за фронтом (x меньше):'));
      }
    }
  }
  st.xtPick = pick;
  st.playing = false;
  st.t = Math.max(0, Math.min(st.tEnd, pick.t));
  syncPlay();
  drawAll();
}

/** Круглые деления оси от 0 до max (4–6 штук). */
function axisTicks(max: number): number[] {
  const raw = max / 5;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((k) => k * mag).find((v) => v >= raw) ?? raw;
  const out: number[] = [];
  for (let v = 0; v <= max + 1e-9; v += step) out.push(+v.toPrecision(10));
  return out;
}

function parseColor(c: string): [number, number, number] {
  const m = c.match(/^#([0-9a-f]{6})$/i);
  if (m) return [parseInt(m[1].slice(0, 2), 16), parseInt(m[1].slice(2, 4), 16), parseInt(m[1].slice(4, 6), 16)];
  const r = c.match(/rgba?\(([^)]+)\)/);
  if (r) return r[1].split(',').slice(0, 3).map((x) => +x) as [number, number, number];
  return [255, 255, 255];
}

function diverging(t: number, base: [number, number, number]): [number, number, number] {
  const a = Math.min(1, Math.abs(t));
  const c = t >= 0 ? [214, 72, 52] : [48, 104, 214];
  return [0, 1, 2].map((k) => Math.round(base[k] + (c[k] - base[k]) * a)) as [number, number, number];
}

let fieldMax = 1;
function computeFieldMax() {
  const r = st.res!;
  let mx = 1e-12;
  const m = normalize(st.material);
  const T0 = traction(m, st.F0)[0];
  for (let f = 0; f < r.frames.length; f += 3)
    for (let i = 0; i < r.N; i += 2) {
      const U = r.frames[f].U;
      const val =
        st.field === 'p'
          ? U[i] - st.F0[0]
          : st.field === 's'
            ? U[2 * r.N + i] - st.F0[2]
            : st.field === 'v1'
              ? U[3 * r.N + i]
              : traction(m, [U[i], U[r.N + i], U[2 * r.N + i]])[0] - T0;
      mx = Math.max(mx, Math.abs(val));
    }
  fieldMax = mx;
}

let heatCache: { key: string; img: HTMLCanvasElement } | null = null;

function drawXt() {
  const canvas = $<HTMLCanvasElement>('s-xt');
  const S = setup(canvas);
  if (!S || !st.res) return;
  const { ctx, w, h } = S;
  const r = st.res;
  const m = { l: 52, r: 14, t: 26, b: 34 };
  const key = `${st.field}|${r.N}|${r.frames.length}|${css('--panel')}|${reqId}`;
  if (!heatCache || heatCache.key !== key) {
    const img = document.createElement('canvas');
    img.width = r.frames.length;
    img.height = r.N;
    const ic = img.getContext('2d')!;
    const data = ic.createImageData(r.frames.length, r.N);
    const base = parseColor(css('--panel'));
    const mm = normalize(st.material);
    for (let f = 0; f < r.frames.length; f++) {
      for (let i = 0; i < r.N; i++) {
        const c = cell(r.frames[f], r.N, i);
        const val = st.field === 's11' ? traction(mm, c.F)[0] - traction(mm, st.F0)[0] : fieldOf(c.F, c.v);
        const [R, G, B] = diverging(val / fieldMax, base);
        const o = 4 * ((r.N - 1 - i) * r.frames.length + f); // время вправо, x вверх
        data.data[o] = R;
        data.data[o + 1] = G;
        data.data[o + 2] = B;
        data.data[o + 3] = 255;
      }
    }
    ic.putImageData(data, 0, 0);
    heatCache = { key, img };
  }
  const win = xtWindow();
  xtGeom = { m, w, h };
  // фон — кусок растра поля, соответствующий окну (фронты сверху рисуются векторно и остаются чёткими)
  const img = heatCache.img;
  const sx = (win.t0 / st.tEnd) * img.width;
  const sw = ((win.t1 - win.t0) / st.tEnd) * img.width;
  const sy = (1 - win.x1 / st.xView) * img.height;
  const sh = ((win.x1 - win.x0) / st.xView) * img.height;
  if (st.fronts || st.hyb) drawExactXt(ctx, m, w, h, win);
  else {
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(img, sx, sy, sw, sh, m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  }
  const T = (t: number) => m.l + ((t - win.t0) / (win.t1 - win.t0)) * (w - m.l - m.r);
  const Xv = (x: number) => h - m.b - ((x - win.x0) / (win.x1 - win.x0)) * (h - m.t - m.b);
  ctx.strokeStyle = css('--border');
  ctx.strokeRect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  ctx.font = '600 12px system-ui, sans-serif';
  ctx.fillStyle = css('--text');
  ctx.fillText(`Диаграмма t–x: ${FIELD_LABELS[st.field]}`, 10, 16);
  ctx.font = '11px system-ui, sans-serif';
  ctx.fillStyle = css('--muted');
  for (const x of rangeTicks(win.x0, win.x1)) ctx.fillText(fmt(x, 3), 20, Xv(x) + 4);
  for (const t of rangeTicks(win.t0, win.t1)) ctx.fillText(fmt(t, 3), T(t) - 12, h - m.b + 18);
  const zoom = Math.min(st.xView / (win.x1 - win.x0), st.tEnd / (win.t1 - win.t0));
  if (zoom > 1.001) {
    ctx.textAlign = 'right';
    ctx.fillText(`увеличение ×${fmt(zoom, 3)} · двойной щелчок — весь вид`, w - m.r, 16);
    ctx.textAlign = 'left';
  }
  ctx.fillText('t c1 / L', w - m.r - 44, h - 4);
  ctx.save();
  ctx.translate(12, (m.t + h - m.b) / 2 + 20);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('x1 / L', 0, 0);
  ctx.restore();
  // фронты (метод волновых фронтов; у гибрида — выделенные ударные волны и тонкие вееры)
  if (segsOf().length && st.showFronts) {
    ctx.save();
    ctx.beginPath();
    ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
    ctx.clip();
    const segs = segsOf();
    const maxStr = Math.max(...segs.map((sg) => sg.strength), 1e-30);
    // слабые фронты рисуются первыми, чтобы сильные были сверху
    const order = [...(segs as Segment[])].sort((a, b) => a.strength - b.strength);
    for (const sg of order) {
      const rel = sg.strength / maxStr;
      if (!st.weakFronts && rel < 1e-4) continue; // по умолчанию самые слабые отражения не рисуются
      ctx.strokeStyle = sg.fam === 'fast' ? css('--front1') : css('--front2');
      if (st.weakFronts) {
        // логарифмическая шкала: от 1e-7 до 1 относительно сильнейшего фронта
        const lg = Math.max(0, Math.min(1, (Math.log10(rel) + 7) / 7));
        ctx.globalAlpha = 0.15 + 0.85 * lg;
        ctx.lineWidth = 0.6 + 1.6 * lg;
        ctx.setLineDash(sg.x1 < sg.x0 ? [3, 2] : []); // фронты, бегущие к границе, — штрихом
      } else {
        ctx.globalAlpha = sg.kind === 'simple' ? 0.55 : 0.95;
        ctx.lineWidth = sg.kind === 'simple' ? 0.8 : Math.max(0.8, Math.min(2.2, 0.8 + 3 * Math.sqrt(rel)));
      }
      ctx.beginPath();
      ctx.moveTo(T(sg.t0), Xv(sg.x0));
      ctx.lineTo(T(sg.t1), Xv(sg.x1));
      ctx.stroke();
    }
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;
    ctx.restore();
  }
  ctx.save();
  ctx.beginPath();
  ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b + 8);
  ctx.clip();
  // импульсы на оси времени (видны, только если граница x = 0 в окне)
  ctx.fillStyle = css('--accent');
  if (win.x0 <= 1e-12)
    for (const p of st.pulses) {
      const x1 = T(Math.min(st.tEnd, p.start));
      const x2 = T(Math.min(st.tEnd, p.start + p.duration));
      ctx.fillRect(x1, h - m.b + 2, Math.max(1.5, x2 - x1), 4);
    }
  // текущее время
  ctx.strokeStyle = css('--text');
  ctx.setLineDash([3, 3]);
  ctx.beginPath();
  ctx.moveTo(T(st.t), m.t);
  ctx.lineTo(T(st.t), h - m.b);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.restore();
  // выбранная точка
  if (st.xtPick) {
    const px = T(st.xtPick.t);
    const py = Xv(st.xtPick.x);
    if (px >= m.l - 2 && px <= w - m.r + 2 && py >= m.t - 2 && py <= h - m.b + 2) {
      ctx.strokeStyle = css('--text');
      ctx.fillStyle = css('--panel');
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(px, py, 5, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();
      ctx.font = '12px system-ui, sans-serif';
      const lw = Math.max(...st.xtPick.lines.map((ln) => ctx.measureText(ln).width)) + 14;
      const lh = 16 * st.xtPick.lines.length + 8;
      let bx = px + 10;
      let by = py - lh - 8;
      if (bx + lw > w - m.r) bx = px - lw - 10;
      if (by < m.t) by = py + 10;
      by = Math.max(m.t, Math.min(by, h - m.b - lh)); // высокая подпись не выходит за поле
      ctx.globalAlpha = 0.94;
      ctx.fillRect(bx, by, lw, lh);
      ctx.globalAlpha = 1;
      ctx.strokeStyle = css('--border');
      ctx.lineWidth = 1;
      ctx.strokeRect(bx, by, lw, lh);
      ctx.fillStyle = css('--text');
      st.xtPick.lines.forEach((ln, i) => ctx.fillText(ln, bx + 7, by + 18 + 16 * i));
    }
  }
  // легенда
  const lx = w - m.r - 150;
  ctx.fillStyle = css('--panel');
  ctx.fillRect(lx - 6, m.t + 6, 150, 34);
  const base = parseColor(css('--panel'));
  for (let i = 0; i < 120; i++) {
    const [R, G, B] = diverging((i / 119) * 2 - 1, base);
    ctx.fillStyle = `rgb(${R},${G},${B})`;
    ctx.fillRect(lx + i, m.t + 10, 1, 10);
  }
  ctx.fillStyle = css('--muted');
  ctx.fillText(`−${fmt(fieldMax, 2)}`, lx, m.t + 34);
  ctx.fillText(`+${fmt(fieldMax, 2)}`, lx + 86, m.t + 34);
}

function setup(canvas: HTMLCanvasElement) {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (!w || !h) return null;
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  }
  const ctx = canvas.getContext('2d')!;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = css('--panel');
  ctx.fillRect(0, 0, w, h);
  return { ctx, w, h };
}

function drawStrip() {
  const canvas = $<HTMLCanvasElement>('s-scene');
  const box3d = $('s-3d');
  canvas.hidden = st.view === '3d';
  box3d.hidden = st.view !== '3d';
  const S = st.view === '3d' ? null : setup(canvas);
  if ((!S && st.view === '2d') || !st.res) return;
  const { ctx, w, h } = S ?? { ctx: null as unknown as CanvasRenderingContext2D, w: 1, h: 1 };
  const r = st.res;
  const fr = frameAt(st.t);
  const N = r.N;
  const left = 70;
  const right = 14;
  const sc = (w - left - right) / st.xView;
  const top = 34;
  const bottom = 48; // место под переключатель вида
  const half = (h - top - bottom) / 2 / sc;
  const oy = top + (h - top - bottom) / 2;
  // перемещения относительно предварительно деформированного состояния: Δu(x) = −∫_x^1 (F − F0) dx'
  let du1: Float64Array;
  let du3: Float64Array;
  if ((st.fronts || st.hyb) && exactFrame && exactFrame.f === fr) {
    ({ du1, du3 } = exactFrame); // точные перемещения по положениям фронтов
  } else {
    du1 = new Float64Array(N + 1);
    du3 = new Float64Array(N + 1);
    for (let i = N - 1; i >= 0; i--) {
      du1[i] = du1[i + 1] - (fr.U[i] - st.F0[0]) * r.dx;
      du3[i] = du3[i + 1] - (fr.U[2 * N + i] - st.F0[2]) * r.dx;
    }
  }
  const K = st.magnify;
  const Kt = st.magnifyT;
  if (st.view === '3d') {
    render3d(fr, du1, du3);
    return;
  }
  const map = (i: number, y: number) => [left + (i * r.dx + K * du1[i]) * sc, oy - (y + Kt * du3[i]) * sc] as const;
  const base = parseColor(css('--panel'));
  for (let i = 0; i < N; i++) {
    const c = cell(fr, N, i);
    const [R, G, B] = diverging(fieldOf(c.F, c.v) / fieldMax, base);
    ctx.fillStyle = `rgb(${R},${G},${B})`;
    const a = map(i, half);
    const b = map(i + 1, half);
    const bb = map(i + 1, -half);
    const aa = map(i, -half);
    ctx.beginPath();
    ctx.moveTo(a[0], a[1]);
    ctx.lineTo(b[0] + 0.5, b[1]);
    ctx.lineTo(bb[0] + 0.5, bb[1]);
    ctx.lineTo(aa[0], aa[1]);
    ctx.closePath();
    ctx.fill();
  }
  // материальные линии x1 = const
  ctx.strokeStyle = css('--grid');
  ctx.lineWidth = 1;
  const every = Math.max(1, Math.round(N / 50));
  for (let i = every; i < N; i += every) {
    const a = map(i, half);
    const b = map(i, -half);
    ctx.beginPath();
    ctx.moveTo(a[0], a[1]);
    ctx.lineTo(b[0], b[1]);
    ctx.stroke();
  }
  for (const y of [-half, -half / 2, 0, half / 2, half]) {
    ctx.beginPath();
    for (let i = 0; i <= N; i += Math.max(1, Math.round(N / 400))) {
      const p = map(i, y);
      if (i === 0) ctx.moveTo(p[0], p[1]);
      else ctx.lineTo(p[0], p[1]);
    }
    ctx.stroke();
  }
  // граница и удар
  const g1 = map(0, half);
  const g2 = map(0, -half);
  ctx.strokeStyle = css('--surface');
  ctx.lineWidth = 3;
  ctx.beginPath();
  ctx.moveTo(g1[0], g1[1]);
  ctx.lineTo(g2[0], g2[1]);
  ctx.stroke();
  const vb = boundaryVelocity(st.pulses, fr.t);
  if (vb !== 0) {
    ctx.strokeStyle = css('--accent');
    ctx.fillStyle = css('--accent');
    ctx.lineWidth = 2;
    for (const yy of [-0.5, 0, 0.5]) {
      const [ex, ey] = map(0, yy * half);
      // удар — стрелка к границе (вглубь), рывок — стрелка от границы (наружу)
      const [tip, tail] = vb > 0 ? [ex - 6, ex - 50] : [ex - 50, ex - 6];
      const sgn = vb > 0 ? 1 : -1;
      ctx.beginPath();
      ctx.moveTo(tail, ey);
      ctx.lineTo(tip, ey);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(tip, ey);
      ctx.lineTo(tip - sgn * 8, ey - 5);
      ctx.lineTo(tip - sgn * 8, ey + 5);
      ctx.closePath();
      ctx.fill();
    }
  }
  ctx.font = '600 12px system-ui, sans-serif';
  ctx.fillStyle = css('--text');
  ctx.fillText('Сечение материала: смещения u1 (→) и u3 (↑); граница x1 = 0 слева', 10, 16);
  ctx.font = '12px system-ui, sans-serif';
  ctx.fillStyle = css('--muted');
  ctx.textAlign = 'right';
  ctx.fillText(`t c1/L = ${fr.t.toFixed(3)} · u1 ×${K}, u3 ×${Kt} · ${vb > 0 ? `удар: d = ${vb}` : vb < 0 ? `рывок наружу: d = ${vb}` : 'граница стоит'}`, w - 10, 16);
  ctx.textAlign = 'left';
}

// ---------------------------------------------------------------- 3D

let view3d: typeof import('./impact3d.ts') | null = null;
let view3dLoading = false;

function render3d(fr: SeriesResult['frames'][number], du1: Float64Array, du3: Float64Array) {
  if (!view3d) {
    if (!view3dLoading) {
      view3dLoading = true;
      import('./impact3d.ts').then((mod) => {
        view3d = mod;
        drawAll();
      });
    }
    return;
  }
  const r = st.res!;
  const N = r.N;
  const step = Math.max(1, Math.round(N / 400));
  const Xs: number[] = [];
  const D: Vec3[] = [];
  // 3D-брусок нормирован на длину показанного участка: X ∈ [0, 1]
  const L = st.xView;
  for (let i = 0; i <= N; i += step) {
    Xs.push((i * r.dx) / L);
    D.push([du1[i] / L, 0, du3[i] / L]);
  }
  if (Xs.at(-1)! < 1) {
    Xs.push(1);
    D.push([du1[N] / L, 0, du3[N] / L]);
  }
  const disp = (X: number): Vec3 => {
    const i = Math.max(0, Math.min(N, Math.round((X * L) / r.dx)));
    return [du1[i] / L, 0, du3[i] / L];
  };
  const base = parseColor(css('--panel'));
  const colorAt = (X: number) => {
    const i = Math.max(0, Math.min(N - 1, Math.floor((X * L) / r.dx)));
    const c = cell(fr, N, i);
    const [R, G, B] = diverging(fieldOf(c.F, c.v) / fieldMax, base);
    return `rgb(${R},${G},${B})`;
  };
  const vb = boundaryVelocity(st.pulses, fr.t);
  view3d.update($('s-3d'), {
    Xs,
    D,
    disp,
    K: st.magnify,
    Kt: st.magnifyT,
    Xmax: 1,
    colorAt,
    fronts: [],
    load: vb > 0 ? [1, 0, 0] : vb < 0 ? [-1, 0, 0] : null,
    loadLabel: vb > 0 ? `удар: d = ${vb}` : vb < 0 ? `рывок наружу: d = ${vb}` : 'граница стоит',
    theme: { bg: css('--panel'), grid: css('--muted'), text: css('--text'), muted: css('--muted'), surface: css('--surface'), accent: css('--accent') },
    caption: `t c1/L = ${fr.t.toFixed(3)} · u1 ×${st.magnify}, u3 ×${st.magnifyT} · цвет: ${FIELD_LABELS[st.field]} · вращение — мышью`,
  });
}

function buildViewToggle() {
  const btns = ['s-view-2d', 's-view-3d'].map((id) => $<HTMLButtonElement>(id));
  const reset = $<HTMLButtonElement>('s-view-reset');
  for (const b of btns)
    b.addEventListener('click', () => {
      st.view = b.dataset.view as typeof st.view;
      for (const x of btns) x.setAttribute('aria-pressed', String(x === b));
      reset.hidden = st.view !== '3d';
      drawAll();
    });
  reset.addEventListener('click', () => view3d?.resetView($('s-3d')));
}

/** Точное решение для первого удара (до момента следующего изменения скорости границы), если применимо. */
function exactFirst(t: number) {
  const changes = [...new Set(st.pulses.flatMap((p) => [p.start, p.start + p.duration]))].filter((x) => x > 0).sort((a, b) => a - b);
  const first = st.pulses.filter((p) => p.start === 0);
  if (!first.length || t <= 0 || (changes.length && t > changes[0])) return null;
  try {
    return solveImpact({ material: normalize(st.material), F0: st.F0, bc: { type: 'velocity', d: boundaryVelocity(st.pulses, 0) } });
  } catch {
    return null;
  }
}

/**
 * Состояние в момент t с непрерывными простыми волнами. Метод фронтов представляет центрированный
 * пучок (в точном решении — непрерывную волну) набором слабых фронтов-лучей с шагом по амплитуде
 * fanStep. Здесь внутри пучка состояние интерполируется линейно: крайние состояния — на крайних
 * лучах, промежуточные — посередине между соседними лучами. Ударные волны остаются разрывами.
 */
function smoothStateFn(f: FrontsFrame, xs: number[]) {
  const runs: { a: number; b: number }[] = [];
  for (let k = 0; k < xs.length; ) {
    let b = k;
    if (f.fan[k] >= 0) while (b + 1 < xs.length && f.fan[b + 1] === f.fan[k]) b++;
    if (b > k) runs.push({ a: k, b });
    k = b + 1;
  }
  const toV3 = (F: number[], v: number[]) => ({ F: [F[0], 0, F[1]] as Vec3, v: [v[0], 0, v[1]] as Vec3 });
  return (x: number) => {
    // пучок, содержащий x
    let lo = 0;
    let hi = runs.length - 1;
    while (lo <= hi) {
      const mid = (lo + hi) >> 1;
      const r = runs[mid];
      if (x < xs[r.a]) hi = mid - 1;
      else if (x > xs[r.b]) lo = mid + 1;
      else {
        // узлы: (xs[a], S[a]), (середины, S[a+1..b]), (xs[b], S[b+1])
        const nx: number[] = [xs[r.a]];
        const ns = [f.states[r.a]];
        for (let k = r.a + 1; k <= r.b; k++) {
          nx.push((xs[k - 1] + xs[k]) / 2);
          ns.push(f.states[k]);
        }
        nx.push(xs[r.b]);
        ns.push(f.states[r.b + 1]);
        let j = 1;
        while (j < nx.length - 1 && nx[j] < x) j++;
        const w = nx[j] > nx[j - 1] ? (x - nx[j - 1]) / (nx[j] - nx[j - 1]) : 1;
        const A = ns[j - 1];
        const B = ns[j];
        return toV3([A.F[0] + w * (B.F[0] - A.F[0]), A.F[1] + w * (B.F[1] - A.F[1])], [A.v[0] + w * (B.v[0] - A.v[0]), A.v[1] + w * (B.v[1] - A.v[1])]);
      }
    }
    let j = 0;
    let h = xs.length;
    while (j < h) {
      const mid = (j + h) >> 1;
      if (xs[mid] < x) j = mid + 1;
      else h = mid;
    }
    const sv = f.states[j];
    return toV3(sv.F, sv.v);
  };
}

const profShared: SharedX = { win: null, full: () => [0, st.xView], redraw: () => renderProfiles() };
let profPlots: ReturnType<typeof createProfilePlot>[] = [];
const PROFILE_ITEMS: [string, (F: Vec3, v: Vec3) => number][] = [
  ['p = u1,1', (F) => F[0]],
  ['s = u3,1', (F) => F[2]],
  ['v1 / c1', (_F, v) => v[0]],
  ['σ11 / (λ+2μ)', (F) => traction(normalize(st.material), F)[0]],
];

/**
 * Профили в текущий момент. Для метода фронтов — точная ступенчатая функция по положениям фронтов
 * в момент t, для КО — значения в ячейках; точное решение §2.3 — функцией. Перерисовка дешёвая,
 * поэтому выполняется в каждом кадре анимации.
 */
function renderProfiles() {
  const box = $('s-profiles');
  if (!profPlots.length) {
    box.innerHTML = '';
    profPlots = PROFILE_ITEMS.map(([title]) => createProfilePlot(box, title, 'x1 / L', profShared));
  }
  if (!st.res) return;
  const r = st.res;
  const t = st.t;
  const ex = exactFirst(t);
  const acc = css('--accent');
  const muted = css('--muted');
  let edges: number[];
  let states: { F: Vec3; v: Vec3 }[];
  let smooth: ((x: number) => { F: Vec3; v: Vec3 }) | null = null;
  if (st.hyb) {
    const f = hybFrameAt(t);
    smooth = hybSampler(t);
    // для автомасштаба — значения в узлах (между узлами профиль линейный)
    edges = Array.from({ length: f.xs.length - 1 }, (_, i) => 0.5 * (f.xs[i] + f.xs[i + 1]));
    states = Array.from({ length: f.xs.length }, (_, i) => ({ F: [st.F0[0] + f.W[4 * i], 0, st.F0[2] + f.W[4 * i + 1]] as Vec3, v: [f.W[4 * i + 2], 0, f.W[4 * i + 3]] as Vec3 }));
  } else if (st.fronts) {
    const { f, xs } = frontsAt(t);
    if (st.smoothFans && f.fan?.some((id, i) => id >= 0 && f.fan[i + 1] === id)) smooth = smoothStateFn(f, xs);
    edges = xs;
    states = f.states.map((sv) => ({ F: [sv.F[0], 0, sv.F[1]] as Vec3, v: [sv.v[0], 0, sv.v[1]] as Vec3 }));
  } else {
    const fr = storedFrameAt(t);
    edges = Array.from({ length: r.N - 1 }, (_, i) => (i + 1) * r.dx);
    states = Array.from({ length: r.N }, (_, i) => cell(fr, r.N, i));
  }
  PROFILE_ITEMS.forEach(([, fn], k) => {
    const label = st.hyb ? 'гибрид (ударные волны + сетка)' : st.fronts ? 'метод фронтов' : 'конечные объёмы';
    const series: PlotSeries[] = [
      smooth
        ? {
            label,
            color: acc,
            fn: (x) => {
              const q = smooth!(x);
              return fn(q.F, q.v);
            },
            steps: { edges, vals: states.map((q) => fn(q.F, q.v)) }, // для автомасштаба (узкие пики)
          }
        : { label, color: acc, steps: { edges, vals: states.map((q) => fn(q.F, q.v)) } },
    ];
    if (ex && t > 0)
      series.push({
        label: 'точное (§2.3)',
        color: muted,
        dash: [5, 4],
        fn: (x) => {
          const q = profileAt(ex, Math.max(0, x) / t);
          return fn(q.F, q.v);
        },
      });
    profPlots[k].draw(series);
  });
  // для браузерного теста
  Object.assign(window, { __profSeries: ex && t > 0 ? 2 : 1, __profWin: profShared.win });
}

function renderStatic() {
  computeFieldMax();
  heatCache = null;
  // график скорости границы
  const box = $('s-boundary');
  box.innerHTML = '';
  const ts: number[] = [];
  const vs: number[] = [];
  for (let i = 0; i <= 600; i++) {
    const t = (st.tEnd * i) / 600;
    ts.push(t);
    vs.push(boundaryVelocity(st.pulses, t));
  }
  box.appendChild(lineChart({ xLabel: 't c1 / L', yLabel: 'd = v*/c1', series: [{ label: 'скорость границы', x: ts, y: vs, color: css('--accent') }], height: 200, width: 520, xRange: [0, st.tEnd] }));
  const ub: number[] = st.res!.frames.map((f) => f.ub);
  box.appendChild(
    lineChart({ xLabel: 't c1 / L', yLabel: 'u1(0, t) / L', series: [{ label: 'перемещение границы', x: st.res!.frames.map((f) => f.t), y: ub, color: css('--front1') }], height: 200, width: 520, xRange: [0, st.tEnd] }),
  );
}

function drawAll() {
  if (!st.res) return;
  drawStrip();
  drawXt();
  renderProfiles();
}

let lastT = 0;
function loop(t: number) {
  const dt = Math.min(0.05, (t - lastT) / 1000 || 0);
  lastT = t;
  if (!st.playing) return;
  st.t += dt * st.speed;
  if (st.t > st.tEnd * 1.1) st.t = 0;
  const saved = st.t;
  st.t = Math.min(st.t, st.tEnd);
  syncPlay();
  if ($('page-series').hidden === false) drawAll();
  st.t = saved;
  requestAnimationFrame(loop);
}

buildControls();
buildViewToggle();
setupXtZoom();
compute();
new ResizeObserver(() => drawAll()).observe($('s-scene'));
new ResizeObserver(() => st.view === '3d' && drawAll()).observe($('s-3d'));
new ResizeObserver(() => drawAll()).observe($('s-xt'));
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
  heatCache = null;
  drawAll();
});
requestAnimationFrame(loop);

(window as unknown as { __series: typeof st }).__series = st;
