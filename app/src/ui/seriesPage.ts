/**
 * Страница «Серия ударов»: несколько нормальных ударов (импульсов скорости границы) с интервалами.
 * Решение нестационарное (волны разных ударов догоняют друг друга) — считается численно (fv1d.ts),
 * а до встречи волн сравнивается с точным автомодельным решением главы 2.
 */
import { type FrontsResult, type Segment, stateAt } from '../core/fronts.ts';
import { type Pulse, type SeriesProblem, type SeriesResult, boundaryVelocity, cell } from '../core/fv1d.ts';
import { profileAt, solveImpact } from '../core/impact.ts';
import { type Material, normalize } from '../core/material.ts';
import { type Vec3, traction } from '../core/material1d.ts';
import { MATERIAL_1, MATERIAL_2 } from '../core/presets.ts';
import { type Series, fmt, lineChart } from './charts.ts';

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
  'Один удар (проверка по §2.3)': [{ start: 0, duration: 1, d: 0.001 }],
};

const matNames = Object.keys(MATERIALS);
const st = {
  matName: matNames[0],
  material: { ...MATERIALS[matNames[0]].material } as Material,
  F0: [...MATERIALS[matNames[0]].F0] as Vec3,
  pulses: PULSE_PRESETS['Три одинаковых удара'].map((p) => ({ ...p })),
  N: 800,
  tEnd: 0.92,
  field: 'p' as FieldKey,
  magnify: 100,
  /** Увеличение поперечных (сдвиговых) перемещений — при нормальном ударе они в сотни раз меньше. */
  magnifyT: 1000,
  view: '2d' as '2d' | '3d',
  playing: true,
  speed: 0.12,
  t: 0.6,
  res: null as SeriesResult | null,
  /** 'fronts' — метод волновых фронтов (как в диссертации), 'fv' — конечные объёмы (для сравнения). */
  method: 'fronts' as 'fronts' | 'fv',
  showFronts: true,
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

  const fsU = fieldset(c, 'Удары по границе');
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
  hint.textContent = 'Каждый удар: с момента t граница движется вглубь со скоростью d·c1 в течение Δt, затем стоит. Время — в единицах L/c1 (L — длина показанного участка).';
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
  msel.add(new Option('конечные объёмы (для сравнения)', 'fv'));
  msel.value = st.method;
  msel.addEventListener('change', () => {
    st.method = msel.value as typeof st.method;
    schedule();
  });
  row(fsN, 'Метод', msel, 's-method');
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
  const nsel = document.createElement('select');
  for (const n of [400, 800, 1600]) nsel.add(new Option(`${n} ячеек${n === 800 ? ' (по умолчанию)' : ''}`, String(n)));
  nsel.value = String(st.N);
  nsel.addEventListener('change', () => {
    st.N = +nsel.value;
    schedule();
  });
  row(fsN, 'Сетка (КО)', nsel, 's-grid');
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
  row(fsN, '× сдвиг', magT, 's-magT');
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
  head.innerHTML = '<span>№</span><span>t начала</span><span>Δt</span><span>d = v/c1</span><span></span>';
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
  st.computing = 0.001; // расчёт запланирован
  timer = window.setTimeout(compute, 250);
}

function problem(): SeriesProblem {
  return { material: normalize(st.material), F0: [...st.F0] as Vec3, pulses: st.pulses.map((p) => ({ ...p })), N: st.N, tEnd: st.tEnd, frames: 240 };
}

async function compute() {
  st.error = '';
  if (!(st.F0[0] < 0)) {
    st.error = 'Предварительная деформация должна быть сжатием: p0 < 0.';
    renderStatus();
    return;
  }
  const id = ++reqId;
  st.computing = 0.001;
  renderStatus();
  const prob = problem();
  const done = (res: SeriesResult, ms: number, fronts: FrontsResult | null = null) => {
    if (id !== reqId) return;
    st.res = res;
    st.fronts = fronts;
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
        else done(msg.res, msg.ms);
      }
      else {
        st.error = msg.message;
        st.computing = 0;
        renderStatus();
      }
    };
    if (st.method === 'fronts') worker.postMessage({ id, method: 'fronts', prob: { material: prob.material, F0: prob.F0, pulses: prob.pulses, tEnd: prob.tEnd, frames: prob.frames } });
    else worker.postMessage({ id, method: 'fv', prob });
  } catch {
    // запасной путь — в основном потоке
    const t0 = performance.now();
    if (st.method === 'fronts') {
      const { trackFronts } = await import('../core/fronts.ts');
      const fr = trackFronts({ material: prob.material, F0: prob.F0, pulses: prob.pulses, tEnd: prob.tEnd, frames: prob.frames });
      done(framesFromFronts(fr), performance.now() - t0, fr);
    } else {
      const { simulate } = await import('../core/fv1d.ts');
      done(simulate(prob), performance.now() - t0);
    }
  }
}

/** Перевести кусочно-постоянное решение метода фронтов в кадры на сетке (для общих графиков). */
function framesFromFronts(fr: FrontsResult): SeriesResult {
  const N = 1000;
  const dx = 1 / N;
  const frames = fr.frames.map((f) => {
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
  const r = st.res!;
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
  st.magnifyT = Math.max(1, Math.min(1e6, Math.round(0.05 / us)));
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
  const r = st.res;
  const fr = st.fronts;
  if (r && fr && !st.computing) {
    add(
      'info',
      `Метод волновых фронтов: ${fr.events} событий (удары, остановки, встречи и отражения фронтов), до ${fr.maxFronts} фронтов одновременно, ${st.ms.toFixed(0)} мс. ` +
        `На каждом фронте точно выполнены условия на разрывах (2.40)–(2.42); простые волны представлены пучками слабых фронтов. ` +
        `Кинематический баланс ${fr.kinematicError.toExponential(1)}. Энергия: работа границы ${fmt(fr.work, 4)} = прирост энергии ${fmt(fr.energyGain, 4)} + диссипация на ударных фронтах (2.45) ${fmt(fr.dissipation, 3)} (расхождение ${(fr.work - fr.energyGain - fr.dissipation).toExponential(1)}).`,
    );
    for (const w of fr.warnings) add('err', w);
  } else if (r && !st.computing) {
    const diss = r.work > 0 ? (r.work - r.energyGain) / r.work : 0;
    add(
      'info',
      `Численное решение: ${r.N} ячеек, ${r.steps} шагов, ${st.ms.toFixed(0)} мс. ` +
        `Кинематический баланс ∫(p − p0)dx + ∫v*dt = ${r.kinematicError.toExponential(1)}. ` +
        `Работа границы ${fmt(r.work, 4)}, прирост энергии среды ${fmt(r.energyGain, 4)}, рассеяно ${(diss * 100).toFixed(2)} % (ударные волны + численная вязкость; убывает с измельчением сетки).`,
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

function frameAt(t: number) {
  const r = st.res!;
  const k = Math.max(0, Math.min(r.frames.length - 1, Math.round((t / st.tEnd) * (r.frames.length - 1))));
  return r.frames[k];
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
let lastProfileT = -1;

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
    img.width = r.N;
    img.height = r.frames.length;
    const ic = img.getContext('2d')!;
    const data = ic.createImageData(r.N, r.frames.length);
    const base = parseColor(css('--panel'));
    const mm = normalize(st.material);
    for (let f = 0; f < r.frames.length; f++) {
      const row = r.frames.length - 1 - f; // время вверх
      for (let i = 0; i < r.N; i++) {
        const c = cell(r.frames[f], r.N, i);
        const val = st.field === 's11' ? traction(mm, c.F)[0] - traction(mm, st.F0)[0] : fieldOf(c.F, c.v);
        const [R, G, B] = diverging(val / fieldMax, base);
        const o = 4 * (row * r.N + i);
        data.data[o] = R;
        data.data[o + 1] = G;
        data.data[o + 2] = B;
        data.data[o + 3] = 255;
      }
    }
    ic.putImageData(data, 0, 0);
    heatCache = { key, img };
  }
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(heatCache.img, m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  const X = (x: number) => m.l + x * (w - m.l - m.r);
  const Y = (t: number) => h - m.b - (t / st.tEnd) * (h - m.t - m.b);
  ctx.strokeStyle = css('--border');
  ctx.strokeRect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  ctx.font = '600 12px system-ui, sans-serif';
  ctx.fillStyle = css('--text');
  ctx.fillText(`Диаграмма x–t: ${FIELD_LABELS[st.field]}`, 10, 16);
  ctx.font = '11px system-ui, sans-serif';
  ctx.fillStyle = css('--muted');
  for (const t of [0, 0.2, 0.4, 0.6, 0.8]) if (t <= st.tEnd) ctx.fillText(t.toFixed(1), 14, Y(t) + 4);
  for (const x of [0, 0.25, 0.5, 0.75, 1]) ctx.fillText(x.toFixed(2), X(x) - 10, h - m.b + 14);
  ctx.fillText('x1 / L', w - m.r - 34, h - 6);
  ctx.save();
  ctx.translate(12, (m.t + h - m.b) / 2 + 30);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('t c1 / L', 0, 0);
  ctx.restore();
  // фронты (метод волновых фронтов)
  if (st.fronts && st.showFronts) {
    ctx.save();
    ctx.beginPath();
    ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
    ctx.clip();
    const segs = st.fronts.segments;
    const maxStr = Math.max(...segs.map((sg) => sg.strength), 1e-30);
    for (const sg of segs as Segment[]) {
      const rel = sg.strength / maxStr;
      if (rel < 1e-4) continue; // самые слабые отражения не загромождают картину
      ctx.strokeStyle = sg.fam === 'fast' ? css('--front1') : css('--front2');
      ctx.globalAlpha = sg.kind === 'simple' ? 0.55 : 0.95;
      ctx.lineWidth = sg.kind === 'simple' ? 0.8 : Math.max(0.8, Math.min(2.2, 0.8 + 3 * Math.sqrt(rel)));
      ctx.beginPath();
      ctx.moveTo(X(sg.x0), Y(sg.t0));
      ctx.lineTo(X(sg.x1), Y(sg.t1));
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
    ctx.restore();
  }
  // импульсы на оси времени
  ctx.fillStyle = css('--accent');
  for (const p of st.pulses) {
    const y1 = Y(Math.min(st.tEnd, p.start));
    const y2 = Y(Math.min(st.tEnd, p.start + p.duration));
    ctx.fillRect(m.l - 6, y2, 4, Math.max(1.5, y1 - y2));
  }
  // текущее время
  ctx.strokeStyle = css('--text');
  ctx.setLineDash([3, 3]);
  ctx.beginPath();
  ctx.moveTo(m.l, Y(st.t));
  ctx.lineTo(w - m.r, Y(st.t));
  ctx.stroke();
  ctx.setLineDash([]);
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
  const sc = (w - left - right) / 1.0;
  const top = 34;
  const bottom = 48; // место под переключатель вида
  const half = (h - top - bottom) / 2 / sc;
  const oy = top + (h - top - bottom) / 2;
  // перемещения относительно предварительно деформированного состояния: Δu(x) = −∫_x^1 (F − F0) dx'
  const du1 = new Float64Array(N + 1);
  const du3 = new Float64Array(N + 1);
  for (let i = N - 1; i >= 0; i--) {
    du1[i] = du1[i + 1] - (fr.U[i] - st.F0[0]) * r.dx;
    du3[i] = du3[i + 1] - (fr.U[2 * N + i] - st.F0[2]) * r.dx;
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
  if (vb > 0) {
    ctx.strokeStyle = css('--accent');
    ctx.fillStyle = css('--accent');
    ctx.lineWidth = 2;
    for (const yy of [-0.5, 0, 0.5]) {
      const [ex, ey] = map(0, yy * half);
      ctx.beginPath();
      ctx.moveTo(ex - 50, ey);
      ctx.lineTo(ex - 6, ey);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(ex - 6, ey);
      ctx.lineTo(ex - 14, ey - 5);
      ctx.lineTo(ex - 14, ey + 5);
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
  ctx.fillText(`t c1/L = ${fr.t.toFixed(3)} · u1 ×${K}, u3 ×${Kt} · ${vb > 0 ? `удар: d = ${vb}` : 'граница стоит'}`, w - 10, 16);
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
  for (let i = 0; i <= N; i += step) {
    Xs.push(i * r.dx);
    D.push([du1[i], 0, du3[i]]);
  }
  if (Xs.at(-1)! < 1) {
    Xs.push(1);
    D.push([du1[N], 0, du3[N]]);
  }
  const disp = (X: number): Vec3 => {
    const i = Math.max(0, Math.min(N, Math.round(X / r.dx)));
    return [du1[i], 0, du3[i]];
  };
  const base = parseColor(css('--panel'));
  const colorAt = (X: number) => {
    const i = Math.max(0, Math.min(N - 1, Math.floor(X / r.dx)));
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
    load: vb > 0 ? [1, 0, 0] : null,
    loadLabel: vb > 0 ? `удар: d = ${vb}` : 'граница стоит',
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

function renderProfiles() {
  const box = $('s-profiles');
  box.innerHTML = '';
  if (!st.res) return;
  const r = st.res;
  const fr = frameAt(st.t);
  const xs: number[] = [];
  const step = Math.max(1, Math.round(r.N / 500));
  for (let i = 0; i < r.N; i += step) xs.push((i + 0.5) * r.dx);
  const m = normalize(st.material);
  const ex = exactFirst(fr.t);
  const items: [string, (F: Vec3, v: Vec3) => number][] = [
    ['p = u1,1', (F) => F[0]],
    ['s = u3,1', (F) => F[2]],
    ['v1 / c1', (_F, v) => v[0]],
    ['σ11 / (λ+2μ)', (F) => traction(m, F)[0]],
  ];
  const acc = css('--accent');
  const muted = css('--muted');
  for (const [title, f] of items) {
    const ys = xs.map((_, j) => {
      const c = cell(fr, r.N, j * step);
      return f(c.F, c.v);
    });
    const series: Series[] = [{ label: st.fronts ? 'метод фронтов' : 'конечные объёмы', x: xs, y: ys, color: acc }];
    if (ex) {
      const ey = xs.map((x) => {
        const q = profileAt(ex, x / fr.t);
        return f(q.F, q.v);
      });
      series.push({ label: 'точное (§2.3)', x: xs, y: ey, color: muted, dash: '5 4' });
    }
    const chart = lineChart({ xLabel: 'x1 / L', yLabel: title, series, allowNegative: true, height: 230, width: 520, xRange: [0, 1] });
    box.appendChild(chart);
  }
}

function renderStatic() {
  computeFieldMax();
  heatCache = null;
  lastProfileT = -1; // новый результат — профили перерисовать
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
  const fr = frameAt(st.t);
  if (fr.t !== lastProfileT) {
    lastProfileT = fr.t;
    renderProfiles();
  }
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
