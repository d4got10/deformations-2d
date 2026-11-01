/**
 * Страница «Удар по границе» (глава 2): нормальный (§2.3) и косой (§2.4) удар по
 * предварительно деформированному полупространству.
 */
import { type ImpactSolution, displacementTable, profileAt, solveImpact } from '../core/impact.ts';
import { type Material, normalize } from '../core/material.ts';
import { type Vec3, traction as tractionModel, traction229 } from '../core/material1d.ts';
import { IMPACT_PRESETS } from '../core/presets.ts';
import { type Series, fmt, lineChart, profileChart } from './charts.ts';

type FieldKey = 'zones' | 'p' | 'g' | 's' | 'shear' | 's11' | 's21' | 's31' | 'speed';
const FIELD_LABELS: Record<FieldKey, string> = {
  zones: 'зоны между волнами',
  p: 'p = u1,1 (продольная деформация)',
  g: 'g = u2,1 (сдвиг в плоскости x1x2)',
  s: 's = u3,1 (сдвиг в плоскости x1x3)',
  shear: '√(g² + s²) (величина сдвига)',
  s11: 'σ11 / (λ+2μ)',
  s21: 'σ21 / (λ+2μ)',
  s31: 'σ31 / (λ+2μ)',
  speed: '|v| / c1 (скорость частиц)',
};

const presetNames = Object.keys(IMPACT_PRESETS);
const first = IMPACT_PRESETS[presetNames[0]];

const st = {
  preset: presetNames[0],
  material: { ...first.material } as Material,
  F0: [...first.F0] as Vec3,
  bcType: first.bc.type as 'velocity' | 'strain',
  d: first.bc.type === 'velocity' ? first.bc.d : 0.001,
  Fb: (first.bc.type === 'strain' ? [...first.bc.F] : [-0.01, 0.005, 0.005]) as Vec3,
  law: 'model' as 'model' | 'printed229',
  magnify: 20,
  field: 'zones' as FieldKey,
  delta: false,
  playing: true,
  speed: 0.15,
  tau: 0.75,
  sweepParam: 'd' as 'd' | 'p' | 'g' | 's',
  view: '2d' as '2d' | '3d',
  sol: null as ImpactSolution | null,
  wTable: null as ((xi: number) => Vec3) | null,
  error: '',
};

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const css = (name: string) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

// ---------------------------------------------------------------- управление

const inputs = new Map<string, { range?: HTMLInputElement; num: HTMLInputElement }>();

function fieldset(parent: HTMLElement, title: string) {
  const fs = document.createElement('fieldset');
  const lg = document.createElement('legend');
  lg.textContent = title;
  fs.appendChild(lg);
  parent.appendChild(fs);
  return fs;
}

function numRow(parent: HTMLElement, key: string, label: string, get: () => number, set: (v: number) => void, range?: [number, number, number], hint?: string, live = true) {
  const row = document.createElement('div');
  row.className = 'row';
  const lab = document.createElement('label');
  lab.textContent = label;
  lab.htmlFor = `i-in-${key}`;
  if (hint) lab.title = hint;
  row.appendChild(lab);
  let r: HTMLInputElement | undefined;
  if (range) {
    r = document.createElement('input');
    r.type = 'range';
    [r.min, r.max, r.step] = range.map(String);
    r.value = String(get());
    r.setAttribute('aria-label', label);
    row.appendChild(r);
  } else row.appendChild(document.createElement('span'));
  const num = document.createElement('input');
  num.type = 'number';
  num.id = `i-in-${key}`;
  num.step = range ? String(range[2]) : 'any';
  num.value = String(get());
  row.appendChild(num);
  const apply = (v: number, fromRange: boolean) => {
    if (!Number.isFinite(v)) return;
    set(v);
    if (!fromRange) r && (r.value = String(v));
    num.value = String(+v.toPrecision(6));
    if (live) recompute();
  };
  r?.addEventListener('input', () => apply(+r!.value, true));
  num.addEventListener('change', () => apply(+num.value, false));
  parent.appendChild(row);
  inputs.set(key, { range: r, num });
}

const customLabel = 'Свои параметры';
let presetSel: HTMLSelectElement;
let bcBox: HTMLElement;
let segBtns: HTMLButtonElement[] = [];

function markCustom() {
  presetSel.value = 'custom';
  st.preset = 'custom';
}

function buildControls() {
  const c = $('i-controls');
  const fsP = fieldset(c, 'Расчётный случай');
  presetSel = document.createElement('select');
  presetSel.id = 'i-preset';
  presetSel.style.width = '100%';
  for (const n of presetNames) presetSel.add(new Option(n, n));
  presetSel.add(new Option(customLabel, 'custom'));
  presetSel.addEventListener('change', () => {
    if (presetSel.value === 'custom') return;
    applyPreset(presetSel.value);
  });
  fsP.appendChild(presetSel);

  const fsM = fieldset(c, 'Упругие модули, 10³ МПа');
  for (const [k, l] of [
    ['lambda', 'λ'],
    ['mu', 'μ'],
    ['nu', 'ν'],
    ['alpha', 'α'],
    ['beta', 'β'],
  ] as [keyof Material, string][])
    numRow(fsM, k, l, () => st.material[k], (v) => ((st.material[k] = v), markCustom()));
  const lawRow = document.createElement('div');
  lawRow.className = 'row wide';
  lawRow.innerHTML = '<label for="i-law" title="Формула напряжений">Формула σ</label>';
  const lawSel = document.createElement('select');
  lawSel.id = 'i-law';
  lawSel.add(new Option('из потенциала (1.9)', 'model'));
  lawSel.add(new Option('как напечатано (2.29)', 'printed229'));
  lawSel.addEventListener('change', () => {
    st.law = lawSel.value as typeof st.law;
    recompute();
  });
  lawRow.appendChild(lawSel);
  fsM.appendChild(lawRow);

  const fsPre = fieldset(c, 'Предварительная деформация (покой)');
  numRow(fsPre, 'p0', 'p0 = u1,1', () => st.F0[0], (v) => ((st.F0[0] = v), markCustom()), [-0.012, -0.00005, 0.00005], 'p0 < 0: среда сжата');
  numRow(fsPre, 'g0', 'g0 = u2,1', () => st.F0[1], (v) => ((st.F0[1] = v), markCustom()), [-0.01, 0.01, 0.0005]);
  numRow(fsPre, 's0', 's0 = u3,1', () => st.F0[2], (v) => ((st.F0[2] = v), markCustom()), [0, 0.01, 0.0005]);

  const fsB = fieldset(c, 'Удар по границе x1 = 0');
  const seg = document.createElement('div');
  seg.className = 'seg';
  for (const [t, l] of [
    ['velocity', 'Нормальный (§2.3)'],
    ['strain', 'Косой (§2.4)'],
  ] as const) {
    const b = document.createElement('button');
    b.textContent = l;
    b.dataset.bc = t;
    b.id = `i-bc-${t}`;
    b.addEventListener('click', () => {
      st.bcType = t;
      st.sweepParam = t === 'velocity' ? 'd' : 'g';
      autoMag = true;
      markCustom();
      syncBc();
      buildSweepControls();
      recompute();
    });
    segBtns.push(b);
    seg.appendChild(b);
  }
  fsB.appendChild(seg);
  bcBox = document.createElement('div');
  fsB.appendChild(bcBox);

  const fsV = fieldset(c, 'Отображение');
  const fr = document.createElement('div');
  fr.className = 'row wide';
  fr.innerHTML = '<label for="i-field">Цвет</label>';
  const fsel = document.createElement('select');
  fsel.id = 'i-field';
  for (const [k, v] of Object.entries(FIELD_LABELS)) fsel.add(new Option(v, k));
  fsel.addEventListener('change', () => {
    st.field = fsel.value as FieldKey;
    draw();
  });
  fr.appendChild(fsel);
  fsV.appendChild(fr);
  numRow(fsV, 'mag', 'Увеличение', () => st.magnify, (v) => ((st.magnify = v), draw()), [1, 2000, 1], 'перемещения умножаются на этот множитель (подбирается автоматически при смене случая)', false);
  numRow(fsV, 'tau', 'Время t/T', () => st.tau, (v) => ((st.tau = v), (st.playing = false), syncPlay(), draw()), [0, 1, 0.005], 'T — момент, когда быстрейший фронт доходит до правого края', false);
  numRow(fsV, 'speed', 'Скорость', () => st.speed, (v) => (st.speed = v), [0.02, 0.5, 0.01], 'скорость анимации', false);
  const chk = document.createElement('label');
  chk.className = 'check';
  const cb = document.createElement('input');
  cb.type = 'checkbox';
  cb.id = 'i-delta';
  cb.addEventListener('change', () => {
    st.delta = cb.checked;
    draw();
  });
  chk.append(cb, 'приращение относительно зоны 0');
  fsV.appendChild(chk);
  const btn = document.createElement('button');
  btn.id = 'i-play';
  btn.className = 'primary';
  btn.addEventListener('click', () => {
    st.playing = !st.playing;
    if (st.playing && st.tau >= 1) st.tau = 0;
    syncPlay();
    if (st.playing) requestAnimationFrame(loop);
  });
  fsV.appendChild(btn);
  const hint = document.createElement('p');
  hint.className = 'hint';
  hint.textContent =
    'Перемещения зависят только от x1 и t, поэтому показаны два сечения: в плоскости x1x2 видны u1 и u2, в плоскости x1x3 — u1 и u3. Горизонтальные материальные линии изламываются на фронтах волн.';
  fsV.appendChild(hint);
  syncPlay();
}

function syncPlay() {
  const b = document.getElementById('i-play');
  if (b) b.textContent = st.playing ? 'Пауза' : 'Пуск';
  const t = inputs.get('tau');
  if (t) {
    t.num.value = st.tau.toFixed(3);
    if (t.range) t.range.value = String(st.tau);
  }
}

function syncBc() {
  for (const b of segBtns) b.setAttribute('aria-pressed', String(b.dataset.bc === st.bcType));
  bcBox.innerHTML = '';
  for (const k of ['d', 'pb', 'gb', 'sb']) inputs.delete(k);
  if (st.bcType === 'velocity') {
    numRow(bcBox, 'd', 'd* = v1*/c1', () => st.d, (v) => ((st.d = v), markCustom()), [0.0002, 0.012, 0.0001], 'скорость движения границы вглубь, в единицах c1');
  } else {
    numRow(bcBox, 'pb', 'p*', () => st.Fb[0], (v) => ((st.Fb[0] = v), markCustom()), [-0.02, -0.0005, 0.0005], 'деформация u1,1 на границе (< p0: дополнительное сжатие)');
    numRow(bcBox, 'gb', 'g*', () => st.Fb[1], (v) => ((st.Fb[1] = v), markCustom()), [-0.012, 0.012, 0.0005], 'сдвиг u2,1 на границе');
    numRow(bcBox, 'sb', 's*', () => st.Fb[2], (v) => ((st.Fb[2] = v), markCustom()), [0, 0.012, 0.0005], 'сдвиг u3,1 на границе');
  }
}

function syncAll() {
  const vals: Record<string, number> = { ...st.material, p0: st.F0[0], g0: st.F0[1], s0: st.F0[2] };
  for (const [k, v] of Object.entries(vals)) {
    const i = inputs.get(k);
    if (!i) continue;
    i.num.value = String(v);
    if (i.range) i.range.value = String(v);
  }
  syncBc();
}

function applyPreset(name: string) {
  const p = IMPACT_PRESETS[name];
  st.preset = name;
  st.material = { ...p.material };
  st.F0 = [...p.F0] as Vec3;
  st.bcType = p.bc.type;
  if (p.bc.type === 'velocity') st.d = p.bc.d;
  else st.Fb = [...p.bc.F] as Vec3;
  st.sweepParam = p.bc.type === 'velocity' ? 'd' : 'g';
  presetSel.value = name;
  autoMag = true;
  syncAll();
  buildSweepControls();
  recompute();
}

// ---------------------------------------------------------------- расчёт

function problem() {
  return {
    material: normalize(st.material),
    F0: [...st.F0] as Vec3,
    bc: st.bcType === 'velocity' ? ({ type: 'velocity', d: st.d } as const) : ({ type: 'strain', F: [...st.Fb] as Vec3 } as const),
    stressLaw: st.law,
  };
}

/** Подобрать увеличение заново (при смене расчётного случая или типа удара). */
let autoMag = true;

function recompute() {
  st.error = '';
  st.sol = null;
  try {
    if (!(st.F0[0] < 0)) throw new Error('Предварительная деформация должна быть сжатием: p0 < 0.');
    st.sol = solveImpact(problem());
    if (!st.sol) st.error = 'Автомодельное решение не найдено (ни одна комбинация ударных и простых волн не сошлась).';
  } catch (e) {
    st.error = (e as Error).message;
  }
  st.wTable = st.sol ? displacementTable(st.sol) : null;
  if (st.sol && autoMag) {
    // Увеличение по умолчанию: самый большой скачок деформации ×K ≈ 0.25
    const z0 = st.sol.zones[0].F;
    const jump = Math.max(...st.sol.zones.flatMap((z) => z.F.map((v, k) => Math.abs(v - z0[k]))), 1e-6);
    st.magnify = Math.max(1, Math.min(2000, Math.round(0.25 / jump)));
    const i = inputs.get('mag');
    if (i) {
      i.num.value = String(st.magnify);
      if (i.range) i.range.value = String(st.magnify);
    }
    autoMag = false;
  }
  renderStatus();
  renderResults();
  renderProfiles();
  scheduleSweep();
  draw();
}

// ---------------------------------------------------------------- вывод

function renderStatus() {
  const box = $('i-status');
  box.innerHTML = '';
  const add = (cls: string, text: string) => {
    const d = document.createElement('div');
    d.className = `msg ${cls}`;
    d.textContent = text;
    box.appendChild(d);
  };
  if (st.preset !== 'custom') add('info', IMPACT_PRESETS[st.preset].note);
  if (st.law === 'printed229') add('warn', 'Используется формула (2.29) в напечатанном виде — она расходится с потенциалом (1.9) в слагаемом с β (см. docs/VERIFICATION.md).');
  if (st.error) add('err', st.error);
  for (const w of st.sol?.warnings ?? []) add('warn', w);
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

const KIND: Record<string, string> = { shock: 'ударная', fan: 'простая (веер)', rotation: 'нейтральная (поворот сдвига)' };
const FAM: Record<string, string> = { fast: 'квазипродольная', slow: 'квазипоперечная', neutral: 'нейтральная' };

function renderResults() {
  const box = $('i-results');
  box.innerHTML = '';
  box.className = 'tables';
  const s = st.sol;
  if (!s) return;
  box.appendChild(
    table(
      ['Волна', 'Тип', 'Семейство', 'ξ = G/c1', 'Диссипация D', 'Лакс: c впереди · G · c позади'],
      s.waves.map((w) => [
        w.label,
        KIND[w.kind],
        FAM[w.family],
        w.kind === 'fan' ? `${w.xiFront.toFixed(4)} … ${w.xiBack.toFixed(4)}` : w.xiFront.toFixed(5),
        w.kind === 'fan' ? '—' : { v: w.dissipation.toExponential(2), bad: w.kind === 'shock' && w.dissipation < 0 },
        w.lax ? { v: `${w.lax.ahead.toFixed(4)} · ${w.xiFront.toFixed(4)} · ${w.lax.behind.toFixed(4)}`, bad: !w.lax.ok } : '—',
      ]),
    ),
  );
  box.appendChild(
    table(
      ['Зона', 'p = u1,1', 'g = u2,1', 's = u3,1', 'd = v1/c1', 'f = v2/c1', 'h = v3/c1', 'σ̃11', 'σ̃21', 'σ̃31'],
      s.zones.map((z) => [z.name, ...[...z.F, ...z.v, ...z.T].map((v) => fmt(v, 4))]),
    ),
  );
}

function sampleXi(s: ImpactSolution, xmax: number, n = 320): number[] {
  const xs: number[] = [];
  for (let i = 0; i <= n; i++) xs.push((xmax * i) / n);
  for (const w of s.waves) {
    xs.push(w.xiFront + 1e-9, w.xiFront - 1e-9, w.xiBack + 1e-9, w.xiBack - 1e-9);
    if (w.profile) for (const q of w.profile) xs.push(q.xi);
  }
  return xs.filter((x) => x >= 0 && x <= xmax).sort((a, b) => a - b);
}

function renderProfiles() {
  const box = $('i-profiles');
  box.innerHTML = '';
  const s = st.sol;
  if (!s) return;
  const xmax = s.waves[0].xiFront * 1.25;
  const xs = sampleXi(s, xmax);
  const pts = xs.map((x) => {
    const q = profileAt(s, x);
    const T = traction(q.F);
    return { F: q.F, v: q.v, T };
  });
  const marks = s.waves.map((w) => ({ x: w.xiFront, label: w.kind === 'fan' ? 'ξ⁺' : w.label, fan: w.kind === 'fan' ? ([w.xiFront, w.xiBack] as [number, number]) : undefined }));
  for (const w of s.waves) if (w.kind === 'fan') marks.push({ x: w.xiBack, label: 'ξ⁻', fan: undefined });
  const color = css('--accent');
  const items: [string, (q: (typeof pts)[number]) => number][] = [
    ['p = u1,1', (q) => q.F[0]],
    ['g = u2,1', (q) => q.F[1]],
    ['s = u3,1', (q) => q.F[2]],
    ['d = v1/c1', (q) => q.v[0]],
    ['f = v2/c1', (q) => q.v[1]],
    ['h = v3/c1', (q) => q.v[2]],
    ['σ11 / (λ+2μ)', (q) => q.T[0]],
    ['σ21 / (λ+2μ)', (q) => q.T[1]],
    ['σ31 / (λ+2μ)', (q) => q.T[2]],
  ];
  for (const [t, f] of items) box.appendChild(profileChart(t, xs, pts.map(f), marks, color));
}

function traction(F: Vec3): Vec3 {
  const s = st.sol!;
  const z = s.zones.find((zz) => zz.F === F);
  if (z) return z.T;
  // внутри веера — пересчёт по закону, использованному в решении
  return lawFor(s)(F);
}

function lawFor(s: ImpactSolution) {
  const m = s.problem.material;
  return s.problem.stressLaw === 'printed229' ? (F: Vec3) => traction229(m, F, 1) : (F: Vec3) => tractionModel(m, F);
}

// ---------------------------------------------------------------- зависимость от воздействия

let sweepTimer = 0;
function scheduleSweep() {
  clearTimeout(sweepTimer);
  sweepTimer = window.setTimeout(renderSweep, 150);
}

function buildSweepControls() {
  const box = $('i-sweep-controls');
  box.innerHTML = '<label for="i-sweep-param">Параметр по оси X:</label>';
  const sel = document.createElement('select');
  sel.id = 'i-sweep-param';
  const opts: [typeof st.sweepParam, string][] =
    st.bcType === 'velocity'
      ? [['d', 'скорость границы d* = v1*/c1']]
      : [
          ['g', 'сдвиг на границе g*'],
          ['s', 'сдвиг на границе s*'],
          ['p', 'продольная деформация на границе p*'],
        ];
  for (const [v, l] of opts) sel.add(new Option(l, v));
  sel.value = st.sweepParam;
  sel.addEventListener('change', () => {
    st.sweepParam = sel.value as typeof st.sweepParam;
    renderSweep();
  });
  box.appendChild(sel);
  const note = document.createElement('span');
  note.className = 'hint';
  note.textContent = 'Точечный пунктир — решение с недопустимой волной; для веера показаны оба фронта.';
  box.appendChild(note);
}

function renderSweep() {
  const box = $('i-sweep');
  box.innerHTML = '';
  if (!st.sol) return;
  const base = problem();
  let values: number[] = [];
  let cur = 0;
  let make: (v: number) => ReturnType<typeof problem>;
  let xLabel = '';
  const lin = (a: number, b: number, n: number) => Array.from({ length: n + 1 }, (_, i) => +(a + ((b - a) * i) / n).toPrecision(8));
  if (st.sweepParam === 'd') {
    values = lin(0.0005, 0.012, 23);
    cur = st.d;
    make = (v) => ({ ...base, bc: { type: 'velocity', d: v } });
    xLabel = 'd* = v1*/c1';
  } else {
    const k = { p: 0, g: 1, s: 2 }[st.sweepParam];
    values = st.sweepParam === 'p' ? lin(-0.02, Math.min(-0.001, st.F0[0] - 0.0005), 24) : lin(st.sweepParam === 's' ? 0.0005 : -0.012, 0.012, 24);
    cur = st.Fb[k];
    make = (v) => {
      const F = [...st.Fb] as Vec3;
      F[k] = v;
      return { ...base, bc: { type: 'strain', F } };
    };
    xLabel = ['p*', 'g*', 's*'][k];
  }
  const rows = values.map((v) => {
    let s: ImpactSolution | null = null;
    try {
      s = solveImpact(make(v));
    } catch {
      s = null;
    }
    return { v, s };
  });
  const ok = rows.filter((r) => r.s);
  const get = (s: ImpactSolution, pred: (w: ImpactSolution['waves'][number]) => boolean, edge: 'xiFront' | 'xiBack') => {
    const w = s.waves.find(pred);
    return w ? w[edge] : NaN;
  };
  const series = (label: string, color: string, pred: (w: ImpactSolution['waves'][number]) => boolean, edge: 'xiFront' | 'xiBack' = 'xiFront', dash?: string): Series => {
    const pts = ok.map((r) => ({ x: r.v, y: get(r.s!, pred, edge), good: r.s!.warnings.length === 0 })).filter((p) => Number.isFinite(p.y));
    return { label, color, dash, x: pts.map((p) => p.x), y: pts.map((p) => p.y), ok: pts.map((p) => p.good) };
  };
  const c1 = css('--front1');
  const c2 = css('--front2');
  const c3 = css('--accent');
  const fast = [series('Σ1 (квазипродольная)', c1, (w) => w.family === 'fast')];
  const shear = [
    series('Σ2 / передний фронт веера', c2, (w) => w.family === 'slow'),
    series('задний фронт веера', c2, (w) => w.family === 'slow' && w.kind === 'fan', 'xiBack', '4 3'),
    series('Σ3 (нейтральная)', c3, (w) => w.family === 'neutral'),
  ].filter((s) => s.x.length > 0);
  const s = st.sol;
  const fastNow = s.waves.find((w) => w.family === 'fast')!.xiFront;
  box.appendChild(
    lineChart({ xLabel, yLabel: 'G1 / c1', series: fast, markers: [{ x: cur, y: fastNow, color: c1, label: 'текущее решение', shape: 'ring' }], vline: cur, allowNegative: true, height: 260, width: 520 }),
  );
  const shearMarkers = s.waves.filter((w) => w.family !== 'fast').map((w) => ({ x: cur, y: w.xiFront, color: w.family === 'neutral' ? c3 : c2, label: 'текущее решение', shape: 'ring' as const }));
  box.appendChild(lineChart({ xLabel, yLabel: 'G / c1 (сдвиговые волны)', series: shear, markers: shearMarkers, vline: cur, allowNegative: true, height: 260, width: 520 }));
}

// ---------------------------------------------------------------- сцена

function fieldAt(s: ImpactSolution, xi: number): number {
  const q = profileAt(s, xi);
  const T = traction(q.F);
  const val = (F: Vec3, v: Vec3, TT: Vec3) => {
    switch (st.field) {
      case 'p':
        return F[0];
      case 'g':
        return F[1];
      case 's':
        return F[2];
      case 'shear':
        return Math.hypot(F[1], F[2]);
      case 's11':
        return TT[0];
      case 's21':
        return TT[1];
      case 's31':
        return TT[2];
      case 'speed':
        return Math.hypot(...v);
      default:
        return 0;
    }
  };
  const x = val(q.F, q.v, T);
  if (!st.delta) return x;
  const z0 = s.zones[0];
  return x - val(z0.F, z0.v, z0.T);
}

/** Номер «участка» для раскраски зон: 2i — зона i, 2i+1 — внутри веера i-й волны. */
function zoneSlot(s: ImpactSolution, xi: number): number {
  for (let i = 0; i < s.waves.length; i++) {
    const w = s.waves[i];
    if (xi > w.xiFront) return 2 * i;
    if (w.kind === 'fan' && xi >= w.xiBack) return 2 * i + 1;
  }
  return 2 * s.waves.length;
}

function draw() {
  drawScene();
  drawXt();
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

function drawScene() {
  const canvas = $<HTMLCanvasElement>('i-scene');
  const box3d = $('i-3d');
  const is3d = st.view === '3d';
  canvas.hidden = is3d;
  box3d.hidden = !is3d;
  const S = is3d ? null : setup(canvas);
  const leg = $('i-legend');
  leg.innerHTML = '';
  if (!is3d && !S) return;
  const { ctx, w, h } = S ?? { ctx: null as unknown as CanvasRenderingContext2D, w: 0, h: 0 };
  const s = st.sol;
  const W = st.wTable;
  if (!s || !W) return;
  const xi1 = s.waves[0].xiFront;
  const tau = Math.max(st.tau, 1e-6);
  const K = st.magnify;
  /** Перемещение (в единицах длины окна) материальной точки X̂ в момент τ. */
  const disp = (X: number): Vec3 => {
    const wv = W((X * xi1) / tau);
    return [(tau / xi1) * wv[0], (tau / xi1) * wv[1], (tau / xi1) * wv[2]];
  };
  const xiOf = (X: number) => (X * xi1) / tau;

  // значения поля для легенды и цветовой шкалы
  let vmax = 0;
  if (st.field !== 'zones') {
    for (let i = 0; i <= 200; i++) vmax = Math.max(vmax, Math.abs(fieldAt(s, (1.3 * xi1 * i) / 200)));
    vmax = Math.max(vmax, Math.abs(fieldAt(s, 0)), Math.abs(fieldAt(s, 2 * xi1)));
    vmax ||= 1;
  }
  const tints = [css('--zone0'), css('--zone1'), css('--zone2'), css('--zone1'), css('--zone2'), css('--zone1'), css('--zone2')];
  const colorAt = (xi: number): string => {
    if (st.field === 'zones') {
      const slot = zoneSlot(s, xi);
      if (slot % 2 === 1) return css('--zone1');
      return tints[slot / 2] ?? tints[1];
    }
    const t = fieldAt(s, xi) / vmax;
    const a = 0.1 + 0.75 * Math.abs(t);
    return t >= 0 ? `rgba(214, 72, 52, ${a})` : `rgba(48, 104, 214, ${a})`;
  };

  const panels: { comp: 1 | 2; title: string }[] = [
    { comp: 1, title: 'Сечение x3 = const: смещения u1 (→) и u2 (↑)' },
    { comp: 2, title: 'Сечение x2 = const: смещения u1 (→) и u3 (↑)' },
  ];
  const gap = 10;
  const ph = (h - gap) / 2;
  const left = 70;
  const right = 14;
  const sc = (w - left - right) / 1.12;
  const Xmax = 1.12;

  const breaks: number[] = [];
  for (const wv of s.waves) breaks.push((wv.xiFront / xi1) * tau, (wv.xiBack / xi1) * tau);
  const Xs: number[] = [];
  for (let i = 0; i <= 420; i++) Xs.push((Xmax * i) / 420);
  for (const b of breaks) if (b > 0 && b < Xmax) Xs.push(b - 1e-9, b + 1e-9);
  for (const wv of s.waves) if (wv.profile) for (const q of wv.profile) Xs.push((q.xi / xi1) * tau);
  Xs.sort((a, b) => a - b);
  const D = Xs.map(disp);

  if (is3d) render3d(s, Xs, D, disp, Xmax, (X) => colorAt(xiOf(X)), xi1, tau);
  else panels.forEach((P, pi) => {
    const top = pi * (ph + gap);
    const titleH = 22;
    const oy = top + titleH + (ph - titleH) / 2;
    const half = (ph - titleH) / 2 / sc - 0.02;
    const map = (X: number, y: number, d: Vec3) => [left + (X + K * d[0]) * sc, oy - (y + K * d[P.comp]) * sc] as const;
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, top, w, ph);
    ctx.clip();
    // полосы
    for (let i = 0; i < Xs.length - 1; i++) {
      const xa = Xs[i];
      const xb = Xs[i + 1];
      if (xb - xa < 1e-12) continue;
      ctx.fillStyle = colorAt(xiOf(0.5 * (xa + xb)));
      const a1 = map(xa, half, D[i]);
      const b1 = map(xb, half, D[i + 1]);
      const b2 = map(xb, -half, D[i + 1]);
      const a2 = map(xa, -half, D[i]);
      ctx.beginPath();
      ctx.moveTo(a1[0], a1[1]);
      ctx.lineTo(b1[0] + 0.6, b1[1]);
      ctx.lineTo(b2[0] + 0.6, b2[1]);
      ctx.lineTo(a2[0], a2[1]);
      ctx.closePath();
      ctx.fill();
    }
    // материальная сетка
    const hg = 0.05;
    ctx.strokeStyle = css('--grid');
    ctx.lineWidth = 1;
    for (let X = hg; X < Xmax; X += hg) {
      const d = disp(X);
      const a = map(X, half, d);
      const b = map(X, -half, d);
      ctx.beginPath();
      ctx.moveTo(a[0], a[1]);
      ctx.lineTo(b[0], b[1]);
      ctx.stroke();
    }
    for (let y = -Math.floor(half / hg) * hg; y <= half + 1e-9; y += hg) {
      ctx.beginPath();
      Xs.forEach((X, i) => {
        const p = map(X, y, D[i]);
        if (i === 0) ctx.moveTo(p[0], p[1]);
        else ctx.lineTo(p[0], p[1]);
      });
      ctx.stroke();
    }
    // фронты
    s.waves.forEach((wv) => {
      const color = wv.family === 'fast' ? css('--front1') : wv.family === 'slow' ? css('--front2') : css('--accent');
      for (const xi of wv.kind === 'fan' ? [wv.xiFront, wv.xiBack] : [wv.xiFront]) {
        const X = (xi / xi1) * tau;
        if (X > Xmax) continue;
        const d = disp(X);
        const a = map(X, half, d);
        const b = map(X, -half, d);
        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        ctx.setLineDash(wv.kind === 'rotation' ? [3, 4] : [7, 5]);
        ctx.beginPath();
        ctx.moveTo(a[0], a[1]);
        ctx.lineTo(b[0], b[1]);
        ctx.stroke();
        ctx.setLineDash([]);
      }
    });
    // граница
    const d0 = disp(0);
    const g1 = map(0, half, d0);
    const g2 = map(0, -half, d0);
    ctx.strokeStyle = css('--surface');
    ctx.lineWidth = 3;
    ctx.beginPath();
    ctx.moveTo(g1[0], g1[1]);
    ctx.lineTo(g2[0], g2[1]);
    ctx.stroke();
    // нагрузка
    const Tb = s.zones.at(-1)!.T;
    let dir: [number, number];
    if (s.problem.bc.type === 'velocity') dir = [1, 0];
    else {
      const tx = -Tb[0];
      const ty = -Tb[P.comp];
      const n = Math.hypot(tx, ty) || 1;
      dir = [tx / n, ty / n];
    }
    ctx.strokeStyle = css('--accent');
    ctx.fillStyle = css('--accent');
    ctx.lineWidth = 2;
    for (let k = -2; k <= 2; k++) {
      const yy = (k / 2.6) * half;
      const [ex, ey] = map(0, yy, d0);
      const L = 40;
      const sx = ex - 6 - dir[0] * L;
      const sy = ey + dir[1] * L;
      const tx = ex - 6;
      ctx.beginPath();
      ctx.moveTo(sx, sy);
      ctx.lineTo(tx, ey);
      ctx.stroke();
      const ang = Math.atan2(ey - sy, tx - sx);
      ctx.beginPath();
      ctx.moveTo(tx, ey);
      ctx.lineTo(tx - 8 * Math.cos(ang - 0.4), ey - 8 * Math.sin(ang - 0.4));
      ctx.lineTo(tx - 8 * Math.cos(ang + 0.4), ey - 8 * Math.sin(ang + 0.4));
      ctx.closePath();
      ctx.fill();
    }
    ctx.restore();
    // подписи
    ctx.font = '600 12px system-ui, sans-serif';
    ctx.fillStyle = css('--text');
    ctx.fillText(P.title, 10, top + 15);
    ctx.font = '600 12px system-ui, sans-serif';
    const used: number[] = [];
    s.waves.forEach((wv) => {
      const X = (wv.xiFront / xi1) * tau;
      if (X > Xmax) return;
      let px = left + (X + K * disp(X)[0]) * sc;
      while (used.some((u) => Math.abs(u - px) < 36)) px -= 36;
      used.push(px);
      ctx.fillStyle = wv.family === 'fast' ? css('--front1') : wv.family === 'slow' ? css('--front2') : css('--accent');
      ctx.fillText(wv.kind === 'fan' ? 'веер' : wv.label, px - 8, top + 34);
    });
  });
  if (!is3d) {
    ctx.font = '12px system-ui, sans-serif';
    ctx.fillStyle = css('--muted');
    ctx.textAlign = 'right';
    ctx.fillText(`перемещения ×${K}, t/T = ${st.tau.toFixed(2)}`, w - 10, 15);
    ctx.textAlign = 'left';
  }

  // легенда
  const lg = (html: string) => {
    const d = document.createElement('div');
    d.innerHTML = html;
    leg.appendChild(d);
  };
  if (st.field === 'zones') {
    lg('<div style="margin-bottom:4px">Зоны (0 — невозмущённая)</div>');
    s.zones.forEach((z, i) => lg(`<span class="sw" style="background:${tints[i] ?? tints[1]}"></span>зона ${z.name}`));
    if (s.waves.some((wv) => wv.kind === 'fan')) lg(`<span class="sw" style="background:${css('--zone1')}"></span>простая волна`);
  } else {
    lg(`<div style="margin-bottom:4px">${st.delta ? 'Δ ' : ''}${FIELD_LABELS[st.field]}</div>`);
    lg(`<span class="sw" style="background:rgba(48,104,214,0.85)"></span>−${fmt(vmax, 3)}`);
    lg(`<span class="sw" style="background:rgba(214,72,52,0.85)"></span>+${fmt(vmax, 3)}`);
  }
}

// ---------------------------------------------------------------- 3D

let view3d: typeof import('./impact3d.ts') | null = null;
let view3dLoading = false;

function render3d(s: ImpactSolution, Xs: number[], D: Vec3[], disp: (X: number) => Vec3, Xmax: number, colorAt: (X: number) => string, xi1: number, tau: number) {
  if (!view3d) {
    if (!view3dLoading) {
      view3dLoading = true;
      import('./impact3d.ts').then((m) => {
        view3d = m;
        draw();
      });
    }
    return;
  }
  const Tb = s.zones.at(-1)!.T;
  const load: Vec3 = s.problem.bc.type === 'velocity' ? [1, 0, 0] : [-Tb[0], -Tb[1], -Tb[2]];
  const fronts: { X: number; color: string; label: string; dashed: boolean }[] = [];
  for (const w of s.waves) {
    const color = w.family === 'fast' ? css('--front1') : w.family === 'slow' ? css('--front2') : css('--accent');
    if (w.kind === 'fan') {
      fronts.push({ X: (w.xiFront / xi1) * tau, color, label: 'веер ξ⁺', dashed: true });
      fronts.push({ X: (w.xiBack / xi1) * tau, color, label: 'ξ⁻', dashed: true });
    } else fronts.push({ X: (w.xiFront / xi1) * tau, color, label: w.label, dashed: w.kind === 'rotation' });
  }
  view3d.update($('i-3d'), {
    Xs,
    D,
    disp,
    K: st.magnify,
    Xmax,
    colorAt,
    fronts,
    load,
    loadLabel: s.problem.bc.type === 'velocity' ? 'v* (нормальный удар)' : 'нагрузка σ*',
    theme: { bg: css('--panel'), grid: css('--grid').startsWith('rgba') ? css('--muted') : css('--grid'), text: css('--text'), muted: css('--muted'), surface: css('--surface'), accent: css('--accent') },
    caption: `перемещения ×${st.magnify}, t/T = ${st.tau.toFixed(2)} · ось бруска выделена · вращение — мышью, масштаб — колесом`,
  });
}

function buildViewToggle() {
  const btns = ['i-view-2d', 'i-view-3d'].map((id) => $<HTMLButtonElement>(id));
  const reset = $<HTMLButtonElement>('i-view-reset');
  for (const b of btns)
    b.addEventListener('click', () => {
      st.view = b.dataset.view as typeof st.view;
      for (const x of btns) x.setAttribute('aria-pressed', String(x === b));
      reset.hidden = st.view !== '3d';
      draw();
    });
  reset.addEventListener('click', () => view3d?.resetView());
}

function drawXt() {
  const canvas = $<HTMLCanvasElement>('i-xt');
  const S = setup(canvas);
  if (!S) return;
  const { ctx, w, h } = S;
  const s = st.sol;
  if (!s) return;
  const xi1 = s.waves[0].xiFront;
  const m = { l: 40, r: 14, t: 44, b: 34 };
  const X = (x: number) => m.l + (x / 1.12) * (w - m.l - m.r);
  const Y = (t: number) => h - m.b - t * (h - m.t - m.b);
  ctx.font = '600 12px system-ui, sans-serif';
  ctx.fillStyle = css('--text');
  ctx.fillText('Диаграмма x–t: волны — лучи x = G t', 10, 16);
  ctx.strokeStyle = css('--border');
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(X(0), Y(0));
  ctx.lineTo(X(1.12), Y(0));
  ctx.moveTo(X(0), Y(0));
  ctx.lineTo(X(0), Y(1));
  ctx.stroke();
  ctx.font = '11px system-ui, sans-serif';
  ctx.fillStyle = css('--muted');
  ctx.fillText('x1', X(1.12) - 12, Y(0) + 16);
  ctx.fillText('t', X(0) - 14, Y(1) + 4);
  ctx.fillText('граница', X(0) + 4, Y(0) + 26);
  // вееры
  for (const wv of s.waves)
    if (wv.kind === 'fan') {
      ctx.fillStyle = css('--front2');
      ctx.globalAlpha = 0.18;
      ctx.beginPath();
      ctx.moveTo(X(0), Y(0));
      ctx.lineTo(X(wv.xiFront / xi1), Y(1));
      ctx.lineTo(X(wv.xiBack / xi1), Y(1));
      ctx.closePath();
      ctx.fill();
      ctx.globalAlpha = 1;
    }
  const used: number[] = [];
  for (const wv of s.waves) {
    const color = wv.family === 'fast' ? css('--front1') : wv.family === 'slow' ? css('--front2') : css('--accent');
    for (const xi of wv.kind === 'fan' ? [wv.xiFront, wv.xiBack] : [wv.xiFront]) {
      ctx.strokeStyle = color;
      ctx.lineWidth = 2;
      ctx.setLineDash(wv.kind === 'rotation' ? [3, 4] : wv.kind === 'fan' ? [] : [7, 5]);
      ctx.beginPath();
      ctx.moveTo(X(0), Y(0));
      ctx.lineTo(X(xi / xi1), Y(1));
      ctx.stroke();
      ctx.setLineDash([]);
    }
    let px = X(wv.xiFront / xi1);
    while (used.some((u) => Math.abs(u - px) < 30)) px -= 30;
    used.push(px);
    ctx.fillStyle = color;
    ctx.font = '600 11px system-ui, sans-serif';
    ctx.fillText(wv.kind === 'fan' ? 'веер' : wv.label, px - 10, Y(1) - 6);
  }
  // зоны: подписи на высоте t = 0.7 между соседними лучами
  ctx.fillStyle = css('--text');
  ctx.font = '600 12px system-ui, sans-serif';
  const tt = 0.7;
  const rayX = (xi: number) => (xi / xi1) * tt;
  s.zones.forEach((z, i) => {
    const right = i === 0 ? Math.min(1.1, rayX(s.waves[0].xiFront) + 0.25) : rayX(s.waves[i - 1].xiBack);
    const leftEdge = i === 0 ? rayX(s.waves[0].xiFront) : i < s.waves.length ? rayX(s.waves[i].xiFront) : 0;
    if ((right - leftEdge) * (w - m.l - m.r) < 14) return; // слишком узкая зона
    ctx.fillText(z.name, X(0.5 * (leftEdge + right)) - 4, Y(tt) + 4);
  });
  // текущее время
  ctx.strokeStyle = css('--muted');
  ctx.setLineDash([2, 3]);
  ctx.beginPath();
  ctx.moveTo(X(0), Y(st.tau));
  ctx.lineTo(X(1.12), Y(st.tau));
  ctx.stroke();
  ctx.setLineDash([]);
}

// ---------------------------------------------------------------- анимация

let lastT = 0;
function loop(t: number) {
  const dt = Math.min(0.05, (t - lastT) / 1000 || 0);
  lastT = t;
  if (!st.playing) return;
  st.tau += dt * st.speed;
  if (st.tau > 1.15) st.tau = 0;
  const shown = Math.min(st.tau, 1);
  const saved = st.tau;
  st.tau = shown;
  syncPlay();
  draw();
  st.tau = saved;
  requestAnimationFrame(loop);
}

buildControls();
buildViewToggle();
syncAll();
buildSweepControls();
recompute();
new ResizeObserver(() => draw()).observe($('i-scene'));
new ResizeObserver(() => st.view === '3d' && draw()).observe($('i-3d'));
new ResizeObserver(() => draw()).observe($('i-xt'));
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => recompute());
requestAnimationFrame(loop);

(window as unknown as { __impact: typeof st }).__impact = st;
