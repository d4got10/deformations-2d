/**
 * Генерирует docs/VERIFICATION_RESULTS.md: сверка расчёта с числами и графиками диссертации.
 * Запуск: npm run report
 */
import { type MovingLoadProblem, type MovingLoadSolution, solveMovingLoad, sweep } from '../src/core/movingLoad.ts';
import { c2OverC1, scaleNonlinear } from '../src/core/material.ts';
import { BASE_CASE, DISSERTATION, LOAD_SWEEP_CASE, MATERIAL_1, MATERIAL_2 } from '../src/core/presets.ts';

const out: string[] = [];
const p = (s = '') => out.push(s);
const f = (x: number, d = 4) => x.toFixed(d);
const range = (from: number, to: number, step: number) => {
  const r: number[] = [];
  const n = Math.round((to - from) / step);
  for (let i = 0; i <= n; i++) r.push(+(from + i * step).toFixed(6));
  return r;
};

p('# Результаты сверки с диссертацией (сгенерировано `npm run report`)');
p();
p(`Дата генерации: ${new Date().toISOString().slice(0, 10)}. Методика — в [VERIFICATION.md](VERIFICATION.md).`);
p();

const cases: [string, typeof MATERIAL_1, { G1: number; G2: number; theta1: number; theta2: number }][] = [
  ['Материал 1 (ν < 0), рис. 21', MATERIAL_1, DISSERTATION.material1],
  ['Материал 2 (ν > 0), рис. 23', MATERIAL_2, DISSERTATION.material2],
];

p('## 1. Базовый случай: M = 1.7, σ̃11* = −0.01, σ̃21* = 0.003, a0 = −0.005, k0 = 0.001');
p();
p('| Случай | Величина | Расчёт | Диссертация | Разница | Линейная теория |');
p('|---|---|---|---|---|---|');
const sols: MovingLoadSolution[] = [];
for (const [name, mat, ref] of cases) {
  const prob: MovingLoadProblem = { material: mat, ...BASE_CASE };
  const s = solveMovingLoad(prob)!;
  const lin = solveMovingLoad({ ...prob, material: scaleNonlinear(mat, 0) })!;
  sols.push(s);
  const rows: [string, number, number, number, number][] = [
    ['θ1, °', s.fronts[0].thetaDeg, ref.theta1, lin.fronts[0].thetaDeg, 2],
    ['θ2, °', s.fronts[1].thetaDeg, ref.theta2, lin.fronts[1].thetaDeg, 2],
    ['G1/c1', s.fronts[0].G, ref.G1, lin.fronts[0].G, 4],
    ['G2/c1', s.fronts[1].G, ref.G2, lin.fronts[1].G, 4],
  ];
  for (const [q, calc, dis, l, d] of rows)
    p(`| ${name} | ${q} | ${f(calc, d)} | ${dis} | ${calc - dis >= 0 ? '+' : ''}${f(calc - dis, d)} | ${f(l, d)} |`);
}
p();
p('Самосогласованность чисел автора: G = M·sin θ.');
for (const [name, , ref] of cases)
  p(`- ${name}: 1.7·sin ${ref.theta1}° = ${f(1.7 * Math.sin((ref.theta1 * Math.PI) / 180))} (в тексте ${ref.G1}); 1.7·sin ${ref.theta2}° = ${f(1.7 * Math.sin((ref.theta2 * Math.PI) / 180))} (в тексте ${ref.G2}).`);
p();

p('## 2. Скачки по зонам (рис. 21, материал 1)');
p();
p('Значения диссертации грубо считаны с графика (точность порядка ±0.0005).');
p();
p('| Величина | Зона 0: расчёт / рис. | Зона I: расчёт / рис. | Зона II: расчёт / рис. |');
p('|---|---|---|---|');
const s1 = sols[0];
const fig = DISSERTATION.fig21;
const val = (z: MovingLoadSolution['zones'][number], key: keyof typeof fig) =>
  key === 's11' ? z.stress.s11 : key === 's21' ? z.stress.s12 : key === 's22' ? z.stress.s22 : z.grad[key];
for (const key of Object.keys(fig) as (keyof typeof fig)[]) {
  const label = { a: 'a = u1,1', b: 'b = u1,2', k: 'k = u2,1', n: 'n = u2,2', s11: 'σ11/(λ+2μ)', s21: 'σ21/(λ+2μ)', s22: 'σ22/(λ+2μ)' }[key];
  p(`| ${label} | ${s1.zones.map((z, i) => `${f(val(z, key), 5)} / ${fig[key][i]}`).join(' | ')} |`);
}
p();

p('## 3. Диссипация на фронтах (термодинамическое условие (3.30))');
p();
p('D отнесена к (λ+2μ)·s·cos θ. Масштаб для сравнения — ξ·|[W]|.');
p();
p('| Случай | Фронт | D | ξ·|[W]| | D / ξ|[W]| | Допустим? |');
p('|---|---|---|---|---|---|');
cases.forEach(([name], ci) => {
  const s = sols[ci];
  s.fronts.forEach((fr, i) => {
    const sc = Math.abs(fr.xi * (s.zones[i].W - s.zones[i + 1].W));
    p(`| ${name} | ${fr.name} | ${fr.dissipation.toExponential(2)} | ${sc.toExponential(2)} | ${(fr.dissipation / sc).toExponential(2)} | ${fr.dissipation >= 0 ? 'да' : '**нет**'} |`);
  });
});
p();

p('## 4. Зависимость углов от числа Маха (рис. 22а, 24а)');
p();
for (const [name, mat, ref] of cases) {
  const base: MovingLoadProblem = { material: mat, ...BASE_CASE };
  const Ms = range(2.5, 1.0, -0.005);
  const sw = sweep(base, Ms, (q, M) => ({ ...q, M }));
  const firstWarn = sw.points.find((pt) => pt.solution.warnings.length > 0);
  p(`### ${name}`);
  p();
  p(`- Решение существует при M ≥ ${sw.points.at(-1)?.value} (шаг 0.005); в диссертации M_min = ${'Mmin' in ref ? (ref as { Mmin: number }).Mmin : '—'}.`);
  p(`- Первое M (при убывании), где нарушено какое-либо ограничение: ${firstWarn ? `${firstWarn.value} — ${firstWarn.solution.warnings.join(' ')}` : 'нет'}.`);
  p();
  p('| M | θ1, ° | θ2, ° | G1/c1 | G2/c1 | Линейная: θ1 / θ2 |');
  p('|---|---|---|---|---|---|');
  const c2 = c2OverC1(mat);
  for (const pt of sw.points.filter((q) => [2.5, 2.2, 2, 1.7, 1.5, 1.4, 1.3, 1.2, 1.11, 1.05, 1.01].includes(q.value))) {
    const fr = pt.solution.fronts;
    const l1 = (Math.asin(1 / pt.value) * 180) / Math.PI;
    const l2 = (Math.asin(c2 / pt.value) * 180) / Math.PI;
    p(`| ${pt.value} | ${f(fr[0].thetaDeg, 2)} | ${f(fr[1].thetaDeg, 2)} | ${f(fr[0].G)} | ${f(fr[1].G)} | ${f(l1, 2)} / ${f(l2, 2)} |`);
  }
  p();
}
p('Значения с рис. 22а (материал 1, считаны с графика и пометок на полях):');
for (const r of DISSERTATION.fig22a) p(`- M = ${r.M}: θ1 ≈ ${r.theta1}°, θ2 ≈ ${r.theta2}°.`);
p();

p('## 5. Зависимость углов от нагрузки при M = 1.4, a0 = −0.003, k0 = 0.001 (рис. 22б,в и 24б,в)');
p();
for (const [name, mat] of cases) {
  const base: MovingLoadProblem = { material: mat, ...LOAD_SWEEP_CASE };
  // б) |σ11*| растёт; для материала 1 σ21* = 0.003 (как в рис. 22), для материала 2 — σ21* = 0.3|σ11*| (текст к рис. 24б)
  const isM2 = mat === MATERIAL_2;
  const swB = sweep(base, range(0.006, 0.07, 0.004), (q, v) => ({ ...q, s11: -v, s21: isM2 ? 0.3 * v : q.s21 }));
  const swC = sweep(base, range(0.0005, 0.007, 0.0005), (q, v) => ({ ...q, s11: -0.01, s21: v }));
  const fmt = (sw: ReturnType<typeof sweep>) =>
    sw.points.map((pt) => `${pt.value}: ${f(pt.solution.fronts[0].thetaDeg, 3)} / ${f(pt.solution.fronts[1].thetaDeg, 3)}`).join('; ');
  const trend = (sw: ReturnType<typeof sweep>, i: 0 | 1) => {
    const a = sw.points[0]?.solution.fronts[i].thetaDeg;
    const b = sw.points.at(-1)?.solution.fronts[i].thetaDeg;
    return a === undefined || b === undefined ? '—' : b > a ? 'растёт' : 'убывает';
  };
  p(`### ${name}`);
  p();
  p(`- б) |σ̃11*| ↑${isM2 ? ' (σ̃21* = 0.3|σ̃11*|)' : ' (σ̃21* = 0.003)'}: θ1 ${trend(swB, 0)}, θ2 ${trend(swB, 1)}. ${swB.failedAt !== null ? `Нет решения при ${swB.failedAt}.` : ''}`);
  p(`  - ${fmt(swB)}`);
  p(`- в) σ̃21* ↑ (σ̃11* = −0.01): θ1 ${trend(swC, 0)}, θ2 ${trend(swC, 1)}. ${swC.failedAt !== null ? `Нет решения при ${swC.failedAt}.` : ''}`);
  p(`  - ${fmt(swC)}`);
  p();
}
p('В диссертации: для материала 1 рост |σ̃11*| уменьшает θ1, θ2, рост σ̃21* увеличивает (рис. 22б,в; θ1 ≈ 46.05…46.25°, θ2 ≈ 23.35…23.6°). Для материала 2 — наоборот (рис. 24б,в).');

console.log(out.join('\n'));
