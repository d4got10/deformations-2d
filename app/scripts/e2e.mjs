/**
 * Сквозная проверка интерфейса: собирает приложение, открывает в Chromium,
 * проверяет отсутствие ошибок в консоли, совпадение показанных углов с расчётным ядром
 * и делает скриншоты. Запуск: npm run build && npm run e2e [папка_для_скриншотов]
 */
import { spawn } from 'node:child_process';
import { mkdirSync } from 'node:fs';
import { createRequire } from 'node:module';

const require = createRequire(import.meta.url);
const { chromium } = require('playwright');

const outDir = process.argv[2] ?? 'e2e-screenshots';
mkdirSync(outDir, { recursive: true });
const port = 4179;
const server = spawn('npx', ['vite', 'preview', '--port', String(port), '--strictPort'], { stdio: 'pipe', detached: true });
await new Promise((resolve, reject) => {
  server.stdout.on('data', (d) => String(d).includes(String(port)) && resolve());
  server.on('exit', reject);
  setTimeout(resolve, 5000);
});

const errors = [];
const fail = (msg) => {
  errors.push(msg);
  console.error('FAIL:', msg);
};
const browser = await chromium.launch();
try {
  for (const scheme of ['light', 'dark']) {
    const page = await browser.newPage({ viewport: { width: 1280, height: 1500 }, colorScheme: scheme });
    page.on('console', (m) => m.type() === 'error' && fail(`console: ${m.text()}`));
    page.on('pageerror', (e) => fail(`pageerror: ${e.message}`));
    await page.goto(`http://localhost:${port}/`);

    // ---------- Глава 2 (вкладка по умолчанию) ----------
    await page.waitForFunction(() => window.__impact?.sol);
    const ch2 = await page.evaluate(() => window.__impact.sol.waves.map((w) => w.xiFront));
    if (Math.abs(ch2[0] - 1.00697) > 1e-4 || Math.abs(ch2[1] - 0.56712) > 1e-4) fail(`§2.3, материал 1: ${ch2}`);
    if ((await page.locator('#i-profiles .chart').count()) !== 9) fail('профилей по ξ должно быть 9');
    await page.click('#i-play');
    await page.fill('#i-in-tau', '0.8');
    await page.dispatchEvent('#i-in-tau', 'change');
    await page.waitForSelector('#i-sweep svg polyline');
    await page.screenshot({ path: `${outDir}/ch2-normal-${scheme}.png`, fullPage: true });
    // ненапряжённая среда: решение главы 2 существует (по волновым кривым)
    for (const [key, val] of [['p0', '0'], ['s0', '0']]) {
      await page.fill(`#i-in-${key}`, val);
      await page.dispatchEvent(`#i-in-${key}`, 'change');
    }
    await page.waitForTimeout(500);
    {
      const z = await page.evaluate(() => ({ err: window.__impact.error, n: window.__impact.sol?.waves.length, xi: window.__impact.sol?.waves.map((w) => w.xiFront) }));
      if (z.err || !z.n || Math.abs(z.xi[0] - 1.00736) > 1e-3) fail(`глава 2, p0 = s0 = 0: ${JSON.stringify(z)}`);
    }
    await page.selectOption('#i-preset', { index: 0 });
    await page.waitForTimeout(300);
    // §2.4.2 с веером
    await page.selectOption('#i-preset', { index: 4 });
    await page.waitForTimeout(300);
    const kinds = await page.evaluate(() => window.__impact.sol.waves.map((w) => w.kind).join(','));
    if (kinds !== 'shock,rotation,fan') fail(`§2.4.2: структура волн ${kinds}`);
    await page.selectOption('#i-field', 'g');
    await page.fill('#i-in-tau', '0.8');
    await page.dispatchEvent('#i-in-tau', 'change');
    await page.waitForTimeout(300);
    await page.screenshot({ path: `${outDir}/ch2-oblique-${scheme}.png`, fullPage: false });
    // 3D-вид: холст WebGL появился и что-то нарисовал
    await page.click('#i-view-3d');
    await page.waitForSelector('#i-3d canvas');
    await page.waitForTimeout(700);
    const colors = await page.evaluate(() => {
      const c = document.querySelector('#i-3d canvas');
      const tmp = document.createElement('canvas');
      tmp.width = 120;
      tmp.height = 80;
      const ctx = tmp.getContext('2d');
      ctx.drawImage(c, 0, 0, 120, 80);
      const d = ctx.getImageData(0, 0, 120, 80).data;
      const set = new Set();
      for (let i = 0; i < d.length; i += 4) set.add((d[i] >> 4) * 256 + (d[i + 1] >> 4) * 16 + (d[i + 2] >> 4));
      return { n: set.size, h: c.clientHeight };
    });
    if (colors.n < 12) fail(`3D-вид почти пустой: ${colors.n} цветов`);
    if (colors.h < 250) fail(`3D-холст слишком низкий: ${colors.h}px`);
    await page.screenshot({ path: `${outDir}/ch2-3d-${scheme}.png`, fullPage: false });
    await page.click('#i-view-2d');
    // переключение типа удара
    await page.click('#i-bc-velocity');
    await page.waitForTimeout(300);
    if (!(await page.evaluate(() => window.__impact.sol?.problem.bc.type === 'velocity'))) fail('переключение на нормальный удар');

    // ---------- Серия ударов ----------
    await page.click('#tab-series');
    await page.waitForFunction(() => window.__series?.res && window.__series.computing === 0, null, { timeout: 60000 });
    const ser = await page.evaluate(() => ({ kin: window.__series.res.kinematicError, n: window.__series.res.frames.length, w: window.__series.res.work, e: window.__series.res.energyGain }));
    if (Math.abs(ser.kin) > 1e-8) fail(`серия ударов: кинематический баланс ${ser.kin}`);
    const ftxt = await page.locator('#s-status').innerText();
    if (!ftxt.includes('Метод волновых фронтов')) fail('по умолчанию должен быть метод волновых фронтов');
    if (await page.evaluate(() => window.__series.fronts.warnings.length)) fail('метод фронтов выдал предупреждения');
    await page.screenshot({ path: `${outDir}/series-fronts-${scheme}.png`, fullPage: false });
    if (ser.n < 100) fail(`серия ударов: мало кадров ${ser.n}`);
    if (!(ser.w >= ser.e)) fail('серия ударов: прирост энергии больше работы границы');
    // масштабирование диаграммы x–t и выбор точки щелчком
    {
      const box = await page.locator('#s-xt').boundingBox();
      const cx = box.x + box.width * 0.4;
      const cy = box.y + box.height * 0.6;
      await page.mouse.move(cx, cy);
      await page.mouse.wheel(0, -600);
      await page.waitForTimeout(150);
      if (!(await page.evaluate(() => window.__series.xtWin))) fail('колесо мыши не приблизило диаграмму x–t');
      await page.mouse.click(cx, cy);
      await page.waitForTimeout(150);
      const pk = await page.evaluate(() => ({ pick: window.__series.xtPick, t: window.__series.t, playing: window.__series.playing }));
      if (!pk.pick) fail('щелчок не выбрал точку на диаграмме');
      else if (Math.abs(pk.t - pk.pick.t) > 1e-9 || pk.playing) fail('щелчок не перевёл время анимации на выбранный момент');
      if (pk.pick && !pk.pick.lines.some((l) => l.includes('встреча') || l.includes('у границы')) && !pk.pick.lines.some((l) => l.includes('p = u1,1')))
        fail('в подписи точки нет значений деформаций');
      await page.locator('#s-xt').screenshot({ path: `${outDir}/xt-zoom-pick-${scheme}.png` });
      await page.mouse.dblclick(cx, cy);
      await page.waitForTimeout(150);
      if (await page.evaluate(() => window.__series.xtWin)) fail('двойной щелчок не вернул весь вид');
    }
    // растяжение: шаблон «сжатие → растяжение → сжатие», слабые фронты
    await page.selectOption('#s-pulse-preset', 'Сжатие → растяжение → сжатие');
    await page.waitForFunction(() => window.__series.computing === 0 && window.__series.fronts && window.__series.pulses.some((p) => p.d < 0), null, { timeout: 60000 });
    if (await page.evaluate(() => window.__series.fronts.warnings.length)) fail('растяжение: метод фронтов выдал предупреждения');
    await page.check('#s-weak-fronts');
    await page.selectOption('#s-field', 's');
    await page.waitForTimeout(300);
    await page.locator('#s-xt').screenshot({ path: `${outDir}/series-tension-weak-${scheme}.png` });
    await page.uncheck('#s-weak-fronts');
    await page.selectOption('#s-field', 'p');
    // материал 2 со слабым предсжатием: сжатие и растяжение бегут с разной скоростью
    await page.selectOption('#s-material', { index: 2 });
    await page.selectOption('#s-pulse-preset', 'Удар, затем рывок (сравнить скорости сжатия и растяжения)');
    await page.waitForFunction(() => window.__series.computing === 0 && window.__series.fronts && window.__series.F0[0] === -0.0005 && window.__series.pulses.length === 2, null, { timeout: 90000 });
    if (await page.evaluate(() => window.__series.fronts.warnings.length)) fail('материал 2, растяжение: предупреждения метода фронтов');
    await page.locator('#s-xt').screenshot({ path: `${outDir}/series-m2-tension-${scheme}.png` });
    // мелкий шаг веера: расчёт идёт с новым шагом и без предупреждений
    await page.selectOption('#s-fan-step', '0.0001');
    await page.waitForFunction(() => window.__series.computing === 0 && window.__series.fronts && window.__series.fanStep === 1e-4, null, { timeout: 120000 });
    await page.waitForTimeout(300);
    if (await page.evaluate(() => window.__series.computing !== 0 || window.__series.fronts.warnings.length)) fail('мелкий шаг веера: расчёт не завершился или с предупреждениями');
    await page.selectOption('#s-fan-step', '0.00025');
    // конечные объёмы считаются до того же момента и на той же длине, что и метод фронтов
    {
      await page.waitForTimeout(300);
      await page.waitForFunction(() => window.__series.computing === 0 && window.__series.fronts && window.__series.fanStep === 2.5e-4, null, { timeout: 120000 });
      const ref = await page.evaluate(() => ({ tEnd: window.__series.tEnd, xView: window.__series.xView }));
      await page.selectOption('#s-method', 'fv');
      await page.waitForFunction(() => window.__series.computing === 0 && window.__series.res && !window.__series.fronts && window.__series.fvMeta, null, { timeout: 180000 });
      const fv = await page.evaluate(() => ({ tEnd: window.__series.tEnd, xView: window.__series.xView, kin: window.__series.res.kinematicError }));
      if (Math.abs(fv.tEnd - ref.tEnd) > 1e-6 * ref.tEnd || Math.abs(fv.xView - ref.xView) > 1e-6 * ref.xView) fail(`КО: время/длина не совпали с методом фронтов: ${JSON.stringify({ fv, ref })}`);
      if (Math.abs(fv.kin) > 1e-12) fail(`КО на длинном отрезке: нарушен кинематический баланс ${fv.kin}`);
      await page.locator('#s-xt').screenshot({ path: `${outDir}/series-fv-long-${scheme}.png` });
      // гибрид: ударные волны выделены, вееры на сетке; без лавины событий и с точными балансами
      await page.selectOption('#s-method', 'hybrid');
      await page.waitForFunction(() => window.__series.computing === 0 && window.__series.hyb, null, { timeout: 180000 });
      const hy = await page.evaluate(() => ({ w: window.__series.hyb.warnings, kin: window.__series.hyb.kinematicError, ev: window.__series.hyb.events, segs: window.__series.hyb.segments.length, work: window.__series.hyb.work, gain: window.__series.hyb.energyGain, diss: window.__series.hyb.shockDissipation }));
      if (hy.w.length) fail(`гибрид: предупреждения ${hy.w}`);
      if (Math.abs(hy.kin) > 1e-7) fail(`гибрид: кинематический баланс ${hy.kin}`);
      if (hy.ev > 500 || !hy.segs) fail(`гибрид: событий ${hy.ev}, отрезков фронтов ${hy.segs}`);
      // остаток — численная диссипация на ударных волнах, возникших на сетке из волн сжатия (не выделяются);
      // здесь (слабое предсжатие, переход через p = 0) ~8 % на 200 ячеек/L и убывает с измельчением
      if (Math.abs(hy.work - hy.gain - hy.diss) > 0.15 * Math.abs(hy.work)) fail(`гибрид: энергетический баланс ${JSON.stringify(hy)}`);
      await page.evaluate(() => {
        const r = document.getElementById('s-time');
        r.value = '0.3';
        r.dispatchEvent(new Event('input'));
      });
      await page.waitForTimeout(300);
      if ((await page.locator('#s-profiles canvas').count()) !== 4) fail('гибрид: профилей должно быть 4');
      await page.locator('#s-xt').screenshot({ path: `${outDir}/series-hybrid-${scheme}.png` });
      await page.locator('#s-profiles').screenshot({ path: `${outDir}/series-hybrid-profiles-${scheme}.png` });
      await page.selectOption('#s-method', 'fronts');
      await page.waitForFunction(() => window.__series.computing === 0 && window.__series.fronts, null, { timeout: 120000 });
      // ненапряжённая среда (p0 = s0 = 0) и предварительное растяжение: оба точных метода считают без ошибок
      for (const idx of [3, 4]) {
        await page.selectOption('#s-material', { index: idx });
        await page.waitForTimeout(300);
        await page.waitForFunction(() => window.__series.computing === 0 && window.__series.fronts, null, { timeout: 120000 });
        const fr = await page.evaluate(() => ({ err: window.__series.error, w: window.__series.fronts.warnings, kin: window.__series.fronts.kinematicError, F0: window.__series.F0 }));
        // у материала 2 около p = 0 метод фронтов уходит в огрублённый режим (пороги 1e-7): баланс ~1e-5
        if (fr.err || fr.w.length || Math.abs(fr.kin) > 1e-4) fail(`фронты, F0 = ${fr.F0}: ${JSON.stringify(fr)}`);
        await page.selectOption('#s-method', 'hybrid');
        await page.waitForFunction(() => window.__series.computing === 0 && window.__series.hyb, null, { timeout: 180000 });
        const h2 = await page.evaluate(() => ({ err: window.__series.error, w: window.__series.hyb.warnings, kin: window.__series.hyb.kinematicError, F0: window.__series.F0 }));
        if (h2.err || h2.w.length || Math.abs(h2.kin) > 1e-10) fail(`гибрид, F0 = ${h2.F0}: ${JSON.stringify(h2)}`);
        if (idx === 4) await page.locator('#s-xt').screenshot({ path: `${outDir}/series-zero-state-${scheme}.png` });
        await page.selectOption('#s-method', 'fronts');
        await page.waitForFunction(() => window.__series.computing === 0 && window.__series.fronts, null, { timeout: 120000 });
      }
    }
    await page.selectOption('#s-material', { index: 0 });
    // один удар: численное решение совпадает с точным (пунктир)
    await page.selectOption('#s-pulse-preset', 'Один удар (проверка по §2.3)');
    await page.waitForFunction(() => window.__series.computing === 0 && window.__series.res && window.__series.pulses.length === 1, null, { timeout: 60000 });
    await page.waitForTimeout(400);
    await page.click('#s-play');
    // момент t = 0.5 — до окончания одиночного удара (там показывается точное решение §2.3)
    await page.evaluate(() => {
      const r = document.getElementById('s-time');
      r.value = String(0.5 / window.__series.tEnd);
      r.dispatchEvent(new Event('input'));
    });
    await page.waitForTimeout(300);
    if ((await page.locator('#s-profiles canvas').count()) !== 4) fail('профилей серии должно быть 4');
    const nser = await page.evaluate(() => window.__profSeries);
    if (nser !== 2) fail(`нет точного решения на профилях: серий ${nser}`);
    // масштаб профилей колесом
    {
      await page.locator('#s-profiles canvas').first().scrollIntoViewIfNeeded();
      const pb = await page.locator('#s-profiles canvas').first().boundingBox();
      await page.mouse.move(pb.x + pb.width * 0.3, pb.y + pb.height * 0.5);
      await page.mouse.wheel(0, -500);
      await page.waitForTimeout(100);
      if (!(await page.evaluate(() => window.__profWin))) fail('колесо не приблизило профили');
      await page.mouse.dblclick(pb.x + pb.width * 0.3, pb.y + pb.height * 0.5);
      await page.waitForTimeout(100);
      if (await page.evaluate(() => window.__profWin)) fail('двойной щелчок не вернул весь вид профилей');
    }
    await page.selectOption('#s-field', 's');
    await page.screenshot({ path: `${outDir}/series-${scheme}.png`, fullPage: true });
    // 3D-вид серии: отдельный экземпляр, 2D-холст скрыт
    await page.click('#s-view-3d');
    await page.waitForSelector('#s-3d canvas');
    await page.waitForTimeout(700);
    const s3 = await page.evaluate(() => {
      const c = document.querySelector('#s-3d canvas');
      const tmp = document.createElement('canvas');
      tmp.width = 120;
      tmp.height = 80;
      const ctx = tmp.getContext('2d');
      ctx.drawImage(c, 0, 0, 120, 80);
      const d = ctx.getImageData(0, 0, 120, 80).data;
      const set = new Set();
      for (let i = 0; i < d.length; i += 4) set.add((d[i] >> 4) * 256 + (d[i + 1] >> 4) * 16 + (d[i + 2] >> 4));
      return { n: set.size, strip: document.getElementById('s-scene').offsetParent !== null };
    });
    if (s3.n < 12) fail(`3D серии почти пустой: ${s3.n}`);
    if (s3.strip) fail('2D-полоса не скрыта в режиме 3D');
    await page.screenshot({ path: `${outDir}/series-3d-${scheme}.png`, fullPage: false });
    await page.click('#s-view-2d');

    // ---------- §3.2 ----------
    await page.click('#tab-load');
    await page.waitForFunction(() => window.__app?.solution);
    const th = await page.evaluate(() => window.__app.solution.fronts.map((f) => f.thetaDeg));
    if (Math.abs(th[0] - 36.366) > 0.01 || Math.abs(th[1] - 19.49) > 0.01) fail(`углы материала 1: ${th}`);
    const rows = await page.locator('#results table').first().locator('tbody tr').count();
    if (rows !== 2) fail(`таблица фронтов: ${rows} строк`);
    const charts = await page.locator('#steps .chart').count();
    if (charts !== 7) fail(`ступенчатых графиков: ${charts}`);
    await page.waitForSelector('#sweep svg polyline');
    await page.screenshot({ path: `${outDir}/material1-${scheme}.png`, fullPage: true });

    // Материал 2 + раскраска по σ22 (приращение)
    await page.selectOption('#preset', { index: 1 });
    await page.selectOption('#field', 's22');
    await page.check('#chk-delta');
    await page.waitForTimeout(300);
    const th2 = await page.evaluate(() => window.__app.solution.fronts.map((f) => f.thetaDeg));
    if (Math.abs(th2[0] - 38.554) > 0.01) fail(`углы материала 2: ${th2}`);
    const warn = await page.locator('#status .warn').count();
    if (warn < 1) fail('для материала 2 должно быть предупреждение о (3.30) на Σ2');
    await page.screenshot({ path: `${outDir}/material2-s22-${scheme}.png`, fullPage: false });

    // Ползунок M и недопустимая область
    await page.fill('#in-M', '1.2');
    await page.dispatchEvent('#in-M', 'change');
    await page.waitForTimeout(300);
    const M = await page.evaluate(() => window.__app.problem.M);
    if (M !== 1.2) fail(`M не обновилось: ${M}`);
    await page.fill('#in-M', '0.9');
    await page.dispatchEvent('#in-M', 'change');
    await page.waitForTimeout(200);
    const err = await page.locator('#status .err').count();
    if (err !== 1) fail('при M < 1 должно быть сообщение об ошибке');
    await page.fill('#in-M', '1.7');
    await page.dispatchEvent('#in-M', 'change');

    // Мобильная ширина: нет горизонтальной прокрутки страницы
    await page.setViewportSize({ width: 390, height: 900 });
    await page.waitForTimeout(200);
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
    if (overflow > 1) fail(`горизонтальная прокрутка на телефоне: ${overflow}px`);
    await page.screenshot({ path: `${outDir}/mobile-${scheme}.png`, fullPage: false });
    await page.close();
  }
} finally {
  await browser.close();
  process.kill(-server.pid); // вся группа процессов: npx и дочерний vite
}
if (errors.length) {
  console.error(`${errors.length} проблем(ы)`);
  process.exit(1);
}
console.log('e2e: OK');
