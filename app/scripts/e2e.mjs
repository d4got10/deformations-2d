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
