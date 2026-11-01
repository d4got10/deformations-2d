/**
 * Отрисовка деформированного полупространства в подвижной системе координат,
 * связанной с передней точкой нагрузки L0: η = x2 − st (по горизонтали), x1 — вниз.
 * В этой системе картина стационарна; материальные линии x2 = const «протекают» влево.
 *
 * Поле перемещений кусочно-линейно и непрерывно, поэтому каждая зона (клин) отображается
 * аффинно и рисуется точно; линии сетки — ломаные с изломами ровно на фронтах.
 */
import type { Grad } from '../core/material.ts';
import { type MovingLoadSolution, displacement, zoneIndex } from '../core/movingLoad.ts';

export type FieldKey = 'none' | 's11' | 's12' | 's22' | 's33' | 'vol' | 'speed';

export const FIELD_LABELS: Record<FieldKey, string> = {
  none: 'нет (только зоны)',
  s11: 'σ11 / (λ+2μ)',
  s12: 'σ12 / (λ+2μ)',
  s22: 'σ22 / (λ+2μ)',
  s33: 'σ33 / (λ+2μ)',
  vol: 'объёмная деформация u1,1 + u2,2',
  speed: 'скорость частиц |v| / s',
};

export interface SceneOptions {
  magnify: number;
  showPrestrain: boolean;
  field: FieldKey;
  /** Отображать приращение поля относительно зоны 0. */
  delta: boolean;
  showGrid: boolean;
  /** Сдвиг материальной сетки (пропорционален s·t), в единицах длины. */
  phase: number;
}

export function fieldValue(sol: MovingLoadSolution, zone: 0 | 1 | 2, key: FieldKey, delta: boolean): number {
  const v = (i: 0 | 1 | 2) => {
    const z = sol.zones[i];
    switch (key) {
      case 's11':
        return z.stress.s11;
      case 's12':
        return z.stress.s12;
      case 's22':
        return z.stress.s22;
      case 's33':
        return z.stress.s33;
      case 'vol':
        return z.grad.a + z.grad.n;
      case 'speed':
        return Math.hypot(z.v[0], z.v[1]);
      default:
        return 0;
    }
  };
  return delta ? v(zone) - v(0) : v(zone);
}

const css = (name: string) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

export function drawScene(canvas: HTMLCanvasElement, sol: MovingLoadSolution | null, opt: SceneOptions): { min: number; max: number } | null {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  }
  const ctx = canvas.getContext('2d')!;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  ctx.fillStyle = css('--panel');
  ctx.fillRect(0, 0, w, h);
  if (!sol) return null;

  // Масштаб: 1 единица длины = sc пикселей; L0 — в 80 % ширины, граница — на 44 px сверху.
  const sc = Math.min(w, 1.6 * h) / 1.25;
  const ox = w * 0.8;
  const oy = 44;
  const etaMin = -ox / sc - 0.2;
  const etaMax = (w - ox) / sc + 0.2;
  const x1Max = (h - oy) / sc + 0.2;
  const g0 = sol.zones[0].grad;
  const K = opt.magnify;

  const disp = (g: Grad, x1: number, eta: number): [number, number] => {
    const u = displacement(g, x1, eta);
    if (opt.showPrestrain) return u;
    const u0 = displacement(g0, x1, eta);
    return [u[0] - u0[0], u[1] - u0[1]];
  };
  /** Отсчётная точка (x1, η) → экран; зона определяется по отсчётной точке. */
  const mapIn = (zone: 0 | 1 | 2, x1: number, eta: number): [number, number] => {
    const [u1, u2] = disp(sol.zones[zone].grad, x1, eta);
    return [ox + (eta + K * u2) * sc, oy + (x1 + K * u1) * sc];
  };
  const map = (x1: number, eta: number) => mapIn(zoneIndex(sol, x1, eta), x1, eta);

  const th1 = Math.atan(sol.fronts[0].xi);
  const th2 = Math.atan(sol.fronts[1].xi);
  const R = 4 * Math.max(-etaMin, x1Max, etaMax);
  const ray = (th: number): [number, number] => [R * Math.sin(th), -R * Math.cos(th)]; // (x1, η)

  // Зоны: полигоны в отсчётной конфигурации (x1, η)
  const polys: [0 | 1 | 2, [number, number][]][] = [
    [0, [[0, 0], ray(th1), [R, -R * Math.cos(th1)], [R, R], [0, R]]],
    [1, [[0, 0], ray(th1), ray(th2)]],
    [2, [[0, 0], ray(th2), [0, -R]]],
  ];
  const vals = [0, 1, 2].map((i) => fieldValue(sol, i as 0 | 1 | 2, opt.field, opt.delta));
  const vmax = Math.max(...vals.map(Math.abs)) || 1;
  const zoneTint = [css('--zone0'), css('--zone1'), css('--zone2')];

  ctx.save();
  ctx.beginPath();
  ctx.rect(0, 0, w, h);
  ctx.clip();
  for (const [zi, poly] of polys) {
    ctx.beginPath();
    poly.forEach(([x1, eta], i) => {
      const [X, Y] = mapIn(zi, x1, eta);
      if (i === 0) ctx.moveTo(X, Y);
      else ctx.lineTo(X, Y);
    });
    ctx.closePath();
    ctx.fillStyle = zoneTint[zi];
    ctx.fill();
    if (opt.field !== 'none') {
      const t = vals[zi] / vmax;
      const alpha = 0.12 + 0.68 * Math.abs(t);
      ctx.fillStyle = t >= 0 ? `rgba(214, 72, 52, ${alpha})` : `rgba(48, 104, 214, ${alpha})`;
      ctx.fill();
    }
  }

  // Материальная сетка
  if (opt.showGrid) {
    const hStep = 26 / sc;
    ctx.strokeStyle = css('--grid');
    ctx.lineWidth = 1;
    const fronts = [sol.fronts[0].xi, sol.fronts[1].xi];
    // вертикальные линии x2 = const: η = j·h − phase
    const p = ((opt.phase % hStep) + hStep) % hStep;
    for (let eta = Math.floor(etaMin / hStep) * hStep - p; eta <= etaMax; eta += hStep) {
      const xs = [0, x1Max];
      if (eta < 0) for (const xi of fronts) xs.push(-eta * xi);
      xs.sort((a, b) => a - b);
      ctx.beginPath();
      xs.filter((x) => x <= x1Max).forEach((x1, i) => {
        const [X, Y] = map(x1, eta);
        if (i === 0) ctx.moveTo(X, Y);
        else ctx.lineTo(X, Y);
      });
      ctx.stroke();
    }
    // горизонтальные линии x1 = const
    for (let x1 = hStep; x1 <= x1Max; x1 += hStep) {
      const es = [etaMin, etaMax, 0];
      for (const xi of fronts) es.push(-x1 / xi);
      es.sort((a, b) => a - b);
      ctx.beginPath();
      es.filter((e) => e >= etaMin && e <= etaMax).forEach((eta, i) => {
        const [X, Y] = map(x1, eta);
        if (i === 0) ctx.moveTo(X, Y);
        else ctx.lineTo(X, Y);
      });
      ctx.stroke();
    }
  }

  // Фронты
  const frontColors = [css('--front1'), css('--front2')];
  [th1, th2].forEach((th, i) => {
    const [x1e, ee] = ray(th);
    const [X0, Y0] = mapIn(i as 0 | 1, 0, 0);
    const [X1, Y1] = mapIn(i as 0 | 1, x1e, ee);
    ctx.strokeStyle = frontColors[i];
    ctx.lineWidth = 2.5;
    ctx.setLineDash([8, 5]);
    ctx.beginPath();
    ctx.moveTo(X0, Y0);
    ctx.lineTo(X1, Y1);
    ctx.stroke();
    ctx.setLineDash([]);
  });

  // Граница x1 = 0 и нагрузка
  ctx.strokeStyle = css('--surface');
  ctx.lineWidth = 2.5;
  ctx.beginPath();
  [etaMin, 0, etaMax].forEach((eta, i) => {
    const [X, Y] = map(0, eta);
    if (i === 0) ctx.moveTo(X, Y);
    else ctx.lineTo(X, Y);
  });
  ctx.stroke();
  const t1 = -sol.problem.s11; // вглубь (вниз)
  const t2 = -sol.problem.s21; // вдоль границы
  const tl = Math.hypot(t1, t2) || 1;
  const dir: [number, number] = [t2 / tl, t1 / tl]; // экранные (X, Y)
  ctx.strokeStyle = css('--accent');
  ctx.fillStyle = css('--accent');
  ctx.lineWidth = 1.8;
  const arrowStep = 30 / sc;
  for (let eta = -arrowStep * 0.5; eta > etaMin; eta -= arrowStep) {
    const [X, Y] = map(0, eta);
    const L = 24;
    const sx = X - dir[0] * L;
    const sy = Y - 3 - dir[1] * L;
    ctx.beginPath();
    ctx.moveTo(sx, sy);
    ctx.lineTo(X, Y - 3);
    ctx.stroke();
    const a = Math.atan2(dir[1], dir[0]);
    ctx.beginPath();
    ctx.moveTo(X, Y - 3);
    ctx.lineTo(X - 7 * Math.cos(a - 0.4), Y - 3 - 7 * Math.sin(a - 0.4));
    ctx.lineTo(X - 7 * Math.cos(a + 0.4), Y - 3 - 7 * Math.sin(a + 0.4));
    ctx.closePath();
    ctx.fill();
  }

  // Подписи
  ctx.font = '600 13px system-ui, sans-serif';
  ctx.fillStyle = css('--text');
  const [Lx, Ly] = map(0, 0);
  ctx.fillText('L0', Lx + 6, Ly - 8);
  ctx.fillText('s →', Lx + 30, Ly - 8);
  const labelAt = (zi: 0 | 1 | 2, th: number, dist: number, text: string, color?: string) => {
    const x1 = dist * Math.sin(th);
    const eta = -dist * Math.cos(th);
    const [X, Y] = mapIn(zi, x1, eta);
    ctx.fillStyle = color ?? css('--text');
    // подписи фронтов — чуть ниже-правее линии, чтобы не перекрывать её
    ctx.fillText(text, color ? X + 10 : X, color ? Y + 14 : Y);
  };
  const D = Math.min(-etaMin, x1Max) * 0.72;
  labelAt(2, th2 * 0.45, D * 0.95, 'II');
  labelAt(1, (th1 + th2) / 2, D, 'I');
  labelAt(0, Math.min(th1 + (Math.PI / 2 - th1) * 0.5, th1 + 0.35), D * 0.9, '0');
  labelAt(0, th1, D * 1.12, `Σ1  θ1 = ${sol.fronts[0].thetaDeg.toFixed(2)}°`, frontColors[0]);
  labelAt(1, th2, D * 1.12, `Σ2  θ2 = ${sol.fronts[1].thetaDeg.toFixed(2)}°`, frontColors[1]);
  ctx.restore();

  ctx.font = '12px system-ui, sans-serif';
  ctx.fillStyle = css('--muted');
  ctx.fillText(`x2 →   перемещения ×${K}${opt.showPrestrain ? ', с предварительной деформацией' : ', относительно зоны 0'}`, 10, 18);
  ctx.fillText('x1 ↓', 10, h - 10);
  return { min: Math.min(...vals), max: Math.max(...vals) };
}
