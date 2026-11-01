/** Фоновый расчёт серии ударов (чтобы интерфейс не подвисал). */
import { type SeriesProblem, simulate } from '../core/fv1d.ts';

self.onmessage = (e: MessageEvent<{ id: number; prob: SeriesProblem }>) => {
  const { id, prob } = e.data;
  try {
    const t0 = performance.now();
    const res = simulate(prob, (f) => self.postMessage({ id, type: 'progress', f }));
    const ms = performance.now() - t0;
    self.postMessage({ id, type: 'done', res, ms }, { transfer: res.frames.map((fr) => fr.U.buffer) });
  } catch (err) {
    self.postMessage({ id, type: 'error', message: (err as Error).message });
  }
};
