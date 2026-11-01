/** Фоновый расчёт серии ударов (чтобы интерфейс не подвисал): метод фронтов или конечные объёмы. */
import { type FrontsProblem, trackFronts } from '../core/fronts.ts';
import { type SeriesProblem, simulate } from '../core/fv1d.ts';

type Msg = { id: number; method: 'fv'; prob: SeriesProblem } | { id: number; method: 'fronts'; prob: FrontsProblem };

self.onmessage = (e: MessageEvent<Msg>) => {
  const { id } = e.data;
  try {
    const t0 = performance.now();
    if (e.data.method === 'fronts') {
      const res = trackFronts(e.data.prob);
      self.postMessage({ id, type: 'done', method: 'fronts', res, ms: performance.now() - t0 });
    } else {
      const res = simulate(e.data.prob, (f) => self.postMessage({ id, type: 'progress', f }));
      self.postMessage({ id, type: 'done', method: 'fv', res, ms: performance.now() - t0 }, { transfer: res.frames.map((fr) => fr.U.buffer) });
    }
  } catch (err) {
    self.postMessage({ id, type: 'error', message: (err as Error).message });
  }
};
