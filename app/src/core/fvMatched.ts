/**
 * Конечные объёмы на том же отрезке времени и длины, что и метод фронтов в режиме «до затухания»:
 * сначала методом фронтов находится момент, после которого фронты больше не встречаются, и длина
 * участка, до которой они успевают дойти; затем на этом участке и до этого момента считаются КО.
 */
import { trackFronts } from './fronts.ts';
import { type SeriesProblem, type SeriesResult, simulate } from './fv1d.ts';

export interface FvMatch {
  fanStep: number;
  /** Число кадров метода фронтов — то же, что в основном расчёте (от него зависит путь решателя). */
  frames: number;
  /** Предел общего числа ячеек (время КО растёт как квадрат числа ячеек). */
  maxCells: number;
}

export interface FvMatchedResult {
  res: SeriesResult;
  tEnd: number;
  length: number;
  /** Сколько ячеек получилось на единицу длины L (меньше заданного, если сработал предел). */
  cellsPerL: number;
  capped: boolean;
}

export function simulateMatched(prob: SeriesProblem, match: FvMatch, onProgress?: (f: number) => void): FvMatchedResult {
  const fr = trackFronts({ material: prob.material, F0: prob.F0, pulses: prob.pulses, tEnd: prob.tEnd, frames: match.frames, untilQuiet: true, fanStep: match.fanStep });
  const length = fr.xView;
  const tEnd = fr.tEnd;
  const want = Math.round(prob.N * length);
  const N = Math.max(prob.N, Math.min(want, match.maxCells));
  const res = simulate({ ...prob, N, length, tEnd }, onProgress);
  return { res, tEnd, length, cellsPerL: N / length, capped: N < want };
}
