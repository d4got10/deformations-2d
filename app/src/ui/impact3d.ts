/**
 * Трёхмерный вид для главы 2: брусок материала x1 ∈ [0, Xmax], |x2|, |x3| ≤ h.
 * Перемещения зависят только от x1, поэтому сечение x1 = const сдвигается как целое
 * на (u1, u2, u3), а продольные материальные линии изгибаются в пространстве: на
 * квазипродольной волне — сжатие, на квазипоперечной — перекос, на нейтральной —
 * поворот направления перекоса. Этот модуль загружается лениво (three.js ~600 КБ).
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { Line2 } from 'three/examples/jsm/lines/Line2.js';
import { LineGeometry } from 'three/examples/jsm/lines/LineGeometry.js';
import { LineMaterial } from 'three/examples/jsm/lines/LineMaterial.js';

type Vec3 = [number, number, number];

export interface Scene3DData {
  /** Отсчётные координаты x1 (в долях окна) с точками разрывов. */
  Xs: number[];
  /** Перемещения в этих точках (уже в единицах окна, без увеличения). */
  D: Vec3[];
  disp: (X: number) => Vec3;
  K: number;
  Xmax: number;
  /** CSS-цвет материала в отсчётной точке X. */
  colorAt: (X: number) => string;
  fronts: { X: number; color: string; label: string; dashed: boolean }[];
  /** Направление нагрузки на границе в осях (x1, x2, x3). */
  load: Vec3;
  loadLabel: string;
  theme: { bg: string; grid: string; text: string; muted: string; surface: string; accent: string };
  caption: string;
}

const H = 0.26; // полуширина бруска
const SHIFT = -0.56; // центрирование по x1

interface State {
  renderer: THREE.WebGLRenderer;
  scene: THREE.Scene;
  camera: THREE.PerspectiveCamera;
  controls: OrbitControls;
  dynamic: THREE.Group;
  overlay: HTMLDivElement;
  labels: { el: HTMLDivElement; pos: THREE.Vector3 }[];
  centerMat: LineMaterial;
  container: HTMLElement;
}

let S: State | null = null;

function init(container: HTMLElement): State {
  const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(window.devicePixelRatio || 1);
  container.appendChild(renderer.domElement);
  renderer.domElement.setAttribute('aria-label', 'Трёхмерный вид бруска материала; вращение — мышью');
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(32, 1, 0.01, 50);
  camera.position.set(0.95, 0.75, 2.05);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.target.set(0, 0, 0);
  scene.add(new THREE.AmbientLight(0xffffff, 1.6));
  const sun = new THREE.DirectionalLight(0xffffff, 1.6);
  sun.position.set(1.5, 2, 2.5);
  scene.add(sun);
  const dynamic = new THREE.Group();
  scene.add(dynamic);
  const overlay = document.createElement('div');
  overlay.className = 'overlay3d';
  container.appendChild(overlay);
  const centerMat = new LineMaterial({ color: 0xffffff, linewidth: 4, depthTest: false, transparent: true });
  const st: State = { renderer, scene, camera, controls, dynamic, overlay, labels: [], centerMat, container };
  controls.addEventListener('change', () => render(st));
  const tick = () => {
    if (!S) return;
    if (controls.update()) render(st);
    requestAnimationFrame(tick);
  };
  requestAnimationFrame(tick);
  return st;
}

function cssColor(str: string, bg: THREE.Color): THREE.Color {
  const m = str.match(/rgba?\(([^)]+)\)/);
  if (m) {
    const [r, g, b, a = '1'] = m[1].split(',').map((x) => x.trim());
    const c = new THREE.Color(+r / 255, +g / 255, +b / 255);
    return bg.clone().lerp(c, +a);
  }
  return new THREE.Color(str);
}

function resize(st: State) {
  const w = st.container.clientWidth;
  const h = st.container.clientHeight;
  if (!w || !h) return;
  const size = st.renderer.getSize(new THREE.Vector2());
  if (size.x !== w || size.y !== h) {
    st.renderer.setSize(w, h, false);
    st.renderer.domElement.style.width = `${w}px`;
    st.renderer.domElement.style.height = `${h}px`;
    st.camera.aspect = w / h;
    st.camera.updateProjectionMatrix();
    st.centerMat.resolution.set(w, h);
  }
}

function render(st: State) {
  resize(st);
  st.renderer.render(st.scene, st.camera);
  const w = st.container.clientWidth;
  const h = st.container.clientHeight;
  for (const l of st.labels) {
    const p = l.pos.clone().project(st.camera);
    const visible = p.z < 1;
    l.el.style.display = visible ? '' : 'none';
    l.el.style.transform = `translate(${((p.x + 1) / 2) * w}px, ${((1 - p.y) / 2) * h}px) translate(-50%, -100%)`;
  }
}

function disposeGroup(g: THREE.Group) {
  for (const c of [...g.children]) {
    g.remove(c);
    const o = c as THREE.Mesh;
    o.geometry?.dispose();
    const mat = o.material as THREE.Material | THREE.Material[] | undefined;
    if (Array.isArray(mat)) mat.forEach((m) => m.dispose());
    else if (mat && mat !== S?.centerMat) mat.dispose();
  }
}

export function update(container: HTMLElement, data: Scene3DData) {
  if (!S || S.container !== container) S = init(container);
  const st = S;
  const bg = new THREE.Color(data.theme.bg);
  st.renderer.setClearColor(bg);
  disposeGroup(st.dynamic);
  st.overlay.innerHTML = '';
  st.labels = [];
  const K = data.K;
  const P = (X: number, y: number, z: number, d: Vec3) => new THREE.Vector3(X + K * d[0] + SHIFT, y + K * d[1], z + K * d[2]);
  // оси three: x ← x1 (вправо), y ← x2 (вверх), z ← x3 (к зрителю)

  // --- боковые грани
  const faces: { fixed: 'y' | 'z'; sign: number }[] = [
    { fixed: 'z', sign: 1 },
    { fixed: 'z', sign: -1 },
    { fixed: 'y', sign: 1 },
    { fixed: 'y', sign: -1 },
  ];
  const n = data.Xs.length;
  // цвет вершины — цвет материала в её точке; точки разрывов сдвоены (X ± 1e-9), граница цвета резкая
  const cols = data.Xs.map((X) => cssColor(data.colorAt(X), bg));
  const faceMat = new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1 });
  for (const f of faces) {
    const pos: number[] = [];
    const col: number[] = [];
    const idx: number[] = [];
    data.Xs.forEach((X, i) => {
      const d = data.D[i];
      const a = f.fixed === 'z' ? P(X, -H, f.sign * H, d) : P(X, f.sign * H, -H, d);
      const b = f.fixed === 'z' ? P(X, H, f.sign * H, d) : P(X, f.sign * H, H, d);
      pos.push(a.x, a.y, a.z, b.x, b.y, b.z);
      const cc = cols[i];
      col.push(cc.r, cc.g, cc.b, cc.r, cc.g, cc.b);
      if (i < n - 1) idx.push(2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 1, 2 * i + 3, 2 * i + 2);
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
    g.setIndex(idx);
    g.computeVertexNormals();
    st.dynamic.add(new THREE.Mesh(g, faceMat));
  }
  // торцы: граница x1 = 0 и дальний торец
  for (const [X, d, c] of [
    [0, data.disp(0), cssColor(data.colorAt(1e-6), bg)],
    [data.Xmax, data.disp(data.Xmax), cssColor(data.colorAt(data.Xmax), bg)],
  ] as [number, Vec3, THREE.Color][]) {
    const g = new THREE.BufferGeometry().setFromPoints([P(X, -H, -H, d), P(X, H, -H, d), P(X, H, H, d), P(X, -H, H, d)]);
    g.setIndex([0, 1, 2, 0, 2, 3]);
    g.computeVertexNormals();
    st.dynamic.add(new THREE.Mesh(g, new THREE.MeshLambertMaterial({ color: c, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1 })));
  }

  // --- материальная сетка на поверхности
  const seg: number[] = [];
  const push = (a: THREE.Vector3, b: THREE.Vector3) => seg.push(a.x, a.y, a.z, b.x, b.y, b.z);
  const hg = 0.065;
  for (let X = 0; X <= data.Xmax + 1e-9; X += hg) {
    const d = data.disp(X);
    const c = [P(X, -H, -H, d), P(X, H, -H, d), P(X, H, H, d), P(X, -H, H, d)];
    for (let k = 0; k < 4; k++) push(c[k], c[(k + 1) % 4]);
  }
  const steps = Math.round((2 * H) / hg);
  for (let k = 0; k <= steps; k++) {
    const t = -H + (2 * H * k) / steps;
    for (const [y, z] of [
      [t, H],
      [t, -H],
      [H, t],
      [-H, t],
    ]) {
      for (let i = 0; i < n - 1; i++) push(P(data.Xs[i], y, z, data.D[i]), P(data.Xs[i + 1], y, z, data.D[i + 1]));
    }
  }
  const gridG = new THREE.BufferGeometry();
  gridG.setAttribute('position', new THREE.Float32BufferAttribute(seg, 3));
  st.dynamic.add(new THREE.LineSegments(gridG, new THREE.LineBasicMaterial({ color: new THREE.Color(data.theme.grid), transparent: true, opacity: 0.55 })));

  // --- ось бруска (видна сквозь грани): её излом показывает поляризацию сдвига
  const cpos: number[] = [];
  data.Xs.forEach((X, i) => {
    const p = P(X, 0, 0, data.D[i]);
    cpos.push(p.x, p.y, p.z);
  });
  const lg = new LineGeometry();
  lg.setPositions(cpos);
  st.centerMat.color = new THREE.Color(data.theme.accent);
  const center = new Line2(lg, st.centerMat);
  center.renderOrder = 10;
  st.dynamic.add(center);

  // --- фронты волн: полупрозрачные плоскости
  for (const fr of data.fronts) {
    if (fr.X > data.Xmax) continue;
    const d = data.disp(fr.X);
    const e = H * 1.25;
    const g = new THREE.BufferGeometry().setFromPoints([P(fr.X, -e, -e, d), P(fr.X, e, -e, d), P(fr.X, e, e, d), P(fr.X, -e, e, d)]);
    g.setIndex([0, 1, 2, 0, 2, 3]);
    const m = new THREE.MeshBasicMaterial({ color: new THREE.Color(fr.color), transparent: true, opacity: fr.dashed ? 0.16 : 0.24, side: THREE.DoubleSide, depthWrite: false });
    const mesh = new THREE.Mesh(g, m);
    mesh.renderOrder = 5;
    st.dynamic.add(mesh);
    const edge = new THREE.LineLoop(new THREE.BufferGeometry().setFromPoints([P(fr.X, -e, -e, d), P(fr.X, e, -e, d), P(fr.X, e, e, d), P(fr.X, -e, e, d)]), new THREE.LineBasicMaterial({ color: new THREE.Color(fr.color) }));
    st.dynamic.add(edge);
    addLabel(st, fr.label, P(fr.X, e, -e, d), fr.color);
  }

  // --- нагрузка на границе
  const dir = new THREE.Vector3(data.load[0], data.load[1], data.load[2]).normalize();
  const d0 = data.disp(0);
  const accent = new THREE.Color(data.theme.accent);
  for (const y of [-0.6, 0, 0.6])
    for (const z of [-0.6, 0, 0.6]) {
      const tip = P(-0.015, y * H, z * H, d0);
      const L = 0.22;
      st.dynamic.add(new THREE.ArrowHelper(dir, tip.clone().addScaledVector(dir, -L), L, accent, 0.05, 0.03));
    }
  addLabel(st, data.loadLabel, P(-0.2, H * 1.15, 0, d0), data.theme.accent);

  // --- оси координат
  const o = new THREE.Vector3(SHIFT - 0.05, -H - 0.12, H + 0.12);
  const axes: [Vec3, string][] = [
    [[1, 0, 0], 'x1'],
    [[0, 1, 0], 'x2'],
    [[0, 0, 1], 'x3'],
  ];
  for (const [a, name] of axes) {
    const v = new THREE.Vector3(...a);
    st.dynamic.add(new THREE.ArrowHelper(v, o, 0.18, new THREE.Color(data.theme.muted), 0.035, 0.02));
    addLabel(st, name, o.clone().addScaledVector(v, 0.24), data.theme.muted);
  }
  const cap = document.createElement('div');
  cap.className = 'caption3d';
  cap.textContent = data.caption;
  st.overlay.appendChild(cap);
  render(st);
}

function addLabel(st: State, text: string, pos: THREE.Vector3, color: string) {
  const el = document.createElement('div');
  el.className = 'label3d';
  el.textContent = text;
  el.style.color = color;
  st.overlay.appendChild(el);
  st.labels.push({ el, pos });
}

/** Сбросить камеру в исходное положение. */
export function resetView() {
  if (!S) return;
  S.camera.position.set(0.95, 0.75, 2.05);
  S.controls.target.set(0, 0, 0);
  S.controls.update();
  render(S);
}
