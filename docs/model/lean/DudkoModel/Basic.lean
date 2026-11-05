import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity

/-!
# Алгебраические соотношения модели (docs/model/model.tex)

Здесь формально проверены те формулы статьи, которые являются алгебраическими тождествами
или следствиями из явно выписанных гипотез. Нумерация разделов совпадает с приложением A статьи.

Что *не* проверяется здесь (и почему): существование и единственность решений нелинейных
задач о распаде разрыва, допустимость ударных волн, сходимость схемы — это утверждения анализа,
для них в статье приведены численные проверки (`app/tests/model.test.ts`).
-/

set_option linter.style.header false
set_option linter.style.longLine false

noncomputable section

namespace Dudko

open Finset

/-! ## A.1. Одномерная деформация: инварианты и вектор напряжений -/

/-- Тензор малых деформаций при `u = u(x₁, t)`: `e₁₁ = p, e₁₂ = g/2, e₁₃ = s/2`. -/
def strain (p g s : ℝ) : Fin 3 → Fin 3 → ℝ :=
  ![![p, g / 2, s / 2], ![g / 2, 0, 0], ![s / 2, 0, 0]]

/-- Произведение матриц 3×3. -/
def mmul (a b : Fin 3 → Fin 3 → ℝ) : Fin 3 → Fin 3 → ℝ := fun i j => ∑ k, a i k * b k j

/-- След. -/
def tr (a : Fin 3 → Fin 3 → ℝ) : ℝ := ∑ i, a i i

theorem J1_eq (p g s : ℝ) : tr (strain p g s) = p := by
  simp [tr, strain, Fin.sum_univ_three]

theorem J2_eq (p g s : ℝ) :
    tr (mmul (strain p g s) (strain p g s)) = p ^ 2 + (g ^ 2 + s ^ 2) / 2 := by
  simp [tr, mmul, strain, Fin.sum_univ_three]
  ring

theorem J3_eq (p g s : ℝ) :
    tr (mmul (mmul (strain p g s) (strain p g s)) (strain p g s))
      = p ^ 3 + 3 / 4 * p * (g ^ 2 + s ^ 2) := by
  simp [tr, mmul, strain, Fin.sum_univ_three]
  ring

/-- Закон (1.10): `σ = A·δ + B·e + C·e²` (коэффициенты A, B, C — функции инвариантов). -/
def stress (A B C p g s : ℝ) (i j : Fin 3) : ℝ :=
  A * (if i = j then 1 else 0) + B * strain p g s i j + C * mmul (strain p g s) (strain p g s) i j

/-- Вектор напряжений на площадке `x₁ = const` — ровно то, что вычисляет `traction` в коде. -/
theorem traction_eq (A B C p g s : ℝ) :
    stress A B C p g s 0 0 = A + B * p + C * (p ^ 2 + (g ^ 2 + s ^ 2) / 4) ∧
    stress A B C p g s 1 0 = (B + C * p) / 2 * g ∧
    stress A B C p g s 2 0 = (B + C * p) / 2 * s := by
  refine ⟨?_, ?_, ?_⟩ <;>
    simp [stress, mmul, strain, Fin.sum_univ_three, -mul_eq_mul_left_iff, -mul_eq_mul_right_iff] <;>
    ring

/-! ## A.2. Условия на разрыве -/

/-- Для `w = (F, v)`, `f(w) = −(v, T)` условие `[f] = G[w]` (покомпонентно) равносильно
`[v] = −G[F]`, `[T] = G²[F]` (условия (2.40)–(2.42) диссертации). -/
theorem rh_iff (dF dv dT G : ℝ) :
    (-dv = G * dF ∧ -dT = G * dv) ↔ (dv = -G * dF ∧ dT = G ^ 2 * dF) := by
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨by linear_combination -h1, by linear_combination -h2 + G * h1⟩
  · rintro ⟨h1, h2⟩
    exact ⟨by linear_combination -h1, by linear_combination -h2 - G * h1⟩

/-- Из `[T] = G²[F]` следуют коллинеарность скачков и формула скорости (МНК в коде). -/
theorem shock_speed (dF₁ dF₂ dT₁ dT₂ G : ℝ) (h₁ : dT₁ = G ^ 2 * dF₁) (h₂ : dT₂ = G ^ 2 * dF₂)
    (hF : dF₁ ^ 2 + dF₂ ^ 2 ≠ 0) :
    dT₁ * dF₂ - dT₂ * dF₁ = 0 ∧
    G ^ 2 = (dT₁ * dF₁ + dT₂ * dF₂) / (dF₁ ^ 2 + dF₂ ^ 2) := by
  refine ⟨by rw [h₁, h₂]; ring, ?_⟩
  rw [eq_div_iff hF, h₁, h₂]
  ring

/-- Для линейного материала `T = A F` скорость разрыва — характеристическая: `A[F] = G²[F]`. -/
theorem linear_shock_is_characteristic (a b d dF₁ dF₂ G : ℝ)
    (h₁ : a * dF₁ + b * dF₂ = G ^ 2 * dF₁) (h₂ : b * dF₁ + d * dF₂ = G ^ 2 * dF₂) :
    (a - G ^ 2) * (d - G ^ 2) - b ^ 2 = 0 ∨ (dF₁ = 0 ∧ dF₂ = 0) := by
  by_cases hz : dF₁ = 0 ∧ dF₂ = 0
  · exact Or.inr hz
  · left
    have e₁ : ((a - G ^ 2) * (d - G ^ 2) - b ^ 2) * dF₁ = 0 := by
      linear_combination (d - G ^ 2) * h₁ - b * h₂
    have e₂ : ((a - G ^ 2) * (d - G ^ 2) - b ^ 2) * dF₂ = 0 := by
      linear_combination (a - G ^ 2) * h₂ - b * h₁
    rcases not_and_or.mp hz with h | h
    · exact (mul_eq_zero.mp e₁).resolve_right h
    · exact (mul_eq_zero.mp e₂).resolve_right h

/-! ## A.3. Энергия и диссипация на ударной волне -/

/-- Диссипация (2.45) в виде, в котором она считается в коде: `D = [T·v] + G[W + |v|²/2]`.
При условиях на разрыве `D = G·([W] − ⟨T⟩·[F])`, `⟨T⟩` — полусумма. -/
theorem dissipation_identity
    (Fm₁ Fm₂ Fp₁ Fp₂ vm₁ vm₂ vp₁ vp₂ Tm₁ Tm₂ Tp₁ Tp₂ Wm Wp G : ℝ)
    (hv₁ : vp₁ - vm₁ = -G * (Fp₁ - Fm₁)) (hv₂ : vp₂ - vm₂ = -G * (Fp₂ - Fm₂))
    (hT₁ : Tp₁ - Tm₁ = G ^ 2 * (Fp₁ - Fm₁)) (hT₂ : Tp₂ - Tm₂ = G ^ 2 * (Fp₂ - Fm₂)) :
    (Tp₁ * vp₁ + Tp₂ * vp₂ - Tm₁ * vm₁ - Tm₂ * vm₂)
        + G * ((Wp + (vp₁ ^ 2 + vp₂ ^ 2) / 2) - (Wm + (vm₁ ^ 2 + vm₂ ^ 2) / 2))
      = G * ((Wp - Wm) - ((Tp₁ + Tm₁) / 2 * (Fp₁ - Fm₁) + (Tp₂ + Tm₂) / 2 * (Fp₂ - Fm₂))) := by
  linear_combination ((Tp₁ + Tm₁) / 2 + G * (vp₁ + vm₁) / 2) * hv₁
    + ((Tp₂ + Tm₂) / 2 + G * (vp₂ + vm₂) / 2) * hv₂
    + (vp₁ + vm₁) / 2 * hT₁ + (vp₂ + vm₂) / 2 * hT₂

/-- Смысл `D`: энергия на отрезке `[a, b]`, содержащем разрыв `x = x₀ + G t` с постоянными
состояниями по бокам (энергия `Em`, `Ep`, мощность напряжений `P = T·v`: `Pm`, `Pp`), за время
`t₁ → t₂` меняется на приток через концы (поток энергии `−T·v`, приток `−Pm − (−Pp)`) минус `D·(t₂ − t₁)`. -/
theorem energy_balance_moving_jump (a b x₀ G t₁ t₂ Em Ep Pm Pp : ℝ) :
    let content := fun t : ℝ => Em * (x₀ + G * t - a) + Ep * (b - x₀ - G * t)
    let D := (Pp - Pm) + G * (Ep - Em)
    content t₂ - content t₁ = (t₂ - t₁) * ((-Pm) - (-Pp)) - (t₂ - t₁) * D := by
  intro content D
  simp only [content, D]
  ring

/-- Для кубического потенциала `W = k₂F²/2 + k₃F³/6` (скалярный случай) `[W] − ⟨T⟩[F] = −k₃[F]³/12`:
диссипация слабой ударной волны кубична по скачку, а её знак задаёт выпуклость (`k₃`). -/
theorem weak_shock_cubic (k₂ k₃ Fm Fp : ℝ) :
    let W := fun F : ℝ => k₂ * F ^ 2 / 2 + k₃ * F ^ 3 / 6
    let T := fun F : ℝ => k₂ * F + k₃ * F ^ 2 / 2
    (W Fp - W Fm) - (T Fp + T Fm) / 2 * (Fp - Fm) = -k₃ * (Fp - Fm) ^ 3 / 12 := by
  intro W T
  simp only [W, T]
  ring

/-! ## A.3a. Ненапряжённое состояние и кусочно-линейный материал -/

/-- Из ненапряжённого состояния: при однородном `T` (`T(εq) = εT(q)`) скорость разрыва
`0 → εq` по формуле (2.12) не зависит от амплитуды `ε` — только от направления `q`. -/
theorem ray_shock_speed (q₁ q₂ t₁ t₂ ε : ℝ) (hε : ε ≠ 0) (hq : q₁ ^ 2 + q₂ ^ 2 ≠ 0) :
    ((ε * t₁) * (ε * q₁) + (ε * t₂) * (ε * q₂)) / ((ε * q₁) ^ 2 + (ε * q₂) ^ 2)
      = (t₁ * q₁ + t₂ * q₂) / (q₁ ^ 2 + q₂ ^ 2) := by
  have h2 : (ε * q₁) ^ 2 + (ε * q₂) ^ 2 ≠ 0 := by
    have : (ε * q₁) ^ 2 + (ε * q₂) ^ 2 = ε ^ 2 * (q₁ ^ 2 + q₂ ^ 2) := by ring
    rw [this]; exact mul_ne_zero (pow_ne_zero 2 hε) hq
  rw [div_eq_div_iff h2 hq]
  ring

/-- Кусочно-линейный материал (`s = 0`): `T = E₋p` при `p < 0`, `T = E₊p` при `p > 0`,
`W = Ep²/2`. Для разрыва, бегущего вглубь, с `p₁` впереди и `p₂` позади:
`[W] − ⟨T⟩[F] = p₁p₂(E₋ − E₊)/2` ([·] = справа − слева, справа — впереди). -/
theorem bimodular_jump (Em Ep p₁ p₂ : ℝ) :
    (Em * p₁ ^ 2 / 2 - Ep * p₂ ^ 2 / 2) - (Em * p₁ + Ep * p₂) / 2 * (p₁ - p₂)
      = p₁ * p₂ * (Em - Ep) / 2 := by
  ring

/-- Разгрузка через `p = 0` одним разрывом (сжатие `p₁ < 0` впереди, растяжение `p₂ > 0` позади)
при `E₋ > E₊` даёт `[W] − ⟨T⟩[F] < 0`, т.е. `D < 0` при `G > 0`: такой разрыв недопустим. -/
theorem bimodular_jump_inadmissible (Em Ep p₁ p₂ G : ℝ) (h₁ : p₁ < 0) (h₂ : 0 < p₂) (hE : Ep < Em)
    (hG : 0 < G) : G * (p₁ * p₂ * (Em - Ep) / 2) < 0 := by
  have hp : p₁ * p₂ < 0 := mul_neg_of_neg_of_pos h₁ h₂
  have : p₁ * p₂ * (Em - Ep) / 2 < 0 := by nlinarith
  nlinarith

/-- Разрыв внутри одной линейной ветви — контактный: `[W] − ⟨T⟩[F] = 0`, `D = 0`. Поэтому
решение — два контактных разрыва (скорости `√E₋`, `√E₊`) с полкой `p = 0` между ними. -/
theorem linear_contact (E p₁ p₂ : ℝ) :
    (E * p₁ ^ 2 / 2 - E * p₂ ^ 2 / 2) - (E * p₁ + E * p₂) / 2 * (p₁ - p₂) = 0 := by
  ring

/-! ## A.4. Акустическая матрица: спектр и собственные векторы -/

section Spectrum

variable (a b d D : ℝ)

/-- `λ± = (a+d)/2 ± D`, `D = √(((a−d)/2)² + b²)` — корни характеристического уравнения. -/
theorem eigen_char (hD : D ^ 2 = ((a - d) / 2) ^ 2 + b ^ 2) :
    (a - ((a + d) / 2 + D)) * (d - ((a + d) / 2 + D)) - b ^ 2 = 0 ∧
    (a - ((a + d) / 2 - D)) * (d - ((a + d) / 2 - D)) - b ^ 2 = 0 := by
  constructor <;> linear_combination hD

/-- Вектор `(b, λ − a)` (так он строится в коде) — собственный. -/
theorem eigen_vec (l : ℝ) (hl : (a - l) * (d - l) - b ^ 2 = 0) :
    a * b + b * (l - a) = l * b ∧ b * b + d * (l - a) = l * (l - a) := by
  constructor
  · ring
  · linear_combination -hl

/-- Собственные векторы двух семейств ортогональны (матрица симметрична). -/
theorem eigen_orth (hD : D ^ 2 = ((a - d) / 2) ^ 2 + b ^ 2) :
    b * b + ((a + d) / 2 + D - a) * ((a + d) / 2 - D - a) = 0 := by
  linear_combination -hD

/-- Строгая гиперболичность: при положительно определённой `A` обе скорости вещественны,
`0 < λ₋ ≤ λ₊`. -/
theorem eigen_pos (hD : D ^ 2 = ((a - d) / 2) ^ 2 + b ^ 2) (hD0 : 0 ≤ D) (ha : 0 < a)
    (hdet : 0 < a * d - b ^ 2) :
    0 < (a + d) / 2 - D ∧ (a + d) / 2 - D ≤ (a + d) / 2 + D := by
  have hd : 0 < d := by nlinarith [sq_nonneg b]
  have hm : 0 < (a + d) / 2 + D := by linarith
  have hprod : ((a + d) / 2 - D) * ((a + d) / 2 + D) = a * d - b ^ 2 := by
    linear_combination -hD
  refine ⟨?_, by linarith⟩
  nlinarith

end Spectrum

/-! ## A.5. Квадратный корень из акустической матрицы -/

section Sqrt

variable (a b d s t : ℝ)

/-- `C = (A + sI)/t`, `s = √det A`, `t = √(tr A + 2s)`: `C² = A`. -/
theorem sqrt_sq (hs : s ^ 2 = a * d - b ^ 2) (ht : t ^ 2 = a + d + 2 * s) (ht0 : t ≠ 0) :
    (a + s) / t * ((a + s) / t) + b / t * (b / t) = a ∧
    (a + s) / t * (b / t) + b / t * ((d + s) / t) = b ∧
    b / t * (b / t) + (d + s) / t * ((d + s) / t) = d := by
  have t2 : t ^ 2 ≠ 0 := pow_ne_zero 2 ht0
  refine ⟨?_, ?_, ?_⟩
  · have e : (a + s) / t * ((a + s) / t) + b / t * (b / t) = ((a + s) ^ 2 + b ^ 2) / t ^ 2 := by
      field_simp
    rw [e, div_eq_iff t2]
    linear_combination hs - a * ht
  · have e : (a + s) / t * (b / t) + b / t * ((d + s) / t) = b * (a + d + 2 * s) / t ^ 2 := by
      field_simp
      ring
    rw [e, div_eq_iff t2]
    linear_combination (-b) * ht
  · have e : b / t * (b / t) + (d + s) / t * ((d + s) / t) = ((d + s) ^ 2 + b ^ 2) / t ^ 2 := by
      field_simp
      ring
    rw [e, div_eq_iff t2]
    linear_combination hs - d * ht

/-- `det C = s = √det A`; при `a > 0, s > 0` матрица `C` положительно определена. -/
theorem sqrt_det (hs : s ^ 2 = a * d - b ^ 2) (ht : t ^ 2 = a + d + 2 * s) (ht0 : t ≠ 0) :
    (a + s) / t * ((d + s) / t) - b / t * (b / t) = s := by
  have t2 : t ^ 2 ≠ 0 := pow_ne_zero 2 ht0
  have e : (a + s) / t * ((d + s) / t) - b / t * (b / t) = ((a + s) * (d + s) - b ^ 2) / t ^ 2 := by
    field_simp
  rw [e, div_eq_iff t2]
  linear_combination -hs - s * ht

/-- `C` коммутирует с `A` (это многочлен от `A`), поэтому `|J| = diag(C, C)` коммутирует с
якобианом потока `J` и `|J|² = J²`. -/
theorem sqrt_comm (t' : ℝ) :
    let c₁₁ := (a + s) * t'
    let c₁₂ := b * t'
    let c₂₂ := (d + s) * t'
    c₁₁ * a + c₁₂ * b = a * c₁₁ + b * c₁₂ ∧
    c₁₁ * b + c₁₂ * d = a * c₁₂ + b * c₂₂ ∧
    c₁₂ * a + c₂₂ * b = b * c₁₁ + d * c₁₂ ∧
    c₁₂ * b + c₂₂ * d = b * c₁₂ + d * c₂₂ := by
  intro c₁₁ c₁₂ c₂₂
  simp only [c₁₁, c₁₂, c₂₂]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> ring

end Sqrt

/-! ## A.6. Якобиан потока, простые волны, характеристический поток -/

/-- Якобиан `J = ∂f/∂w` для `f = −(v, T(F))`: `J(e_F, e_v) = −(e_v, A e_F)`. Вектор
`e = (r, −σ r)` с `A r = σ² r` — собственный, `J e = σ e` (направление простой волны в коде). -/
theorem simple_wave_dir (a b d r₁ r₂ σ : ℝ)
    (h₁ : a * r₁ + b * r₂ = σ ^ 2 * r₁) (h₂ : b * r₁ + d * r₂ = σ ^ 2 * r₂) :
    -(-σ * r₁) = σ * r₁ ∧ -(-σ * r₂) = σ * r₂ ∧
    -(a * r₁ + b * r₂) = σ * (-σ * r₁) ∧ -(b * r₁ + d * r₂) = σ * (-σ * r₂) := by
  refine ⟨by ring, by ring, ?_, ?_⟩
  · linear_combination -h₁
  · linear_combination -h₂

/-- `J² = diag(A, A)`: дважды применённый `J` даёт `(A e_F, A e_v)`. -/
theorem jacobian_sq (a b d F₁ F₂ v₁ v₂ : ℝ) :
    let J := fun (x : ℝ × ℝ × ℝ × ℝ) =>
      ((-x.2.2.1), (-x.2.2.2), -(a * x.1 + b * x.2.1), -(b * x.1 + d * x.2.1))
    J (J (F₁, F₂, v₁, v₂)) = (a * F₁ + b * F₂, b * F₁ + d * F₂, a * v₁ + b * v₂, b * v₁ + d * v₂) := by
  intro J
  simp only [J, Prod.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> ring

/-- Характеристический поток `f* = ½(f(L) + f(R)) − ½|J|(R − L)` согласован: `f*(w, w) = f(w)`. -/
theorem flux_consistent (v₁ T₁ c₁₁ c₁₂ F₁ F₂ : ℝ) :
    -(v₁ + v₁) / 2 - (c₁₁ * (F₁ - F₁) + c₁₂ * (F₂ - F₂)) / 2 = -v₁ ∧
    -(T₁ + T₁) / 2 - (c₁₁ * (v₁ - v₁) + c₁₂ * (v₁ - v₁)) / 2 = -T₁ := by
  constructor <;> ring

/-- Фиктивная ячейка у границы `x = 0` (`F_g = F`, `v_g = 2v* − v`): кинематическая часть
потока через границу равна точно `−v*` — граница сдвигается со скоростью `v*`. -/
theorem ghost_flux (F₁ F₂ v₁ v₂ d c₁₁ c₁₂ c₂₂ : ℝ) :
    -((2 * d - v₁) + v₁) / 2 - (c₁₁ * (F₁ - F₁) + c₁₂ * (F₂ - F₂)) / 2 = -d ∧
    -((0 - v₂) + v₂) / 2 - (c₁₂ * (F₁ - F₁) + c₂₂ * (F₂ - F₂)) / 2 = 0 := by
  constructor <;> ring

/-! ## A.7. Инвариант Римана у границы и граничное состояние -/

section Boundary

variable (c₁₁ c₁₂ c₂₂ : ℝ)

/-- Линейная акустика `F_t = v_x`, `v_t = A F_x` с `A = C²`: величина `q = v + C F`
переносится со скоростью `−C`, т.е. `q_t − C q_x = 0` (покомпонентно). -/
theorem riemann_invariant (Fx₁ Fx₂ vx₁ vx₂ : ℝ) :
    let a := c₁₁ * c₁₁ + c₁₂ * c₁₂
    let b := c₁₁ * c₁₂ + c₁₂ * c₂₂
    let d := c₁₂ * c₁₂ + c₂₂ * c₂₂
    -- F_t = v_x, v_t = A F_x
    let Ft₁ := vx₁
    let Ft₂ := vx₂
    let vt₁ := a * Fx₁ + b * Fx₂
    let vt₂ := b * Fx₁ + d * Fx₂
    (vt₁ + c₁₁ * Ft₁ + c₁₂ * Ft₂) - (c₁₁ * (vx₁ + c₁₁ * Fx₁ + c₁₂ * Fx₂) + c₁₂ * (vx₂ + c₁₂ * Fx₁ + c₂₂ * Fx₂)) = 0 ∧
    (vt₂ + c₁₂ * Ft₁ + c₂₂ * Ft₂) - (c₁₂ * (vx₁ + c₁₁ * Fx₁ + c₁₂ * Fx₂) + c₂₂ * (vx₂ + c₁₂ * Fx₁ + c₂₂ * Fx₂)) = 0 := by
  intro a b d Ft₁ Ft₂ vt₁ vt₂
  simp only [a, b, d, Ft₁, Ft₂, vt₁, vt₂]
  constructor <;> ring

/-- Граничное состояние: из `v_b + C F_b = v + C F` следует `F_b = F − C⁻¹ (v_b − v)`.
(В коде `hybrid.ts` стоял знак «+»; ошибка найдена при выводе этой формулы и исправлена.) -/
theorem boundary_state (F₁ F₂ v₁ v₂ Fb₁ Fb₂ vb₁ vb₂ : ℝ)
    (hdet : c₁₁ * c₂₂ - c₁₂ ^ 2 ≠ 0)
    (h₁ : vb₁ + (c₁₁ * Fb₁ + c₁₂ * Fb₂) = v₁ + (c₁₁ * F₁ + c₁₂ * F₂))
    (h₂ : vb₂ + (c₁₂ * Fb₁ + c₂₂ * Fb₂) = v₂ + (c₁₂ * F₁ + c₂₂ * F₂)) :
    Fb₁ = F₁ - (c₂₂ * (vb₁ - v₁) - c₁₂ * (vb₂ - v₂)) / (c₁₁ * c₂₂ - c₁₂ ^ 2) ∧
    Fb₂ = F₂ - (-c₁₂ * (vb₁ - v₁) + c₁₁ * (vb₂ - v₂)) / (c₁₁ * c₂₂ - c₁₂ ^ 2) := by
  have k₁ : c₂₂ * (vb₁ - v₁) - c₁₂ * (vb₂ - v₂) = (F₁ - Fb₁) * (c₁₁ * c₂₂ - c₁₂ ^ 2) := by
    linear_combination c₂₂ * h₁ - c₁₂ * h₂
  have k₂ : -c₁₂ * (vb₁ - v₁) + c₁₁ * (vb₂ - v₂) = (F₂ - Fb₂) * (c₁₁ * c₂₂ - c₁₂ ^ 2) := by
    linear_combination -c₁₂ * h₁ + c₁₁ * h₂
  constructor
  · rw [k₁, mul_div_assoc, div_self hdet]; ring
  · rw [k₂, mul_div_assoc, div_self hdet]; ring

end Boundary

/-! ## A.8. Линеаризованные задачи: распад разрыва и краевая задача -/

section Linearized

/- Ортонормированный базис собственных векторов: `r_f = (x, y)`, `r_s = (−y, x)`. -/
variable (x y cf cs : ℝ)

/-- Линеаризованная задача о распаде разрыва (начальное приближение Ньютона в коде):
`ΔF = α₁ r_f + α₂ r_s − β₂ r_s − β₁ r_f`, `Δv = c_f α₁ r_f + c_s α₂ r_s + c_s β₂ r_s + c_f β₁ r_f`.
Явное решение через проекции на `r_f`, `r_s`. -/
theorem lin_riemann_solution (ΔF₁ ΔF₂ Δv₁ Δv₂ : ℝ) (hn : x ^ 2 + y ^ 2 = 1) (hf : cf ≠ 0) (hs : cs ≠ 0) :
    let Pf := fun u₁ u₂ : ℝ => u₁ * x + u₂ * y
    let Ps := fun u₁ u₂ : ℝ => -u₁ * y + u₂ * x
    let α₁ := (Pf ΔF₁ ΔF₂ + Pf Δv₁ Δv₂ / cf) / 2
    let β₁ := (Pf Δv₁ Δv₂ / cf - Pf ΔF₁ ΔF₂) / 2
    let α₂ := (Ps ΔF₁ ΔF₂ + Ps Δv₁ Δv₂ / cs) / 2
    let β₂ := (Ps Δv₁ Δv₂ / cs - Ps ΔF₁ ΔF₂) / 2
    α₁ * x + α₂ * (-y) - β₂ * (-y) - β₁ * x = ΔF₁ ∧
    α₁ * y + α₂ * x - β₂ * x - β₁ * y = ΔF₂ ∧
    cf * α₁ * x + cs * α₂ * (-y) + cs * β₂ * (-y) + cf * β₁ * x = Δv₁ ∧
    cf * α₁ * y + cs * α₂ * x + cs * β₂ * x + cf * β₁ * y = Δv₂ := by
  intro Pf Ps α₁ β₁ α₂ β₂
  simp only [Pf, Ps, α₁, β₁, α₂, β₂]
  refine ⟨?_, ?_, ?_, ?_⟩
  · field_simp
    linear_combination 2 * cf * cs * ΔF₁ * hn
  · field_simp
    linear_combination 2 * cf * cs * ΔF₂ * hn
  · field_simp
    linear_combination 2 * Δv₁ * hn
  · field_simp
    linear_combination 2 * Δv₂ * hn

/-- Единственность: однородная линеаризованная задача имеет только нулевое решение
(волны не рождаются, если состояния совпадают). -/
theorem lin_riemann_unique (α₁ α₂ β₁ β₂ : ℝ) (hn : x ^ 2 + y ^ 2 = 1) (hf : cf ≠ 0) (hs : cs ≠ 0)
    (e₁ : α₁ * x + α₂ * (-y) - β₂ * (-y) - β₁ * x = 0)
    (e₂ : α₁ * y + α₂ * x - β₂ * x - β₁ * y = 0)
    (e₃ : cf * α₁ * x + cs * α₂ * (-y) + cs * β₂ * (-y) + cf * β₁ * x = 0)
    (e₄ : cf * α₁ * y + cs * α₂ * x + cs * β₂ * x + cf * β₁ * y = 0) :
    α₁ = 0 ∧ α₂ = 0 ∧ β₁ = 0 ∧ β₂ = 0 := by
  have d₁ : α₁ - β₁ = 0 := by linear_combination x * e₁ + y * e₂ - (α₁ - β₁) * hn
  have d₂ : α₂ - β₂ = 0 := by linear_combination -y * e₁ + x * e₂ - (α₂ - β₂) * hn
  have s₁ : cf * (α₁ + β₁) = 0 := by linear_combination x * e₃ + y * e₄ - cf * (α₁ + β₁) * hn
  have s₂ : cs * (α₂ + β₂) = 0 := by linear_combination -y * e₃ + x * e₄ - cs * (α₂ + β₂) * hn
  have p₁ : α₁ + β₁ = 0 := (mul_eq_zero.mp s₁).resolve_left hf
  have p₂ : α₂ + β₂ = 0 := (mul_eq_zero.mp s₂).resolve_left hs
  refine ⟨by linarith, by linarith, by linarith, by linarith⟩

/-- Линеаризованная краевая задача §2.3: у границы задана скорость, `Δv = −c_f α r_f − c_s β r_s`
⇒ `α = −(Δv·r_f)/c_f`, `β = −(Δv·r_s)/c_s`. -/
theorem lin_boundary_solution (Δv₁ Δv₂ : ℝ) (hn : x ^ 2 + y ^ 2 = 1) (hf : cf ≠ 0) (hs : cs ≠ 0) :
    let a₁ := -(Δv₁ * x + Δv₂ * y) / cf
    let a₂ := -(-Δv₁ * y + Δv₂ * x) / cs
    (-cf) * a₁ * x - cs * a₂ * (-y) = Δv₁ ∧ (-cf) * a₁ * y - cs * a₂ * x = Δv₂ := by
  intro a₁ a₂
  simp only [a₁, a₂]
  constructor
  · field_simp
    linear_combination (Δv₁) * hn
  · field_simp
    linear_combination (Δv₂) * hn

end Linearized

end Dudko

end
