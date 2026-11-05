# Математическая модель (статья)

- `model.tex`, `model.pdf` — описание модели, классификация столкновений волн, разбор варианта
  «материал 2, слабое предсжатие, удар — рывок — удар», проверка формул.
- `lean/` — формальная проверка алгебраических соотношений статьи (Lean 4 + Mathlib), 33 теоремы
  в `lean/DudkoModel/Basic.lean`.
- `figures/` — рисунки (`make_figures.py` строит их из `case.json` и `uniaxial.json`, которые пишут `app/scripts/uniaxial-case.ts` и
  `app/scripts/model-case.ts`; сам `case.json` в репозиторий не кладётся — 4 МБ).

## Сборка

```bash
cd app && npm test                                   # численные проверки (tests/model.test.ts и др.)
cd app && node scripts/model-case.ts ../docs/model/figures/case.json
python3 docs/model/figures/make_figures.py           # нужны numpy, matplotlib
cd docs/model && latexmk -pdf model.tex              # texlive + cyrillic
cd docs/model/lean && lake build                     # Lean 4.36.0-rc2, Mathlib v4.36.0-rc2
```

`lake build` без кэша Mathlib компилирует нужную часть библиотеки (~2200 модулей, ~15 мин на 4 ядрах).
