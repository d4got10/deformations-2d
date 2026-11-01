# deformations-2d

Визуализация плоских деформаций разномодульной упругой среды по мотивам диссертации
О.В. Дудко «Динамика разномодульной изотропной упругой среды при ударных воздействиях» (1998).

Сейчас реализована задача §3.2: постоянная нагрузка движется со сверхсейсмической скоростью по границе
предварительно деформированного полупространства (модель Мясникова–Олейникова, две ударные волны).

## Запуск

```bash
cd app
npm install
npm run dev
```

Откройте адрес, который напечатает Vite (обычно http://localhost:5173).

## Публикация

Сайт: https://d4got10.github.io/deformations-2d/ (ветка `gh-pages`).
Обновить после изменений:

```bash
bash app/scripts/deploy.sh
```

## Проверка

```bash
cd app
npm test
npm run report
npm run build && npm run e2e
```

- [docs/PLAN.md](docs/PLAN.md) — план и структура кода
- [docs/VERIFICATION.md](docs/VERIFICATION.md) — методика проверки и открытые вопросы к тексту диссертации
- [docs/VERIFICATION_RESULTS.md](docs/VERIFICATION_RESULTS.md) — сгенерированные таблицы сверки
- `sources/disser.pdf` — текст диссертации (скан)
