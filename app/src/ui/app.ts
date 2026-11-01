/** Вкладки: глава 2 (удар по границе) и §3.2 (подвижная нагрузка). */
import './impactPage.ts';
import './main.ts';
import './seriesPage.ts';

const tabs = Array.from(document.querySelectorAll<HTMLButtonElement>('.tabs [role=tab]'));

function show(page: string) {
  for (const t of tabs) {
    const on = t.dataset.page === page;
    t.setAttribute('aria-selected', String(on));
    document.getElementById(`page-${t.dataset.page}`)!.hidden = !on;
  }
  try {
    history.replaceState(null, '', `#${page}`);
  } catch {
    /* предпросмотр без истории */
  }
  window.dispatchEvent(new Event('resize'));
}

for (const t of tabs) t.addEventListener('click', () => show(t.dataset.page!));
show(['#load', '#series'].includes(location.hash) ? location.hash.slice(1) : 'impact');
