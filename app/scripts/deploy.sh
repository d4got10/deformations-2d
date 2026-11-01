#!/usr/bin/env bash
# Собирает приложение и публикует app/dist в ветку gh-pages (GitHub Pages).
# Запуск из корня репозитория: bash app/scripts/deploy.sh
set -euo pipefail
root=$(git rev-parse --show-toplevel)
rev=$(git -C "$root" rev-parse --short HEAD)
(cd "$root/app" && npm run build)
tmp=$(mktemp -d)
git -C "$root" fetch origin gh-pages
git -C "$root" worktree add "$tmp" FETCH_HEAD --detach
find "$tmp" -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +
cp -r "$root/app/dist/." "$tmp/"
touch "$tmp/.nojekyll"
git -C "$tmp" add -A
# Дата коммита — на 2 часа раньше текущей: если часы машины спешат, коммит «из будущего»
# может заблокировать публикацию GitHub Pages (так было в облачной среде).
past=$(date -u -d '-2 hours' '+%Y-%m-%dT%H:%M:%S+0000' 2>/dev/null || date -u '+%Y-%m-%dT%H:%M:%S+0000')
GIT_AUTHOR_DATE="$past" GIT_COMMITTER_DATE="$past" git -C "$tmp" commit -m "Deploy web app @ $rev" || echo "Нет изменений"
git -C "$tmp" push origin HEAD:gh-pages
git -C "$root" worktree remove --force "$tmp"
