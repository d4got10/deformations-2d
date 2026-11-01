#!/usr/bin/env bash
# Собирает приложение и публикует app/dist в ветку gh-pages (GitHub Pages).
# Запуск из корня репозитория: bash app/scripts/deploy.sh
set -euo pipefail
root=$(git rev-parse --show-toplevel)
rev=$(git -C "$root" rev-parse --short HEAD)
(cd "$root/app" && npm run build)
tmp=$(mktemp -d)
git -C "$root" fetch origin gh-pages
git -C "$root" worktree add "$tmp" origin/gh-pages --detach
find "$tmp" -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +
cp -r "$root/app/dist/." "$tmp/"
touch "$tmp/.nojekyll"
git -C "$tmp" add -A
git -C "$tmp" commit -m "Deploy web app @ $rev" || echo "Нет изменений"
git -C "$tmp" push origin HEAD:gh-pages
git -C "$root" worktree remove --force "$tmp"
