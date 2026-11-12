#!/usr/bin/env bash
# Build web-sized copies of the photos in _photos-src/ (git-ignored originals).
#
#   _photos-src/<album>/<name>.jpg  ->  images/photos/<album>/<name>.jpg        (lightbox, 2000px)
#                                   ->  images/photos/<album>/thumbs/<name>.jpg (grid, 1200px)
#   _photos-src/headshot.jpg        ->  images/headshot.jpg
#
# Folders starting with "_" (e.g. _photos-src/_unpublished/) are skipped.
#
# Metadata (including GPS) is stripped. Existing outputs are skipped unless
# the source is newer, so it is safe to rerun after adding photos.
# Requires ImageMagick 7 (`brew install imagemagick`).
#
# Usage: scripts/optimize-photos.sh
set -euo pipefail
cd "$(dirname "$0")/.."

SRC=_photos-src
OUT=images/photos
FULL_PX=2000
THUMB_PX=1200

resize() { # src dest max_px quality
  if [[ ! -f "$2" || "$1" -nt "$2" ]]; then
    mkdir -p "$(dirname "$2")"
    magick "$1" -auto-orient -resize "${3}x${3}>" -strip \
      -interlace Plane -sampling-factor 4:2:0 -quality "$4" "$2"
    echo "  $2 ($(du -h "$2" | cut -f1))"
  fi
}

shopt -s nullglob nocaseglob
for src in "$SRC"/*/*.{jpg,jpeg,png,tif,tiff}; do
  album=$(basename "$(dirname "$src")")
  [[ "$album" == _* ]] && continue   # e.g. _photos-src/_unpublished/
  name=$(basename "${src%.*}").jpg
  resize "$src" "$OUT/$album/$name"        "$FULL_PX"  82
  resize "$src" "$OUT/$album/thumbs/$name" "$THUMB_PX" 78
done

[[ -f "$SRC/headshot.jpg" ]] && resize "$SRC/headshot.jpg" images/headshot.jpg 1000 85

echo "Photos: $(du -sh "$OUT" | cut -f1) total"
