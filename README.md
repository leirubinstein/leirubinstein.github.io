# [Leilanie Rubinstein's Personal Website](https://leirubinstein.github.io/)

Built with [Quarto](https://quarto.org). The rendered site lives in `docs/` and is served by GitHub Pages.

```bash
quarto preview      # live local preview
quarto render       # rebuild docs/ before committing
```

## Layout

| Path | What it is |
|---|---|
| `index.qmd` | Home / about page |
| `cv.qmd` | Web CV; the PDF lives in `cv/` |
| `blog.qmd` | Blog listing (+ RSS feed at `blog.xml`) |
| `posts/<date-slug>/` | One folder per post. Posts use `freeze: true`, so their results are cached in `_freeze/` and rendering does not re-run the analyses. Raw `data/` folders are git-ignored. |
| `photos.qmd` | Photo gallery. Photos are listed in its front matter and drawn by the `gallery` shortcode. |
| `images/` | Site images: logo, favicon, headshot, and web-sized photos in `images/photos/<album>/` |
| `styles/` | `light.scss` / `dark.scss` themes and `site.css` (CV + gallery layout) |
| `_extensions/` | Font Awesome and the local `gallery` shortcode |
| `scripts/optimize-photos.sh` | Builds web-sized photos and thumbnails from `_photos-src/` |
| `docs/` | Rendered output. Don't edit by hand. |

## Adding photos

Full-resolution originals stay out of git (in `_photos-src/`, which is ignored) to keep the repo small.

1. Put the original in `_photos-src/<album>/` (`nature` or `portraits`; a new folder makes a new album).
2. Run `scripts/optimize-photos.sh` (needs ImageMagick: `brew install imagemagick`). It writes a 2000px copy and a 1200px thumbnail to `images/photos/<album>/`, stripping metadata such as GPS location.
3. Add the file to the album's list in the front matter of `photos.qmd`, with an optional caption:

   ```yaml
   gallery:
     nature:
       - file: DSC_1234.jpg
         caption: "Where it was taken"
   ```

   A new album also needs a `{{< gallery album-name >}}` line in the page body.

## Updating the CV

Replace `cv/Rubinstein_Leilanie_CV.pdf` (keep the filename so links don't break), update `cv.qmd` to match, and re-render.
