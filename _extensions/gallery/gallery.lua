-- {{< gallery ALBUM >}}
--
-- Renders a masonry photo grid with a lightbox from the page's front matter:
--
--   gallery:
--     ALBUM:
--       - file: DSC_1234.jpg
--         caption: "Optional caption"
--
-- Expects web-sized images at   images/photos/ALBUM/FILE
-- and thumbnails at            images/photos/ALBUM/thumbs/FILE
-- (both produced by scripts/optimize-photos.sh). The lightbox (GLightbox) is
-- bundled with this extension and added to any page that uses the shortcode.

local function escape(s)
  return (s:gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;"):gsub('"', "&quot;"))
end

return {
  gallery = function(args, kwargs, meta)
    local album = pandoc.utils.stringify(args[1] or "")
    local photos = meta.gallery and meta.gallery[album]
    if album == "" or not photos then
      quarto.log.warning("gallery: no photos listed under `gallery: " .. album .. "` in front matter")
      return pandoc.Null()
    end

    quarto.doc.add_html_dependency({
      name = "glightbox",
      version = "3.2.0",
      scripts = { "glightbox/glightbox.min.js", "glightbox/gallery-init.js" },
      stylesheets = { "glightbox/glightbox.min.css" },
    })

    local dir = "images/photos/" .. album .. "/"
    local default_alt = pandoc.utils.stringify(kwargs["alt"] or "")
    local out = { '<div class="photo-grid column-screen-inset">' }

    for _, p in ipairs(photos) do
      local file = pandoc.utils.stringify(p.file)
      local caption = p.caption and pandoc.utils.stringify(p.caption) or ""
      local alt = escape(caption ~= "" and caption or default_alt)
      local overlay = caption ~= ""
        and ('<span class="photo-caption">' .. escape(caption) .. "</span>") or ""
      table.insert(out, string.format(
        '<a class="photo" href="%s%s" data-gallery="%s"%s>' ..
        '<img src="%sthumbs/%s" alt="%s" loading="lazy" decoding="async">%s</a>',
        dir, file, album,
        caption ~= "" and (' data-description="' .. escape(caption) .. '"') or "",
        dir, file, alt, overlay))
    end

    table.insert(out, "</div>")
    return pandoc.RawBlock("html", table.concat(out, "\n"))
  end
}
