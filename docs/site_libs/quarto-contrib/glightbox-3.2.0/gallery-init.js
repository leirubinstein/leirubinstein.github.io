// Open gallery photos in a lightbox, grouped by album (data-gallery).
document.addEventListener("DOMContentLoaded", function () {
  if (window.GLightbox && document.querySelector(".photo-grid .photo")) {
    GLightbox({ selector: ".photo-grid .photo", loop: true, openEffect: "zoom", closeEffect: "zoom" });
  }
});
