(TeX-add-style-hook
 "sec-dust-bow-intro"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-add-symbols
    "PaperI"
    "hii")
   (LaTeX-add-labels
    "sec:introduction"
    "fig:3-types-bow"
    "sec:different-types-bow"
    "fig:zones-v-n-plane"
    "sec:imperf-coupl-betw"
    "fig:decouple-v-n-plane"
    "fig:decouple-v40-versus-n"
    "sec:case-inside-out"))
 :latex)

