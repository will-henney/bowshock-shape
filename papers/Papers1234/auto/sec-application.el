(TeX-add-style-hook
 "sec-application"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-add-symbols
    '("C" 1))
   (LaTeX-add-labels
    "sec:application"
    "eq:nprop"
    "eq:ngen"
    "sec:methodology"
    "sec:results"
    "fig:radii-measures-example"
    "fig:char-radii-obs"
    "fig:conic-xi"
    "fig:pressure"
    "tab:arc-fits"))
 :latex)

