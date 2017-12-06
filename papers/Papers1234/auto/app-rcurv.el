(TeX-add-style-hook
 "app-rcurv"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "app:rc-analytic"
    "eq:R-exp"
    "eq:AR1"
    "eq:AR2"
    "eq:th1th-app"
    "eq:th1th-small"
    "eq:r-small-theta"
    "eq:app-gamma"))
 :latex)

