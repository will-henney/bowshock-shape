(TeX-add-style-hook
 "app-intermediate"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:radius-curvature"
    "eq:Rcurv-general"
    "eq:Rcurv-polar"
    "sec:plane-sky-projection"
    "eq:Rotation-matrix-y"
    "eq:xunit-prime"
    "eq:yunit-prime"
    "eq:zunit-prime"
    "eq:Rotation-matrix-x"))
 :latex)

