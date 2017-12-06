(TeX-add-style-hook
 "sec-intro"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-add-symbols
    "Mach")
   (LaTeX-add-labels
    "sec:intro"
    "fig:bow-terminology"
    "fig:2-winds"
    "fig:crw-schema"
    "eq:stagnation-radius"
    "eq:beta-definition"
    "sec:plan-alat-bow"
    "fig:characteristic-radii"
    "eq:radius-curvature"
    "eq:taylor-R-theta"
    "eq:planitude"
    "eq:alatude"))
 :latex)

