(TeX-add-style-hook
 "app-shape-parameters"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:thin-shell-shapes"
    "sec:ancantoid-planitude"
    "eq:planitude-from-2nd-derivative"
    "eq:taylor-expansion-implicit"
    "eq:C-k-beta"
    "eq:taylor-R-over-D"
    "eq:again-R0-over-D"
    "eq:final-second-derivative"
    "eq:final-planitude"
    "sec:ancantoid-alatude"
    "eq:Lambda-from-theta-1-90"
    "eq:theta-1-90-implicit"
    "eq:Lambda-beta-xi-theta-1-90"
    "eq:xi-k"
    "eq:theta-1-90-Taylor"
    "eq:theta-1-90-approx"
    "eq:Lambda-approx"
    "sec:cantoid-wilkinoid-shapes"
    "eq:cantoid-Pi-Lambda"
    "eq:wilkinoid-Pi-Lambda"
    "sec:asympt-open-angle"
    "eq:ancantoid-theta-inf"
    "eq:cantoid-theta-inf"))
 :latex)

