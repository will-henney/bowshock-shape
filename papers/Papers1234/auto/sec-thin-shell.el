(TeX-add-style-hook
 "sec-thin-shell"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-add-symbols
    "thC"
    "CRW"
    "head"
    "tail")
   (LaTeX-add-labels
    "sec:crw-scenario"
    "eq:wilkinoid-R-theta"
    "sec:ancantoid"
    "fig:anisotropic-arrows"
    "eq:ancantoid-density"
    "fig:r-beta"
    "eq:Rmom"
    "eq:ancantoid-momenta"
    "eq:ancantoid-mass-loss"
    "eq:ancantoid-I-integral"
    "eq:ancantoid-momenta-outer"
    "fig:ancantoid-Pi-lambda-true"
    "eq:ancantoid-momentum-ratio"
    "eq:crw-angles"
    "eq:ancantoid-theta-theta1-implicit"
    "fig:ancantoid-angles"
    "fig:ancantoid-departure"
    "fig:xyprime"
    "fig:xyprime-ancantoid"
    "sec:true-cantoids-ancantoids"
    "fn:discontinuity"
    "eq:thetaQ-head"
    "eq:thetaQ-tail"
    "fig:thin-shell-R90-Rc"
    "sec:proj-shap-cant"
    "fig:convergence-cantoid-wilkinoid"))
 :latex)

