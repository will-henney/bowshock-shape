(TeX-add-style-hook
 "sec-projection"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:projection"
    "sec:ref-frames"
    "eq:body-frame"
    "eq:Trans"
    "fig:projection-pos"
    "sec:unit-vectors-normal"
    "fig:unitvec"
    "eq:uvec-phi0"
    "eq:alpha"
    "eq:omega"
    "eq:nprime"
    "eq:tprime"
    "sec:tangent-line"
    "eq:tangent-line-condition"
    "eq:tanphi"
    "eq:tangential"
    "eq:R-prime-theta-prime"
    "eq:thetapar"
    "eq:Rpar"
    "eq:R90prime"
    "eq:th90"
    "eq:projected-radius-curvature"
    "sec:line-sight-veloc"
    "eq:vlos"))
 :latex)

