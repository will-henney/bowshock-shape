(TeX-add-style-hook
 "sec-conic"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-add-symbols
    "Q"
    "thetaQ"
    "fQi")
   (LaTeX-add-labels
    "fig:conics"
    "fig:conics-family"
    "sec:conic"
    "eq:conic-x0"
    "eq:par-xy"
    "eq:sin-sinh-etc"
    "eq:Tq"
    "eq:thetaQ"
    "fig:conic-departure"
    "fig:quadric-projection"
    "fig:quadric-projection-continued"
    "eq:Rc-conic"
    "eq:conic-R90"
    "eq:R0-from-Q-a-x0"
    "eq:Pi-from-Q-a-x0"
    "eq:Lambda-from-Q-a-x0"
    "eq:Tq-from-Pi-Lambda"
    "sec:parab-depart-funct"
    "eq:confocal-polar"
    "eq:departure-function"
    "eq:xyz-XYZ"
    "eq:quadric-XYZ"
    "eq:phit-quadric"
    "eq:conic-projected-XY"
    "eq:conic-projected-XY-conic"
    "eq:a-prime"
    "eq:b-prime"
    "eq:fQi-factor"
    "eq:Q-prime"
    "eq:XYZ-xyz-prime"
    "eq:x0-prime"
    "fig:projected-R90-Rc-snapshots"
    "eq:R0-prime"
    "eq:Pi-prime"
    "eq:Lambda-prime"))
 :latex)

