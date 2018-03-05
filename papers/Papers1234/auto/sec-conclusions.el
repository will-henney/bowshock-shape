(TeX-add-style-hook
 "sec-conclusions"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:more-realistic-bow"
    "fig:meyer-trace"
    "fig:sim-depart"
    "fig:sim-xyp"
    "fig:sim-Pi-Lambda"
    "fig:sim-R0-prime"
    "fn:meyer-correction"
    "fig:sim-histograms"
    "fig:000-400-fit"
    "fig:000-400-planitude-alatude"
    "fig:000-400-Delta-theta"
    "tab:m42-fit"
    "sec:obs"
    "eq:alatude-average-and-asym"
    "sec:conc"))
 :latex)

