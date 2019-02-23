(TeX-add-style-hook
 "sec-dust-bow-discuss"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-add-symbols
    "pc")
   (LaTeX-add-labels
    "sec:summary-discussion"
    "eq:tau-empirical"
    "P1"
    "P2"
    "P3"
    "eq:eta-shell"
    "tab:observations"
    "fig:All-sources-eta-tau"
    "fig:mass-loss-vs-luminosity"
    "sec:case-inside-out"
    "sec:conclusions"))
 :latex)

