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
    "fig:All-sources-eta-tau"
    "sec:energy-trapp-vers"
    "eq:tau-empirical"
    "P1"
    "P2"
    "P3"
    "eq:eta-shell"
    "sec:eta-tau-diagnostic"
    "tab:observations"
    "sec:rand-syst-uncert"
    "sec:syst-uncert-due"
    "sec:spec-treatm-part"
    "sec:cand-radi-supp"
    "sec:stellar-wind-mass"
    "fig:mass-loss-vs-luminosity"
    "sec:case-inside-out"
    "sec:conclusions"))
 :latex)

