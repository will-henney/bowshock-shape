(TeX-add-style-hook
 "app-dust-ifront-trap"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:furth-deta-ioniz"
    "eq:shell-recombination-rate"
    "eq:shell-ionizing-flux"
    "eq:shocked-shell-column"
    "fn:temperature-dependence"
    "eq:ifront-trap-x-cubed-taustar"
    "eq:ifront-trap-taustar-bow-shock"
    "eq:ion-tau-gas"
    "eq:ion-tau-gas-expanded"
    "eq:outer-shell-ionization-balance"
    "eq:sigma-vs-tau"))
 :latex)

