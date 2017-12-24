(TeX-add-style-hook
 "sec-dust-bow-intro"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-add-symbols
    "PaperI"
    "hii"
    "alphaB"
    "shell")
   (LaTeX-add-labels
    "sec:introduction"
    "fig:3-types-bow"
    "sec:different-types-bow"
    "fig:zones-v-n-plane"
    "sec:strong-gas-grain"
    "eq:wind-efficiency"
    "eq:rad-press-balance-thick"
    "eq:Rstar"
    "eq:rad-accel"
    "eq:rad-poten"
    "eq:rad:R0"
    "eq:rad-press-balance-tau"
    "eq:tau-thin"
    "eq:tau-star"
    "eq:rad-full-x"
    "eq:x-cases"
    "tab:stars"
    "eq:stellar-parameters"
    "eq:wind-eta-typical"
    "eq:fiducial-typical"
    "fn:meyer-velocities-too-low"
    "sec:trapp-ioniz-front"
    "eq:shell-recombination-rate"
    "eq:shell-ionizing-flux"
    "eq:shocked-shell-column"
    "fn:temperature-dependence"
    "eq:ifront-trap-x-cubed-taustar"
    "eq:ifront-trap-taustar-bow-shock"
    "eq:ion-tau-gas"
    "eq:ion-tau-gas-expanded"
    "eq:outer-shell-ionization-balance"
    "eq:sigma-vs-tau"
    "sec:radi-cool-lengths"
    "sec:imperf-coupl-betw"
    "fig:decouple-v-n-plane"
    "fig:decouple-v40-versus-n"
    "sec:case-inside-out"))
 :latex)

