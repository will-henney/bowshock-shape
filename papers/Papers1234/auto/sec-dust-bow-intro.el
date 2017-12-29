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
    "shell"
    "M"
    "cool")
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
    "eq:Rstar-typical"
    "eq:taustar-typical"
    "fn:meyer-velocities-too-low"
    "fig:M-supergiant"
    "sec:trapp-ioniz-front"
    "fig:B-supergiant"
    "eq:shell-recombination-rate"
    "eq:shell-ionizing-flux"
    "eq:shocked-shell-column"
    "eq:isothermal-shell-density"
    "fn:temperature-dependence"
    "eq:ifront-trap-x-cubed-taustar"
    "eq:ifront-trap-density-RBS"
    "eq:ifront-trap-taustar-RBW"
    "sec:radi-cool-lengths"
    "eq:shock-n-jump"
    "eq:shock-T-jump"
    "eq:shock-v-jump"
    "eq:dcool"
    "eq:cooling-coefficient"
    "eq:heating-coefficient"
    "eq:strong-cooling-h0"
    "sec:imperf-coupl-betw"
    "fig:decouple-v-n-plane"
    "fig:decouple-v40-versus-n"
    "sec:case-inside-out"))
 :latex)

