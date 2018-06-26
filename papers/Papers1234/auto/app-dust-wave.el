(TeX-add-style-hook
 "app-dust-wave"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:gas-free-bow"
    "fig:dust-trajectories"
    "eq:dust-rad-force"
    "fn:Qp"
    "eq:dust-r0"
    "eq:dust-r-theta"
    "sec:dust-parallel"
    "eq:dust-r-in"
    "sec:dust-divergent"
    "eq:dust-divergent-r-in"
    "sec:tight-magn-coupl"
    "fig:inertia-thB0"
    "eq:perpendicular-drift"
    "eq:vdrift-over-vinfinity"
    "sec:parall-magn-field"
    "eq:parallel-energy-balance"
    "eq:parallel-radiation-potential"
    "eq:thB-0-shape"
    "eq:thB-0-density"
    "sec:perp-magn-field"
    "fig:inertia-thB90"
    "eq:ode-perp-bfield"))
 :latex)

