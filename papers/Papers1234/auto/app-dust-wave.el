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
   (TeX-add-symbols
    "Larmor")
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
    "eq:perpendicular-drift"
    "eq:vdrift-over-vinfinity"
    "eq:parallel-energy-balance"
    "eq:parallel-radiation-potential"
    "eq:thB-0-shape"
    "fig:dust-coupling-1d"
    "fig:dust-wave-coupling"
    "sec:bow-wave-drag"
    "eq:dust-fdrag"
    "eq:dust-wdrift"
    "eq:dust-alpha"
    "sec:bow-wave-with"
    "fig:inertia-thB90"
    "fig:magnetic-dust-waves"
    "fig:projected-magnetic-dust-waves-15"
    "fig:projected-magnetic-dust-waves-45"
    "fig:projected-magnetic-dust-waves-75"
    "fig:stream-3d"
    "sec:resolv-larm-radi"
    "sec:dust-applicability"
    "sec:dust-wave-apparent"
    "fig:dragoid-xy-prime"
    "fig:dragoid-div-xy-prime"
    "fig:dragoid-Rc-R90"
    "sec:shape-bow-wave"))
 :latex)

