(TeX-add-style-hook
 "app-dust-equations"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:equat-moti-grains"
    "eq:s-velocity"
    "eq:ds79"
    "eq:G0"
    "eq:G2"
    "fig:divergent-dragoids"
    "eq:dust-XY"
    "eq:dust-UV"
    "eq:dust-motion"
    "eq:dust-gas-velocities"
    "eq:dust-R1"
    "eq:dust-gas-density"))
 :latex)

