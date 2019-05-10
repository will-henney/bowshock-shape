(TeX-add-style-hook
 "app-p-values"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "figs/mipsgal-summary-stats-table-body")
   (LaTeX-add-labels
    "sec:distr-p-values"
    "fig:histo-p-values"
    "tab:big-p"
    "sec:perturbed-bows"
    "fig:perturb-shapes"
    "fig:perturb-xy-prime"
    "eq:standing-wave"
    "eq:fractional-phase"))
 :latex)

