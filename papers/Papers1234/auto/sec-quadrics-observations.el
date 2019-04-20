(TeX-add-style-hook
 "sec-quadrics-observations"
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
   (LaTeX-add-labels
    "sec:mid-infrared-arcs"
    "sec:comp-with-observ"
    "fig:mipsgal-examples"
    "sec:autom-trac-fitt"
    "step:R0"
    "step:Rc"
    "step:R90"
    "step:Pi-Lambda"
    "sec:subj-eval-fit"
    "fig:mipsgal-shapes"
    "sec:ob-shapes"
    "sec:corr-size"
    "fig:mipsgal-pairplot"
    "fig:mipsgal-uncorrelated"
    "sec:corr-shape"
    "fig:mipsgal-correlated"
    "fig:mipsgal-boxplot"
    "sec:far-infrared-arcs"
    "fig:herschel-arc-fits"
    "fig:herschel-arc-fits-poor"
    "fig:herschel-compare-mipsgal"
    "sec:stat-emiss-line"
    "fig:ll-arcs"
    "fig:ll-compare-mipsgal"))
 :latex)

