(TeX-add-style-hook
 "quadrics-bowshock"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("mnras" "useAMS" "usenatbib" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("newtxmath" "varvw" "smallerops")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "sec-intro"
    "sec-projection"
    "sec-conic"
    "sec-thin-shell"
    "sec-conclusions"
    "app-intermediate"
    "app-parabola"
    "app-shape-parameters"
    "app-rcurv-empirical"
    "mnras"
    "mnras10"
    "graphicx"
    "dblfloatfix"
    "microtype"
    "xcolor"
    "fixltx2e"
    "booktabs"
    "siunitx"
    "color"
    "enumerate"
    "pdflscape"
    "rotating"
    "nicefrac"
    "dsfont"
    "bbm"
    "fontenc"
    "inputenc"
    "newtxtext"
    "newtxmath"
    "hyperref"
    "xr-hyper"
    "etoolbox"
    "bm"
    "aastex-compat")
   (TeX-add-symbols
    '("TODO" 1)
    '("uvec" 1)
    '("Abs" 1)
    '("abs" 1)
    "hmmax"
    "bmmax"
    "hii"
    "AddressCRyA"
    "ecc"
    "w"
    "C"
    "T")
   (LaTeX-add-labels
    "firstpage"
    "lastpage")
   (LaTeX-add-environments
    "Vector")
   (LaTeX-add-bibliographies
    "bowshocks-biblio"))
 :latex)

