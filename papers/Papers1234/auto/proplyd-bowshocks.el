(TeX-add-style-hook
 "proplyd-bowshocks"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("mnras" "useAMS" "usenatbib")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("newtxmath" "varg")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "sec-intro-proplyd"
    "sec-application"
    "sec-conclusions-proplyd"
    "amsmath"
    "mnras"
    "mnras10"
    "newtxmath"
    "newtxtext"
    "graphicx"
    "microtype"
    "xcolor"
    "fixltx2e"
    "booktabs"
    "hyperref"
    "siunitx"
    "color"
    "bm"
    "aastex-compat")
   (TeX-add-symbols
    "hmmax"
    "bmmax"
    "AddressCRyA"
    "thC"
    "CRW")
   (LaTeX-add-bibliographies
    "bowshocks-biblio"))
 :latex)

