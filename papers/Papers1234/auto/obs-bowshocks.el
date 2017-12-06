(TeX-add-style-hook
 "obs-bowshocks"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("mnras" "useAMS" "usenatbib" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("newtxmath" "varg")))
   (TeX-run-style-hooks
    "latex2e"
    "sec-obs-intro"
    "sec-quadrics-observations"
    "sec-obs-conclusions"
    "app-p-values"
    "mnras"
    "mnras10"
    "newtxmath"
    "newtxtext"
    "graphicx"
    "microtype"
    "xcolor"
    "fixltx2e"
    "booktabs"
    "siunitx"
    "color"
    "enumerate"
    "pdflscape"
    "rotating"
    "hyperref"
    "musixtex"
    "fontawesome"
    "staves"
    "etoolbox"
    "bm"
    "aastex-compat")
   (TeX-add-symbols
    "textfermata"
    "hmmax"
    "bmmax"
    "hii"
    "AddressCRyA"
    "extractline")
   (LaTeX-add-bibliographies
    "bowshocks-biblio"))
 :latex)

