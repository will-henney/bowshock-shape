(TeX-add-style-hook
 "dusty-bow-wave"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("mnras" "useAMS" "usenatbib" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("newtxmath" "varvw" "smallerops")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "sec-dust-bow-intro"
    "app-dust-wave"
    "sec-standing-wave"
    "sec-dust-case-studies"
    "sec-dust-bow-discuss"
    "app-dust-equations"
    "mnras"
    "mnras10"
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
    "xr-hyper"
    "hyperref"
    "fontenc"
    "inputenc"
    "newtxtext"
    "newtxmath"
    "etoolbox"
    "bm"
    "aastex-compat")
   (TeX-add-symbols
    '("uvec" 1)
    '("TODO" 1)
    '("Abs" 1)
    '("abs" 1)
    "hmmax"
    "bmmax"
    "AddressCRyA"
    "sgn"
    "Sin"
    "Cos"
    "Cot"
    "GammaFunc"
    "w"
    "C"
    "T"
    "Qp"
    "grain"
    "xsec"
    "frad"
    "thm"
    "drag"
    "gas"
    "drift"
    "sound"
    "soundspeed")
   (LaTeX-add-labels
    "firstpage"
    "lastpage")
   (LaTeX-add-bibliographies
    "bowshocks-biblio"))
 :latex)

