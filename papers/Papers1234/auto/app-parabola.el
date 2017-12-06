(TeX-add-style-hook
 "app-parabola"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "app:parabola"
    "eq:parabola-xy"
    "eq:parabola-xy-prime-phi"
    "eq:parabola-xy-prime-final"
    "eq:parabola-R0-prime"
    "eq:parabola-xy-all-primes"
    "eq:parabola-Pi-prime"
    "eq:parabola-t-prime"
    "eq:parabola-Lambda-prime"))
 :latex)

