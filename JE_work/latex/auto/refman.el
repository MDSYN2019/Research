(TeX-add-style-hook
 "refman"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("book" "twoside")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("adjustbox" "export") ("inputenc" "utf8") ("wasysym" "nointegrals") ("xcolor" "table") ("fontenc" "T1") ("helvet" "scaled=.90") ("tocloft" "titles") ("hyperref" "pdftex" "pagebackref=true" "ps2pdf")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "annotated"
    "classJarzynskiFreeEnergy"
    "book"
    "bk10"
    "fixltx2e"
    "calc"
    "doxygen"
    "adjustbox"
    "graphicx"
    "inputenc"
    "makeidx"
    "multicol"
    "multirow"
    "textcomp"
    "wasysym"
    "xcolor"
    "fontenc"
    "helvet"
    "courier"
    "amssymb"
    "sectsty"
    "geometry"
    "fancyhdr"
    "natbib"
    "tocloft"
    "ifpdf"
    "hyperref"
    "caption")
   (TeX-add-symbols
    "clearemptydoublepage"))
 :latex)

