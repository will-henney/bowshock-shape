# Paper Guide

This folder contains the neccesary elements to produce my first scientific publication in an
astronomical journal. Here briefly describe each one of them.

- Elements with extensions *.tex and *.bib

Actually the written paper. The main document is *proplyd-bowshocks.tex*, while the
documents whose name begin with *sec* are the corresponding sections of the article.
*bowshocks-biblio.bib* is our references library.

- Elements with extension *.pdf

The outputs: The file *proplyd-bowshocks.pdf* is the actual paper, while the rest of them
are figures.

- Elements with extension *.py

Some python scripts which creates some figures for the article. 
  - **conic_fig.py** produces the outputs **conic1.pdf** (default) and **conic2.pdf* *(option --fig 1)
  - **shape_beta.py** produces the outputs **r-beta.pdf** (default) and **ch-radii.pdf** (option --fig 1)
  - The rest of the figures are made or edited via Inkscape
  - **bowshocks-toolbox.py** is a library of functions that will help me to do a set of other figures.
  - **space-param-tab.py** produces a json file which could be as follows:
    - **RcVsR0.json** gives a dictionary which contains $R_c$, $R_0$ and $i$ for a sample of $\beta$ values
    - **Rc-R0Vsbeta.json** gives a dictionary which contains $R_c, $R_0$ and $\beta$ for a sample of inclinations
  - **json-plot.py** makes a plot using the json files mentioned above. The choices of axis are made via 
    command line 
