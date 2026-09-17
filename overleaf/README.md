# overleaf/ — the manuscript, and only the manuscript

Twenty files, 19 MB. This is everything `main.tex` needs to compile and nothing
else. Upload this folder to Overleaf, or push it to its own small GitHub
repository and sync Overleaf with that.

    main.tex
    rilabRxiv.cls      the document class, which sets \bibliographystyle{genetics}
    genetics.bst
    SoLD.bib
    logo.png
    fig/main/          the 6 main figures, at the paths \includegraphics already uses
    fig/supp/          the 9 supplementary figures

All fifteen `\includegraphics` paths in `main.tex` were checked against this
folder and every one resolves. Nothing in `main.tex` uses `\input` or
`\include`, so there are no other dependencies.

## Why this exists

Overleaf's GitHub sync was failing on the full repository, first refusing a
filename and then returning a bare `GithubApiError`. The cause is size, not any
one file. `data/LD_mapped` alone holds 582 files and 863 MB; `scripts/` adds
about 190 more; `table/`, `fig/`, `output/`, `docs/`, `www/`, `tables/`,
`new_new_figs/`, `new_new_tables/`, `figures_old/` and `rsconnect/` add
hundreds more on top. Overleaf caps a project at 2,000 files and its GitHub
sync gives up well before a repository of that size finishes transferring.

None of it belongs in Overleaf anyway. LaTeX reads a `.tex`, a `.cls`, a `.bst`,
a `.bib` and some images. It cannot read an `.xlsx`, an `.R` script or a `.csv`
of candidate genes, so syncing them buys nothing and costs the whole project.

## What is deliberately absent

The ten supplementary workbooks in `final/table/workbooks/`. LaTeX cannot
include a spreadsheet, and the manuscript cites Supplementary Tables S1 to S10
by name rather than by path. They are journal submission files, uploaded through
the submission portal alongside the compiled PDF, not compile inputs.

## Keeping it current

`main.tex` here is a copy. When you edit it in Overleaf, that copy is the live
one and the root `main.tex` falls behind — decide which is authoritative and do
not edit both. The figures are copies of `final/fig/`, so after rebuilding any
figure, re-run `scripts/new_new_script/99_assemble_final.R` and then copy
`final/fig/` over `overleaf/fig/` again.
