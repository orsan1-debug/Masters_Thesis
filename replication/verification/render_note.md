`thesis_render.pdf` was made on 15 Oct 2026. `verify_render.R` ran on the Mac, which has no TeX,
so Quarto stopped after writing the LaTeX, kept in `render_tex/`. That LaTeX was compiled with
XeLaTeX (TeX Live 2022) in a Linux virtual machine on the same Mac, with the Mac's Times New
Roman. The stored PDF was compiled elsewhere, with a later TeX Live (its xdvipdfmx is dated
2026-03-17, the virtual machine's 2021-06-09). In a few places a line breaks at another point,
so the comparison goes by numbers and words rather than layout.
