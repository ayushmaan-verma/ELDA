.RECIPEPREFIX = >

TEX = elda_booklet.tex
OUTDIR = build
PDF = $(OUTDIR)/elda_booklet.pdf

.PHONY: pdf clean distclean

pdf:
> mkdir -p $(OUTDIR)
> latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=$(OUTDIR) $(TEX)

clean:
> latexmk -c -outdir=$(OUTDIR) $(TEX)

distclean: clean
> rm -f $(PDF)
