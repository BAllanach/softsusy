docdir = ${datadir}/doc/${PACKAGE}

doc_DATA = nmssmManual.pdf \
           rpcManual.pdf \
           rpvManual.pdf \
           neutManual.pdf

EXTRA_DIST = \
	$(doc_DATA) \
	rpcManual.tex \
	rpvManual.tex \
	neutManual.tex \
	nmssmManual.tex

doc: $(doc_DATA)

%.pdf: %.tex
	latex $<
	latex $<
	dvips -f $(*F) > $(*F).ps
	ps2pdf $(*F).ps > $(*F).pdf

clean-local:
	-rm -f *.log *.spl *.aux *.dvi *.ps *.toc *.out

distclean-local: clean-local
	-rm -f $(doc_DATA)
