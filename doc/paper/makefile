all: percolator.pdf diff.pdf coverletter.pdf

percolator.pdf: percolator.tex percolator.bib img/*

diff.tex: percolator.tex
	git show v1:./$? > /tmp/$?
	latexdiff /tmp/$? $? > diff.tex
	
%.pdf: %.tex
	pdflatex -shell-escape $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		pdflatex -shell-escape $*.tex ; \
	fi
	pdflatex -shell-escape $*.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *~
