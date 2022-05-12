all: percolator.pdf #diff.pdf coverletter.pdf coverletter2.pdf

percolator.pdf: percolator.tex # percolator.bib img/*

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

%-eps-converted-to.pdf: %.eps
	repstopdf --outfile=$@ $?

clean_img: clean
	rm -f img/*converted-to.pdf

clean:
	rm -f *.aux *.bbl *.blg *.log *~
