LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
	Rplots.pdf\
	*.dvi\
	*.map\
        *.out\
	RforProteomics.tex

all:
	make vignette
	make clean

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

vignette:
	make stub
	"$(R_HOME)/bin/R" CMD Sweave RforProteomics.Rnw

stub:
	cp vigsrc/RforProteomics-stub.Rnw RforProteomics.Rnw 

build: 
	## Only locally
	rm -f $(LATEXFILES)
	cp vigsrc/RforProteomics.Rnw .
	if [ -d 'vigsrc/cache' ]; then \
		ln -sf vigsrc/cache .; \
		ln -sf vigsrc/figure .; \
	fi
	if [ -f 'vigsrc/F063721.dat-mztab.txt' ]; then \
		ln -sf vigsrc/P00924.fasta .; \
		ln -sf vigsrc/F063721.dat-mztab.txt .; \
		ln -sf vigsrc/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML .; \
	fi
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('RforProteomics.Rnw');"
	bibtex RforProteomics
	"$(R_HOME)/bin/R" CMD pdflatex RforProteomics.tex
	scp RforProteomics.pdf lg390@proteome.sysbiol.cam.ac.uk:/var/www/htdocs/lgatto/RforProteomics/.
	cp RforProteomics.pdf vigsrc/.
	make stub


clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -f RforProteomics.R
	rm -f PRIDE_Exp_Complete_Ac_22134.xml
	rm -f F063721.dat-mztab.txt TMT*mzXML* P00924.fasta

allclean:
	make clean
	rm -f cache figure 

