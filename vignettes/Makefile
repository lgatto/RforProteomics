patched=0

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
ifeq (${patched},1)
R_HOME=/opt/Rpatched
else 
R_HOME=	$(shell R RHOME)
endif

vignette:
	make clean
	make stub
	"$(R_HOME)/bin/R" CMD Sweave RforProteomics.Rnw
	cp src/RforProteomics.pdf .
	cp src/RforProteomics.Rnw .

build: 
	rm -f $(LATEXFILES)
	cp src/RforProteomics.Rnw .
	if [ -d 'src/cache' ]; then \
		ln -sf src/cache .; \
		ln -sf src/figure .; \
	fi
	if [ -f 'src/F063721.dat-mztab.txt' ]; then \
		ln -sf src/P00924.fasta .; \
		ln -sf src/F063721.dat-mztab.txt .; \
		ln -sf src/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML .; \
	fi
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('RforProteomics.Rnw');"
	bibtex RforProteomics
	"$(R_HOME)/bin/R" CMD pdflatex RforProteomics.tex
	cp RforProteomics.pdf src/.
	make stub

stub:
	cp src/RforProteomics-stub.Rnw RforProteomics.Rnw 

clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -f RforProteomics.R
	rm -f PRIDE_Exp_Complete_Ac_22134.xml
	rm -f F063721.dat-mztab.txt TMT*mzXML* P00924.fasta

allclean:
	make clean
	rm -fr cache figure 
