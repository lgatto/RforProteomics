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
	make vignettes
	make clean

setvars:
ifeq (${R_HOME},)
## R_HOME=	$(shell R RHOME)
R_HOME=/opt/Rpatched/lib/R
endif

vignettes: 
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('RforProteomics.Rnw');"
	## bibtex RforProteomics.tex
	"$(R_HOME)/bin/R" CMD pdflatex RforProteomics.tex

clean:
	rm -f $(LATEXFILES)
	rm -f *~
