LATEXFILES = *.aux\
	*.bbl\
	*.blg\
	*.ilg\
	*.log\
	*.nlo\
	*.nls\
	*.toc\
	*.aux\
	*.dvi\
	*.map\
	*.out\
	*.dep\
	*.dpth

all: 
	make vig
	make viz
	make poster
	make clean


setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif


vig:
	"$(R_HOME)/bin/R" CMD pdflatex RforProteomics.Rnw
	cp vigsrc/RforProteomics.pdf .

poster:
	"$(R_HOME)/bin/R" CMD pdflatex Bioc2013_poster.Rnw
	cp vigsrc/Bioc2013_poster.pdf .

viz: 
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2html('ProtViz.Rmd');"

.PHONY: clean allclean 

cache:
	mv cache/ vigsrc/.
	mv figure/ vigsrc/.
	mv TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML P00924.fasta F063721.dat-mztab.txt vigsrc/.

stub:
	cp vigsrc/RforProteomics-stub.Rnw ./RforProteomics.Rnw 


vig2: 
	## Only locally, using default R-devel
	rm -f $(LATEXFILES)
	cp vigsrc/RforProteomics.Rnw .
	cp vigsrc/RforProteomics.bib .
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
	"$(R_HOME)/bin/R" CMD pdflatex RforProteomics.tex
	cp RforProteomics.pdf vigsrc/.
	## scp RforProteomics.pdf lg390@proteome.sysbiol.cam.ac.uk:/var/www/htdocs/lgatto/RforProteomics/.
	cp vigsrc/RforProteomics-stub.Rnw RforProteomics.Rnw
	rm RforProteomics.bib RforProteomics.tex

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -f *~
	rm -f */*~
	rm -f PRIDE_Exp_Complete_Ac_22134.xml
	rm -f F063721.dat-mztab.txt TMT*mzXML* P00924.fasta
	rm -f PXD000001_mztab.txt
	rm -rf output*xml
	rm -f RforProteomics.tex 
	rm -f ProtViz.md 

allclean: clean 
	rm -fr cache 
	find -path "./figure/*" ! -name shinyMA.png -delete

srcclean: 
	rm -fr vigsrc/cache vigsrc/figure 
	rm -f vigsrc/F063721.dat-mztab.txt
	rm -f vigsrc/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
	rm -f vigsrc/P00924.fasta
