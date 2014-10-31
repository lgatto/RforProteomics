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
	make RforProteomics.pdf
	make poster
	make RProtVis.html
	# make clean


setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif


RforProteomics.pdf: 
	"$(R_HOME)/bin/R" CMD pdflatex RforProteomics.Rnw
	cp vigsrc/RforProteomics.R .
	cp vigsrc/RforProteomics.pdf .

poster:
	cp vigsrc/Bioc2013_poster.pdf ./Bioc2013poster.pdf 
	cp vigsrc/HUPO2014_rbioc_poster.pdf ./HUPO2014poster.pdf

RProtVis.html: RProtVis.Rmd
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2html('RProtVis.Rmd'); purl('RProtVis.Rmd');"

ProteomicsWorkflow.html: ProteomicsWorkflow.Rmd
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2html('ProteomicsWorkflow.Rmd'); purl('ProteomicsWorkflow.Rmd');"


viz2: 
	if [ -d 'vigsrc/cache' ]; then  \
		ln -sf vigsrc/cache .;  \
		ln -sf vigsrc/figure .; \
		ln -sf vigsrc/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML .; \
	fi
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2html('RProtVis.Rmd');"


.PHONY: clean allclean 

cache:
	mv cache/ vigsrc/.
	mv figure/ vigsrc/.
	mv TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML P00924.fasta F063721.dat-mztab.txt vigsrc/.

stub:
	cp vigsrc/RforProteomics-stub.Rnw ./RforProteomics.Rnw 


vig2: 
	## Only locally, using default R-devel
	# make allclean
	cp vigsrc/RforProteomics.Rnw .
	cp vigsrc/RforProteomics.bib .
	ln -s vigsrc/cache .
	ln -s vigsrc/figure .
	ln -s vigsrc/P00924.fasta .
	ln -s vigsrc/F063721.dat-mztab.txt .
	ln -s vigsrc/PXD000001_mztab.txt
	ln -s vigsrc/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML .
	"$(R_HOME)/bin/R" CMD Sweave --engine=knitr::knitr --pdf RforProteomics.Rnw
	bibtex RforProteomics
	"$(R_HOME)/bin/R" CMD pdflatex RforProteomics.tex
	"$(R_HOME)/bin/R" CMD pdflatex RforProteomics.tex
	"$(R_HOME)/bin/R" --vanilla -e "knitr::purl('RforProteomics.Rnw');"
	cp RforProteomics.pdf vigsrc/.
	mv RforProteomics.R vigsrc/.
	## scp RforProteomics.pdf lg390@proteome.sysbiol.cam.ac.uk:/var/www/htdocs/lgatto/RforProteomics/.
	cp vigsrc/RforProteomics-stub.Rnw RforProteomics.Rnw
	rm RforProteomics.bib RforProteomics.tex

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -f Bioc2013_poster.pdf HUPO2014_poster.pdf
	rm -f *~
	rm -f .R
	rm -f */*~
	rm -f PRIDE_Exp_Complete_Ac_22134.xml
	rm -f F063721.dat-mztab.txt TMT*mzXML* P00924.fasta
	rm -f PXD000001_mztab.txt
	rm -rf output*xml
	rm -f RforProteomics.tex 
	rm -f RProtVis.md RProtVis.R
	rm -f Bioc2013_poster.tex HUPO2014_poster.tex
	rm -f Bioc2013_poster.R HUPO2014_poster.R
	rm -f ProteomicsWorkflow.md ProteomicsWorkflow.html 

allclean: 
	make clean 
	rm -fr cache 
	rm -fr figure
	rm -f RforProteomics.R RProtVis.R ProteomicsWorkflow.R

srcclean: 
	rm -fr vigsrc/cache vigsrc/figure
	rm -f vigsrc/F063721.dat-mztab.txt
	rm -f vigsrc/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
	rm -f vigsrc/P00924.fasta
