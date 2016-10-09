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

R = /opt/Rpatched/bin/R

all: 
	make RforProteomics.pdf
	make poster

RforProteomics.pdf: 
	cp vigsrc/RforProteomics.R .
	cp vigsrc/RforProteomics.pdf .

vis:
	cp vigsrc/RProtVis.html .
	cp vigsrc/RProtVis.R .

poster:
	cp vigsrc/Bioc2013_poster.pdf ./Bioc2013poster.pdf 
	cp vigsrc/HUPO2014_rbioc_poster.pdf ./HUPO2014poster.pdf

.PHONY: clean allclean stub vig2 vis vis2

stub:
	cp vigsrc/RforProteomics-stub.Rnw ./RforProteomics.Rnw 
	cp vigsrc/RProtVis-stub.Rmd ./RProtVis.Rmd


vig2: 
	## Only locally, using default R-devel
	# make allclean
	cp vigsrc/RforProteomics.Rnw .
	cp vigsrc/RforProteomics.bib .
	# ln -s vigsrc/cache .
	# ln -s vigsrc/figure .
	# ln -s vigsrc/F063721.dat-mztab.txt .
	# ln -s vigsrc/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML .
	# ln -s vigsrc/c_elegans_A_3_1_21Apr10_Draco_10-03-04.mzML vigsrc/c_elegans_A_3_3_21Apr10_Draco_10-03-04.mzML   \
	# 	vigsrc/c_elegans_B_2_1_21Apr10_Draco_10-03-05.mzML vigsrc/c_elegans_B_2_3_21Apr10_Draco_10-03-05.mzML \
	# 	vigsrc/c_elegans_C_1_1_21Apr10_Draco_10-03-06.mzML vigsrc/c_elegans_C_1_3_21Apr10_Draco_10-03-06.mzML \
	# 	vigsrc/c_elegans_D_1_1_21Apr10_Draco_10-03-07.mzML vigsrc/c_elegans_D_1_3_21Apr10_Draco_10-03-07.mzML \
	# 	vigsrc/c_elegans_E_3_1_21Apr10_Draco_10-03-04.mzML vigsrc/c_elegans_E_3_3_21Apr10_Draco_10-03-04.mzML \
	# 	vigsrc/c_elegans_fwd_rev.fasta .
	${R} CMD Sweave --engine=knitr::knitr --pdf RforProteomics.Rnw
	bibtex RforProteomics
	${R} CMD pdflatex RforProteomics.tex
	${R} CMD pdflatex RforProteomics.tex
	${R} --vanilla -e "knitr::purl('RforProteomics.Rnw');"
	cp RforProteomics.pdf vigsrc/.
	mv RforProteomics.R vigsrc/.
	cp vigsrc/RforProteomics-stub.Rnw RforProteomics.Rnw
	rm RforProteomics.bib RforProteomics.tex

dataclean:
	rm -f F063721.dat-mztab.txt TMT*mzXML* P00924.fasta TMT*mzML
	rm -f PXD000001_mztab.txt
	rm -f output*xml
	rm -f P00924.fasta
	rm -f c_elegans*

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -f Bioc2013_poster.pdf HUPO2014_poster.pdf
	rm -f *~
	rm -f .R
	rm -f */*~
	rm -f PRIDE_Exp_Complete_Ac_22134.xml
	rm -f RforProteomics.tex 
	rm -f RProtVis.md RProtVis.R
	rm -f Bioc2013_poster.tex HUPO2014_poster.tex
	rm -f Bioc2013_poster.R HUPO2014_poster.R
	rm -f ProteomicsWorkflow.md ProteomicsWorkflow.html
	rm -rf vigsrc/auto
	make dataclean

allclean: 
	make clean 
	rm -fr cache 
	rm -fr .Rcache
	rm -fr RProtVis_cache 
	rm -fr figure
	rm -fr RProtVis_files 
	rm -f RforProteomics.R RProtVis.R ProteomicsWorkflow.R

srcclean: 
	rm -fr vigsrc/cache vigsrc/figure
	rm -f vigsrc/F063721.dat-mztab.txt
	rm -f vigsrc/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
	rm -f vigsrc/P00924.fasta
