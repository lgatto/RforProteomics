R_HOME = "$(shell R RHOME)"

RforProteomics.html: RforProteomics.Rmd
	"$(R_HOME)/bin/Rscript" -e 'rmarkdown::render("RforProteomics.Rmd")'

pkgsdownloads.html: pkgsdownloads.Rmd
	"$(R_HOME)/bin/Rscript" -e 'rmarkdown::render("pkgsdownloads.Rmd")'

clean: 
	rm -fr auto cache figure
	rm -f output.20*
	rm -f P00924.fasta PXD000001_mztab.txt
	rm -f TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML
	rm -f TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
	rm -rf RProtVis_cache RforProteomics_cache biocprot_cache

.PHONY: clean
