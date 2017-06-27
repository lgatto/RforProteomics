RforProteomics.html: RforProteomics.Rmd
	Rscript -e 'rmarkdown::render("RforProteomics.Rmd")'

clean: 
	rm -fr auto cache figure
	rm -f output.20*
	rm -f P00924.fasta PXD000001_mztab.txt
	rm -rf RforProteomics_cache RforProteomics_files 
	rm -f TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML
	rm -f TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
	rm -rf RProtVis_files

.PHONY: clean
