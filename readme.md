# Ularcirc
An R package that provides analysis and visualisation of canonical and backsplice junctions.
Takes output provided by the STAR aligner as well as CIRI2 and circExplorer2 output and enables circRNA downstream analysis.

Author and maintainer: David Humphreys (d.humphreys  at      victorchang  dot   edu   dot    au)

Ularcirc manuscript now available through [Nucleic Acids Research](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz718/5552786). 

# Installation
You can install Ularcirc using the 'devtools' package.  

    > install.packages("devtools")
    > library(devtools)
    > devtools::install_github("VCCRI/Ularcirc", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))

Ularcirc can annotate circRNA with overlapping gene information. This is obtained from available 
bioconductor databases. Use the following command to identify what databases to download:

    > library("Ularcirc")
    > all_dbs <- Compatible_Annotation_DBs() # This will return all compatible databases
    > mmu_dbs <- Compatible_Annotation_DBs(search_term = 'mm10') # returns mm10 compatible databases
    > # Lets see what is stored in mmu_dbs
    > mmu_dbs
   annotation     genome                         txdb                                
16 "org.Mm.eg.db" "BSgenome.Mmusculus.UCSC.mm10" "TxDb.Mmusculus.UCSC.mm10.ensGene"  
17 "org.Mm.eg.db" "BSgenome.Mmusculus.UCSC.mm10" "TxDb.Mmusculus.UCSC.mm10.knownGene"
    
    > # Now lets download all of the above databases
    > if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")   # Make sure R is looking at bioconductor repository
    > BiocManager::install(c(mmu_dbs))
    
    
To start Ularcirc shiny app

    > library('Ularcirc')
	  > Ularcirc()

	
# Documentation
Please refer to vignette within R. Additionally  there are a number of screen casts that highlights how to get going with Ularcirc.

##Screen casts

Please click [this link to view a ~5 minute screen cast](https://youtu.be/96rcxlh_aLA) that walks through a simple circRNA analysis using Ularcirc.

The following link demonstrates how to [upload and recover sequence information from BSJ and FSJ](https://youtu.be/txWAI-LJCVw)


##  Features

* Friendly user interface
* Circular RNA detection independent of gene annotation.
* Provides visualisation of forward AND backsplice junctions
* Recover predicted circRNA sequence
* Recover sequence of backsplice junctions and forward splice junctions
* Support both single-read and paired-end sequencing (paired end prefered).
* Detect miRNA binding sites
* detect putative open reading frame of circRNA
