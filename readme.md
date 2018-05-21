# Ularcirc
An R package that provides analysis and visualisation of canonical and backsplice junctions.
Requires output provided by the STAR aligner

Author and maintainer: David Humphreys (d.humphreys  at      victorchang  dot   edu   dot    au)

Pre-print for [Ularcirc now on bioarchives](https://www.biorxiv.org/content/early/2018/05/15/318436). 

# Installation
You can install Ularcirc using the 'devtools' package

    > install.packages("devtools")
    > library(devtools)
    > devtools::install_github("VCCRI/Ularcirc")

Ularcirc can annotate circRNA with overlapping gene information. This is obtained from available 
bioconductor databases. Use the following command to identify what databases to download:

    > library("Ularcirc")
    > all_dbs <- Compatible_Annotation_DBs()
    > # List database IDs
    > names(all_dbs)
    >
    > # Select a database and display the commands needed to install
    > # Use noquote to correctly format output
    > noquote(all_dbs[['Hsapiens.UCSC.hg38']])
    
    
To start Ularcirc shiny app

    > library('Ularcirc')
	  > Ularcirc()

	
# Documentation
Please refer to pdf manual under the docs folder. This document is currently being built so please check back for regular updates.

Please click [this link to view a 12 minute screen cast](https://www.youtube.com/watch?v=BQkCurUF2XQ) that walks through a simple circRNA analysis using Ularcirc.


##  Features

* Friendly user interface
* Circular RNA detection independent of gene annotation.
* Provides visualisation of forward AND backsplice junctions
* Recover predicted circRNA sequence
* Recover sequence of backsplice junctions and forward splice junctions
* Support both single-read and paired-end sequencing (paired end prefered).
* Detect miRNA binding sites
* detect putative open reading frame of circRNA
