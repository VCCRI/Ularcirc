library(shiny)
library(data.table)
library(DT)
library(gsubfn)			# Use strapplyc
library(Biostrings)
require(GenomicFeatures)
library(Sushi)
library(moments)

options(shiny.maxRequestSize=900*1024^2)  # Set upper limit at 900MB

#############################################################################################
# ----- function which create the button into the table
## Code copied from http://stackoverflow.com/questions/33540389/r-shiny-delete-button-into-data-table-doesnt-work-well
shinyInput <- function(FUN, len, id, ...) {
  inputs <- len
  for (i in seq(len)) {
    inputs[i] <- as.character(FUN(paste0(id, len[i]), ...))
  }
  inputs
}

############################################################3
## Filter_by_Data_Set
##
##  Will select all records that have corresponding ID associated in DataSet column of data.table All_junctions
##
Filter_by_Data_Set <- function(fileID, All_junctions)
{
  if (fileID[1] != -1)
  {
    for (i in 1:length(fileID))
    {
      if (i == 1)
        temp <- All_junctions[(DataSet==fileID[i]),]
      else
        temp <- rbind (temp,All_junctions[(DataSet==fileID[i]),])
    }
    All_junctions <- temp
  }
  return(All_junctions)
}

########################
##
## ege circJunctions(m379,"chr9",46240800,46243500)  # Apoa4
circJunctions<-function(All_junctions, chrom, chromstart, chromend, fileID = c(-1))     #Need to consider strand!!!
{	cat (paste("\nRequest junctions on chrom",chrom,":",chromstart,"-",chromend))
	# 1 = chromDonor, 2 = startDonor, 3 = strandDonor, 4 = chromAcceptor, 5=startAcceptor, 6=strandAcceptor
	if (length(All_junctions) == 0)
	{
		return (NULL)
	}

	# First select data from files that user has selected
  withProgress(message="Finding junctions to display", # style="old",
               value=0, {
    incProgress(1/3, detail = paste("Short listing all junctions for sample"))
  	if (fileID[1] != -1)
	  {
		  for (i in 1:length(fileID))
		  {
			  if (i == 1)
				  temp <- All_junctions[(DataSet==fileID[i]),]
			  else
				  temp <- rbind (temp,All_junctions[(DataSet==fileID[i]),])
		  }
		  All_junctions <- temp
	  }
    incProgress(1/3, detail = paste("Short listing junctions for gene"))

	# Select junctions within the range defined by chrom, chromstart and chromend.
	## STILL NEED TO MAKE DATA STRANDED !! ie accurately accept illumina or Proton data.
    chromstart <- as.numeric(chromstart)
    chromend <- as.numeric(chromend)
    junc = All_junctions[chromDonor == chrom & chromAcceptor == chrom &
							startDonor >= chromstart & startDonor <= chromend &
							startAcceptor >= chromstart & startAcceptor <= chromend & strandDonor == strandAcceptor,]

  })

	if (nrow(junc) == 0)	# No data !!
	{ return (-1)	}

    BS_juncs <- junc[(strandDonor=='-' & (startAcceptor > startDonor))  | (strandDonor=='+' & (startDonor > startAcceptor) ),]     # Select the backsplice junctions, this dataframe will be to return to user
	  tmp <- junc[,length(BSjuncName),by="chromDonor,startDonor,chromAcceptor,startAcceptor"]
	  setnames(tmp,5,c("score"))
	  junc.bed<-data.table(junc[,.(chromDonor,startDonor,startDonor,chromAcceptor,startAcceptor,startAcceptor, BSjuncName)], 1, junc[,.(strandDonor,strandAcceptor)],500)

    setnames(junc.bed, 1:11,c("chrom1","start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1",  "strand2","JunctionType"))     # Assign names

    junc.uniques.bed <- unique(junc.bed)	            # junc.uniques.bed is the dataframe used to create loop graph
	  junc.uniques.bed$score <-  tmp[chromDonor == junc.uniques.bed$chrom1 & startDonor == junc.uniques.bed$start1 & chromAcceptor== junc.uniques.bed$chrom2 & startAcceptor== junc.uniques.bed$start2,.(score)]
	  junc.uniques.bed$JunctionType[junc.uniques.bed[,(strand1=='-' & start2>start1) | (strand1=='+' & start1>start2)]] <- 1    # backsplice

    ## Build a data set for histogram analysis
	  junc.bed<-data.frame(junc.bed)
    bed_dataset <- junc.bed[,c(1,2,3)]
    a <- junc.bed[,c(4,5,6)]
    colnames(bed_dataset) -> colnames(a)
    dim(rbind(bed_dataset,a))          # [1] 174560      3
    rbind(bed_dataset,a) -> bed_dataset
    bed_dataset <- cbind(bed_dataset,apply(bed_dataset[,c(1,2,3)],1,paste,collapse="_"))

    tmp<-table(bed_dataset[,4])
    unique.data.frame(bed_dataset) -> bed_dataset
    gsub(pattern=" ",replacement="",x=names(tmp)) -> names(tmp)
    gsub(pattern=" ",replacement="",x=bed_dataset[,4])-> bed_dataset[,4]

    cbind(bed_dataset,tmp[bed_dataset[,4]]) -> bed_dataset
    bed_dataset[,c(1,2,3,5)] -> bed_dataset
    colnames(bed_dataset) <- c("chrom","start","end","score")

 	  cat(paste("\njunc.uniques.bed is",dim(junc.uniques.bed),"\nbed_dataset is",dim(bed_dataset),
				"\nBS_juncs is",dim(BS_juncs),  "\ncanonical_juncs (all input junctions) is",dim(junc)))

    return(list(uniques.bed=junc.uniques.bed, bed=bed_dataset,
                BS_juncs=BS_juncs, # All backsplice junctions after filtering for input file.
                canonical_juncs=junc))
}

##############################################################################################
## circFigure_template1
##
## This function draws genomic graphs from the sushi package. If zoom_coords is NULL then function will display
##     3 genomic graphs (i) transcripts of a gene (ii) junctions graphs for circRNA and (iii) linear RNA.
##    IF zoom_coords contains coordinates then the function will display 4 genomic graphs.
##
##
##
circFigure_template1 <-function(GeneObject, chrom, chromstart, chromend, zoom_coords, JunctionOption, Junction_highlight=list(BSjunc=NULL, Canonical=NULL))
{
  junc <- GeneObject$Junctions
  if (length(junc) == 0)
	{	cat(paste("length of zero (junc) in circFigure_template1"))
		return (NULL)
	}

  a <- GeneObject$Transcript
  display_transcript <- data.frame(chrom=a$chrom, start=as.numeric(a$start), stop=as.numeric(a$stop), gene=a$gene, score=as.numeric(a$score),
                                        strand=as.numeric(a$strand), type=a$type)
  if (display_transcript$strand == 0)
  { display_transcript$strand <- -1 }


	bedJunctions <- junc$uniques.bed
	if (JunctionOption > 0)	# Junction type filter is set therefore extract the junction type to view
	{
	  # Plot Gene Transripts
	  pg = plotGenes(display_transcript,chrom,chromstart,chromend ,
	                 #                types = display_transcript$type,
	                 #                colorby=log10(display_transcript$score+0.001),
	                 #                colorbycol= SushiColors(5),colorbyrange=c(0,1.0),
	                 labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
	  labelgenome( chrom, chromstart,chromend,n=3,scale="Mb")


	  ## Plot linear junctions
	  bedJunctions <- GeneObject$Transcript_Canonical_juncs

	  color_to_graph <- rep(1,nrow(bedJunctions))
	  if (! is.null(Junction_highlight$Canonical$Chr))
	  {
	    BSjunc_idx <- bedJunctions[(chrom1 == as.character(Junction_highlight$Canonical$Chr)) &
	                     (start1 == as.numeric(Junction_highlight$Canonical$Start)) &
	                 (end2 == as.numeric(Junction_highlight$Canonical$End)),which = TRUE]
	    color_to_graph[BSjunc_idx] <- 2      # Colour user selected junction; 2= red; 3 = light green ; 5 = light blue
	  }

	  pbpe = plotBedpe(bedJunctions, chrom, chromstart, chromend, heights = bedJunctions$score, plottype="loops", color = color_to_graph)
	  labelgenome(chrom, chromstart,chromend,n=3,scale="Mb")
	  axis(side=2,las=2,tcl=.2)
	  mtext("SequencingDepth ",side=2,line=3,cex=.75,font=1.5)
	  mtext("Linear junctions",side=3,line=0,cex=.75,font=1.5)

	  ### Plot backsplice junctions
	  bedJunctions <- junc$uniques.bed[JunctionType==JunctionOption,]	}
    ## set line colours
	  color_to_graph <- rep(1,nrow(bedJunctions))  # Default colour = black
	  if (! is.null(Junction_highlight$BSjunc))
	  {
	    BSjunc_idx <- which(bedJunctions$name == Junction_highlight$BSjunc)
	    color_to_graph[BSjunc_idx] <- 5      # Colour user selected junction light blue
	  }

    if (max(bedJunctions$score) < 5)
    {      pbpe = plotBedpe(bedJunctions, chrom, chromstart, chromend, heights = bedJunctions$score, plottype="loops", ymax = 5, color = color_to_graph)  }
    else   # I did try setting ymax to max(bedJunctions$score but did not work for some cases?!?
    {       pbpe = plotBedpe(bedJunctions, chrom, chromstart, chromend, heights = bedJunctions$score, plottype="loops", color = color_to_graph)     }


    labelgenome(chrom, chromstart,chromend,n=3,scale="Mb")
    axis(side=2,las=2,tcl=.2)
    mtext("SequencingDepth ",side=2,line=3,cex=.75,font=1.5)  # increasing value of line moves y-axis label from axis
    mtext("Backsplice junctions",side=3,line=0,cex=.75,font=1.5)


 #    zoomregion1 = c(46241505,46241540)     # To center around 46241520
 #    zoomregion2 = c(46242335,46242370)     # To center around 46242348
 #    zoomregion3 = c(46242405,46242445)     # To center around 46242431

#     zoomsregion(zoomregion1,extend=c(0.01,0.13),wideextend=0.05,offsets=c(0,0.580))
#     zoomsregion(zoomregion2,extend=c(0.01,0.13),wideextend=0.05, offsets=c(0.580,0))
#     zoomsregion(zoomregion3,extend=c(0.01,0.13),wideextend=0.05, offsets=c(0.580,0))

     plotJunctionFrequency<-function(signal, chrom, start, stop)
     {
          signal[, 1] = as.character(signal[, 1])
          signaltrack = signal[which(signal[, 1] == chrom & ((signal[, 2] > start & signal[, 2] < stop) | (signal[, 3] > start & signal[, 3] < stop))), (2:4)]
          plot(signaltrack$start,signaltrack$score,type="h",lwd=10,xlim=c(start,stop))
     }

 #    plotJunctionFrequency(junc$bed,"chr9",46241512,46241528)
 #    plotJunctionFrequency(junc$bed,"chr9",46242340,46242356)
 #    plotJunctionFrequency(junc$bed,"chr9",46242424,46242439)

}

#################################################################################
### Extracts genomic coordinates of a Gene, draws image and returns exon boundaries as dataframe object
##
##
## BS_Junctions            <- raw data containing backsplice junctions
## Canonical_Junctions    <- raw data containing canonical splice junctions
##
##Prepare_Gene_Object<-function(GeneName, BS_Junctions, transcript_reference, File_idx = c(-1), Canonical_Junctions)
Prepare_Gene_Object <- function(GeneName, BS_Junctions, GeneList, File_idx = c(-1), Canonical_Junctions, Genome_Coords = NULL)
{
  if (is.null(GeneName ))
  {  return(NULL)  }
  if ((GeneName == "Novel") || (GeneName =="Unknown"))
  { warning("No annotated gene name")
    return(NULL)
  }

  a <- try(select(GeneList$Annotation_Library, keys = GeneName, columns=c("ENTREZID", "SYMBOL", "ENSEMBL"),keytype="SYMBOL"),silent=TRUE)
  lookupID <- a$ENTREZID   # Default lookup
  if (length(grep(pattern="Error", x=a)))
  {
    showModal(modalDialog(title="Cannot retrieve genomic information for this gene",
            "Please check that you have loaded correct annotation database.
             Nothing will be displayed.",easyClose=TRUE,footer=NULL))
    return ( NULL )
  }


  if (length(which(keys(GeneList$transcript_reference) == a$ENTREZID) ) == 0)
  {
    if (length(which(keys(GeneList$transcript_reference) == a$ENSEMBL) ) == 0)
    {
      ExampleGeneIDs <- head(keys(GeneList$transcript_reference))
      showModal(modalDialog(title="Cannot retrieve or convert genomic information for this gene",
                            "For some reason this gene does not have an entry. If you are seeing this for
                            multiple entries please contact let developer know so we can look into it.
                            Nothing will be displayed.",easyClose=TRUE,footer=NULL))


      return ( NULL )

    }
    else
      lookupID <- a$ENSEMBL

  }

	b <- select(GeneList$transcript_reference, keys = lookupID, columns=c('GENEID', 'TXCHROM', 'EXONSTART',  'EXONEND','TXID', 'EXONSTRAND'),keytype="GENEID")
	b <- b[,c("TXCHROM","EXONSTART","EXONEND","TXID","EXONSTRAND")]
	names(b) <- c('chrom', 'start', 'stop', 'gene','strand')
	b$score <- 0

	if (names(table(b$strand)) == "-") #(b$strand[1]== "-")
	{ b$strand <- 0 }
	else
	{ b$strand <- 1 }
	strand <- b$strand
	b$type ="exon"
#		c <- merge(a, b, 'ENTREZID')
	Transcript <- data.table(b[,c('chrom','start','stop','gene','score','strand','type')])

	# Prepare coordinates that will flank transcript to set boundaries for plots
	chrom <- as.character(Transcript[1,.(chrom)])
	chromstart = as.numeric(min(Transcript[,.(start)]))-50000
	chromend = as.numeric(max(Transcript[,.(stop)]))+500000

#	if (ActivePanel == "Genome_View")
#	{
#	  chrom <- Genome_Coords$chrom
#	  chromstart <- Genome_Coords$chromstart
#	  chromend <- Genome_Coords$chromend
#	  strand <- 0
#	  if (Genome_Coords$chromstrand == "+")
#	    strand <- 1
#	}

	### Filter all circRNA junctions that meet criteria
	Junctions <- circJunctions(BS_Junctions ,chrom, chromstart, chromend,fileID=File_idx)

	### Filter canonical junctions within transcript
	if (strand == 0)  # negative strand as assigned by reference transcript
	   strand = 2
	else             # Positive strand as assigned by reference transcript
	   strand = 1

	Transcript_Canonical_juncs = Canonical_Junctions[chrom1 == chrom & start1 > chromstart &
	                                                   end2 < chromend &   strand1 == strand & DataSet %in% File_idx,]

	return (list(Transcript=Transcript, Junctions=Junctions, Transcript_Canonical_juncs= Transcript_Canonical_juncs))
}


####################################################################################
## Draw_Transcript_Exons
##
##   Takes in a gene name and draws transcript. Currently has issues drawing arrows
## for negative strand. Draws positive strand arrows correctly. Also there is potential
## issue with the case of the gene name (ie not checked)
##
## GeneObject is a list with at least 3 elements:
##    Transcript                  (parental transcript coords)
##    Junctions                   (Junc.bed, )
##    Transcript_Canonical_juncs (data frame containing canonical junctions)
##  Zoom_coords is either NULL or user defined coordinates
##
Draw_Transcript_Exons<-function(GeneObject, JunctionOption, Zoom_coords, GenomeDisplay=FALSE,
                                Genome_Coords=list(chrom=NULL, chromstart=NULL, chromend=NULL, chromstrand=NULL),
                                Junction_highlight = list(BSjunc=NULL, Canonical=NULL))
{
	if (length(GeneObject) == 0)
	{ return (NULL)	}

  chrom <- {}
  chromstart <- {}
  chromend <- {}

  if (GenomeDisplay) # User requesting to view a specific region of genome
  {
    if ( (is.null(Genome_Coords$chrom)) || (is.null(Genome_Coords$chromstart)) || (is.null(Genome_Coords$chromend)) )
    {  return(NULL)  }
    chrom <- Genome_Coords$chrom
    chromstart <- as.numeric(Genome_Coords$chromstart)
    chromend <- as.numeric(Genome_Coords$chromend)
  }


	if (! GenomeDisplay)# Prepare genome coordinates that will flank transcript to set boundaries for plots
	{ chrom <- as.character(GeneObject$Transcript[1,.(chrom)])
	  chromstart = as.numeric(min(GeneObject$Transcript[1,.(start)]))-15
	  chromend = as.numeric(max(GeneObject$Transcript[nrow(GeneObject$Transcript),.(stop)]))+15
	}
	if (! is.null(Zoom_coords))
	{ #if ((chromstart-550 < Zoom_coords[1]) && (chromend+550 > Zoom_coords[2]))
	  {   chromstart  <- Zoom_coords[1]
	      chromend    <- Zoom_coords[2]
	   }
	}

	if (length(GeneObject$Junctions) == 4)	# Have some circRNAs to display
	{   circFigure_template1(GeneObject = GeneObject,chrom = chrom,chromstart=chromstart, chromend=chromend, JunctionOption=JunctionOption, Junction_highlight=Junction_highlight )
	}

	if (length(GeneObject$Junctions) == 1)	# No circRNAs to display. Just show transcript
	{ if (length(GeneObject$Transcript_Canonical_juncs) == 13 ) # We have canonical junctions to view but no Backsplice.
  	{
	    # Make a dummy backsplice junction object
	    GeneObject$Junctions <- list()
	    GeneObject$Junctions$uniques.bed <- GeneObject$Transcript_Canonical_juncs[1,]
	    GeneObject$Junctions$uniques.bed$score <- 0
	    GeneObject$Junctions$uniques.bed$JunctionType <- 1
      circFigure_template1(GeneObject = GeneObject,chrom = chrom,chromstart=chromstart, chromend=chromend, JunctionOption=JunctionOption )

	  }
	  if (length(GeneObject$Transcript_Canonical_juncs) != 13 )
	  { plotGenes(GeneObject$Transcript, chrom, chromstart, chromend, maxrows=50,height=0.4,plotgenetype="box")	}

	  } # if (length(GeneObject$Junctions) == 1)
}

########################################################################################################
## Gene_Transcript_Features
##
##  This function returns a numeric array containing number of exons for every transcript. Transcript
##     names are assigned to every element of the array.
##
##  Future ideas:
##      - Identify coding vs non-coding transcripts?
##      - Intron analysis ??
##      - Average junction read count
##      - If multiple transcripts are expressed: estimate transcript expression levels
##      - Novel exon usage
##
##   INPUTS:
##  Gene_Symbol is the gene symbol ID  eg "SLC8A1"
##
##
## Written by David Humphreys 13th Jan 2016
##
##  Note I think I have doubled up on GeneList and the gene transcript present in GeneObject. LEaving for now.. it works
Gene_Transcript_Features <- function (GeneList, Gene_Symbol, GeneObject)
{
  a <- select(GeneList$Annotation_Library, keys = Gene_Symbol, columns=c("ENTREZID", "SYMBOL"),keytype="SYMBOL")
  b <- select(GeneList$transcript_reference, keys = a$ENTREZID, columns=c('GENEID', 'TXNAME', 'EXONRANK'),keytype="GENEID")
  Num_Transcripts <- length(unique(b$TXNAME))
  Num_Exons_Per_Transcript <- {}
  for(i in 1:length(unique(b$TXNAME)))
  {     Num_Exons_Per_Transcript <- c(Num_Exons_Per_Transcript , length(which(b$TXNAME == unique(b$TXNAME)[i])))          }
  Num_Exons_Per_Transcript <- as.numeric(Num_Exons_Per_Transcript)
  names(Num_Exons_Per_Transcript) <- unique(b$TXNAME)


##  start <- min(GeneObject$Transcript$start)
##  end <- max(GeneObject$Transcript$stop)
##  chrom <- GeneObject$Transcript$chrom[1]
##  strand <- GeneObject$Transcript$strand[1]
  LinearJunctions <- GeneObject$Transcript_Canonical_juncs    # Note that this is already subsetted prior to calling this function.
  if (length(LinearJunctions) > 0)
    t_idx <- order(x=LinearJunctions$score, decreasing=TRUE)[seq(from=1, to = max(Num_Exons_Per_Transcript))]
  else
    return (NULL)

##  debugme("Gene_Transcript_Features at line 335", LinearJunctions, Num_Exons_Per_Transcript)
  #a <- order(x = m47_Rbm47Juncs$UniqueMapped, decreasing = TRUE )[seq(from=1, to=Answer$exon_number -1)]
  Av_junc_count <-  ave(LinearJunctions$score[t_idx])[1]

  return(list(Num_Exons_Per_Transcript= Num_Exons_Per_Transcript, Av_junc_count=Av_junc_count))
}

########################################################################################################
## Abundant_circRNAs
##
##    This function will determine the abundance of circRNA junctions compared to linear junctions
##
Abundant_circRNAs<- function()
{

}

#########################################################################################
## BS_Distribution
##
## Currently have just pasted code below. Need to work through. Ideally would like to pass a backsplice
##   junction unique key
##
##  This function uses calles from gsubfn library (strapplyc)
Fragment_Alignment_Distribution<-function(BSJ_table)
{
	if (nrow(BSJ_table) == 0)
	{	return (NULL)	}		# Dataset has not been loaded.

#	setnames(junction,1:15,c("chromDonor","startDonor","strandDonor",
#								"chromAcceptor","startAcceptor","strandAcceptor",
#								"JuncType", "RepeatLength_L", "RepeatLength_R",
#								"ReadName","FirstBase_1stSeq","CIGAR_1stSeg",
#								"FirstBase_2ndSeq","CIGAR_2ndSeg", "BSjuncName"))

	## Below is what tmp looks like (paired end data):
#		  row.names BS_Jnc_A     BS_Jnc_B       FragA        FragB      CIGAR_FragA    CIGAR_FragB
#	1     8818      76881840     76857903     76881769     76857904     71M80S        71S63M17S
#	2     17973     76881840     76857903     76881723     76857904     1S117M33S     118S33M
#	3     34504     76881840     76857903     76881756     76857904     84M67S        84S67M
#	4     41226     76881840     76857903     76881759     76857904     2S81M68S      83S68M
#	5     45541     76881840     76857903     76881769     76857904     71M80S        71S80M
#	6     51752     76881840     76857903     76881789     76857904     51M100S       51S87M13S

# The match in column "CIGAR_FragA" is the length of sequence 5' of BS_Jnc_A (ie add Match to column "FragA" = column "BS_Jnc_A"
# The match in column "CIGAR_FragB" is the length of sequence 3' of BS_Jnc_B
# Therefore length of backsplice junction is (BS_Junc_A - match length in CIGAR_FragA) to (BS_Junc_B + match length in CIGAR_FragB)
# Will provide average read length, and range of read lengths across backsplice junction
# Will be counting nucleotide frequency on each side of backsplice junction for EVERY read <== this will be assembled into a data frame (BED format with two count columns)
#

# The following examples from Apoa4 cause the following lines of code to fail. Think need to add ALL M values for each CIGAR.
# Example 1 20M684N80M-5p77M25S		76S23M2S				HWI-D00119:42:H0U5HADXX:1:1210:7553:28529
# Example 2 100M-6p15M87S			14S73M684N13M1S			HWI-D00119:42:H0U5HADXX:1:1116:6120:18060
# Example 3 76S23M2S				98M-77p77M25S			HWI-D00119:42:H0U5HADXX:1:1210:16462:37981
#
# 		- 684N  specifies the intron
#		- Examples listed above wrap around. The -5p indicates a wrap.
#		- NEED to add all M values from each CIGAR.
#			. Need to split strings that have a p value. Ignore right hand side of string. Left hand side of string:
#			. Need to incorporate intron (ie N) scores.

	# Extract First segment matches from CIGAR strings
	FirstSeg_Match <- lapply(strapplyc(as.character(BSJ_table[,c(CIGAR_1stSeg)]),pattern="([-0-9]+)M"), FUN = function(x) {sum(as.numeric(x))})

	# Adjust for negative padding values eg 150M-63p114M36S
	FirstSeg_Negative_Padding <- lapply(strapplyc(as.character(BSJ_table[,c(CIGAR_1stSeg)]),pattern="([-0-9]+)p"), FUN = function(x) {sum(as.numeric(x))})
  FirstSeg_NoMatch <- lapply(strapplyc(as.character(BSJ_table[,c(CIGAR_1stSeg)]),pattern="([-0-9]+)N"),FUN=function(x) { x[duplicated(x)]})
	# Calculate the sum of matched positions. This is the end position of fragment.
	FirstSeg_position <- mapply(sum,FirstSeg_Match, FirstSeg_Negative_Padding)   #, as.numeric(tmp[,c(FirstBase_1stSeq)]))


  # Repeat for 	other segment
	SecondSeg_Match <- lapply(strapplyc(as.character(BSJ_table[,c(CIGAR_2ndSeg)]),pattern="([-0-9]+)M"), FUN = function(x) {sum(as.numeric(x))})
	SecondSeg_Negative_Padding <- lapply(strapplyc(as.character(BSJ_table[,c(CIGAR_2ndSeg)]),pattern="([-0-9]+)p"), FUN = function(x) {sum(as.numeric(x))})
	SecondSeg_NoMatch <- lapply(strapplyc(as.character(BSJ_table[,c(CIGAR_2ndSeg)]),pattern="([-0-9]+)N"),FUN=function(x) { x[duplicated(x)]})
	SecondSeg_position <- mapply(sum,SecondSeg_Match, SecondSeg_Negative_Padding)   #, as.numeric(tmp[,c(FirstBase_1stSeq)]))

	# Adjust for any negative values that arise due to short fragments
	#     eg 25S98M6914N27M-7010p69M6914N56M1906N24M1S
	# Resolve by adding the duplicated intronic regions eg 6914
	for (i in 1: length(FirstSeg_position))
	{
	  if (FirstSeg_position[i] < 0)
	  {   if (length(FirstSeg_NoMatch[[i]]) > 0)
	      { FirstSeg_position[i] <- FirstSeg_position[i] + as.numeric(FirstSeg_NoMatch[[i]]) }
	  }
	  if (SecondSeg_position[i] < 0)
	  {   if (length(SecondSeg_NoMatch[[i]]) > 0)
	       {  SecondSeg_position[i] <- SecondSeg_position[i] + as.numeric(SecondSeg_NoMatch[[i]]) }
	  }
	}

  # Calculate the range of where alignments fall
	if ((length(FirstSeg_position[FirstSeg_position > 0]) == 0) || (length(SecondSeg_position[SecondSeg_position > 0]) ==0))
	{
#	  browser()
	  a <- "Need to identify why lengths are incorrect"
	}
	FirstSeg_position <- FirstSeg_position[FirstSeg_position > 0]   # Remove -ve values as have incorrectly calculated these
	FirstSeg_Range <- range(FirstSeg_position)
	if (FirstSeg_Range[2] > 250)  # Assuming most library fragment lengths are approximately 250nt. Will get larger fragments, resize back to 250
	{ FirstSeg_Range[2] <- 250 }

	SecondSeg_position <- SecondSeg_position[SecondSeg_position > 0] # Remove -ve values as have incorrectly calculated these
	SecondSeg_Range <- range(SecondSeg_position)
	if (SecondSeg_Range[2] > 250)  # Same as above
	{ SecondSeg_Range[2] <- 250 }

	# Calculate duplication frequency ratio
	FirstSeg_DFR <- 5
	FirstSeg_DFR <- 5

	FirstSeg_Dups  <- table(duplicated(FirstSeg_position))
	SecondSeg_Dups <- table(duplicated(SecondSeg_position))

	if (! is.na(FirstSeg_Dups["FALSE"]))
	{  FirstSeg_DFR <- length(FirstSeg_position) / FirstSeg_Dups["FALSE"] }

	if (! is.na(SecondSeg_Dups["FALSE"]))
	{  SecondSeg_DFR <- length(SecondSeg_position) / SecondSeg_Dups["FALSE"] }

	######## Calculate RAD score #########
  FirstSeg_RAD <- 	1 - ((FirstSeg_Range[2] - FirstSeg_Range[1]) / 250)   # * FirstSeg_DFR
  SecondSeg_RAD <- 	1 - ((SecondSeg_Range[2] - SecondSeg_Range[1]) / 250) # * SecondSeg_DFR

  FirstSeg_RAD <- skewness(FirstSeg_position)  # kurtosis(FirstSeg_position)
  SecondSeg_RAD <- skewness(SecondSeg_position) # kurtosis(SecondSeg_position)

  return (list(FirstSeg_RAD = FirstSeg_RAD, SecondSeg_RAD = SecondSeg_RAD, FirstSeg_position = FirstSeg_position, SecondSeg_position = SecondSeg_position))

	# Now make a BED file data frame. For TTN example column "BS_Jnc_B" defines backsplice junction point.
	# Will extend 5' 3' by largest range
	BSJ <- as.numeric(gsub(pattern=".*_",replacement = "",x =junctionID))

	# Create a nucleotide frequency histogram
	BS_bed_dataset<- rep.int(0,times= range(Col4Match)[2] * 2 +1)
	names(BS_bed_dataset) <- seq(from=(BSJ-range(Col4Match)[2]), to=(BSJ+range(Col4Match)[2]))
	for(i in 1:length(Col2Match))
	{
	    idx <- names(BS_bed_dataset)>(BSJ - Col2Match[i] )  & names(BS_bed_dataset) < BSJ
	    BS_bed_dataset[idx] <- BS_bed_dataset[idx] + 1
	    idx <- names(BS_bed_dataset)>(BSJ  )  & names(BS_bed_dataset) < (BSJ + Col4Match[i])
	    BS_bed_dataset[idx] <- BS_bed_dataset[idx] + 1
	}
	plot(BS_bed_dataset)

	# Following code will create a read based representation of data, data stored in bedpe format:
	# chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 samplenumber

#	BS_bedpe_dataset <- data.frame()

#	for(i in 1:length(Col2Match))
#	{
#		BS_bedpe_dataset <- rbind(BS_bedpe_dataset, c(BSJ - Col2Match[i], BSJ, BSJ, BSJ + Col4Match[i] ))       # Start stop  Start stop     <== linearised backsplice junction
#	}


}

########################################################################################
## LoadJunctionData
##
## multi_option: This currently will collate and merge different data sets (i.e. union).
##					Would like to implement:
##						- Difference (junctions that are distinct between data sets).
##						- Intersection (junctions that are common between data sets).
## 						- Confident junctions (only accept junctions that are known to exist)
## 						- Keep a link that identifies what file
##
##  This function will be required to read in different file types and return the appropriate data lists.
##
LoadJunctionData <- function(filename, ChromFilter, StrandFilter, Genomic_Distance, Junction_abundance, RAD_filter, MultipleDataOption,CanonicalJuncs, SubSelect)
{	ncolumns <- 1

	# AddColumns function populates a column with a unique identifier for backsplice junctions plus a column for input file ID
	AddColumns <- function(junc, File_Index)
	{	setnames(junc,1:14,c("chromDonor","startDonor","strandDonor",
								"chromAcceptor","startAcceptor","strandAcceptor",
								"JuncType", "RepeatLength_L", "RepeatLength_R",
								"ReadName","FirstBase_1stSeq","CIGAR_1stSeg",
								"FirstBase_2ndSeq","CIGAR_2ndSeg"))
		junc$BSjuncName <- paste(junc$chromDonor,junc$startDonor,junc$chromAcceptor, junc$startAcceptor,sep="_")
		junc$DataSet <-  File_Index	# This adds an index to identify the file
		return(junc)
	}
	# name; size; datapath; type
	AllData <- {}
	SummarisedData <- {}
	Canonical_AllData <- {}
	ReadsPerGene_Data <- {}
	DataType <- c("Unknown")

	# Identify which files are chimeric and which files are linear.
	# Filenames should conform to following format :
	#                 <SampleID>.Chimeric.out.junction
	#                 <SampleID>.SJ.out.tab
	#
#	IDs <- unique(gsub(pattern="\\..*?$",replacement="",x=filename$name))  # Extract ID to associate with samples
	IDs <- unique(gsub(pattern="\\.Chimeric.out.junction$",replacement="",x=filename$name))
	IDs <- unique(gsub(pattern="\\.SJ.out.tab$",replacement="",x=IDs))

withProgress(message="Importing data. This could take a few minutes", value=0, {

	for (i in 1: length(IDs))
	{ incProgress(1/(length(IDs)*2), detail = paste("Loading BSJ for file number ", i))
	  Chimeric_Idx <- grep(pattern=paste(IDs[i],".Chimeric.out.junction",sep=""), x=filename$name )
	  if (length(Chimeric_Idx) > 0)
	  {
	    data_set <- fread(filename$datapath[Chimeric_Idx[1]], sep="\t")
	    total_rows <- nrow(data_set)
	    n_UniqueJunctions <- length(table(data_set$BSjuncName))
	    if (ncol(data_set) != 14)
	    {   ### Need to set an alert if ncol was not what was expected !!!

	    }
	    if (ncol(data_set) == 14)           # Chimeric Junction data file from STAR aligner
	    {	setnames(data_set,1:14,c("chromDonor","startDonor","strandDonor", "chromAcceptor","startAcceptor","strandAcceptor",
	                           "JuncType", "RepeatLength_L", "RepeatLength_R",  "ReadName","FirstBase_1stSeq","CIGAR_1stSeg",
	                           "FirstBase_2ndSeq","CIGAR_2ndSeg"))
	      data_set$BSjuncName <- paste(data_set$chromDonor,data_set$startDonor,data_set$chromAcceptor, data_set$startAcceptor,sep="_")
	      data_set$DataSet <-  i	# This adds an index to identify the file
	      DataLists <- FilterChimeric(All_junctions = data_set, ChromFilter=ChromFilter, StrandFilter=StrandFilter, Genomic_Distance=Genomic_Distance, Junction_abundance=Junction_abundance, RAD_filter=RAD_filter, CanonicalJuncs=CanonicalJuncs)
######################################### Can possibly delete this line
  	    DataLists$SummaryData <- SelectUniqueJunctions(DataLists$RawData)   # May not need this line
#########################################

  		  DataType <- c("BackSplice")
	      n_UniqueJunctions_PostFilt <- length(table(data_set$BSjuncName))

	      # Merge Dataset to master table
	      if (! is.null(AllData))
	      {  # We should have 15 columns of data at this stage. Merge everything.
  	      if ((ncol(DataLists$RawData) > 1) && (ncol(AllData) == ncol(DataLists$RawData)) )
	        {
	          AllData<-rbind(AllData,DataLists$RawData)
	          SummarisedData <- rbind(SummarisedData, DataLists$SummaryData)
  	      }
	      }
	      if (is.null(AllData))
	      {
	        SummarisedData <- DataLists$SummaryData
	        AllData <- DataLists$RawData
	      }
	    } # if (ncol(data_set) == 14)
	  } # if (Chimeric_Idx > 0)
	  else   # There is no chimeric data imported. Make a fake entry
	  {  data_set <-data.table(chromDonor="chr1",startDonor=1,strandDonor="+", chromAcceptor="chr1",startAcceptor=10,strandAcceptor="+",
	                           JuncType=0, RepeatLength_L=0, RepeatLength_R=0,  ReadName="fake",FirstBase_1stSeq=0,CIGAR_1stSeg="0M",
	                           FirstBase_2ndSeq=0,CIGAR_2ndSeg="0M")
	    data_set$BSjuncName <- paste(data_set$chromDonor,data_set$startDonor,data_set$chromAcceptor, data_set$startAcceptor,sep="_")
	    data_set$DataSet <-  1
	    DataLists <- FilterChimeric(All_junctions = data_set, ChromFilter=ChromFilter, StrandFilter=StrandFilter, Genomic_Distance=Genomic_Distance, Junction_abundance=Junction_abundance, RAD_filter=RAD_filter, CanonicalJuncs=CanonicalJuncs)
	    SummarisedData <- DataLists$SummaryData
	    AllData <- DataLists$RawData
	    total_rows=1
	    n_UniqueJunctions=1
	    n_UniqueJunctions_PostFilt=0
	  }

	  incProgress(1/(length(IDs)*2), detail = paste("Loading FSJ for file number ", i))
	  Linear_Idx <- grep(pattern=paste(IDs[i],".SJ.out.tab",sep=""), x=filename$name )
	  if (length(Linear_Idx) > 0)
	  {
	    data_set <- fread(filename$datapath[Linear_Idx[1]], sep="\t")
	    if (ncol(data_set) == 9)            # Linear Junction data file from STAR aligner
	    {
	      DataType <- c("Canonical")
	      ## Canonical junctions originally look like this:
	      # c("chrom1","start1","end1","strand","intron_motif", "annotated", "unique_counts", "multi_mapped_counts", "overhang")
	      ## Need to convert it to the following
	      # chrom1   start1     end1 chrom2   start2     end2   name score  strand1 strand2
	      data_set <- data_set[,.(V1,V2,V2, V1,V3,V3,   V1,  V7, V4, V4,V8 ,V9)]
	      setnames(data_set,1:12,c("chrom1","start1","end1","chrom2","start2","end2",   "name","score", "strand1", "strand2", "multimappers", "overhang"))
	    #  Canonical_SummarisedData <- data_set
	      data_set$DataSet <-  i	   # This adds an index to identify the file

	      # Merge data
	      if (! is.null(Canonical_AllData))
	      {  Canonical_AllData <- rbind(Canonical_AllData, data_set)  }
	      else
	      {  Canonical_AllData <- data_set         }

	    }
	  } # 	  if (length(Linear_Idx) > 0)
	  ReadsPerGene_Idx <- grep(pattern=paste(IDs[i],".ReadsPerGene.out.tab",sep=""), x=filename$name )
	  if (length(ReadsPerGene_Idx) > 0)
	  { data_set <- {}
	    data_set[[i]] <- fread(filename$datapath[ReadsPerGene_Idx[1]], sep="\t")
	    setnames(data_set[[i]],1:5, c("geneName","unstranded","Fstrand","Rstrand"))
	    data_set$DataSet <- i
	    if (! is.null(ReadsPerGene_Data))
	    {    ReadsPerGene_Data <- rbind(ReadsPerGene_Data, data_set) }
	    else
	    {    ReadsPerGene_Data <- data_set    }
	  }
	} # 	for (i in 1: length(IDs))

}) # withProgress


	return(list(Junctions=AllData, SummarisedData=SummarisedData,
	            Original_Junction_Numbers=total_rows, # the number of rows in the last file read. This is pointless if multiple files are read in.
	            Original_n_unique_junctions=n_UniqueJunctions,  # the number of columns in the last file read. This is pointless if multiple files are read in.
	            Original_Postfilt_n_unique_junctions=n_UniqueJunctions_PostFilt,   # This is identical to n_UniqueJunctions !! probably can delete.
	            DataType = DataType,                             # "Backsplice"  or "Canonical" so far  <= This is now redundant !!!
	            SampleIDs = IDs,
	            Canonical_AllData = Canonical_AllData,
	            ReadsPerGene_Data = ReadsPerGene_Data
	))

}


#########################################################################################################################
## FilterChimeric
##
## * Filter_options
##		- chromosomes: Equal, Any
##		- Genomic distance :assuming chromosomes are equal then startDonor and startAcceptor to be within a defined range
##		- strands are the same
##		- abundance : minimum frequency of reads
##		- RAD_filter
##
## Have not implemented this function yet. Need to implement
##
FilterChimeric <- function(All_junctions, ChromFilter, StrandFilter, Genomic_Distance, Junction_abundance, RAD_filter, CanonicalJuncs=TRUE, fileID= c(-1), ChrM_Filter=FALSE, InvertReads = FALSE, SummaryNumber = 50)
{

withProgress(message="Calculating RAD scores", value=0, {
  incProgress(1/4, detail = paste("Selecting junctions from selected data sets "))   # 1 of 4


	if (fileID[1] != -1)  # Select data from files that user has selected
	{	All_junctions <- Filter_by_Data_Set(fileID, All_junctions)	}

  incProgress(1/4, detail = paste("Applying chromosome filters"))           # 2 of 4
	if (ChromFilter == TRUE)
	{	All_junctions = All_junctions[chromDonor == chromAcceptor,]		}

  if (ChrM_Filter == TRUE)
  { All_junctions = All_junctions[chromDonor != "chrM",]	}

  incProgress(1/4, detail = paste("Applying strand filter"))	              # 3 of 4
	if (StrandFilter == TRUE)
	{	All_junctions = All_junctions[strandDonor == strandAcceptor,]	}

  if (InvertReads == TRUE)
  {
    pos_strand_idx <- All_junctions$strandDonor == "+"
    All_junctions$strandDonor[pos_strand_idx] = "-"
    All_junctions$strandDonor[! pos_strand_idx] = "+"
    tmp <- All_junctions$startDonor[pos_strand_idx]
    All_junctions$startDonor[pos_strand_idx]  <- All_junctions$startAcceptor[pos_strand_idx]
    All_junctions$startAcceptor[pos_strand_idx] <- tmp
    }

  incProgress(1/4, detail = paste("Applying genomic distance filter"))        # 4 of 4
	if ((ChromFilter == TRUE) && (StrandFilter == TRUE)) # Apply genomic distance filter
	{	# Classical BSjunctions are defined as:
		# (strandDonor=='-' & (startAcceptor > startDonor))  | (strandDonor=='+' & (startDonor > startAcceptor) )

		BS_junctions = All_junctions[((strandDonor == "+" ) &
				((startDonor - startAcceptor) > Genomic_Distance[1]) & ((startDonor - startAcceptor) < Genomic_Distance[2])) |
						((strandDonor == "-" ) &
				((startAcceptor - startDonor) > Genomic_Distance[1]) & ((startAcceptor - startDonor) < Genomic_Distance[2])),]
		BS_junctions$type = "bs"

		if (CanonicalJuncs == TRUE)
		{	# Merge both canonical and BS junctions together
			Canonical_junctions =  All_junctions[((strandDonor == "-" ) &
				((startDonor - startAcceptor) > Genomic_Distance[1]) & ((startDonor - startAcceptor) < Genomic_Distance[2])) |
						((strandDonor == "+" ) &
				((startAcceptor - startDonor) > Genomic_Distance[1]) & ((startAcceptor - startDonor) < Genomic_Distance[2])),]
			Canonical_junctions$type = "c"
			All_junctions <- rbindlist(list(BS_junctions, Canonical_junctions))
		}
		if (CanonicalJuncs != TRUE)
			All_junctions <- BS_junctions
	}
  else  # Will only come in here if searching for F-circRNA or fusion reads?
  { BS_junctions <- All_junctions
    BS_junctions <- "Any"

  }
}) # withProgress(message="Calculating RAD scores", value=0, {

  filterlist =list(BSjuncName=NULL, SortDir="Descending", IndexNumber=1, DisplayNumber=SummaryNumber, DisplayRAD_score= FALSE)
  return(list(RawData=All_junctions,  SummaryData= SelectUniqueJunctions(All_junctions, filterlist = filterlist)	))
}

################################################################################################
## SelectUniqueJunctions
##
## This is going to be the workhorse for displaying collated junctions from the input data. It will return
##    selected rows of data (annotated) that will enable enhanced browsing of raw data on the fly.
##
##		Filter options: Junction abundance. Sort
##    library_strand = "Opposing strand" or "Same strand" or "Unstranded"
SelectUniqueJunctions <- function(All_junctions, filterlist = list(BSjuncName=NULL, SortDir="Descending", IndexNumber=1, DisplayNumber=10, DisplayRAD_score= FALSE), library_strand = "Opposing strand")
{
  filterlist$DisplayNumber <- as.numeric(filterlist$DisplayNumber)
  if (! is.null(filterlist$BSjuncName))  # Just return entries that relate this this filter.
  {
    toReturn <- All_junctions[All_junctions$BSjuncName == filterlist$BSjuncName]
    return(toReturn)
  }

	junc_count <- table(All_junctions$BSjuncName)				# Will use this to quantify the number of the backsplice junction

	junc_number <- length(junc_count)
	if (junc_number < filterlist$DisplayNumber)
	{	filterlist$DisplayNumber <- junc_number	}

	# junc_count <- junc_count[junc_count>Junction_abundance]	# This selects only junctions with a minimum number of counts

	# Making a table with following columns:
	#		"BSjuncName" (eg chr9_46241521_chr9_46242431), "type" (canonical/BS), "JuncType" (0 1 2), Count, PCR duplicate score.
	# 		NOT IMPLEMENTED: annotation,

	First_CIGAR <- {}
	Second_CIGAR <- {}


	# Sort table
	if (filterlist$SortDir == "Descending")
	{	junc_count <- junc_count[names(sort(junc_count, decreasing = TRUE))]	}
	else
	{	junc_count <- junc_count[names(sort(junc_count, decreasing = FALSE))]	}

	# The following loop calculates PCR duplicate frequency of backsplice junction.
	First_CIGAR <- {}
	Second_CIGAR <- {}
	UniqueJunctions <- {}

withProgress(message="Calculating distributions", value=0, {
  Max_Display <- filterlist$IndexNumber + filterlist$DisplayNumber

	for(i in filterlist$IndexNumber : Max_Display)
	{
	  incProgress(1/Max_Display, detail = paste("Updating Row ",i))
		# Check to make sure the loop increment has not jumped past the max numbers of entries
		if (i > length(junc_count))
			break
		BS_Junc_ID <- names(junc_count)[i]

		OneJunctionReads <- (All_junctions[All_junctions$BSjuncName == BS_Junc_ID,])

		if (length(grep(pattern = 'type', x = colnames(OneJunctionReads))))
		{    temp <- OneJunctionReads[1,.(BSjuncName, type, JuncType, strandDonor)]	  }
		else  # On very rare occations when applying NO filter for BSJ will enter this block
		{
		    temp <- OneJunctionReads[1,.(BSjuncName, JuncType, strandDonor)]
		    temp$type <- "unknown"

		}

		additional_counts <- 0
		if (library_strand == "Unstranded")  # Need to append data from palindrome ID
		{  # Only keep IDs that have a larger coordinate in first position of ID
		  decode_juncID <- strsplit(BS_Junc_ID  , "_")
		  if (as.numeric(decode_juncID[[1]][4]) > as.numeric(decode_juncID[[1]][2]))
		    next

		  # Extract matching coordinate on other strand and add value
		  matching_coord_ID <- paste(decode_juncID[[1]][1],decode_juncID[[1]][4],decode_juncID[[1]][1],decode_juncID[[1]][2],sep = "_")

      matching_coord_entries <- (All_junctions[All_junctions$BSjuncName == matching_coord_ID,])
      if (dim(matching_coord_entries)[1])
      {  		  OneJunctionReads <- rbind(All_junctions[All_junctions$BSjuncName == matching_coord_ID,], matching_coord_entries)
		          #temp <- rbind(temp, OneJunctionReads[1,.(BSjuncName, type, JuncType, strandDonor)])
      }

		  if (! is.na(junc_count[matching_coord_ID]))  # Great! Found the junctions from opposing reads
		  { additional_counts <- junc_count[matching_coord_ID] }
      else # Sometimes opposing reads are offset by a couple of nucleotides. Lets search for it:
      {
        for(j in -5:5)
        {
          matching_coord_ID <- paste(decode_juncID[[1]][1],as.numeric(decode_juncID[[1]][4])+j,decode_juncID[[1]][1],as.numeric(decode_juncID[[1]][2])+j,sep = "_")
          if (! is.na(junc_count[matching_coord_ID]))  # Great! Found the junctions from opposing reads
          { additional_counts <- junc_count[matching_coord_ID]
            break
          }
        }
      }

		}
		temp$Freq <- junc_count[i] + additional_counts

		if (filterlist$DisplayRAD_score)
		{  FAD <- Fragment_Alignment_Distribution(BSJ_table = OneJunctionReads)
			 temp$PCR_dup_A <- round(FAD$FirstSeg_RAD,2)
		   temp$PCR_dup_B <- round(FAD$SecondSeg_RAD, 2)
		}
##### TypeII_III should move to above if statement?
		temp$TypeII_TypeIII <- round(length(grep(pattern = "M.*M",x = OneJunctionReads$CIGAR_1stSeg))/length(OneJunctionReads$CIGAR_1stSeg),2)

		if (length(UniqueJunctions)[1] == 0)
		{	UniqueJunctions <- temp	}
		else
		{	UniqueJunctions <- try(rbind(UniqueJunctions, temp)  )		}
	}
}) # withProgress(message="Calculating RAD scores", value=0, {

	return(UniqueJunctions)
}

#######################################################################################################
## BS_Junc_details
##
## JuncCoords : Junction coordinates in following format:
##				chr1_33667149_chr1_33668857
##
## This function will provide following information regarding the splice junction.
##		- Backsplice junction genomic distance (less than 200nt is classed as a small circRNA)
##
## Would like to display the following information:
			#	- Number of supporting reads (PCR duplicates, Type I II, III and IV reads)
			#	- Junction type, i.e. backsplice/canonical/other
			# 	- If junction is flanked by splicing acceptor/donor?
			#	- Is there evidence of splicing between junction?
			#	- Is there evidence of alternative junction (perhaps providing links to load these)?
			#	- Has this junction been reported before (ie high confidence junction)?
			#   - List nearby repeat regions?
			# 	- What proportion of reads came from each input data file.
BS_Junc_details <- function(JuncCoords)
{
	BS_details <- strsplit(JuncCoords,split = "_")
	output_BS_details <- c("Error, Junction incorrect format. Please report", JuncCoords)
	if (length(BS_details[[1]]) == 4)
	{
		# First determine if a backsplice junction or canonical
		JuncType <- c("Canonical junction")
		JuncA <-as.numeric(BS_details[[1]][4])
		JuncB <-as.numeric(BS_details[[1]][2])
		if (JuncA > JuncB)
		{	JuncType <- c("Backsplice junction")  }

		# Now determine genomic distance that makes the junction
		BS_Junction_Width <- JuncA - JuncB
		output_BS_details <- c('')
		if (BS_Junction_Width < 200)
		{	output_BS_details<- paste("Small ")	}
		output_BS_details<- paste(JuncType, " . ", output_BS_details, "circle of width ",BS_Junction_Width)
	}
	return(output_BS_details)
}

###################################################################################################
##
##  Annotate_BS_Junc
##
##
## OnlyAnnotated:   Ony return entries that are annotated with gene name. Can be a reasonable filter
##
##
Annotate_BS_Junc<- function(DataSet, GeneList, MaxDisplay = 15, Library_Strand = "", OnlyAnnotated="FALSE")
{
	if (nrow(DataSet) < MaxDisplay)
		MaxDisplay <- nrow(DataSet)

#	t_GR <- transcripts(GeneList$transcript_reference)
	g_GR <- genes(GeneList$transcript_reference)
	# select(org.Hs.eg.db, keys = '1', columns=c("ENTREZID", "SYMBOL","OMIM"),keytype="ENTREZID")
	#   ENTREZID SYMBOL   OMIM
	# 1        1   A1BG 138670

	withProgress(message="Annotating table", value=0, {

	  for(i in 1:MaxDisplay)
	  {
		  BSjuncDetails <- strsplit(DataSet$BSjuncName[i], split = "_")
		  strand <- DataSet$strandDonor[i]
		  if (Library_Strand == "Same Strand")
		  {  if (length(strand) == 0)
		    {
		       # browser()
		        cat("Catch error")
		    }
		      if (strand == "-")
		        strand <- "+"
		     else
		        strand <- "-"
		  }
		  else if (Library_Strand == "Unstranded")
		  { strand <- '*'  }

		  bs_junc_gr <- GRanges(seqnames=BSjuncDetails[[1]][1], ranges = as.numeric(min(BSjuncDetails[[1]][c(2,4)]),min(BSjuncDetails[[1]][c(2,4)])),strand = strand)
		  t_start <- findOverlaps(invertStrand(bs_junc_gr),g_GR, type=c("within"))
		  bs_junc_gr <- GRanges(seqnames=BSjuncDetails[[1]][1], ranges = as.numeric(max(BSjuncDetails[[1]][c(2,4)]),max(BSjuncDetails[[1]][c(2,4)])),strand = strand)
		  t_end <- findOverlaps(invertStrand(bs_junc_gr),g_GR, type=c("within"))

		  entrezID <- c("Novel")
	#	  browser()
		  if ((length(t_start) > 0) && (length(t_end) > 0))
		  {	 entrezID_start <- g_GR[subjectHits(t_start)]$gene_id
         entrezID_end <- g_GR[subjectHits(t_end)]$gene_id
         entrezID <- intersect(entrezID_start, entrezID_end)
         entrezID <- try(select(GeneList$Annotation_Library, keys = entrezID, columns=c("SYMBOL"),keytype="ENTREZID")[,'SYMBOL'],silent=TRUE) # Convert to Symbol

         if(length(grep(pattern = "Error in ", x = entrezID)))  # This is start of error message when lookup is not linked
         {  entrezID <- intersect(entrezID_start, entrezID_end)
            entrezID <- try(select(GeneList$Annotation_Library, keys = entrezID, columns=c("SYMBOL"),keytype="ENSEMBL"),silent=TRUE)
            entrezID <- entrezID$SYMBOL
            if(length(grep(pattern = "Error in ", x = entrezID)))
            {    entrezID <- c("Unknown") }
         }
		  }
      DataSet$Gene[i] <- paste(unique(entrezID),collapse=",")
      incProgress(1/MaxDisplay, detail = paste("Updating Row ",i))
	  }  # for(i in 1:MaxDisplay)

	}) # withProgress

	return(DataSet)
}

######################################################################################################
##
## extractGenomeSequence
##
##        Need to pass library prep so that can extract correct sequence.
##
##   August 2016
extractGenomeSequence <- function(chr, start, end, strand, GeneList)
{
  genome <- GeneList$Genome
  BS_junc_Seq<-getSeq(genome,chr,start=start,end=end,strand=strand)
  return(BS_junc_Seq)
}

######################################################################################################
## Grab_BS_Junc_Sequence
##
##
## SelectUniqueJunct_value should be a dataframe with columns names startDonor, strandDonor, startAcceptor.
##        Currently only extracts information from the first row, but in theory could be rewritten to work on multiple entries
##
## NOTE: This was originally designed purely for backsplice junction analysis and because of this the column names are from
## STAR junction output tables. It is now also used for canonical splice junctions.
##
Grab_BS_Junc_Sequence <- function(SelectUniqueJunct_Value, GeneList)
{
  toDisplay <- c("Error, No coordinate data to extract")
  if (nrow(SelectUniqueJunct_Value) >= 1)
  { # Initially assume we have backsplice junciton
    Seq_length <- 125


    Donor <- list()
    Acceptor <- list()
    # Following if statements work for illumina Truseq data. Need to identify if this works for same stranded libraries
    TranscriptStrand <- SelectUniqueJunct_Value$strandDonor[1]
    Transcriptchrom <- SelectUniqueJunct_Value$chromDonor[1]   # This should be same for donor and acceptor

    if ((TranscriptStrand == '-') && (SelectUniqueJunct_Value$type[1] == 'bs'))  # Need to swap donor and acceptor for BS junctions only
    {
      tmp <- SelectUniqueJunct_Value$startDonor[1]
      SelectUniqueJunct_Value$startDonor[1] <- SelectUniqueJunct_Value$startAcceptor[1]
      SelectUniqueJunct_Value$startAcceptor[1] <- tmp
    }


    if (TranscriptStrand == '+')
    { Donor$start <- SelectUniqueJunct_Value$startDonor[1] - Seq_length -1
      Donor$end   <- SelectUniqueJunct_Value$startDonor[1]  - 1    # +1 ensures start in exon as STAR passes coordinates of intron donor sequence
    }
    if (TranscriptStrand == '-')
    {  Donor$end   <- SelectUniqueJunct_Value$startDonor[1] -1
       Donor$start <- SelectUniqueJunct_Value$startDonor[1]  - Seq_length -1
    }

    if (TranscriptStrand == '+')
    { Acceptor$start <- SelectUniqueJunct_Value$startAcceptor[1] +  1    # -1 ensures start in exon as STAR passes coordinates of intron donor sequence
      Acceptor$end   <- SelectUniqueJunct_Value$startAcceptor[1] + Seq_length +1
    }
    if (TranscriptStrand == '-')
    {  Acceptor$end   <- SelectUniqueJunct_Value$startAcceptor[1] + Seq_length + 1
        Acceptor$start <- SelectUniqueJunct_Value$startAcceptor[1] +1
    }

    # Need to test for kit type . i.e. Same strand or opposing.
    if (TranscriptStrand =="+")
      TranscriptStrand <- "-"
    else
      TranscriptStrand <- "+"

    DonorSequence <- extractGenomeSequence(Transcriptchrom, Donor$start, Donor$end, TranscriptStrand, GeneList = GeneList)
    AcceptorSequence <- extractGenomeSequence(Transcriptchrom, Acceptor$start, Acceptor$end, TranscriptStrand, GeneList = GeneList)


    if  (SelectUniqueJunct_Value$type[1] == 'c') # Canonical junction
    { toDisplay <- paste(DonorSequence, AcceptorSequence, sep=".") # Positive strand
      if (TranscriptStrand =="-")
        toDisplay <- paste(AcceptorSequence, DonorSequence, sep=".") # Negative strand
    }

    if (SelectUniqueJunct_Value$type[1] == 'bs')  # backsplice junction
    { toDisplay <- paste(DonorSequence,AcceptorSequence,sep = ".") # Positive strand  <== currently back to front
      if (TranscriptStrand =="-")
        toDisplay <- paste(AcceptorSequence, DonorSequence, sep=".") # Negative strand
    }
  }
  toDisplay <- gsub(pattern=" ",replacement="", x=toDisplay)              # Clean up space characters

  return(toDisplay)

}

#############################################################################################
## miR_Analysis
##
## This function loads miRBase, extracts the species specific miRNAs and then searches for seed patterns
## within sequence
##
## Written by David Humphreys, May 2017
##
miR_binding_site_Analysis <- function(Sequence_to_examine, species_code, seed_length=6)
{
  if (length(seed_length) == 0)
  {  seed_length <- 6 }

  library("mirbase.db")
  species_code <- paste(species_code,"-",sep="")
  ## Need to put checks in regarding species code
          ########### Extracting HUMAN stem loop miRNA  sequences
          x <- mirbaseSEQUENCE
          # Get the microRNA identifiers that are mapped to a SEQUENCE
          mapped_keys <- mappedkeys(x)    # This is a vector of miRNA names
          hsa_mapped_keys <- grep(pattern = species_code,x = mapped_keys)
          stem_loops <- as.list(x[hsa_mapped_keys ])          # Convert human miRNAs stem loop sequences to a list

          ############# Extract start stop coordinates for seed regions
          x <- mirbaseMATURE
          mapped_keys <- mappedkeys(x)
          hsa_mapped_keys <- grep(pattern = species_code,x = mapped_keys)

          # Get the MATURE for ALL elements of xx
          mature_miRs <- mget(mapped_keys[hsa_mapped_keys], x)   # This is a mirnaMATURE class
          b <- lapply(mature_miRs,matureFrom)      # Start positions.   Replace matureFrom with matureTo for end positions
          seed_start <- lapply(b, FUN=function(x) { x + 1})

          # Function to Extract seed sequence
          build_seed_list <- function(stem_loops, seed_start, seed_length)
          {     seedlist <- substr(x = stem_loops, start = as.numeric(seed_start[1]), stop = as.numeric(seed_start[1])+ as.numeric( seed_length)-1)

                if (length(seed_start) > 1)
                {   names(seedlist) <- paste(names(stem_loops), "-5p",sep="")
                    seedlist <- c(seedlist ,substr(x = stem_loops, start = as.numeric(seed_start[2]), stop = as.numeric(seed_start[2])+as.numeric(seed_length) -1) )
                    names(seedlist)[2] <- paste(names(stem_loops), "-3p",sep="")
                }
                return(seedlist)
          }

          allseeds <- unlist( mapply(build_seed_list, stem_loops, seed_start, seed_length) )
          all_unique_seeds <- unique(allseeds)  #unique destroys names
          seed_df <- as.data.frame(allseeds)
          Sequence_to_examine <- reverseComplement(Sequence_to_examine)
          Sequence_to_examine <- gsub(pattern = "T",replacement = "U",x = Sequence_to_examine)
          ## Find seed matches in sequence
          SeedMatchResult <-  mapply(FUN = function(x,y)
                    {   seed_regex <- paste("(?=",x,")",sep="")
                        return( gregexpr(seed_regex ,y, perl=TRUE) )
                    } , all_unique_seeds,as.character(Sequence_to_examine))

          return(list(SeedMatchResult=SeedMatchResult, Total_miR_number=length(all_unique_seeds), miR_Seed_Lookup=seed_df))
}

################################################################################
##  Annotate_with_av_FSJ_coverage
#
##
## File_idx is a list of file indexes
################################################################################
Annotate_with_av_FSJ_coverage <- function(toDisplay_df, GeneList, File_idx, BS_Junctions, Canonical_Junctions, LibraryStrandType)
{
  AllGenes <- unique(toDisplay_df$Gene)
  Gene_FSJ_coverage <- matrix(ncol=length(File_idx), nrow= length(AllGenes), data=0) # rep(0,length(AllGenes))
  row.names(Gene_FSJ_coverage) <- AllGenes
  toDisplay_df$BSJ_vs_FSJ <- 0
  Num_Columns <- ncol(toDisplay_df)   # Record number of columns in data table.
  canonical_flanking_Junctions <- matrix(data = 0,nrow = length(toDisplay_df$Gene),ncol = (Num_Columns-3))
  Internal_FSJ_of_circRNA <- matrix(data = 0,nrow = length(toDisplay_df$Gene),ncol = (Num_Columns-3))    # To record FSJ counts that are internal of BSJ
  External_FSJ_of_circRNA <- Internal_FSJ_of_circRNA      # To record FSJ counts that are outside BSJ of circRNA

withProgress(message="Annotating with FSJ coverage", value=0, {

  for( i in 1: length(AllGenes))
  { incProgress(1/length(AllGenes), detail = paste("Updating ",i, " of ", length(AllGenes) ))
    GeneName <- AllGenes[i]
    if ((GeneName == "Novel") || (GeneName == "Unknown"))
      next
    if (length(grep(pattern=",",x = GeneName)) > 0)   # multiple possible genes. Skip.
      next

    GeneIdx <- which(toDisplay_df$Gene == GeneName)   # This stored index of mutiple or single gene entries
    cat(paste("\n Gene is ",GeneName,". Indexes of:",GeneIdx))
    if (GeneName == '')
    { next }
    for (g in 1:length(GeneIdx))
    {
      for (j in 1: length(File_idx))
      {   current_BSJ_coords <- range(as.numeric(unlist(strsplit(x = toDisplay_df$BSjuncName[GeneIdx[g]],split = "_"))[c(2,4)]))   # Grab genomic coordinates of BSJ
          PGO<-Prepare_Gene_Object(GeneName, BS_Junctions = BS_Junctions,
                               GeneList= GeneList, File_idx = File_idx[[j]],
                               Canonical_Junctions = Canonical_Junctions)

          if (dim(PGO$Transcript_Canonical_juncs[start1< current_BSJ_coords[1] & end2 > current_BSJ_coords[2],])[1] > 0)
          { # A FSJ that may have resulted from BSJ formation. At this stage save all/any counts to report.
            canonical_flanking_Junctions[GeneIdx[g],j] <- canonical_flanking_Junctions[GeneIdx[g],j] + sum(PGO$Transcript_Canonical_juncs[start1< current_BSJ_coords[1] & end2 > current_BSJ_coords[2],.(score)])
          }
          if (dim(PGO$Transcript_Canonical_juncs[start1>= current_BSJ_coords[1] & end2 <= current_BSJ_coords[2],])[1] > 0)
          { # Record FSJ counts that are located inside BSJ coordinate.
            Internal_FSJ_of_circRNA[GeneIdx[g],j] <- round(mean(unlist(PGO$Transcript_Canonical_juncs[start1> current_BSJ_coords[1] & end2 < current_BSJ_coords[2],.(score)])),1)
          }
          if (dim(PGO$Transcript_Canonical_juncs[start1< current_BSJ_coords[1] | end2 > current_BSJ_coords[2],])[1] > 0)
          { # Record FSJ counts that are located inside BSJ coordinate.
            External_FSJ_of_circRNA[GeneIdx[g],j] <- round(mean(unlist(PGO$Transcript_Canonical_juncs[start1< current_BSJ_coords[1] | end2 > current_BSJ_coords[2],.(score)])),1)
          }

          GeneFeatures <- Gene_Transcript_Features(GeneList=GeneList, Gene_Symbol=GeneName, GeneObject=PGO)
          if (is.null(GeneFeatures))
          {
            Gene_FSJ_coverage[GeneName,j] <-  -1    }
          else
          {  Gene_FSJ_coverage[GeneName,j] <- GeneFeatures$Av_junc_count }
      }
    } # for (g in 1:length(GeneIdx))
    DT_idx <- which(toDisplay_df$Gene == GeneName)

#    if ((! is.na(Gene_FSJ_coverage[GeneName])) && (Gene_FSJ_coverage[GeneName] > 0))   # Be surprised if this condition was FALSE
    { if (length(which(colnames(toDisplay_df) == "Freq") >0))
      {   toDisplay_df$BSJ_vs_FSJ[DT_idx] <- round(unlist(toDisplay_df[DT_idx,.(Freq)]/Gene_FSJ_coverage[GeneName,1]),2)   }
      else # Have group data
      {
        toDisplay_df <- as.data.frame(toDisplay_df)
        for (j in 1:(Num_Columns-3))
        {
           if (is.null(toDisplay_df[DT_idx,Num_Columns+j]))
           {   toDisplay_df[,Num_Columns+j] <- 0  }

           toDisplay_df[DT_idx,Num_Columns+j] <- round(unlist(toDisplay_df[DT_idx,j+1]/Gene_FSJ_coverage[GeneName,j]),2)

        } # for (j in 1:(Num_Columns-2))
        toDisplay_df <- as.data.table(toDisplay_df)

      }
    }

#    toDisplay_df$FSJ_av[DT_idx] <- round(Gene_FSJ_coverage[GeneName],1)
  } # for( i in 1: length(AllGenes))
})
  colnames(Internal_FSJ_of_circRNA) <- paste("Internal",1:ncol(Internal_FSJ_of_circRNA),sep="_")
  colnames(External_FSJ_of_circRNA) <- paste("External",1:ncol(Internal_FSJ_of_circRNA),sep="_")
  colnames(canonical_flanking_Junctions) <- paste("Canonical",1:ncol(canonical_flanking_Junctions),sep="_")
  toDisplay_df <- cbind(toDisplay_df,canonical_flanking_Junctions, Internal_FSJ_of_circRNA, External_FSJ_of_circRNA)
  return(toDisplay_df)
}
################################################################################
## Predicted_CircRNA_Sequence
##
## Extracts all exon sequences that are within BSJ coordinate and concatenates them
##
Predicted_CircRNA_Sequence <- function(circRNA_exons, genelist)
{
  circRNA_exon_lengths <-by(circRNA_exons,circRNA_exons$gene,identity )  # This makes a list of all transcripts
  if (length(circRNA_exon_lengths) ==0)
    return(NULL)
  a <- lapply(X = circRNA_exon_lengths,FUN = function(x) { return(abs(sum(x$start-x$stop)))})
  a_idx <- order(unlist(a), decreasing = TRUE)
  largest_transcript_ID <- names(a[a_idx])[1]
  transcript_ID_idx <- which(names(circRNA_exon_lengths) == largest_transcript_ID)  # Get Index of transcript ID name
  transcript_ID_idx <- which(circRNA_exons$gene == names(circRNA_exon_lengths)[transcript_ID_idx])  # Get Row index(es) corresponding to transcript
  Exons_of_Interest <- circRNA_exons[transcript_ID_idx,]

  circRNA_Sequence <- ''
  FSJs <- c(1)# This will contain start positions for ALL Forward splice junctions
  for (i in 1:length(transcript_ID_idx))   # Need to stitch together multiple exons
  { transcript_strand <- "+"
    if (Exons_of_Interest$strand[i] == 0)
    {   transcript_strand <- "-"  }
    tmp <- as.character(extractGenomeSequence(Exons_of_Interest$chrom[i],Exons_of_Interest$start[i],
                                            Exons_of_Interest$stop[i], transcript_strand,
                                            GeneList = genelist)  )
    FSJs <- c(FSJs, FSJs[i] + nchar(tmp))
    circRNA_Sequence <- DNAString(paste(circRNA_Sequence,tmp,sep="",collapse = ""))
  }

  circRNA_Sequence <- DNAString(circRNA_Sequence)
  return(circRNA_Sequence)
}


debugme<-function(aa, bb, cc)
{	x <- 1 + 1
	y <- 1+ 1
}

######################################################################################################################################
######################################################################################################################################
###########                       SHINY  SERVER                                              #########################################
######################################################################################################################################
######################################################################################################################################

shinyServer(function(input, output, session)
{	# input$JunctionFile will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.

debug(debugme)

##  showReactLog()

# Global variables
  extdata_path <- as.character(DataPath())    # The data path used to save and load project data from

  volumes <- c('Ularcirc'= extdata_path)      ## R.home()    or getVolumes()
  shinyDirChoose(input, 'dir', roots = volumes, session=session, restrictions=system.file(package='base'))
  dir <- reactive(input$dir)


  Ularcirc_data = reactiveValues(
    GenePanelLoaded = TRUE,
    Groupings = list(),                    # This holds the group information
    SelectedGene_from_BSJ_Table = NULL,    # The Junction selected in BSJ table. USed to update gene list.
    Current_SelectedGene = NULL,           # The current gene name of focus. Selected from Pull down list or by clicking table row of data
    Current_Selected_BS_Junction = NULL,   # The Back splice junction of focus
    Current_Selected_BS_Junction_RAWData = NULL,    # A datatable of raw data for the Current_Selected_BS_Junction
    SelectedTranscript = NULL,             # The transcript selected in Transcript table on gene view page
    SelectedCanonical = list(Chr= NULL, Start=NULL, End= NULL, strand=NULL),  # User selected canonical junction
    BackSpliceJunctionCountTable = NULL,   # Backsplice junction table
    CanonicalJunctionCountTable = NULL,    # Canonical Junction Count table
    GenomeCanonicalJunctionCountTable = NULL,    # Canonical Junction Count table for Genome tab menu
    selected_circRNA_stats = list(miRNA_BS_Sites = NULL, ORFs= NULL),
    Genome_Coordinates = list(chrom=NULL, chromstart=NULL, chromend=NULL, chromstrand=NULL),
    SelectedGenome_FSJ = list(Chr=NULL, Start=NULL, End=NULL, Strand=NULL),   # This is for FSJ selected in "Genome" tab
    ProjectData  = list()                  # This contains all raw data of a project file.
  )

	captionText <- reactiveValues ()

	output$FileNameDataTable <- renderTable({
		# input$JunctionFile is made up of: name, size, type, data path

		if (is.null(input$JunctionFile))
			return(NULL)

		inFile<-as.data.frame(input$JunctionFile)
		DataSet <- Ularcirc_data$ProjectData # m379()  #Junctions=AllData, Column_Numbers=total_rows)))
		ojn <- DataSet$Original_Junction_Numbers
		onuj <- DataSet$Original_n_unique_junctions
		opfnuj <- DataSet$Original_Postfilt_n_unique_junctions

		fuj <- length(table(DataSet$Junctions$BSjuncName))	# filtered unique junctions

		captionText$output <- paste("Total number of junctions =",nrow(DataSet$Junctions),
					"<br>Final number of unique junctions",fuj," ",date())
		cbind(inFile[,1:2], Input_Junctions=ojn, Unique_Junctions=onuj, filtered_Unique_Junctions=opfnuj)#,
			#		Filtered_Junction_Entries=nrow(DataSet$Junctions), Filtered_Unique_Junctions= fuj)
		}) # output$FileNameDataTable

	output$FileNameDataTableDetails <- renderText({
			if (is.null(input$JunctionFile))
			{	return(paste("No data file uploaded yet therefore nothing to display"))   }

			DataSet <- Ularcirc_data$ProjectData # m379()
			jn <- nrow(DataSet$Junctions)
			fujn <- length(table(DataSet$Junctions$BSjuncName))	# filtered unique junctions

			c("\ng","\n","After filtering:",jn," number of junctions and ", fujn," unique junctions.\n-",
				length(which(DataSet$Junctions$type == 'c')) ," Canonical junctions\n -",
				length(which(DataSet$Junctions$type == 'bs')) ," Backsplice junctions\n -", sep="")
	})

	output$DisplaySelectedFileNames <- renderDataTable({
	#		debugme(input$JunctionFile, input$SelectedFiles, input)
			SelectedFile <- input$JunctionFile  ### UP TO HERE
			cbind(Selected=, input$SelectedFiles)
	})


	## Following two assignments coordinate between UI.R and Server.R the uploading of a file.
	output$fileUploaded <- reactive({
	  return(!is.null(Ularcirc_data$ProjectData)) # (m379()))
	  })
	outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)


	output$Filename <- renderText({ 	# Would like to change this to a table
		inFile <- input$JunctionFile
		inFile$name		# name is name of input file; datapath: contains path to tempory uploaded data file
		})


	 output$InputFiles <- renderUI ({   # This draws a checkboxInput of all uploaded data sets (files) for UI
	  inFile = Ularcirc_data$ProjectData$SampleIDs  #m379()$SampleIDs
		if (is.null(inFile))
			return(NULL)
		else
		  # Would like to make following column have a colour... cannot get to work currently
		  w <- tags$hr(style="border-color:black")
			w <- paste(column(12,checkboxGroupInput('SelectedFiles','', inFile,inFile),br(),w))#, style="background-color:#4d3a7d;"))
			HTML(w)
	})


  IdentifyDataSets <- function()
  {   idx <- {}
      if (! is.null(input$SelectedFiles))
	    { idx <- which( Ularcirc_data$ProjectData$SampleIDs == input$SelectedFiles) } # m379()$SampleIDs == input$SelectedFiles)  }
	    Selected_DataSets <- warning("No data loaded or selected so nothing to display",
	                               "\nPlease navigate back to PROJECT tab and load a data set")
	    if (length(idx) > 0)
	    {	Selected_DataSets <- paste("\nDisplaying data from ", Ularcirc_data$ProjectData$SampleIDs[idx])}   # m379()$SampleIDs[idx])}

	    Selected_DataSets
  }

	output$ShowDataSets_on_GeneView <- renderText({   IdentifyDataSets()	})
	output$ShowDataSets_on_Genome_View <- renderText({   IdentifyDataSets()	})
	output$ShowDataSets_on_DataView <- renderText({   IdentifyDataSets()	})
	output$ShowDataSets_on_JunctionView <- renderText({   IdentifyDataSets()	})

	m379 <- observeEvent(input$JunctionFile,{
	  #extdata_path <- DataPath()  # function from Global.R
		inFile <- input$JunctionFile

		cat(paste("\nLoading data", inFile))

		Ularcirc_data$ProjectData <- LoadJunctionData(filename = inFile,
			ChromFilter = input$ChromosomeFilter, StrandFilter = input$StrandFilter,
			Genomic_Distance = input$GenomicDistance,
			Junction_abundance = input$Junction_abundance, RAD_filter = input$RAD_filter,
			MultipleDataOption = input$MultipleDataOption, CanonicalJuncs=input$CanonicalJuncs,
			SubSelect= SubsetInputData)

#		return(Ularcirc_data$ProjectData)
	})

	  observe({     # The following code ensure the annotate with gene name does NOT get selected when no genome selected
	    temp <- input$Annotate_with_GeneName
	    if (input$Annotation_lib == "NO_ANNOTATION")
	    {
	      updateCheckboxInput(session, 'Annotate_with_GeneName', value = FALSE)
	    }

	     })

	# This function will build table based on pooling of selected individual data sets.
	# This does NOT pool grouped samples,
	PartialPooledDataSet <- observeEvent(input$Annotate_Option_Submit_Button, { # reactive({

	  Identify_poor_RAD <- function(x) ## Function identifies which table entries have non-complying RAD score
	  {   idx_to_remove <- {}
	  x <- as.numeric(x)
	  idx_to_remove <- c(which(is.na(x)), which(x < input$RAD_filter[1]), which(x>input$RAD_filter[2]))
	  return(idx_to_remove)
	  }

	  if (! input$Annotate_Option_Submit_Button)# If annotation button has not been pressed do nothing.
	  { return(data.frame(c(ACTION_REQUIRED="Please select annotation method on left hand tab"))); }

	  inFile = Ularcirc_data$ProjectData$SampleIDs  #m379()$SampleIDs      # This contains all possible input files (samples)

	  toDisplay <- list(RAW={}, CPM={}, TOTAL_COUNTS = {})
	  if (length(grep(pattern = "Selected sample", x = input$Annotation_Options)) > 0)
	  {
	    inFile_idx <- input$SelectedFiles
	    idx <- list()

	    if (! is.null(inFile_idx))
	    { idx[[1]] <- which( Ularcirc_data$ProjectData$SampleIDs == inFile_idx) } #m379()$SampleIDs == inFile_idx)  }

	    SubsettedData <- Ularcirc_data$ProjectData  # m379()   # By default collect ALL data
	    if (length(idx[[1]]) == 0)
	    {	cat(paste("\nNo data set selected from following samples ", Ularcirc_data$ProjectData$SampleIDs))
	      Ularcirc_data$PartialPooledDataSet <- data.frame(c(ERROR="No data set selected. "))
	      return(-1)
	    }
	    else
	    {	# Should change FilterChimeric to SelectUniqueJunctions
	      BSJ_junctions <- Filter_by_Data_Set(fileID=idx[[1]], All_junctions = Ularcirc_data$ProjectData$Junctions)

        filterlist <- list(BSjuncName=NULL, SortDir="Descending", IndexNumber=1,
                           DisplayNumber=input$MAX_BS_juncs_to_annotate,
                           DisplayRAD_score= input$Display_RAD_Score)
	      SubsettedData <-  SelectUniqueJunctions(All_junctions=BSJ_junctions, filterlist = filterlist, input$LibraryStrandType)

	    }

	    if (nrow(SubsettedData) > 0)
	    {
	        toDisplay$RAW <- SubsettedData
	        toDisplay$TOTAL_COUNTS <- sum(SubsettedData$Freq)

	        if ((input$Annotate_with_GeneName) && (input$Annotation_lib != "NO_ANNOTATION"))
	        { # Annotate all entries with Gene name
	          toDisplay$RAW <- Annotate_BS_Junc(DataSet=toDisplay$RAW, GeneList = GeneList(), MaxDisplay = nrow(toDisplay$RAW), input$LibraryStrandType)
	        }
	        else
	        {  toDisplay$RAW$Gene <- '' }

	        if (input$Display_RAD_Score)
	        {  idx_to_remove<- Identify_poor_RAD(toDisplay$RAW$TypeII_TypeIII)
	           if (length(idx_to_remove) > 0)
	           { toDisplay$RAW <- toDisplay$RAW[-1*idx_to_remove,] }
	        }
	        toDisplay$RAW <- toDisplay$RAW[,.(Gene, Freq,BSjuncName,TypeII_TypeIII, strandDonor)]

	        if (input$Percent_of_Parent)  # If requested annotate with parental transcript abundance
	        { # Get list of all genes and build table
	          toDisplay$RAW <- Annotate_with_av_FSJ_coverage(toDisplay_df=toDisplay$RAW , GeneList=GeneList(), File_idx=idx,
	                                        BS_Junctions=Ularcirc_data$ProjectData$Junctions,
	                                        Canonical_Junctions= Ularcirc_data$ProjectData$Canonical_AllData,
	                                        LibraryStrandType = input$LibraryStrandType)
	        }

	        toDisplay$CPM <- toDisplay$RAW
	        toDisplay$CPM$Freq <- round((toDisplay$RAW$Freq / toDisplay$TOTAL_COUNTS * 1000000),2)
	        Ularcirc_data$PartialPooledDataSet <- toDisplay
	    }
	    else
	      Ularcirc_data$PartialPooledDataSet <- data.frame(c(NO_DATA="After filtering no data returned"))
	  }


	  if ((length(grep(pattern = "Grouped analysis", x = input$Annotation_Options)) > 0)   # Prepare comparison table
	      && (! is.null(PrepareGroupOptions())) )   # Check there are groups defined. There should be always at least one.
	  {
	    AllGroupIDs <- paste("Group",seq(from=1, to=input$Number_BiologicalSamples, by=1),sep="_")
	    inFile = Ularcirc_data$ProjectData$SampleIDs # m379()$SampleIDs      # This contains all possible input files (samples)
	    a<- length(Groupings)
	    SubsettedData <- list()    # This holds the top requested number of BSJ
	    BSJ_junctions <- list()    # This hold ALL BSJs
	    idx <- list()              # Hold all file indexes for each group
	    data_set_idx <- 0          # This is the list index to both SubsettedData and BSJ_junctons

withProgress(message="Calculating BSJ : ", value=0, {
	    for(i in 1:a)   # This loop collects all data. Need to keep a copy of everything so in next loop can collate easily
	    {
	      incProgress(1/(a), detail = paste("Sample ",i))
	      SampleIDs_for_current_group <- Groupings[[i]]

	      if (! is.null(inFile))
	      {	data_set_idx <- data_set_idx + 1
	        for(j in 1:length(SampleIDs_for_current_group))     # Extract file IDs
	        {  if (j == 1)
	           { idx[[data_set_idx]] <- c(which( SampleIDs_for_current_group[j] == inFile))  }
	           else
	           { idx[[data_set_idx]] <- c(idx[[data_set_idx]], which( SampleIDs_for_current_group[j] == inFile))   }
	        }
	      }
	      else
	      {  # browser()
	         stop("inFile is NULL, unexpected and unrecoverable error")
	      }

	      if (length(idx[[data_set_idx]]) > 0)
	      {

	        BSJ_junctions[[data_set_idx]] <- Filter_by_Data_Set(fileID=idx[[data_set_idx]], All_junctions = Ularcirc_data$ProjectData$Junctions)

	        filterlist <- list(BSjuncName=NULL, SortDir="Descending", IndexNumber=1,
	                           DisplayNumber=input$MAX_BS_juncs_to_annotate,
	                           DisplayRAD_score= input$Display_RAD_Score)
	   #     browser()
#	        SelectUniqueJunctions(All_junctions=BSJ_junctions, filterlist = filterlist, input$LibraryStrandType)
	        SubsettedData[[data_set_idx]] <-  SelectUniqueJunctions(All_junctions=BSJ_junctions[[data_set_idx]], filterlist = filterlist, library_strand = input$LibraryStrandType) # This is limited by number BSJ to display on screen

	      }
	    } # for(i in 1:a)
})  # withProgress(message="Calculating BSJ : ", value=0, {


withProgress(message="Fixing blank BSJ : ", value=0, {
	    # This loop collates data by assembling a table and filling in the blanks (NAs) where possible
	    toDisplay$TOTAL_COUNTS <- {}
	    toDisplay$RAD_Score <- {}        # Type II vs type III junctions
	    for ( i in 1:data_set_idx)
	    { OneDataSet <- list()

	      incProgress(1/data_set_idx, detail = paste("Sample ",i))
  	    TotalCounts <- nrow(BSJ_junctions[[i]])   #SubsettedData[[i]]$Freq)
  	    toDisplay$TOTAL_COUNTS <- as.numeric(c(toDisplay$TOTAL_COUNTS, TotalCounts))

  	    OneDataSet$Counts <- data.table(BSjuncName=SubsettedData[[i]]$BSjuncName, CPM=round(x = SubsettedData[[i]]$Freq,digits = 0))
  	    OneDataSet$RAD_Score  <- data.table(BSjuncName=SubsettedData[[i]]$BSjuncName, RAD= SubsettedData[[i]]$TypeII_TypeIII)

  	    if (nrow(OneDataSet$Counts) > 0)  	      # Merge if have extracted data
  	    {
  	      if (! is.null(toDisplay$RAW))
  	      {
  	        # MERGE counts
  	        toDisplay$RAW <-  merge(toDisplay$RAW, OneDataSet$Counts, by="BSjuncName", all=TRUE)
  	        toDisplay$RAW <- as.matrix(toDisplay$RAW)
  	        colnames(toDisplay$RAW)[i+1] <- AllGroupIDs[i]
  	        # Merge RAD scores
  	        toDisplay$RAD_Score <-  merge(toDisplay$RAD_Score, OneDataSet$RAD_Score, by="BSjuncName", all=TRUE)
  	        toDisplay$RAD_Score <- as.matrix(toDisplay$RAD_Score)
  	        colnames(toDisplay$RAD_Score)[i+1] <- paste(AllGroupIDs[i],"_II_III",sep="")

  	        for (j in 1:i)    #data_set_idx)
  	        { # Identify which BSJ in other data sets that do not have a count i.e. are currently annotated as NA.
  	          NA_idx <- which(is.na(toDisplay$RAW[,j+1]))
  	          NA_IDs <- toDisplay$RAW[NA_idx,c("BSjuncName")]
  	          if (length(NA_idx) == 0)  # Great no NA values, no more work to do for this sample.
  	          {  next  }

  	          countMatches <- function(x,RawData)  # This function is used to count BSJ annotated as NA in current sample
  	          {  return(table(RawData == x)) }

  	          b <- sapply(X = NA_IDs, FUN = countMatches, RawData = BSJ_junctions[[j]]$BSjuncName)

  	          if (typeof(b) == "integer")
  	          {   b<- t(b) }
  	          else
  	          { b <- lapply(b, FUN = function(X) { if (is.na(X["TRUE"])) {X["TRUE"] <- 0}; return(X) })
  	            b<- do.call(rbind, b)
  	          }

              if (is.null(dim(b)))      # Could not identify circRNA in current sample. Set them to zero
  	          { b <- 0  }
  	          else
  	          {   if (dim(b)[1] == 1)
  	              {
  	                if (length(grep(pattern = "TRUE",x = colnames(b))) == 0)
  	                {   b <- 0   }
  	                else
  	                {   b <- b[,"TRUE"] }
  	              }
  	              else
      	          {  b <- b[,"TRUE"]  }
  	          }
  	  	      toDisplay$RAW[NA_idx,j+1] <- round(x = b,digits = 0)   # Assign count
  	        }
  	        toDisplay$RAW <- as.data.table(toDisplay$RAW)
  	      }
  	      else # is.null(toDisplay$RAW)     #  First time through loop
  	      { toDisplay$RAW <- OneDataSet$Counts
  	        colnames(toDisplay$RAW) <- c("BSjuncName", AllGroupIDs[i])
  	        toDisplay$RAD_Score <-  OneDataSet$RAD_Score
  	        colnames(toDisplay$RAD_Score) <- c("BSjuncName", paste(AllGroupIDs[i],"_II_III",sep=""))
  	      }
  	    }
  	    else  # No data for some reason. This scenario should not happen.
  	    {
  	      toDisplay$New_____Sample <- 0
  	      colnames(toDisplay$RAW) <- c(colnames(toDisplay$RAW), AllGroupIDs[i])
  	    }
	    } # for ( i in 1:dataset_idx)
	    # Order table for most variable differences.
})   # withProgress(message="Fixing blank BSJ : ", value=0, {

#	    most_variable <- apply(toDisplay$RAW[,-1],1,var)
#	    Top_variable<-toDisplay$RAW[order(most_variable ,decreasing = T ),]
#	    colnames(Top_variable)<-colnames(toDisplay$RAW)
	    # Need to construct a properly constructed data.table:


	    Top_variable <- toDisplay$RAW
	    toDisplay$RAW <- data.table(Top_variable[,1])
	    for(i in 2:ncol(Top_variable))
	    { toDisplay$RAW <-  cbind(toDisplay$RAW, as.numeric(as.matrix(Top_variable)[,i]))  }
	    colnames(toDisplay$RAW) <- colnames(Top_variable)

	    toDisplay$CPM <- toDisplay$RAW
	    toDisplay$CPM[,2:ncol(toDisplay$RAW)] <- round(toDisplay$RAW[,2:ncol(toDisplay$RAW)]/toDisplay$TOTAL_COUNTS*1000000,digits = 0)

	    if (length(input$LibraryStrandType) == 0)
	    { # browser()
	      }

	    if ((input$Annotate_with_GeneName) && (input$Annotation_lib != "NO_ANNOTATION"))
	    { toDisplay$RAW <- Annotate_BS_Junc(DataSet=toDisplay$RAW, GeneList = GeneList(), MaxDisplay = nrow(toDisplay$RAW), input$LibraryStrandType)
	    }
	    else
	    {  toDisplay$RAW$Gene <- '' }

	    if (input$Percent_of_Parent)   	    ## Add parental transcript abundance annotation here
	    {
	        toDisplay$RAW <- Annotate_with_av_FSJ_coverage(toDisplay_df=toDisplay$RAW , GeneList=GeneList(), File_idx=idx,
	                                                       BS_Junctions=Ularcirc_data$ProjectData$Junctions,
	                                                       Canonical_Junctions= Ularcirc_data$ProjectData$Canonical_AllData,
	                                                       LibraryStrandType = input$LibraryStrandType)
	    }

	    ## Following code will identify and remove BSJ that don't meet RAD filter.
	    if (input$Display_RAD_Score)
	    {
  	    idx_to_remove <- {}
  	    for(i in 2 : ncol(toDisplay$RAD_Score))
  	    { if (i == 2)
  	      { idx_to_remove<- Identify_poor_RAD(toDisplay$RAD_Score[,i]) }
  	      else
  	      { idx_to_remove<- intersect(idx_to_remove, Identify_poor_RAD(toDisplay$RAD_Score[,i]))  }
  	    }

  	    if (length(idx_to_remove) > 0)
  	    {
  	      toDisplay$RAD_Score <- toDisplay$RAD_Score[-1*idx_to_remove,]
  	      toDisplay$RAW <- toDisplay$RAW[-1*idx_to_remove,]
  	      toDisplay$CPM <- toDisplay$CPM[-1*idx_to_remove,]
  	    }
  	  }
	    # Append RAD scores onto table
	    toDisplay$RAW <- merge(toDisplay$RAW, toDisplay$RAD_Score, by="BSjuncName", all=TRUE)


	    toDisplay$CPM$Gene <- toDisplay$RAW$Gene
	    # redblueColor<-c(colorRampPalette(c("blue","white"))(35),colorRampPalette(c("white","red"))(35));
	    # heatmap(as.matrix(Top_variable ),col=redblueColor,scale = "none")
	    Ularcirc_data$PartialDataSet <- toDisplay
	  }  # length(grep(pattern = "Grouped analysis", x = input$Annotation_Options)) > 0)
})


  GeneList <- eventReactive( input$LoadTxDb,   # GeneList is called when transcript database is selected
	{ require(GenomicFeatures)
	  cat(paste("\n\nLoading species transcriptome coordinates", date()))
	  Genome_lib <- paste("BSgenome.",input$Species_Genome,sep="")
	  require(as.character(Genome_lib),character.only = TRUE)
	  Genome<-get(Genome_lib)     # Reference object by character string

	  # Determine if using custom or installed gene model
	  TxDb <- {}
	  if (length(grep(pattern = "*.sqlite", x = input$TxDb)) > 0)
	  {
	    TxDb <- loadDb(file=input$TxDb)
	  }
	  else
	  {
	    require(as.character(input$TxDb),character.only = TRUE)
	    TxDb <- get(input$TxDb)       # Reference object by character string
	  }


	  require(as.character(input$Annotation_lib),character.only = TRUE)
	  Annotation_Library <- get(input$Annotation_lib)
	  GL <- as.character(keys(Annotation_Library, "SYMBOL"))

		##### Setup pulldown menu of genes. Default is first item in list ########
		cat(paste("\n  Displaying list of ", length(GL)," genes built", date()))
		updateSelectizeInput(session,inputId="GeneListDisplay", label="Select a gene from this list", choices = GL, selected = GL[1], server=TRUE)

		return(list(transcript_reference=TxDb, Genome=Genome, Annotation_Library = Annotation_Library, GeneList=GL))
	})


  output$List_Loaded_TxDB <- renderText({
    TxDb_information <- c('No Transcript database loaded')
    if (! is.null(GeneList()))
    {
      TxDb_information <- c (input$Species_Genome)
    }

    TxDb_information
  })


	circRNA_Subset <- reactive (  # This code block will prepare data sets surrounding a gene feature.
	  { if (is.null(input$SelectedFiles))
	  { return(NULL)}
	    idx <- seq(from=1, to = length(Ularcirc_data$ProjectData$SampleIDs)) # m379()$SampleIDs))     # Default is to select everything. PROBABLY BETTER TO RETURN A NULL IF NO DATA SELECTED AND DEAL WITH THIS ELSEWHERE??
	    if (! is.null(input$SelectedFiles))
	    { idx <- which( Ularcirc_data$ProjectData$SampleIDs == input$SelectedFiles) } #m379()$SampleIDs == input$SelectedFiles)  }

	  #  if (is.null(input$GeneList))
	   # { cat(paste("\nGene list not currently defined"))
	   #   return(NULL)
	   #}

		  cat(paste("\nRequesting Gene Object for", Ularcirc_data$Current_SelectedGene,date()))

		  PGO<-Prepare_Gene_Object(Ularcirc_data$Current_SelectedGene, BS_Junctions = Ularcirc_data$ProjectData$Junctions,#m379()$Junctions,
		                                      GeneList= GeneList(), File_idx = idx,
		                                      Canonical_Junctions = Ularcirc_data$ProjectData$Canonical_AllData) #m379()$Canonical_AllData,
		                                      # Genome_Coords = Ularcirc_data$Genome_Coordinates)
		  if (is.null(PGO))
		  {  return(PGO) }

		  Ularcirc_data$CanonicalJunctionCountTable <- PGO$Transcript_Canonical_juncs[,.(chrom1,start1,end2,strand1,score,DataSet,multimappers,overhang)]

		  if (is.double(PGO$Junctions))  #i.e. has a value of -1, means there is no data
		    Ularcirc_data$BackSpliceJunctionCountTable <- data.table(Status="No data points")
		  else
		    Ularcirc_data$BackSpliceJunctionCountTable <- PGO$Junctions$uniques.bed[,.(chrom1, start1,end2,strand1,score,name)]
		    # Ularcirc_data$BackSpliceJunctionCountTable <- PGO$Junctions$uniques.bed
		  cat(paste("\n-Prepared Gene Object",date()))
		  return(PGO)    # 	PGO =  (list(Transcript=Transcript, Junctions=Junc.bed, Transcript_Canonical_juncs= Transcript_Canonical_juncs))
	})




  # The output$view depends on both the databaseInput reactive
  # expression and input$obs, so will be re-executed whenever
  # input$dataset or input$obs is changed.
	output$caption <- renderText({ "ExonTable" })

	output$TranscriptTable <- renderDataTable ({
	  Transcript_Stats <- table(circRNA_Subset()$Transcript$gene)
	  Transcript_DT <- data.table(TranscriptID=names(Transcript_Stats), Number_Of_Exons = as.numeric(Transcript_Stats))
	  datatable(Transcript_DT, selection='single',options = list(lengthMenu = c(5, 10, 50), pageLength = 5))
	})

	output$ExonTable <- renderDataTable ({
	    if (is.null(Ularcirc_data$SelectedTranscript))
	    {   return(NULL) }
	    Row_idx <- which(circRNA_Subset()$Transcript$gene == Ularcirc_data$SelectedTranscript)
	    Ularcirc_data$CurrentExonTable <- data.table(circRNA_Subset()$Transcript[Row_idx,])
			datatable(Ularcirc_data$CurrentExonTable, selection='single', options = list(lengthMenu = c(5, 10, 50), pageLength = 5))
			})


	output$DisplayJunctionCountTable<- renderDataTable({   # DisplayAllJunctions <- renderDataTable({
	    if (! input$Annotate_Option_Submit_Button)# If annotation button has not been pressed do nothing.
	    { return(data.frame(c(ACTION_REQUIRED="Please select annotation method on left hand tab"))); }

	    toDisplay <- data.frame(c(ERROR="Cannot build data sets. Please check files have been selected and try again"))

	    if (length(grep(pattern = "Selected sample", x = input$Annotation_Options)) > 0)
	    {
	      if (input$Normalisation == "Raw counts")
  	      toDisplay <- Ularcirc_data$PartialPooledDataSet$RAW
	      else
  	      toDisplay <- Ularcirc_data$PartialPooledDataSet$CPM
	    }

	    if ((length(grep(pattern = "Grouped analysis", x = input$Annotation_Options)) > 0)   # Prepare comparison table
	      && (! is.null(PrepareGroupOptions())) )   # Check there are groups defined. There should be always at least one.
	    {
	      if (input$Normalisation == "Raw counts")
	        toDisplay <- Ularcirc_data$PartialDataSet$RAW
	      else # ,"CPM")
	        toDisplay <- Ularcirc_data$PartialDataSet$CPM
	    }

	    if (is.null(toDisplay))
	    { toDisplay <- data.table(ACTION_REQUIRED="Please select annotation action button in left hand panel")
	    }

	    datatable(toDisplay, selection = 'single', options = list(lengthMenu = c(10,50,500,5000), pageLength = 15))
			#Annotate_BS_Junc(DataSet=toDisplay, GeneList = GeneList(), MaxDisplay = nrow(toDisplay))
	})

	PrepareGroupOptions <- reactive({
	  w <- {}
	  if ((input$Number_BiologicalSamples > 1) && (length(grep(pattern = "Normalised counts", x = input$Annotation_Options)) > 0))
	  {
	    AllGroupIDs <- paste("Group",seq(from=1, to=input$Number_BiologicalSamples, by=1))
	    w<- paste(selectizeInput("GroupCompareA",label="Comparison group A",choices =  AllGroupIDs))
	    w<- paste(w, selectizeInput("GroupCompareB",label="Comparison group B",choices =  AllGroupIDs))
	  }
	  HTML(w)
	})

	output$TwoGroupCompareChoices <- renderUI({   # This will display selectize menu in UI.R to allow user to select groups to compare.
      PrepareGroupOptions()
	})

	output$BS_Junction_Count_Table <- renderDataTable({
	  	if (is.null(Ularcirc_data$BackSpliceJunctionCountTable))
	  	{ return(NULL)}
#Display_DT <- circRNA_Subset()$Junctions$uniques.bed[,.(name,score,chrom1,start1,chrom2,start2,strand1,strand2)]
			datatable(Ularcirc_data$BackSpliceJunctionCountTable, selection='single', options = list(lengthMenu = c(5, 10, 50), pageLength = 5))
  })

	output$CanonicalJunctionCountTable <- renderDataTable({
	      datatable(Ularcirc_data$CanonicalJunctionCountTable, selection='single', options = list(lengthMenu = c(5, 10, 50), pageLength = 5))
	  })

	output$GenomeCanonicalJunctionCountTable <- renderDataTable({
	  datatable(Ularcirc_data$GenomeCanonicalJunctionCountTable, selection='single', options = list(lengthMenu = c(5, 10, 50), pageLength = 5))
	})

  	## Observe if junction button is pressed.
	observeEvent(input$select_button, {
			if (is.null(input$select_button))
				return(NULL)

	#			debugme(input$select_buttton)
	#		output$Specific_BS_Junc <- strsplit(input$select_button, "_")[[1]]
			tmp <- strsplit(input$select_button, "@")
	#			debugme(tmp)
			cat("\nYou selected ", input$select_button, " and the junction is ", tmp[[1]][2])
			# It works, below is an example of the output:
			# button@chr1_33667149_chr1_33668857@Mon Dec 07 2015 17:33:42 GMT+1100 (AUS Eastern Daylight Time)
	})

	## Never got annything out of following. Leaving here just in case I want to look at it again at a later date
	#	eventReactive(input$select_button, {
		#	cat("\nPressed ", input$StrandFilter)
#				cat("\nPressed ", input$input_button)
#	})


	Selected_Junction <- reactive({
			if (is.null(input$select_button))
				return (NULL)
			tmp <- strsplit(input$select_button, "@")
			tmp[[1]][2]
	})


	output$Specific_BS_Junc <- renderText({

	  #			if (is.null(input$select_button))
	  #				return (NULL)
	  #			debugme(Selected_Junction(), input$select_button)
	  #			tmp <- strsplit(input$select_button, "@")
	  #			tmp[[1]][2]

	  #			Selected_Junction()
	  toDisplay <- paste("Displayed junction is: ",Ularcirc_data$Current_Selected_BS_Junction)
	  # Need to display following information
	  #     Canonical /backsplice junction
	  #     Typical Splice junction acceptor
	  #     Type I/II/III/IV reads supporting junction?
	  #     Sequence ove backsplice junction?

	  #temp$startDonor       # can use this to extract sequence
 	  #temp$startAcceptor     # use this to extract sequence
 	  #temp$strandDonor   # use this to extract sequence
 	  # JuncType    # use this to define junction

	  toDisplay

	})


	output$Specific_BS_Junc_details <- renderText({

	  BS_Junc_details(JuncCoords = Selected_Junction())
	})

	output$DisplayCanonical_sequence <- renderText({  #renderUI ({
	  # Obtain current canonical junction and obtain genomic sequence
	  if (is.null(Ularcirc_data$SelectedCanonical$Chr))
	  {  return (NULL)}
	  strandDonor <- "+"
	  if (Ularcirc_data$SelectedCanonical$Strand == 1)
	  {   strandDonor <- "-" }

	  Canonical_Junc_Entry <- data.frame(chromDonor=as.character(Ularcirc_data$SelectedCanonical$Chr),
	                                     startDonor=as.numeric(Ularcirc_data$SelectedCanonical$Start),
	                                     startAcceptor=as.numeric(Ularcirc_data$SelectedCanonical$End),
	                                     strandDonor=strandDonor,
	                                     type="c" )
    toDisplay <- Grab_BS_Junc_Sequence(Canonical_Junc_Entry, GeneList = GeneList())

	  HTML(toDisplay)
	  # Set flag so that can display graphic?
	})

	output$DisplayBS_sequence <- renderUI ({
	  UniqueJunctions <- Ularcirc_data$Current_Selected_BS_Junction_RAWData

	  toDisplay <- Grab_BS_Junc_Sequence(UniqueJunctions, GeneList = GeneList())
	  toDisplay <- HTML(paste('<p><strong>JUNCTION SEQUENCE: </strong></p>  <p> The full stop represents the junction point::
                    <pre   style=display: block;padding: 9.5px;margin: 0 0 10px;font-size: 13px;line-height: 1.42857143;
	                  color: #333;word-break: break-all;word-wrap: break-word;
	                  background-color: rgb(203, 219, 248);border: 1px solid #ccc;border-radius: 4px;>',toDisplay,'</pre></p>'))
	})

	output$Predicted_circRNA_Sequence <- renderUI ({
	  circRNA_Sequence <- Predicted_CircRNA_Sequence(circRNA_exons = Ularcirc_data$circRNA_exons, genelist = GeneList())
	  circRNA_Sequence <- HTML(paste('<p><strong>Predicted circRNA SEQUENCE: </strong></p>  <p>
                            <pre   style=display: block;padding: 9.5px;margin: 0 0 10px;font-size: 13px;line-height: 1.42857143;
	                          color: #333;word-break: break-all;word-wrap: break-word;
	                          background-color: rgb(203, 219, 248);border: 1px solid #ccc;border-radius: 4px;>',circRNA_Sequence,'</pre></p>'))
	})

	output$Predicted_Genomic_Junction_Sequence <- renderUI ({
	  FSJ_info <- Ularcirc_data$SelectedGenome_FSJ
	  strandDonor <- "+"
	  if (FSJ_info$Strand == 1)
	  {   strandDonor <- "-" }
	  Canonical_Junc_Entry <- data.frame(chromDonor=as.character(FSJ_info$Chr),
	                                     startDonor=as.numeric(FSJ_info$Start),
	                                     startAcceptor=as.numeric(FSJ_info$End),
	                                     strandDonor=strandDonor,
	                                     type="c" )
	  Genomic_FSJ_Sequence <- Grab_BS_Junc_Sequence(Canonical_Junc_Entry, GeneList = GeneList())

	  Genomic_FSJ_Sequence <- HTML(paste('<p><strong>Predicted Genomic FSJ SEQUENCE: </strong></p>  <p>
                            <pre   style=display: block;padding: 9.5px;margin: 0 0 10px;font-size: 13px;line-height: 1.42857143;
	                          color: #333;word-break: break-all;word-wrap: break-word;
	                          background-color: rgb(203, 219, 248);border: 1px solid #ccc;border-radius: 4px;>',Genomic_FSJ_Sequence,'</pre></p>'))
	})

	Fastq_Generate <- observeEvent(input$PE_Fastq_Request, #eventReactive(input$Update_Genome_Position,
	   {	# This prepares a list of genome coordinates as submitted by user
	      circRNA_Sequence <- Predicted_CircRNA_Sequence(circRNA_exons = Ularcirc_data$circRNA_exons, genelist = GeneList())
	      circ_length <- nchar(circRNA_Sequence)
	      circRNA_Sequence <- paste(circRNA_Sequence, circRNA_Sequence, sep = "")

	      FragmentLength <- input$FragSize #  300
	      ReadLength <- input$ReadLength   # 100
	      if (FragmentLength < ReadLength)
	      {  FragmentLength <- ReadLength + 1 }

	      if (circ_length > 350)
	      {   # First generate 3 x Type III reads. Back splice junction will be kept within 20nt from any read end
	          Read_One <- {};  Read_Two <- {};
	          typeII_III_IV_Offset <- c(280, 180, 80)  # These are the offset to generate the appropriate alignment types. Assuming fragment size of 300
	          typeII_III_IV_Offset <- c(FragmentLength-20, FragmentLength-ReadLength-80, ReadLength-20)
	          TypeIV_Label <- c("TypeIV_80","TypeIV_60", "TypeIV_40","TypeIV_20")
	          if (typeII_III_IV_Offset[2] < ReadLength)  # May get a situation where type IV reads cannot be made. Just set to zero.
	          {
	            typeII_III_IV_Offset[2] <- 0
	            TypeIV_Label <- paste(TypeIV_Label,"FAIL",sep="_")
	          }
	          for (j in 1:3)
  	        { start_pos <- circ_length - as.numeric(typeII_III_IV_Offset[j])
  	          for (i in 1:4)
  	          { if ((j == 2) && (  (typeII_III_IV_Offset[j] + 20 * (i-1)) > (FragmentLength-ReadLength) ))  # See if type IV reads turn into type III reads
  	            {   TypeIV_Label[i] <- paste(TypeIV_Label[i],"FAIL",sep="_") }

  	            Read_One <- c(Read_One, substr(x = circRNA_Sequence, start = start_pos, stop = start_pos + ReadLength))
  	            Read_Two <- c(Read_Two, substr(x = circRNA_Sequence, start = start_pos+FragmentLength-ReadLength, stop = start_pos + FragmentLength))
  	            start_pos <- start_pos + 20
  	          }
	          }
	      }
        names(Read_One) <- c("TypeIII_80","TypeIII_60", "TypeIII_40","TypeIII_20", TypeIV_Label,"TypeII_80","TypeII_60", "TypeII_40", "TypeII_20")
        names(Read_Two) <- c("TypeIII_80","TypeIII_60", "TypeIII_40","TypeIII_20", TypeIV_Label,"TypeII_80","TypeII_60", "TypeII_40", "TypeII_20")
        Read_One <- DNAStringSet(x=Read_One)
        Read_Two <- reverseComplement(DNAStringSet(x=Read_Two))
        browser()
        a<- "You now have the chance to save this output"
        a<- "quick before it is too late"
        # writeXStringSet( Read_One,"test_R1.fastq.gz",compress = TRUE, format="fastq")
        # writeXStringSet( Read_Two,"test_R2.fastq.gz",compress = TRUE, format="fastq")

	   })


	output$DisplayBS_sequence_details <- renderUI ({
	  if (is.null(Ularcirc_data$Current_Selected_BS_Junction))
	  {   return(NULL) }

	  withProgress(message="BE PATIENT: Assembling backsplice junction information", value=0, {
	    UniqueJunctions <- Ularcirc_data$Current_Selected_BS_Junction_RAWData
  	  toDisplay <- "No entries, check input junction"
  	  if (nrow(UniqueJunctions) > 1)
  	  {
  	    GeneName <- Annotate_BS_Junc(DataSet=UniqueJunctions, GeneList = GeneList(), MaxDisplay = 1, input$LibraryStrandType)
  	    GeneName <- GeneName$Gene[1]          # Grab gene name
  	    JuncType <- 'Unknown'
  	    if (UniqueJunctions$type[1] == 'c')    { JuncType <- 'Canonical'}
  	    if (UniqueJunctions$type[1] == 'bs')   { JuncType <- 'Backsplice'}

  	    incProgress(1/3, detail = paste("Extracting genomic sequence"))
  	    idx <- which( Ularcirc_data$ProjectData$SampleIDs == input$SelectedFiles) #m379()$SampleIDs == input$SelectedFiles)

  	    incProgress(1/3, detail = paste("Extracting annotation data"))
  	    # Calling Prepare_Gene_Object will get all data in correct format for Gene_Transcript_Features
  	    BS_Junctions <- Ularcirc_data$ProjectData$Junctions
  	    ## Thinking about load FAD score so can plot data?
  	    # FAD <- Fragment_Alignment_Distribution(All_junctions[All_junctions$BSjuncName == BS_Junc_ID,])  Ularcirc_data$ProjectData$Junctions

  	    PGO<-Prepare_Gene_Object(GeneName, BS_Junctions = BS_Junctions,GeneList= GeneList(), File_idx = idx,
  	                             Canonical_Junctions = Ularcirc_data$ProjectData$Canonical_AllData)

  	    incProgress(1/3, detail = paste("Extracting junction data on parental gene"))
  	    GeneFeatures <- Gene_Transcript_Features(GeneList=GeneList(), Gene_Symbol=GeneName, GeneObject=PGO)
  	    BS_Junc_idx <- {}

  	    if (length(PGO$Junctions) == 4)  # Should be a list with four dataframe entries
  	    {  BS_Junc_idx <- which(PGO$Junctions$uniques.bed$name==Ularcirc_data$Current_Selected_BS_Junction) }
  	    if ( (length(PGO$Junctions) != 4) || (length(BS_Junc_idx) == 0))   # No BSJ recovered. Attempt a direct retrieval
  	    {  PGO$Junctions <- BS_Junctions[BSjuncName == Ularcirc_data$Current_Selected_BS_Junction,]
  	    #& DataSet == idx
  	      ## Need to flag that I have pooled counts for ALL samples.
  	    }


  	    BS_Junc_Count <- PGO$Junctions$uniques.bed$score[BS_Junc_idx]
  	    BS_junc_proportion <- BS_Junc_Count/GeneFeatures$Av_junc_count

  	    Abundant_BS_alert <- c('')

  	    ### Identify which exons are included within circRNA
  	    Exons_tmp <- {}

  	    if ( ((PGO$Transcript$strand[1] == 1) && (input$LibraryStrandType == "Opposing strand"))  # Positive strand.    opposing strand library   This is Illumina Tru seq
  	         || ((PGO$Transcript$strand[1] == 0) && (input$LibraryStrandType == "Same Strand"))  )  # Negative strand same strand library.
  	    { # negative strand genes
  	      Exons_tmp <- PGO$Transcript[stop > as.numeric(UniqueJunctions[1,.(startDonor)])  & start < as.numeric(UniqueJunctions[1,.(startAcceptor)]),]
  	      true_candidates_stop <- which(abs(Exons_tmp$start - as.numeric(UniqueJunctions[1,.(startDonor)])) == 1)
  	      true_candidates_start <- which(abs(Exons_tmp$stop - as.numeric(UniqueJunctions[1,.(startAcceptor)])) == 1)
  	    }
  	    else if (input$LibraryStrandType != "Unstranded")   # positive strand genes
  	    { Exons_tmp <- PGO$Transcript[stop < as.numeric(UniqueJunctions[1,.(startDonor)])  & start > as.numeric(UniqueJunctions[1,.(startAcceptor)]),]
    	    true_candidates_stop <- which(abs(Exons_tmp$stop - as.numeric(UniqueJunctions[1,.(startDonor)])) == 1)
    	    true_candidates_start <- which(abs(Exons_tmp$start - as.numeric(UniqueJunctions[1,.(startAcceptor)])) == 1)
  	    }
  	    else
  	    {
  	     # browser()
  	      print("\nHave not implemented unstranded data yet.....")
  	    }

        ## Occasionally may have a transcript that has multiple exons within circRNA boundaries but no exon that finishes at BSJ. We want to get rid of these.
  	    ## Working example is PTK2, transcript "7894".

  	    true_candidate_IDs <- intersect(Exons_tmp$gene[true_candidates_start], Exons_tmp$gene[true_candidates_stop])
  	    tmp_idx <- Exons_tmp$gene %in% true_candidate_IDs
  	    Ularcirc_data$circRNA_exons <- Exons_tmp[tmp_idx,]


  	    ### Now to work out maximum length by sifting through tx entries and adding up exon lengths
  	    circRNA_exon_lengths <-by(Ularcirc_data$circRNA_exons,Ularcirc_data$circRNA_exons$gene,identity )  # This makes a list of all transcripts
  	    circRNA_Size <- max(range(lapply(X = circRNA_exon_lengths,FUN = function(x) { return(abs(sum(x$start-x$stop)))}))) +1

  	    if ((! is.na(GeneFeatures$Av_junc_count)) && (length(BS_junc_proportion) > 0))
  	    { if (BS_junc_proportion > 0.1)
  	       Abundant_BS_alert <- c(' (ABUNDANT circRNA relative to parental transcript)')
  	    }
        toDisplay <- HTML(paste('<p><strong>JUNCTION: </strong>',Ularcirc_data$Current_Selected_BS_Junction,'</p>',
                  '<p><strong>JUNCTION TYPE: </strong>',JuncType,'</p>',

                  '<p ><pre  style=display: block;padding: 9.5px;margin: 0 0 10px;font-size: 13px;line-height: 1.42857143;
                            color: #333;word-break: break-all;word-wrap: break-word;
                            background-color: rgb(203, 219, 248);border: 1px solid #ccc;border-radius: 4px;>
                            <strong>PARENTAL GENE: </strong>',GeneName,'<br>',
                            '<strong>Number of transcript isoforms:</strong>',length(GeneFeatures$Num_Exons_Per_Transcript),'<br>',
                            '<strong>Avg linear junction count:</strong>',round(GeneFeatures$Av_junc_count,2),'</pre></p>',

                  '<p><strong>Backsplice junction count:</strong>',BS_Junc_Count,' <span style="color:red;"> ', Abundant_BS_alert,'</span></p>',
                  '<p><strong>Type II vs Type III reads:</strong>',round(Ularcirc_data$BSJ_TypeI_vs_TypeII_Ratio,1),'</p>',
                  '<p><strong>CircRNA length:</strong>',circRNA_Size, '</p>'))
  	  }
	  })  # withProgress
	  HTML(toDisplay)
	})

	output$Plot_RAD_Histogram <- renderPlot({
	  layout(matrix(c(#1,1,1,
	    2,2,2,           # canonical junctions
	    1,1,1           #  Transcripts
	  ), 2, 3, byrow = TRUE))      # 2 rows, three columns. First seg  | Second Seg  | total length
	  par(mar=c(3,2,1,1))
#browser()
	  RAD_data_set <- Fragment_Alignment_Distribution(Ularcirc_data$Current_Selected_BS_Junction_RAWData)
	  a <- 1
	  a<- 1
	  hist(RAD_data_set$FirstSeg_position,breaks = 50,main = "First Segment distribution")
	  hist(RAD_data_set$SecondSeg_position,breaks = 50, main="Second segment distribution")
	  # Would like a third histogram to show length distribution

	})

	output$miRNA_Options <- renderUI({
	  w <- paste(selectizeInput("miRNA_Seed_Length",label="miRNA seed length",choices = 6:15, selected=7))
	  w<- paste(w, selectizeInput("miRNA_Multiples",label="minimum number of miRNA multiples",choices = 1:15, selected=2))
	  HTML(w)
	})

	update_miRNA_Sites <- reactive({
	  if (is.null(input$miRNA_Seed_Length))
	  {  return(NULL) }
	  circRNA_Sequence <- Predicted_CircRNA_Sequence(circRNA_exons = Ularcirc_data$circRNA_exons, genelist = GeneList())
	  Binding_Sites <- miR_binding_site_Analysis(circRNA_Sequence, "hsa", seed_length=as.numeric(input$miRNA_Seed_Length))   # SeedMatchResult=SeedMatchResult, Total_miR_number=length(allseeds),  miR_Seed_Lookup
    return(list(circRNA_Sequence=circRNA_Sequence, Binding_Sites=Binding_Sites)	)
	})

	output$circRNA_Sequence_Analysis_miRNA <- renderPlot({
	  par(mar=c(2, 2, 2, 2));
	  plot(c(1,900), c(1,900), type="n", axes=FALSE, xlab="", ylab="", main="");
	  draw.arc(xc=400,yc=400,r=400,w1=270,w2=630,col="blue", lwd=6)    # Draws circRNA

	      BS <- update_miRNA_Sites()
	      Binding_Sites <- BS$Binding_Sites
	      circRNA_Sequence <- BS$circRNA_Sequence

	      if(is.null(Binding_Sites))
	      {   Explanation <- paste("0 miRNA sites (seed = ",input$miRNA_Seed_Length,")")
	          text(400,400,Explanation)
	          return()
	      }
	      Binding_Sites$SeedMatchResult <- lapply(Binding_Sites$SeedMatchResult, FUN= function(x) { if (length(x) >= as.numeric(input$miRNA_Multiples)) return(x) })

    	  hist_data <- unlist(lapply(X = Binding_Sites$SeedMatchResult,FUN = function(x) {as.numeric(x)}))
    	  Ularcirc_data$selected_circRNA_stats$miRNA_BS_Sites <- hist_data
    	  if (length(which(hist_data > 0)))
    	  {
    	    miR_match_idx <- which(hist_data > 0)
    	    hist_data <- table(sapply(hist_data[miR_match_idx],length))
    	    unmatched_miR <- Binding_Sites$Total_miR_number - sum(hist_data)
    	  }
    	  else
    	  { hist_data <- c("0"=length(hist_data))  }
    	  # barplot(hist_data, xlab="# binding sites", ylab="Numer of miRNA", main = paste(Ularcirc_data$SelectedGene_from_BSJ_Table, " circRNA", sep=""))

    	  Binding_Report <- unlist(Binding_Sites$SeedMatchResult)
    	  Binding_Report <- Binding_Report[Binding_Report > 0]
    	  if (length(Binding_Report) == 0)
    	  {
    	    Explanation <- paste("0 miRNA sites (multiple = ",input$miRNA_Multiples,")")
    	    text(400,400,Explanation)
    	    return()
    	  }

    	  names(Binding_Report) <- substr(names(Binding_Report),1,as.numeric(input$miRNA_Seed_Length))

    	  # Extract IDs
    	  miR_ID_idx <- list()
    	  miR_IDs <- list()
    	  for(i in 1:length(Binding_Report))
    	  {   miR_ID_idx[[i]] <- which(Binding_Sites$miR_Seed_Lookup == names(Binding_Report)[i])
    	       miR_IDs[[i]] <-  row.names(Binding_Sites$miR_Seed_Lookup)[miR_ID_idx[[i]]]
    	  }


    	  # Calculate positions to draw in circle
    	  miR_positions <- round(Binding_Report/nchar(as.character(circRNA_Sequence)) *360,0)  # Coordinates
    	  miR_length <- round(25/nchar(as.character(circRNA_Sequence)) *360,0)

##	  browser()
    	  for (i in 1:length(miR_positions))
    	  {  if (miR_positions[i] > -1)
    	       draw.arc(xc=400,yc=400,r=370,w1= miR_positions[i], w2 = miR_positions[i] + miR_length, col = "black",lwd=4, draw_tick=FALSE, Internal_Labels=miR_IDs[[i]][1])	  # Draws miRNA binding site

    	  }
    	  draw.text.w(xc=400,yc=400,r=430, w=270, n=paste(Ularcirc_data$SelectedGene_from_BSJ_Table, " circRNA", sep=""), col = "blue")

	})

	output$circRNA_Sequence_Analysis_ORF <- renderPlot({
  	 if (input$circRNA_Sequence_Analysis == "Open reading frame analysis")
	  { # browser()
  	  circRNA_Sequence <- Predicted_CircRNA_Sequence(circRNA_exons = Ularcirc_data$circRNA_exons, genelist = GeneList())
	    tmp <- as.character(circRNA_Sequence)
	    ORF_lengths <- ''
	    for (i in 1:3)
	    {
	       concat_circRNA <- DNAString(paste(substr(tmp,i,nchar(tmp)),tmp,sep=""))
	       concat_circRNA <- translate(concat_circRNA)
	       ORFs <- strapplyc(as.character(concat_circRNA), pattern="M.*?\\*"  )
	       ORF_lengths <- c(ORF_lengths,unlist(lapply(ORFs[[1]],nchar)))
	    }
	    ORF_lengths <- as.numeric(ORF_lengths)
	    ORF_lengths <- ORF_lengths[! is.na(ORF_lengths)]
	    hist_data <- data.frame("50-100nt"= length(which((ORF_lengths > 50) & (ORF_lengths < 100))),
	                            "100-200nt" = length(which((ORF_lengths >= 100) & (ORF_lengths <= 200))),
	                            ">200nt" = length(which(ORF_lengths > 200 )) )

	    # Set up parameters for plotting ORFs
	    w1 <- 270;	# TO DO: change this to reflect correct starting position
	    w2 <- w1 +  ( max(ORF_lengths)* 3/ nchar(as.character(circRNA_Sequence)) * 360)
	    r <- 400
	 #   browser()
	    par(mar=c(2, 2, 2, 2));
	    plot(c(1,900), c(1,900), type="n", axes=FALSE, xlab="", ylab="", main="");
	    draw.arc(xc=400,yc=400,r=400,w1=270,w2=630,col="blue", lwd=6)    # Draws circRNA
      draw.arc(xc=400,yc=400,r=370,w1=270, w2 =  w2 ,col = "black",lwd=4, draw_tick=TRUE)	  # Draws largest ORF
      draw.text.w(xc=400,yc=400,r=430, w=270, n=paste(Ularcirc_data$SelectedGene_from_BSJ_Table, " circRNA", sep=""), col = "blue")
      text(350,400, labels=paste("Largest ORF is ", max(ORF_lengths)," amino acids"))
	    #barplot(as.matrix(hist_data), xlab="ORF length", ylab="Number", main= paste(Ularcirc_data$SelectedGene_from_BSJ_Table, " circRNA", sep=""))
	  }

	  # Need to extract circRNA sequence
	  # Then perform either miRNA or ORF analysis and prepare circos plot
	})

	output$circRNA_Sequence_Analysis_Table <- renderDataTable({
#	  browser()

	  a <- Ularcirc_data$selected_circRNA_stats$miRNA_BS_Sites

	  return(NULL)

	})

	Genome_Coords <- observeEvent(input$Update_Genome_Position, #eventReactive(input$Update_Genome_Position,
	  {	# This prepares a list of genome coordinates as submitted by user
	      Ularcirc_data$Genome_Coordinates <- list(chrom=input$GenomeChrom_Input, chromstart=input$GenomeStart_Input,
	          chromend=input$GenomeEnd_Input, chromstrand=input$GenomeStrand_Input)
	})


	Previous_m379 <- observeEvent(input$LoadProjectRequest,
	 { LoadStatus <- c("No project to load")

  	 withProgress(message="Loading projet data. Be patient", value=0, {
	     incProgress(1/2, detail = paste("Warning:: status bar cannot increment"))
  	    if (exists("ProjectGroupings"))   # Remove existing project Grouping before loading in new data
  	    { remove(ProjectGroupings) }
  	    if (exists("meta_data"))
  	    { remove(meta_data)}

    	#  extdata_path <- DataPath()  # function from Global.R
	      ProjectFileName <- paste(extdata_path,"/",input$LoadExistingProject, ".RData", sep="")
	      load( file= ProjectFileName)                  # This will load object called DataSet
	      incProgress(1/2, detail = paste("Done !! "))
	      Ularcirc_data$ProjectData <- DataSet
	      if (exists("ProjectGroupings"))
	      {  # browser()
	         Groupings <- ProjectGroupings
	      }

	      Ularcirc_data$ProjectData$ProjectFileName <- ProjectFileName         # Record the filename of the project in case need to update and re-save

	      # Display notes on current project that will help remind them of the setting to use with this data set
	      if (exists("meta_data"))
	      { blurb <- paste("<b>Project name:</b> ", ProjectFileName, "<br/>", sep="")
  	      blurb <- paste(blurb, "<b>TxDb:</b> ", meta_data$ProjectSpecies , "<br/>", sep="")
	        blurb <- paste(blurb, "<b>Library type:</b> ",meta_data$LibraryStrandType, "<br/>", sep="")
	        ProjNotes <- gsub(pattern = "\n",replacement = "<br/>",x = meta_data$ProjectNotes)
	        blurb <- HTML(paste(blurb, "<b>Additional notes:</b><br/> ",ProjNotes, "<br/>", sep=""))
	        showModal(modalDialog(title="Associated meta data",blurb,easyClose=TRUE,footer=NULL))
	      }

  	 }) # withProgress
	})



	ntext <- eventReactive(input$SaveProjectRequest, {

	  # Need to do following:
	  #   - to confirm overwriting of existing project
	  #   - amount of file space available/used
	  #   - Update the "load" list to see project just saved
	  SaveStatus <- c("No file name provided cannot save")
	  if (length(input$NewProject_Filename) > 0)
	  {
	 #   extdata_path <- as.character(DataPath())    # DataPath function is in Global.R
	    ProjectFileName <- paste(extdata_path,"/",input$NewProject_Filename, ".RData", sep="")
	    if (file.exists(ProjectFileName))
	    {   SaveStatus <- paste("Project already exists, please modify filename")}
	    else
	    {   withProgress(message="Saving project data, please wait..", value=0, {
	          SaveStatus <- paste("Project saved as", ProjectFileName, sep= " ")
	          DataSet <- Ularcirc_data$ProjectData
	          meta_data <- list()
	          ProjectGroupings <- Groupings
	          meta_data$ProjectNotes <- input$ProjectNotes
	          meta_data$ProjectSpecies <- input$Species_Genome    # This will provide species name and genome release
	          meta_data$LibraryStrandType <- input$LibraryStrandType
    	      save(DataSet, ProjectGroupings, meta_data, file= ProjectFileName)
	        })
  	  }
	  }
    SaveStatus

    # Below are the objects returned from LoadJunctionData  which is what m379 calls.
    #  return(list(Junctions=AllData, SummarisedData=SummarisedData,
    #              Original_Junction_Numbers=total_rows, # the number of rows in the last file read. This is pointless if multiple files are read in.
    #             Original_n_unique_junctions=n_UniqueJunctions,  # the number of columns in the last file read. This is pointless if multiple files are read in.
    #            Original_Postfilt_n_unique_junctions=n_UniqueJunctions_PostFilt,   # This is identical to n_UniqueJunctions !! probably can delete.
    #           DataType = DataType                             # "Backsplice"  or "Canonical" so far


	})

	output$Save_Load_Status <- renderText({ntext()})


  output$JunctionTableOther <- renderDataTable({
    if (! is.null(Ularcirc_data$Current_Selected_BS_Junction_RAWData))
    {
      dt <- data.table(Ularcirc_data$Current_Selected_BS_Junction_RAWData)
      datatable(dt, selection='none', options = list(lengthMenu = c(5, 10, 50), pageLength = 5))
    }
  })

  output$circRNA_Read_Distribution <- renderPlot(
    {  # This produces a static plot of the expected TypeI read to other read types for circRNA
      readLength <- 100
      Fragment_Length <- 300
      circRNA_size_Range <- 300:3000
      plot(x = circRNA_size_Range, y=(circRNA_size_Range-Fragment_Length)/circRNA_size_Range*100, xlab="circRNA size", ylab="Proportion of Type I reads (%)", ylim=c(1,100))

    })

  output$Predicted_Read_Distribution <- renderDataTable({

    input$FragmentSize
    input$ReadLength
    input$CircRNA_Size

    TypeI_proportion <- (input$CircRNA_Size - input$FragmentSize)/input$CircRNA_Size

    # Scale BSJ input accordingly
    ReadNumber <- (input$ReadNumber/(1-TypeI_proportion))

    TypeI <- round(TypeI_proportion * ReadNumber, 0)

    if (input$ReadType == "Paired")
    {    TypeIV <- (input$FragmentSize - (input$ReadLength * 2))/ input$FragmentSize  }
    else
    {    TypeIV <- (input$FragmentSize - input$ReadLength) /input$FragmentSize  }
    if (TypeIV < 0) # Read length is greater than fragment length. Therefore entire fragment is covered.
    {    TypeIV <- 0 }

    Remaining_Reads <- ReadNumber - TypeI
    # Now work out where these reads would be discovered
    TypeIV <- round(TypeIV * Remaining_Reads, 0)

    Remaining_Reads <- Remaining_Reads - TypeIV
    TypeII <- round(Remaining_Reads/2,0)
    TypeIII <- TypeII

#    toDisplay <- HTML(paste('<p style=color:red><strong>TypeI'



    dt<-data.table(TypeI, TypeII, TypeIII, TypeIV)
    colnames(dt)[4] <- c('TypeIV*')

    row.names(dt) <- c("Number of reads")
    datatable(dt, TypeIV,selection='none', options = list(lengthMenu = c(1), pageLength = 1, dom = 't')) %>%
      formatStyle('TypeI',  color = 'red', backgroundColor = 'orange', fontWeight = 'bold') %>%
      formatStyle('TypeIV*',  color = 'red', backgroundColor = 'orange', fontWeight = 'bold')


  })


  #renderTable({
  #  head(datasetInput(), n = input$obs)
#	if (length(circRNA_Subset()[[2]]) == 4)
#	{	 head(circRNA_Subset()[[2]][[1]], n = input$obs)		} })

	# Set what junctions user wishes to view
	JunctionOption <- reactive ({
		JO <- 0 # Default is to view ALL junctions
		if (input$JunctionType == "Backsplice")
		{	JO <- 1	}
		if (input$JunctionType == "Alternative Canonical")
		{	JO <- 500	}
		JO
	})


    output$Gene_transcript_plot <- renderPlot({
      GeneObject <- circRNA_Subset()  # Will need to separate gene transcript table from data
      if (length(GeneObject) > 0)
      {

        chrom <- as.character(GeneObject$Transcript[1,.(chrom)])
        chromstart = as.numeric(GeneObject$Transcript[1,.(start)])-15
        chromend = as.numeric(GeneObject$Transcript[nrow(GeneObject$Transcript),.(stop)])+15
        plotGenes(GeneObject$Transcript, chrom, chromstart, chromend, maxrows=50,height=0.4,plotgenetype="box")
      }
    })

  	output$distPlot <- renderPlot ({    ## This will draw circRNA and linear junction abundance graphs.
  	  if ((is.null(Ularcirc_data$Current_SelectedGene)) && (is.null(input$DisplayJunctionCountTable_row_last_clicked)))
  	  { return(NULL)}

  	  if (! Ularcirc_data$GenePanelLoaded)
  	  { return(NULL) }

		  cat(paste("\nStarting to renter transcript plot", date()))
  		layout(matrix(c(#1,1,1,
						4,4,4,
						3,3,3,           ## backsplice junctions
						2,2,2,           # canonical junctions
						1,1,1           #  Transcripts
	#					5,6,7), 5, 3, byrow = TRUE))     # 4 rows of figures, first two rows contain one major image, while third and final row contains three items (will fill two using zoom feature of sushi)
           	), 4, 3, byrow = TRUE))      # 4 rows, three columns
	  	par(mar=c(3,4,1,1))


  		Zoom_coords <- View_Gene_At_Coords()
  		circs <- circRNA_Subset()
	  	DTE <- Draw_Transcript_Exons(circs, JunctionOption(), Zoom_coords, Junction_highlight=list(BSjunc=Ularcirc_data$Current_Selected_BS_Junction, Canonical = Ularcirc_data$SelectedCanonical))
		  cat(paste("\n-Completed rendering transcript plot", date()))
		  DTE
	  })

  	output$genomePlot <- renderPlot ({    ## This will draw circRNA and linear junction abundance graphs on user defined genome coordinates

### Ideas:  GENERATE a list of genomic coorindinates that contain significant FSJ reads.
  	  ### Convert TxDb as follows:      g<- genes(GeneList()$transcript_reference)
  	  ### Invert                       inverse <- gaps(g)
  	  ### Convert raw FSJ to Granges object       FSJ <- makeGrangesFromDataFrame( Ularcirc_data$ProjectData$Canonical_AllData )
  	  ### Intersect raw FSJ with inverse TxDb object      intersect(FSJ, inverse)
  	  ### Display any regions that have counts > some threshold
  	  ###

  	  g<- genes(GeneList()$transcript_reference)
  	  inverse <- gaps(g)
  	  FSJ <- makeGRangesFromDataFrame( Ularcirc_data$ProjectData$Canonical_AllData,seqnames.field = "chrom1",start.field = "start1", end.field = "end2", keep.extra.columns = TRUE )

  	  if (is.null(Ularcirc_data$Genome_Coordinates$chrom))
  	  { return(NULL)  }

  	  idx <- seq(from=1, to = length(Ularcirc_data$ProjectData$SampleIDs)) # m379()$SampleIDs))     # Default is to select everything. PROBABLY BETTER TO RETURN A NULL IF NO DATA SELECTED AND DEAL WITH THIS ELSEWHERE??
  	  if (! is.null(input$SelectedFiles))
  	  { idx <- which( Ularcirc_data$ProjectData$SampleIDs == input$SelectedFiles) }

  	  cat(paste("\nStarting to render GENOME plot", date()))

  	  layout(matrix(c(#1,1,1,
  	    4,4,4,
  	    3,3,3,           ## backsplice junctions
  	    2,2,2,           # canonical junctions
  	    1,1,1           #  Transcripts
  	  ), 4, 3, byrow = TRUE))      # 4 rows, three columns
  	  par(mar=c(3,4,1,1))
 # 	  Zoom_coords <- View_Gene_At_Coords()

  	  strand <- 1  # +ve strand
  	  if (Ularcirc_data$Genome_Coordinates$chromstrand == "-")
  	    strand <- 2

  	  BS_Junctions <- Ularcirc_data$ProjectData$Junctions
  	  BS_Junctions <- circJunctions(BS_Junctions ,Ularcirc_data$Genome_Coordinates$chrom,
  	                             Ularcirc_data$Genome_Coordinates$chromstart, Ularcirc_data$Genome_Coordinates$chromend,fileID=idx)

  	  Canonical_Junctions <- Ularcirc_data$ProjectData$Canonical_AllData
  	  Canonical_Junctions = Canonical_Junctions[chrom1 == Ularcirc_data$Genome_Coordinates$chrom &
  	                                              start1 > as.numeric(Ularcirc_data$Genome_Coordinates$chromstart) &
  	                                              end2 < as.numeric(Ularcirc_data$Genome_Coordinates$chromend) &
  	                                              strand1 == strand & DataSet==idx,]
  	  Ularcirc_data$GenomeCanonicalJunctionCountTable <- Canonical_Junctions    # Save data in case want to display

  	  if (strand == 2) { strand <- -1 }   # Re-format strand for Transcript object
#  	  Transcript <- data.table(b[,c('chrom','start','stop','gene','score','strand','type')])
  	  Transcript <- data.table(chrom=Ularcirc_data$Genome_Coordinates$chrom ,
  	                           start=Ularcirc_data$Genome_Coordinates$chromstart ,
  	                           stop=Ularcirc_data$Genome_Coordinates$chromend,gene ="" ,score=0,
  	                           strand=strand , type='exon')

  	  GeneObject <- list(Transcript= Transcript, Junctions = BS_Junctions, Transcript_Canonical_juncs = Canonical_Junctions)


  	  DTE <- Draw_Transcript_Exons(GeneObject, JunctionOption(), Zoom_coords=NULL, GenomeDisplay=TRUE, Genome_Coords=Ularcirc_data$Genome_Coordinates)
  	  cat(paste("\n-Completed rendering genome plot", date()))
  	  DTE
  	})


  	output$Display_Gene_Zoom_Coords <- renderUI ({
  	  # Draw slider input to define regions to display
  	  Gene_Transcripts = circRNA_Subset()
  	  Gene_min <- min(Gene_Transcripts$Transcript$start)
  	  Gene_max <- max(Gene_Transcripts$Transcript$stop)
  	  if (is.null(Gene_Transcripts))
  	    return(NULL)
  	  else
  	  { 	HTML(paste(sliderInput("Gene_Zoom", "Define region within gene to view:",min = Gene_min-zoomOffset, max = Gene_max+zoomOffset, value = c(Gene_min-15,Gene_max+15)),
  	                actionButton("Navigate_Around_Gene","Navigate")) )
  	   }
  	})


  	output$DisplayGroupings <- renderUI ({   # This displays group IDs and asigned sample IDs
  	  if (is.null(input$Number_BiologicalSamples))
  	  { return() }

  	  SampleIDs <- {}
      for (i in 1:input$Number_BiologicalSamples)
      { lookupID <- paste("Group",i,sep="_")
  	    SampleIDs <- c(SampleIDs, input[[lookupID]])
      }
  	  w <- ""
  	  for(i in 1:input$Number_BiologicalSamples)
  	  { #groupIDs <- paste("Group",SampleIDs[i])
  	    w <-  paste(w,column(3,checkboxGroupInput('SamplesAssignedToGroup', SampleIDs[i], CreateBiologicalGroupings()[[i]])))
  	  }
  	  w <- paste(w,selectizeInput("AssignedGroupID",label="Assign checked samples to group:",choices = SampleIDs))
  	  w <- paste(w, actionButton(inputId = "AssignToGroup",label = "Assign"))
  	  HTML(w)
  	})


  	output$DisplayGroupNumber <- renderUI ({   # This prepares group IDs for sidebar, accessed in UI.R
 ## This is currently not in use. Was thinking of using this to replace line in UI.R so that can  have more control when loading data.
  	  w <-	sliderInput("Number_BiologicalSamples", "Number of biological treatments in data set:",min = 1, max = 10, value = 1)
      HTML(w)
  	})

  	output$DisplayGroupNames <- renderUI ({   # This prepares group IDs for sidebar, accessed in UI.R
      if (is.null(input$Number_BiologicalSamples))
      { return() }

  	  SampleIDs <- seq(from=1, to=input$Number_BiologicalSamples, by=1)
  	  w <- ""
  	  for(i in 1:input$Number_BiologicalSamples)
  	  { #inputId, label, value
  	    groupIDs <- paste("Group",SampleIDs[i],sep = "_")
  	    w <- paste(w,textInput(inputId=groupIDs, label=groupIDs, value=groupIDs))
  	  }
  	  HTML(w)
  	})


  	output$Display_Species_Options <- renderUI ({   # This prepares species annotation and genome selection options, accessed in UI.R
  	  if (is.null(input$Annotation_lib))
  	  { return() }
  	  # First check not running in un-annotated mode
  	  if (input$Annotation_lib == "NO_ANNOTATION")
  	  {
  	    w <- paste(selectizeInput("Species_Genome",label="Genome release",choices = "No_Annotation", selected="No_Annotation"))
  	    w <- paste(p('Operating Ularcirc without annotation. To change this ensure you have downloaded the appropriate databases and selected them in above pulldown menu'))
  	  }
  	  else
  	  { # Need to keep species initials: eg org.Hs.eg.db  to Hs
  	    withProgress(message="Searching and populating installed compatible annotation libraries.", value=0, {

  	    Species_Initials <- gsub(pattern = ".eg.db", replacement = "" ,x = gsub(pattern = "org.", replacement = "", x = input$Annotation_lib))
  	    Species_Genome <- List_Species_Files$GenomeOptions[grep(pattern = Species_Initials, x = List_Species_Files$GenomeOptions)]
        w <- ""
  	    w <- paste(w,selectizeInput("Species_Genome",label="Genome release",choices = Species_Genome, selected=Species_Genome[1]))
  	    }) # withProgress
  	  }
      HTML(w)
  	})

  	output$Display_Txdb_Options <- renderUI ({   # This prepares transcriptome options which are dependent on Display_Species_Options, accessed in UI.R
  	  if (is.null(input$Species_Genome) || (input$Annotation_lib == "NO_ANNOTATION"))
  	  { return() }
    	# Need to keep species Genome information: eg HSapiens.UCSC.Hg38  to Hg38
  	  Species_Genome <- strsplit(input$Species_Genome, split = "\\.")[[1]]
	    Species_Genome <- Species_Genome[length(Species_Genome)]
	    Species_Transcripts <- List_Species_Files$TxDbOptions[grep(pattern = Species_Genome, x = names(List_Species_Files$TxDbOptions))]
	    w <- paste(selectizeInput("TxDb",label="Transcriptome database",choices = Species_Transcripts, selected=Species_Transcripts[1]))
	    w <- paste(w, actionButton("LoadTxDb","LOAD TRANSCRIPT DATABASE"))
  	  HTML(w)
  	})


  	Groupings <- list()

  	CreateBiologicalGroupings <- reactive( {   ## This code will assemble list of samples
  	  ## Currently struggling to see how to not overwrite existing list

  	  inFile = Ularcirc_data$ProjectData$SampleIDs #m379()$SampleIDs
  	  Existing_Group_number <- length(Groupings)
  	  # Following if statements look after group size changes
  	  if ((is.null(input$Number_BiologicalSamples)) || (is.null(inFile)))
  	  { return() }
  	  if (Existing_Group_number > input$Number_BiologicalSamples)  # Make sure any entries in groups that are to be trimmed are assigned elsewhere
  	  {
  	    for(i in Existing_Group_number:(input$Number_BiologicalSamples+1))
  	    { Groupings[[i]] <<- c(Groupings[[1]], Groupings[[i]])  # Copy samples from list that is about to be deleted
  	      Groupings[i] <<- NULL
  	    }
  	  }
  	  if (Existing_Group_number < input$Number_BiologicalSamples)  # Add extra group(s)
  	  {
  	    SampleIDs <- {} # input$groupIDs[Existing_Group_number:input$Number_BiologicalSamples]
  	    if (Existing_Group_number==0)
  	    {  Existing_Group_number <- 1  }
  	    GroupLookupIDs <- paste("Group",seq(from=Existing_Group_number, to=input$Number_BiologicalSamples, by=1),sep="_") # Create unique numeric group IDs
  	    input_names <- names(input)[which(names(input)== GroupLookupIDs)] # Indexes to group


  	   if (length(input_names) == 0)     # This will catch if group names are not assigned in "input" yet.
  	   { SampleIDs <- GroupLookupIDs  }
  	   else
  	   {   for (i in 1:length(input_names))
  	       { SampleIDs<- c(SampleIDs, input[[ input_names[i] ]])     }
  	    }
  	    Groupings <<- c(Groupings, vector("list",(input$Number_BiologicalSamples-Existing_Group_number)))  # Add some blank list elements

  	    if (length(SampleIDs) == length(Groupings))
  	    {        names(Groupings) <<- SampleIDs    }

  	  }

  	  # Following if statements look after re-assignment of samples to groups
  	  if (input$Number_BiologicalSamples == 1)
  	  {   Groupings[[1]] <<-  inFile
  	      names(Groupings) <<- input$Group_1
  	  } # Assign all samples to group 1 as only 1 group
  	  else if ( input$AssignToGroup )      # If Assign to group action Button has been pressed start reassignment
  	  {  Samples_to_Reassign <- input$SamplesAssignedToGroup
  	     GroupID <- input$AssignedGroupID
  	     AllGroupIDs <- input$groupIDs
  	     if (length(AllGroupIDs) == 0)
  	     {    AllGroupIDs <- paste("Group",seq(from=1, to=input$Number_BiologicalSamples, by=1), sep="_")   }
  	     GroupIdx <- which(GroupID == AllGroupIDs)
         RemoveEntry<- function(X, SampleIDs)   	     # This function will be used to remove entries that are about to be reassigned
         {
            if (length(SampleIDs) > 0)
            {  for(i in 1:length(SampleIDs))
               {  SampleIdx <- which(X == SampleIDs[i])
                  if (length(SampleIdx) > 0)
                  {   X <- X[SampleIdx* -1]  }
               }
            }
           return(X)
         }

         Groupings <<- lapply(Groupings, FUN=RemoveEntry, Samples_to_Reassign)
         Groupings[[GroupIdx]] <<- c(Groupings[[GroupIdx]],Samples_to_Reassign)


  	  }

  	  if (length(Groupings) == input$Number_BiologicalSamples)
  	  {
  	     if (! is.null(input$groupIDs))
  	     {   names(Groupings) <<- input$groupIDs   }
  	     else
  	     {  # Following if statement is a hack to correct "Groupings" list names. Have yet to identify how list names are not being generated.
  	       if ((length(which(Groupings=="")) > 0) && (length(Groupings) == input$Number_BiologicalSamples))
  	         names(Groupings) <- paste("Group",seq(from=1, to=input$Number_BiologicalSamples, by=1), sep="_")
  	     }
    	   return(Groupings)   # The list that assigns all samples to a group ID
  	  }
  	  else
  	     return (NULL)
  	})

  	View_Gene_At_Coords <- eventReactive({  # This code is executed when navigate submit button is pressed
  	  input$Navigate_Around_Gene
  	  input$DisplayAllJunctions_row_last_clicked
  	  input$GeneListDisplay
  	  Ularcirc_data$Current_SelectedGene
  	}, {  # First check if new gene was selected. IF so update genomic coordinates
  	  #browser()
  	  Gene_Transcripts = circRNA_Subset()
  	  Gene_min <- min(Gene_Transcripts$Transcript$start) -501
  	  Gene_max <- max(Gene_Transcripts$Transcript$stop) + 501
  	  Gene_Coords <- as.numeric(c(Gene_min, Gene_max))

  	  if (! is.null(input$Gene_Zoom))
  	  {  # This check just ensures zoom coordinates are applied at right time.
  	    if ((input$Gene_Zoom[1] > (Gene_min - zoomOffset))  &&
  	        (input$Gene_Zoom[1] < (Gene_max + zoomOffset)) &&
  	        (input$Gene_Zoom[2] > (Gene_min - zoomOffset)) &&
  	        (input$Gene_Zoom[2] < (Gene_max + zoomOffset)))
  	    {
  	         Gene_Coords <- input$Gene_Zoom   # Navigate button selected, select user defined coords
  	    }
  	  }
  	  Gene_Coords
  	})  # This returns user requested coordinates




  	Selected_Junction_Row <- observeEvent(input$DisplayJunctionCountTable_row_last_clicked,
  	 {   # This function will identify the selected values from a table row
  	    SelectedRow <- input$DisplayJunctionCountTable_row_last_clicked

  	    if (length(grep(pattern = "Selected sample", x = input$Annotation_Options)) > 0)
  	    {
  	      Ularcirc_data$SelectedGene_from_BSJ_Table  <- Ularcirc_data$PartialPooledDataSet$RAW$Gene[SelectedRow]
  	      Ularcirc_data$Current_SelectedGene         <- Ularcirc_data$PartialPooledDataSet$RAW$Gene[SelectedRow]
  	      Ularcirc_data$Current_Selected_BS_Junction <- Ularcirc_data$PartialPooledDataSet$RAW$BSjuncName[SelectedRow]
  	    }
  	    else
  	    {
  	      Ularcirc_data$SelectedGene_from_BSJ_Table  <- Ularcirc_data$PartialDataSet$RAW$Gene[SelectedRow]
  	      Ularcirc_data$Current_SelectedGene         <- Ularcirc_data$PartialDataSet$RAW$Gene[SelectedRow]
  	      Ularcirc_data$Current_Selected_BS_Junction <- Ularcirc_data$PartialDataSet$RAW$BSjuncName[SelectedRow]
  	    }
  	    updateSelectizeInput(session,inputId="GeneListDisplay", choices=GeneList()$GeneList, selected=Ularcirc_data$Current_SelectedGene, server=TRUE)
  	 })

  	Selected_Table_Row <- observeEvent(input$TranscriptTable_row_last_clicked,
  	   {  # This function updates what transcript was last selected
  	     SelectedRow <- input$TranscriptTable_row_last_clicked
  	     Ularcirc_data$SelectedTranscript <- names(table(circRNA_Subset()$Transcript$gene))[SelectedRow]
  	   })

  	Selected_Table_Row <- observeEvent(input$ExonTable_row_last_clicked,
  	   {  # This function is not functional as yet.... what do I do with selected EXONS?

  	     SelectedRow <- input$ExonTable_row_last_clicked
  	     Row_idx <- which(circRNA_Subset()$Transcript$gene == Ularcirc_data$SelectedTranscript)
  	 #    Ularcirc_data$SelectedCanonical$Start <- Ularcirc_data$CurrentExonTable[SelectedRow,.(start)]
  	  #   Ularcirc_data$SelectedCanonical$End  <-  Ularcirc_data$CurrentExonTable[SelectedRow,.(stop)]
  	   #  Ularcirc_data$SelectedCanonical$Chr  <- Ularcirc_data$CurrentExonTable[SelectedRow,.(chrom)]
  	   })

  	Selected_Table_Row <- observeEvent(input$CanonicalJunctionCountTable_row_last_clicked,
  	   {  # This function records coordinates of a selected canonical junction so that it can be highlighted.

  	     SelectedRow <- input$CanonicalJunctionCountTable_row_last_clicked
  	     Ularcirc_data$SelectedCanonical$Start <- Ularcirc_data$CanonicalJunctionCountTable[SelectedRow,.(start1)]
  	     Ularcirc_data$SelectedCanonical$End  <-  Ularcirc_data$CanonicalJunctionCountTable[SelectedRow,.(end2)]
  	     Ularcirc_data$SelectedCanonical$Chr  <- Ularcirc_data$CanonicalJunctionCountTable[SelectedRow,.(chrom1)]
  	     Ularcirc_data$SelectedCanonical$Strand <- Ularcirc_data$CanonicalJunctionCountTable[SelectedRow,.(strand1)]
  	   })

  	Selected_Table_Row <- observeEvent(input$BS_Junction_Count_Table_row_last_clicked,
  	   {
  	     SelectedRow <- input$BS_Junction_Count_Table_row_last_clicked
  	     Ularcirc_data$Current_Selected_BS_Junction <- Ularcirc_data$BackSpliceJunctionCountTable$name[SelectedRow]

  	   })

  	Selected_Table_Row <- observeEvent(input$GenomeCanonicalJunctionCountTable_row_last_clicked,
  	   { SelectedRow <- input$GenomeCanonicalJunctionCountTable_row_last_clicked
  	     Ularcirc_data$SelectedGenome_FSJ$Start  <- Ularcirc_data$GenomeCanonicalJunctionCountTable[SelectedRow,.(start1)]
  	     Ularcirc_data$SelectedGenome_FSJ$End    <-  Ularcirc_data$GenomeCanonicalJunctionCountTable[SelectedRow,.(end2)]
  	     Ularcirc_data$SelectedGenome_FSJ$Chr    <- Ularcirc_data$GenomeCanonicalJunctionCountTable[SelectedRow,.(chrom1)]
  	     Ularcirc_data$SelectedGenome_FSJ$Strand <- Ularcirc_data$GenomeCanonicalJunctionCountTable[SelectedRow,.(strand1)]

  	   })

  	Select_Tab_Update <- observeEvent(input$PanelSelect,
  	   {  # This pre-loads required data when specific panels are selected

  	     if (input$PanelSelect == "Junction_View")
  	     {
  	        if (is.null(Ularcirc_data$Current_Selected_BS_Junction))
  	        { Ularcirc_data$Current_Selected_BS_Junction_RAWData <- data.frame(STATUS="ACTION REQUIRED", ACTION="Select junction from Gene view tab")}
  	        else  # Pre-Assemble RAW counts for junction
  	        { ## Need to filter for selected data set
  	          idx <- seq(from=1, to = length(Ularcirc_data$ProjectData$SampleIDs)) # m379()$SampleIDs))     # Default is to select everything. PROBABLY BETTER TO RETURN A NULL IF NO DATA SELECTED AND DEAL WITH THIS ELSEWHERE??
  	          if (! is.null(input$SelectedFiles))
  	          { idx <- which( Ularcirc_data$ProjectData$SampleIDs == input$SelectedFiles) } #m379()$SampleIDs == input$SelectedFiles)  }

  	          SelectedDataSets_Junctions <- Filter_by_Data_Set(idx , Ularcirc_data$ProjectData$Junctions)  # Short list to selected data set
  	          Ularcirc_data$Current_Selected_BS_Junction_RAWData <- SelectUniqueJunctions(SelectedDataSets_Junctions,
  	                                                                                 list(BSjuncName=Ularcirc_data$Current_Selected_BS_Junction, SortDir="Descending", IndexNumber=1, DisplayNumber=5))

  	          # Calculate where BSJ is detected in read. This should be close to 50% for paired end reads, and 0 % for single end reads
  	          temp <- Ularcirc_data$Current_Selected_BS_Junction_RAWData

  	          Ularcirc_data$BSJ_TypeI_vs_TypeII_Ratio <- length(grep(pattern = "M.*M",x = temp$CIGAR_1stSeg))/length(temp$CIGAR_1stSeg)
  	        }
  	     }

  	   })


  	Selected_Gene_Name <- observeEvent( input$Update_Gene_of_Interest,
  	   {  # This function waits for action button to be pressed from Gene view tab.
  	      # Updates current gene if a new gene is selected from list.
  	      # Also re-setes the flag which indicates page has been built and now allows other reactives to proceed
  	      Ularcirc_data$GenePanelLoaded <- TRUE
  	      Ularcirc_data$Current_SelectedGene <- input$GeneListDisplay
  	   })


  	UpdateProjectWD <- observeEvent(input$UpdateProjectWorkingDirectory,
  	     { # This function will re-write working directory if update button selected

  	       if (dir.exists(input$Project_WD))
  	       {
  	         extdata_path <- paste("DataPath = ", input$Project_WD, sep="")
  	         Ularcirc_extdata_path <- system.file("extdata",package = "Ularcirc")
  	         New_WD <- c(paste("# Updated ",date()),extdata_path)
  	         write(x = New_WD, file = paste(Ularcirc_extdata_path,"extdata.txt", sep="/"))
  	       }
  	       else
  	       {
  	         showNotification("That is not an existing directory. Please try again", type = "error")
  	       }
  	     })

})
