
##########################################################################################
## Return_Splice_Juncs
##
## This function takes a gene symbol and returns some basic properties of that gene. Returns
## the start and end genomic coordinates, the strand, chromosome and number of exons.
##
Return_Splice_Juncs <- function(transcript_reference, GeneName)
{
  strand<-as.numeric(strsplit(as.character(transcript_reference[name2==GeneName,strand]),split = ",")=='+')
  starts<-(strsplit(as.character(transcript_reference[name2==GeneName,exonStarts]),split = ","))
  stops<-(strsplit(as.character(transcript_reference[name2==GeneName,exonEnds]),split = ","))
  exons<- as.character(transcript_reference[name2==GeneName,exonStarts])
  chrom<-strsplit(as.character(transcript_reference[name2==GeneName,chrom]),split = ",")
  
  exon_number <- max(sapply(regmatches(exons, gregexpr(",", exons)), length)) 
  tmp <- do.call(rbind,starts)
  # Initially assume gene is on +ve strand
  transcript_start  <- min(tmp[,1])    
  transcript_stop <- max(tmp[,ncol(tmp)])
  if   (strand == 0)   # -ve strand
  {    transcript_start  <- max(tmp[,ncol(tmp)])  
  transcript_stop  <-  min(tmp[,1])
  strand <- 2                          # To be compatible with STAR output
  }
  chrom <-do.call(rbind,chrom)[1,1]     # Make matrix from a list
  
  return(list(exon_number=exon_number, chrom=chrom, transcript_start = transcript_start, transcript_stop = transcript_stop, strand=strand, GeneName=GeneName ))
}


Obtain_Key_Gene_Coordinates <- function(Gene_Symbol,  GeneList) ## This function is redundant. Look at Gene_Transcript_Features() function in server.R
{
  #  Genomic range object can be used to extract start and stop coordinates
##  g_GR <- genes(GeneList$transcript_reference)
##  entrezID <- g_GR[subjectHits(t)]$gene_id   # Can use this to look up geneID
##  entrezID <- select(GeneList$Annotation_Library, keys = entrezID, columns=c("SYMBOL"),keytype="ENTREZID")[,'SYMBOL'] # Convert to Symbol
  
  
  # Alternatively use the annotation database. This may be easier as can shift between gene and transcript information.
  a <- select(GeneList$Annotation_Library, keys = Gene_Symbol, columns=c("ENTREZID", "SYMBOL"),keytype="SYMBOL")
  b <- select(GeneList$transcript_reference, keys = a$ENTREZID, columns=c('GENEID', 'TXNAME', 'EXONRANK'),keytype="GENEID")
  Num_Transcripts <- length(unique(b$TXNAME))
  Num_Exons_Per_Transcript <- {}
  for(i in 1:length(unique(b$TXNAME)))
  {     Num_Exons_Per_Transcript <- c(Num_Exons_Per_Transcript , length(which(b$TXNAME == unique(b$TXNAME)[i])))          }
  Num_Exons_Per_Transcript <- as.numeric(Num_Exons_Per_Transcript)
  names(Num_Exons_Per_Transcript) <- unique(b$TXNAME)
  
  
}


JunctionReport <- function(JunctionData, transcript_reference, Transcript_Features, Stringency_filter_threshold= 4, Stringency_threshold_minimum_threshold = 5)
{
  # WARNING assuming JunctionData is in correct format. We will label columns 
  JunctionData <- lapply(JunctionData , FUN=function(x) {  setnames(x,1:9,c("chrom","start","stop", "strand", "intron_motif", "Annotated","UniqueMapped","Multimapped","Overhang")) } )
  
  Collated_Junctions_Count <- data.frame()
  for(i in 1:length(JunctionData))
  {
    if (Transcript_Features$strand == 2)     # negative strand
    {     TranscriptJunctions_existing_in_RawData   <-  JunctionData[[i]][chrom==Transcript_Features$chrom & strand==Transcript_Features$strand & start > Transcript_Features$transcript_stop & stop < Transcript_Features$transcript_start,]     }
    else     # positive strand
    {    Transcript_Features$strand <-  1
    TranscriptJunctions_existing_in_RawData   <-  JunctionData[[i]][chrom==Transcript_Features$chrom & strand==Transcript_Features$strand & start > Transcript_Features$transcript_start & stop < Transcript_Features$transcript_stop,]   
    }
    if (nrow(TranscriptJunctions_existing_in_RawData) == 0)
    {     return(NULL)     }
    a <- order(x = TranscriptJunctions_existing_in_RawData$UniqueMapped, decreasing = TRUE )[seq(from=1, to=Transcript_Features$exon_number -1)]
    stringency_filter <- ave(x = TranscriptJunctions_existing_in_RawData$UniqueMapped[a])[1] / Stringency_filter_threshold
    
    #          cat(paste("\nstringency_filter is: ",stringency_filter, " i is:",i)) 
    if (is.na(stringency_filter))
    {     return(NULL)     }
    
    if (stringency_filter  < Stringency_threshold_minimum_threshold )
    {     return(NULL)     }          # Exit the analysis as we are near the level of noise
    
    Junction_Result <-   table(TranscriptJunctions_existing_in_RawData$UniqueMapped > stringency_filter)  
    Junctions_Above_Threshold <- Junction_Result["TRUE"]
    
    if (i == 1)
    {          Collated_Junctions_Count <- data.frame(Junctions_Above_Threshold)     }
    else
    {               Collated_Junctions_Count <- cbind(Collated_Junctions_Count, Junctions_Above_Threshold)  }
  }
  colnames(Collated_Junctions_Count) <- names(JunctionData)
  row.names(Collated_Junctions_Count) <- Transcript_Features$GeneName
  return(Collated_Junctions_Count)
}
