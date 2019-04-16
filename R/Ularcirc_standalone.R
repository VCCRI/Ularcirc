
circRNA_seq_example <- "GGAAGAGGAAGAACGTCTGAGAAATAAAATTCGAGCTGATCATGAGAAGGCCTTGGAAGAAGCAAAAGAAAAATTAAGAAAGTCAAGAGAGGAAATTCGAGCAGAAATTCAGACAGAGAAAAATAAGGTAGTCCAAGAAATGAAGATAAAAGAGAACAAGCCACTGCCACCAGTCCCTATTCCCAACCTTGTAGGAATACGTGGTGGAGACCCAGAAGATAATGACATAAGAGAGAAAAGGGAAAAAATTAAAGAGATGATGAAACATGCTTGGGATAACTATAGGACATATGGGTGGGGACATAATGAACTCAGACCTATTGCAAGGAAAGGACACTCCCCTAACATATTTGGAAGTTCACAAATGGGTGCTACCATAGTAGATGCTTTGGATACCCTTTATATCATGGGACTTCATGATGAATTCCTAGATGGGCAAAGATGGATTGAAGACAACCTTGATTTCAGTGTGAATTCAGAGGTGTCTGTGTTTGAAGTCAACATTCGATTTATTGGAGGCCTACTTGCAGCATATTACCTATCAGGAGAGGAG"




####
#' BSJ_to_circRNA_sequence
#'
#' Takes one BSJ coordinate and generates a predicted circular RNA sequence.
#' @param BSJ : BSJ coordinate
#' @param geneID : The gene ID that the BSJ aligns to. Not essential as this can
#' be identified from the BSJ coordinate, however there is a dramatic speed improvemnt
#' if this can be passed to the function.
#' @param genome : Is the length f the library fragment
#' @param TxDb : The sequence read length
#' @return Returns a DNAstring object.
#' @examples
#'
#' library('Ularcirc')
#' TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' BSJ  <- 'chr2_40430305_chr2_40428472'   # chr2_40430302_chr2_40428472
#' annotationLibrary <- org.Hs.eg.db
#' BSJ_to_circRNA_sequence(BSJ, "SLC8A1", genome,TxDb, annotationLibrary)
#' # chr11_33286412_chr11_33287512
#'
#' TxDb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#' genome_hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' @export
BSJ_to_circRNA_sequence <- function(BSJ, geneID=NULL, genome, TxDb, annotationLibrary)
{
  lookupID <- {}
  BSJ_donor <- {}
  BSJ_acceptor <- {}

  # Need to distinguish the following formats
  # chr10:100923974-100926020:+
  # chr11_33286412_chr11_33287512

  if (length(gregexpr("_",BSJ)[[1]]) == 3) # Ularcirc format
  {
    BSjuncDetails <- strsplit(BSJ, split = "_")
    BSJ_donor <- as.numeric(min(BSjuncDetails[[1]][c(2,4)]))
    BSJ_acceptor <- as.numeric(max(BSjuncDetails[[1]][c(2,4)]))
  }
  else if (length(gregexpr(":",BSJ)[[1]]) == 2) # generic format
  {
    BSjuncDetails <- strsplit(BSJ, split = ":")
    coordinates <- unlist(strsplit(BSjuncDetails[[1]][2], split="-"))
    BSJ_donor <- as.numeric(min(coordinates))
    BSJ_acceptor <- as.numeric(max(coordinates))
  }
  else
  {
    warning("BSJ not in the correct format. Did not detect separating characters.")
    return(NULL)
  }



  if (is.null(geneID))  # Identify potential gene coordinate
  {
    g_GR <- genes(TxDb)

    bs_junc_gr <- GRanges(seqnames=BSjuncDetails[[1]][1], ranges = as.numeric(min(BSjuncDetails[[1]][c(2,4)]),min(BSjuncDetails[[1]][c(2,4)])),strand = strand)
    t_start <- findOverlaps(invertStrand(bs_junc_gr),g_GR, type=c("within"))
    bs_junc_gr <- GRanges(seqnames=BSjuncDetails[[1]][1], ranges = as.numeric(max(BSjuncDetails[[1]][c(2,4)]),max(BSjuncDetails[[1]][c(2,4)])),strand = strand)
    t_end <- findOverlaps(invertStrand(bs_junc_gr),g_GR, type=c("within"))

    entrezID <- c("Novel")
    #	  browser()
    if ((length(t_start) > 0) && (length(t_end) > 0))
    {	entrezID_start <- g_GR[subjectHits(t_start)]$gene_id
    entrezID_end <- g_GR[subjectHits(t_end)]$gene_id
    entrezID <- intersect(entrezID_start, entrezID_end)
    entrezID <- try(select(GeneList$Annotation_Library, keys = entrezID, columns=c("SYMBOL"),keytype="ENTREZID")[,'SYMBOL'],silent=TRUE) # Convert to Symbol

    if(length(grep(pattern = "Error in ", x = entrezID)))  # This is start of error message when lookup is not linked
    {   entrezID <- intersect(entrezID_start, entrezID_end)
    entrezID <- try(select(GeneList$Annotation_Library, keys = entrezID, columns=c("SYMBOL"),keytype="ENSEMBL"),silent=TRUE)
    entrezID <- entrezID$SYMBOL
    if(length(grep(pattern = "Error in ", x = entrezID)))
    {    entrezID <- c("Unknown") }
    }
    }
    DataSet$Gene[i] <- paste(unique(entrezID),collapse=",")

    ################## following code taken from Prepare_Gene_Object
    a <- try(select(GeneList$Annotation_Library, keys = GeneName, columns=c("ENTREZID", "SYMBOL", "ENSEMBL"),keytype="SYMBOL"),silent=TRUE)
    lookupID <- a$ENTREZID
  }

  if (! is.null(geneID))  # Gene name was provided. Trust this and obtain entrezID
  {
    a <- try(select(annotationLibrary, keys = geneID, columns=c("ENTREZID", "SYMBOL", "ENSEMBL"),keytype="SYMBOL"),silent=TRUE)
    if(length(grep(pattern = "Error", x = a)))
    {   # We cannot continue
      return(-1)
    }
    lookupID <- a$ENTREZID   # Default lookup
  }

  if (! is.null(lookupID))
  {
    # Create exon table
    b <- select(TxDb, keys = lookupID, columns=c('GENEID', 'TXCHROM', 'EXONSTART',  'EXONEND','TXID', 'EXONSTRAND'),keytype="GENEID")
    b <- b[,c("TXCHROM","EXONSTART","EXONEND","TXID","EXONSTRAND")]
    names(b) <- c('chrom', 'start', 'stop', 'gene','strand')
    b <- as.data.table(b)

    # Short list exons
#    BSJ_donor <- as.numeric(min(BSjuncDetails[[1]][c(2,4)]))
#    BSJ_acceptor <- as.numeric(max(BSjuncDetails[[1]][c(2,4)]))
    possible_exons <- b[start >= BSJ_donor & stop <=  BSJ_acceptor,]

    # select candidates. In some situations coordinates may be entered in 0 or 1 base.

    true_candidates_stop <- which(abs(possible_exons$start - BSJ_donor) == 1 |
                                    abs(possible_exons$start - BSJ_donor) == 0)
    true_candidates_start <- which(abs(possible_exons$stop - BSJ_acceptor) == 1 |
                                     abs(possible_exons$stop - BSJ_acceptor) == 0 )

    true_candidate_IDs <- intersect(possible_exons$gene[true_candidates_start],
                                    possible_exons$gene[true_candidates_stop])
    exon_idx <- possible_exons$gene %in% true_candidate_IDs
    circRNA_exons <- possible_exons[exon_idx,]
    if (nrow(circRNA_exons) < 1)
    { warning("No exons exist within these coordinates")
      return(NULL)
    }

    ### Now to work out maximum length by sifting through tx entries and adding up exon lengths
    circRNA_exon_lengths <-by(circRNA_exons,circRNA_exons$gene,identity )  # This makes a list of all transcripts

    circRNA_sizes <- lapply(X = circRNA_exon_lengths,FUN = function(x) { return(abs(sum(x$start-x$stop)))})
    circRNA_sizes_idx <- order(unlist(circRNA_sizes), decreasing = TRUE)
    largest_transcript_ID <- names(circRNA_sizes[circRNA_sizes_idx])[1]
    transcript_ID_idx <- which(names(circRNA_exon_lengths) == largest_transcript_ID)  # Get Index of transcript ID name
    transcript_ID_idx <- which(circRNA_exons$gene == names(circRNA_exon_lengths)[transcript_ID_idx])  # Get Row index(es) corresponding to transcript
    Exons_of_Interest <- circRNA_exons[transcript_ID_idx,]
    Exons_of_Interest <- Exons_of_Interest[ ! duplicated(Exons_of_Interest)]  # Sometimes there is duplicated records.

    circRNA_Sequence <- ''
    FSJs <- c(1)# This will contain start positions for ALL Forward splice junctions
    for (i in 1:nrow(Exons_of_Interest))   # Need to stitch together multiple exons
    { tmp <- as.character(getSeq(genome,Exons_of_Interest$chrom[i],
                                 start=Exons_of_Interest$start[i],
                                 end=Exons_of_Interest$stop[i],
                                 strand = Exons_of_Interest$strand[i])  )
    FSJs <- c(FSJs, FSJs[i] + nchar(tmp))
    circRNA_Sequence <- paste(circRNA_Sequence,tmp,sep="",collapse = "")
    }

    circRNA_Sequence <- DNAString(circRNA_Sequence)
    return(circRNA_Sequence)
  }

  return(NULL)

}



####
#' BSJ_Fastq_Generate
#'
#' Takes a circRNA predicted sequence and generates synthetic sequence reads that
#' @param circRNA_Sequence : Linear sequence of a circRNA. i.e. the backspice junction
#'               is the first and last base of this sequence
#' @param fragmentLength : Is the length f the library fragment
#' @param readLength : The sequence read length
#' @return Returns a list DNAstring sets.
#' @examples
#'
#' library('Ularcirc')
#'
#'
#' # Obtain a circRNA sequence
#' circRNA_Sequence <- circRNA_seq_example
#' fastqReads <- BSJ_Fastq_Generate(circRNA_Sequence, fragmentLength=300, readLength=100, variations = 4, headerID='')
#' writeXStringSet( fastqReads$Read_One,"circRNA_Sample_R1.fastq.gz",compress = TRUE, format="fastq")
#' writeXStringSet( fastqReads$Read_Two,"circRNA_Sample_R2.fastq.gz",compress = TRUE, format="fastq")
#'
#' @export
BSJ_Fastq_Generate <- function(circRNA_Sequence, fragmentLength=300, readLength=100, variations = 4, headerID='')
{
  if (variations < 1)
  { warning("Number of fragment variations must be 1 or more. Resetting to 1")
    variations <- 1
  }
  # Variations: Number of read variations prepared for each read type.
  #
  circ_length <- nchar(circRNA_Sequence)
  if (fragmentLength > circ_length)
  {
    warning("Fragment length is larger than circRNA length. Please check input sequence. Returning NULL.")
    return(NULL)
  }

  circRNA_Sequence <- paste(circRNA_Sequence, circRNA_Sequence, sep = "")
  typeII_III_offset_step_size <- round(readLength/(variations + 1))
  typeIV_gap_size <- fragmentLength - (2 * readLength)
  typeIV_offset_step_size <- round(typeIV_gap_size / (variations + 1))



  if (fragmentLength <= readLength)
  { warning("fragment length is shorter than or equal to read length.
            Increasing fragment length to read length + 1")
    fragmentLength <- readLength + 1
  }

  # Calculate offset to generate the appropriate alignment types.
  typeII_offset  <- readLength - typeII_III_offset_step_size
  typeIII_offset <- fragmentLength - typeII_III_offset_step_size
  typeIV_offset  <- readLength + typeIV_gap_size - typeIV_offset_step_size
  common_label <- paste("_F",fragmentLength,"_R",readLength, sep="")

  Read_One <- {};  Read_Two <- {};
  # Generate Type II reads
  start_pos <- circ_length - typeII_offset
  for (i in 1:variations)
  { Read_One <- c(Read_One, substr(x = circRNA_Sequence, start = start_pos, stop = start_pos + readLength))
  Read_Two <- c(Read_Two, substr(x = circRNA_Sequence, start = start_pos+fragmentLength-readLength, stop = start_pos + fragmentLength))
  names(Read_One)[i] <- paste(headerID,"_typeII_",start_pos, common_label,sep="")
  start_pos <- start_pos + typeII_III_offset_step_size
  }

  # Generate Type III reads
  start_pos <- circ_length - typeIII_offset
  for (i in 1:variations)
  { Read_One <- c(Read_One, substr(x = circRNA_Sequence, start = start_pos, stop = start_pos + readLength))
  Read_Two <- c(Read_Two, substr(x = circRNA_Sequence, start = start_pos+fragmentLength-readLength, stop = start_pos + fragmentLength))
  names(Read_One)[length(Read_One)] <- paste(headerID,"_typeIII_",start_pos,common_label,sep="")
  start_pos <- start_pos + typeII_III_offset_step_size
  }


  if (fragmentLength > (readLength*2))
  {  # Generate typeIV reads

    if (typeIV_offset_step_size < 1)  # Ensure we have a step increment.
    {   typeIV_offset_step_size <- 1      }
    else if (round(typeIV_gap_size/typeIV_offset_step_size) < variations)
    { warning("Requested number of fastq are more than what is possible, downsizing number of generated entries")
      variations <- round(typeIV_gap_size/typeIV_offset_step_size)
    }

    start_pos <- circ_length - typeIV_offset
    for (i in 1:variations)
    { Read_One <- c(Read_One, substr(x = circRNA_Sequence, start = start_pos, stop = start_pos + readLength))
    Read_Two <- c(Read_Two, substr(x = circRNA_Sequence, start = start_pos+fragmentLength-readLength, stop = start_pos + fragmentLength))
    names(Read_One)[length(Read_One)] <- paste(headerID,"_typeIV_",start_pos,common_label,sep="")
    start_pos <- start_pos + typeIV_offset_step_size
    }
  }

  names(Read_Two) <- names(Read_One)
  Read_One <- Biostrings::DNAStringSet(x=Read_One)
  Read_Two <- Biostrings::reverseComplement(Biostrings::DNAStringSet(x=Read_Two))
  return(list(readOne=Read_One, readTwo=Read_Two))
}







