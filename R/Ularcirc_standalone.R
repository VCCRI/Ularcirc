## Could make a function that returns a couple of circRNA sequences? Could be a resource
circRNA_seq_example <- "GGAAGAGGAAGAACGTCTGAGAAATAAAATTCGAGCTGATCATGAGAAGGCCTTGGAAGAAGCAAAAGAAAAATTAAGAAAGTCAAGAGAGGAAATTCGAGCAGAAATTCAGACAGAGAAAAATAAGGTAGTCCAAGAAATGAAGATAAAAGAGAACAAGCCACTGCCACCAGTCCCTATTCCCAACCTTGTAGGAATACGTGGTGGAGACCCAGAAGATAATGACATAAGAGAGAAAAGGGAAAAAATTAAAGAGATGATGAAACATGCTTGGGATAACTATAGGACATATGGGTGGGGACATAATGAACTCAGACCTATTGCAAGGAAAGGACACTCCCCTAACATATTTGGAAGTTCACAAATGGGTGCTACCATAGTAGATGCTTTGGATACCCTTTATATCATGGGACTTCATGATGAATTCCTAGATGGGCAAAGATGGATTGAAGACAACCTTGATTTCAGTGTGAATTCAGAGGTGTCTGTGTTTGAAGTCAACATTCGATTTATTGGAGGCCTACTTGCAGCATATTACCTATCAGGAGAGGAG"




####
#' bsj_to_circRNA_sequence
#'
#' Takes one BSJ coordinate and generates a predicted circular RNA sequence.
#' @param BSJ : BSJ coordinate in the format of chr_coordinate_chr_coorindate OR chr:coordinate-coorindate:strand.
#' @param geneID : The gene ID that the BSJ aligns to. Not essential as this can
#' be identified from the BSJ coordinate, however time performance of function improved
#' if this information can be provided.
#' @param genome : Is the length f the library fragment
#' @param TxDb : The sequence read length
#' @param annotationLibrary : annotation database. See details for example.
#' @param reduce_candidates : IF multiple exon entries align to a single BSJ then either return
#' longest entry (TRUE) or all entries (FALSE)
#' @param shiny : If TRUE then will setup shiny progress bars. Default is FALSE where a standard
#'  text progress bar is used.
#' @return Returns a DNAstring object.
#' @details
#' Backsplice junction coordinates are typically reported as a character string. Two formats
#' are recognised, ":" delimited (eg circExplorer, CIRI) or "_" delimited (Ularcirc). The
#' BSJ genomic coordinates are compared against the supplied gene model and exonic sequences
#' from matching splice junctions are concatenated. This means the BSJ is the first and last
#' nucleotide of the returned sequence. The current implementation will automatically check 0 or
#' 1 base coordinates and any match is returned.
#'
#' In some cases one BSJ will match multiple exon combinations. The default setting is to return
#' the longest sequence. Alternatively all possibilities can be returned by setting
#' reduce_candidates to FALSE. BSJ candidates that align to multiple exon combinations are
#' added to duplicated list.
#'
#' BSJ that do not align to any canonical junctions are returned as failed.
#'
#' @examples
#'
#' library('Ularcirc')
#' TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' annotationLibrary <- org.Hs.eg.db::org.Hs.eg.db
#'
#' # Define BSJ. Following two formats are accepted
#' BSJ <- 'chr2:40430305-4 0428472:-'       # SLC8A1
#' BSJ  <- 'chr2_40430305_chr2_40428472'   # SLC8A1
#'
#' circRNA_sequence <- bsj_to_circRNA_sequence(BSJ, "SLC8A1", genome,TxDb, annotationLibrary)
#'
#' # You can also retrieve sequence without passing gene annotation - but this is slower
#' # circRNA_sequence <- bsj_to_circRNA_sequence(BSJ, NULL, genome,TxDb, annotationLibrary)
#'
#' TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' # EXAMPLE1 (3 fail and 2 will produce sequences)
#' BSJ <- c("chr14_99465814_chr14_99458278","chr22_20933778_chr22_20934245",
#'          "chr12_120155720_chr12_120154969", "chr4_143543508_chr4_143543973",
#'          "chr10_7285955_chr10_7276891")
#' GeneIDs <- c("SMARCA5","MSLN","RNF138","KIAA0368","CRKL")
#' circRNA_sequence <- bsj_to_circRNA_sequence(BSJ, GeneIDs, genome,TxDb, annotationLibrary)
#'
#' # Returns a list with three items:
#' # (1) "identified" is a list of DNA strings from BSJ that aligned to FSJ coordinates of the gene model
#' # (2) "failed" is a character object of BSJ that did not align to FSJ coordinates of gene model. Each entry is
#' # named with gene ID.
#' # (3) "duplicates" (not implemented yet) identifies which BSJ returned multiple sequences
#' @export
#'
# load(file="c:/TEMP/TestData.RData")  # BSJ and GeneID in circExplorer2 format
#  bsj_to_circRNA_sequence(BSJ[1:100], as.character(GeneID[1:100]), genome=genome_hg38, TxDb=TxDb_hg38,annotationLibrary = annotationLibrary)
#  Does not work on unannotated data.
#
#
#  CURRENTLY ANNOTATION FOR MULTIPLE ENTRIES DOES NOT WORK
bsj_to_circRNA_sequence <- function(BSJ, geneID=NULL, genome, TxDb,
                                    annotationLibrary, reduce_candidates = TRUE, shiny=FALSE)
{
  lookupID <- {}
  BSJ_donor <- {}
  BSJ_acceptor <- {}

  # Need to distinguish the following formats
  # chr10:100923974-100926020:+
  # chr11_33286412_chr11_33287512

  if (length(gregexpr("_",BSJ)[[1]]) == 3) # Ularcirc format  eg chr14_99465814_chr14_99458278
  {
    BSjuncDetails <- strsplit(BSJ, split = "_")
    BSJ_donor <- as.numeric(unlist(lapply(BSjuncDetails, FUN = function(x){min(x[c(2,4)])})))
    BSJ_acceptor <- as.numeric(unlist(lapply(BSjuncDetails, FUN = function(x){max(x[c(2,4)])})))
  }
  else if (length(gregexpr(":",BSJ)[[1]]) == 2) # generic format eg chr10:100923974-100926020:+
  {
    BSjuncDetails <- strsplit(BSJ, split = ":")
    coordinates <- lapply(BSjuncDetails,FUN = function(x) { strsplit(x[2],split="-")})
    BSJ_donor <- as.numeric(lapply(coordinates, FUN = function(x) { min(unlist(x)) }))
    BSJ_acceptor <- as.numeric(lapply(coordinates, FUN = function(x) { max(unlist(x)) }))
  }
  else
  {
    warning("BSJ not in the correct format. Did not detect separating characters.")
    return(NULL)
  }
  if ((length(BSJ_donor) == 0) || (length(BSJ_acceptor) ==0 ))
  {
    warning("Could not extract coordinates from BSJ (please check). Aborting")
    return(NULL)
  }


  if (is.null(geneID))  # Identify potential gene coordinate
  { warning("This function currently requires annotated data. Please ensure you use parameter geneID")
    return(NULL)
    g_GR <- GenomicFeatures::genes(TxDb)
    strand <- "*"
    bs_junc_gr <- GenomicRanges::GRanges(seqnames=BSjuncDetails[[1]][1], ranges = as.numeric(min(BSjuncDetails[[1]][c(2,4)]),min(BSjuncDetails[[1]][c(2,4)])),strand = strand)
    t_start <- GenomicAlignments::findOverlaps(BiocGenerics::invertStrand(bs_junc_gr),g_GR, type=c("within"))
    bs_junc_gr <- GenomicRanges::GRanges(seqnames=BSjuncDetails[[1]][1], ranges = as.numeric(max(BSjuncDetails[[1]][c(2,4)]),max(BSjuncDetails[[1]][c(2,4)])),strand = strand)
    t_end <- GenomicAlignments::findOverlaps(BiocGenerics::invertStrand(bs_junc_gr),g_GR, type=c("within"))

    entrezID <- c("Novel")
    #	  browser()
    if ((length(t_start) > 0) && (length(t_end) > 0))
    {	entrezID_start <- g_GR[S4Vectors::subjectHits(t_start)]$gene_id
      entrezID_end <- g_GR[S4Vectors::subjectHits(t_end)]$gene_id
      entrezID <- intersect(entrezID_start, entrezID_end)
      lookupID <- entrezID
    }

  }


  if (! is.null(geneID))  # Gene name was provided. Trust this and obtain entrezID
  { message("Extracting entrez lookup ID")
    if (shiny)
    { shiny::incProgress(1/4, detail = paste("Extracting entrez lookup ID")) }

    a <- try(AnnotationDbi::select(annotationLibrary, keys = geneID, columns=c("ENTREZID", "SYMBOL", "ENSEMBL"),keytype="SYMBOL"),silent=TRUE)
    if(length(grep(pattern = "Error", x = a)))
    {   # cannot continue
      warning(paste(a,"Error obtaining annotation information.","\n"))
      return(-1)
    }
    lookupID <- a$ENTREZID   # Default lookup
  }

  if (length(BSJ) != length(geneID))
  {
    warning("Length of BSJ is not same as length of geneID, please check.")
    return(NULL)
  }
  else
  { names(BSJ) <- geneID }


  genes_to_consider <- length(geneID)
  if(genes_to_consider > 1)
    reduce_candidates <- FALSE

  if (! is.null(lookupID))
  { message("Extracting exon coordinates")
    if (shiny)
    { shiny::incProgress(1/4, detail = paste("Extracting exon coordinates")) }

    # Create exon table
    b <- try(AnnotationDbi::select(TxDb, keys = lookupID, columns=c('GENEID', 'TXCHROM', 'EXONSTART',  'EXONEND','TXID', 'EXONSTRAND'),keytype="GENEID"),silent=TRUE)
    if(length(grep(pattern = "Error", x = b)))
    {   # cannot continue
      warning(paste(b,"Error obtaining annotation information.","\n"))
      return(-1)
    }

    b <- b[,c("GENEID","TXCHROM","EXONSTART","EXONEND","TXID","EXONSTRAND")]
    names(b) <- c('gene', 'chrom', 'start', 'stop', 'txid','strand')
    b <- data.table::as.data.table(b)

    # Need to work through all entries (list) to extract possible exon boundaries
    if (shiny)
    { shiny::incProgress(1/4, detail = paste("Short listing exons")) }

    message("short listing exons")
    pb <- R.utils::ProgressBar(max=length(BSJ_donor))
    R.utils::reset(pb)
    circRNA_exons <- {}
    failed_BSJ <- {}
    for (i in seq_along(1:length(BSJ_donor)))
    {
      idx <- b$start >= BSJ_donor[i] & b$stop <= BSJ_acceptor[i]
      possible_exons <- b[idx,]

      # select candidates. In some situations coordinates may be entered in 0 or 1 base.
      true_candidates_stop <- which(abs(possible_exons$start - BSJ_donor[i]) == 1 |
                                      abs(possible_exons$start - BSJ_donor[i]) == 0)
      true_candidates_start <- which(abs(possible_exons$stop - BSJ_acceptor[i]) == 1 |
                                       abs(possible_exons$stop - BSJ_acceptor[i]) == 0 )

      true_candidate_IDs <- intersect(possible_exons$gene[true_candidates_start],
                                      possible_exons$gene[true_candidates_stop])
      exon_idx <- possible_exons$gene %in% true_candidate_IDs

      if (length(exon_idx) == 0)
      { failed_BSJ <- c(failed_BSJ, BSJ[i]) }

      if (length(circRNA_exons) == 0)
      { circRNA_exons <- possible_exons[exon_idx,]}
      else
      {
        circRNA_exons <- rbind(circRNA_exons, possible_exons[exon_idx,])
      }
      R.utils::increase(pb)
    }  # for (i in seq_along(1:length(BSJ_donor)))

    if (nrow(circRNA_exons) == 0)
    {
      warning("No exon junctions aligned with BSJ coorindates, therfore no sequences recovered")
      return(NULL)
    }
    ### Now to work out maximum length by sifting through tx entries and adding up exon lengths
    circRNA_exon_lengths <-by(circRNA_exons,circRNA_exons$gene,identity )  # This makes a list of all transcripts
    circRNA_sizes <- lapply(X = circRNA_exon_lengths,FUN = function(x) {  return(abs(sum(x$start-x$stop)))})
    circRNA_sizes_idx <- order(unlist(circRNA_sizes), decreasing = TRUE)


    # Sometimes a single BSJ can align to multiple exons within a transcript.
    # In these situations will choose the longest entry as the candidate circRNA entry.
    Exons_of_Interest <- circRNA_exons  # default is to select all entries
    duplicated_BSJ <- {}
    if (reduce_candidates)
    { largest_transcript_ID <- names(circRNA_sizes[circRNA_sizes_idx])[1]
      transcript_ID_idx <- which(names(circRNA_exon_lengths) == largest_transcript_ID)  # Get Index of transcript ID name
      transcript_ID_idx <- which(circRNA_exons$gene == names(circRNA_exon_lengths)[transcript_ID_idx])  # Get Row index(es) corresponding to transcript
      Exons_of_Interest <- circRNA_exons[transcript_ID_idx,]
      Exons_of_Interest <- Exons_of_Interest[ ! duplicated(Exons_of_Interest)]  # Sometimes there is duplicated records.
    }
    else # Will return all candidates
    { # would be nice to record how many possible duplicate entries
    #  all_entries <- table(names(circRNA_sizes))


    }

    message("Extracting sequence")
    if (shiny)
    { shiny::incProgress(1/4, detail = paste("Extracting sequence")) }
    else
    {
      pb <- R.utils::ProgressBar(max=length(circRNA_exon_lengths))
      R.utils::reset(pb)
    }

     all_seq <- Biostrings::getSeq(genome,Exons_of_Interest$chrom,start=Exons_of_Interest$start, end=Exons_of_Interest$stop, strand=Exons_of_Interest$strand)
     names(all_seq) <- Exons_of_Interest$gene
     exons_to_stitch <- table(Exons_of_Interest$gene)
     temp <- by(data.frame(as.character(all_seq), names(all_seq)), names(all_seq), identity)
     all_circRNA_seq <- lapply(temp, FUN= function(x) {  R.utils::increase(pb);
          paste(as.character(x[,1]) ,sep="",collapse = "")})

#    message("Extracting sequence - method 2")
#    R.utils::reset(pb)
#    circRNA_exon_lengths <- by(circRNA_exons,Exons_of_Interest$gene,identity )
#    circRNA_sequence <- lapply(X = circRNA_exon_lengths, FUN= function(x) {     R.utils::increase(pb); sequence_from_exon_coords(genome, x) })

    return(list(identified=all_circRNA_seq, failed=failed_BSJ, duplicates=duplicated_BSJ))
  } # if (! is.null(lookupID))

  warning("Cannot find or match gene ID")
  return(NULL)

}

#################################################################333
#' sequence_from_exon_coords
#'
#' @param genome : genome object
#' @param exon_df : data frame of exons. Must have column with names "chrom", "start",
#'    "stop", "strand"
#'
sequence_from_exon_coords <- function(genome, exon_df)
{
  seq <- ''
  if (nrow(exon_df) < 1)
  { return(seq) }


  tmp <- as.character(Biostrings::getSeq(genome,exon_df$chrom,
                                           start=exon_df$start,
                                           end=exon_df$stop,
                                           strand = exon_df$strand)  )
  seq <- paste(tmp,sep="",collapse = "")

  return(seq)
}


####
#' bsj_fastq_generate
#'
#' Takes a circRNA predicted sequence and generates synthetic short sequence reads
#' @param circRNA_Sequence : Linear sequence of a circRNA. i.e. the backspice junction
#'               is the first and last base of this sequence
#' @param fragmentLength : Is the length the library fragment
#' @param readLength : The sequence read length
#' @param variations : Number of sequences returned for each read type. Note each
#' sequence variation will start at a unique location (where possible)
#' @param headerID : Character identifier that will be incorporated into sequence header
#' @return Returns a list of two DNAstring sets labelled "read1" and "read2" which correspond
#' to forward and reverse read pairs.
#'
#' @examples
#'
#' library('Ularcirc')
#'
#'
#' # Generate a 500nt sequence containing A and which is flanked with GG and CC.
#' circRNA_Sequence <- paste(rep('A',500),collapse='')
#' circRNA_Sequence <- paste('GG',circRNA_Sequence, 'CC', sep='')
#' # The GG and CC ends of sequence represent ends of linear exons that are circularised.
#' # Therefore the backsplice junction (BSJ) is GGCC.
#' # Generate reads that alternate over this BSJ
#'
#' fastqReads <- bsj_fastq_generate(circRNA_Sequence, fragmentLength=300, readLength=100,
#'                variations = 4,   # Four type I , II, III, and IV reads generated
#'                headerID='circRNA_example')  # Identifier incorporated in name of each sequence
#' # The following will indicate 12 sequences are present in each list entry
#' length(fastqReads$read1)
#' length(fastqReads$read2)
#'
#' # Can create fastq file as follows
#' Biostrings::writeXStringSet( fastqReads$read1,"circRNA_Sample_R1.fastq.gz",
#'                               compress = TRUE, format="fastq")
#' Biostrings::writeXStringSet( fastqReads$read2,"circRNA_Sample_R2.fastq.gz",
#'                               compress = TRUE, format="fastq")
#' @import Biostrings
#'
#' @export
bsj_fastq_generate <- function(circRNA_Sequence, fragmentLength=300, readLength=100, variations = 4, headerID='')
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

  read_one <- {};  read_two <- {};
  # Generate Type II reads
  start_pos <- circ_length - typeII_offset
  for (i in 1:variations)
  { read_one <- c(read_one, substr(x = circRNA_Sequence, start = start_pos, stop = start_pos + readLength))
  read_two <- c(read_two, substr(x = circRNA_Sequence, start = start_pos+fragmentLength-readLength, stop = start_pos + fragmentLength))
  names(read_one)[i] <- paste(headerID,"_typeII_",start_pos, common_label,sep="")
  start_pos <- start_pos + typeII_III_offset_step_size
  }

  # Generate Type III reads
  start_pos <- circ_length - typeIII_offset
  for (i in 1:variations)
  { read_one <- c(read_one, substr(x = circRNA_Sequence, start = start_pos, stop = start_pos + readLength))
    read_two <- c(read_two, substr(x = circRNA_Sequence, start = start_pos+fragmentLength-readLength, stop = start_pos + fragmentLength))
    names(read_one)[length(read_one)] <- paste(headerID,"_typeIII_",start_pos,common_label,sep="")
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
    { read_one <- c(read_one, substr(x = circRNA_Sequence, start = start_pos, stop = start_pos + readLength))
    read_two <- c(read_two, substr(x = circRNA_Sequence, start = start_pos+fragmentLength-readLength, stop = start_pos + fragmentLength))
    names(read_one)[length(read_one)] <- paste(headerID,"_typeIV_",start_pos,common_label,sep="")
    start_pos <- start_pos + typeIV_offset_step_size
    }
  }

  names(read_two) <- names(read_one)
  read_one <- Biostrings::DNAStringSet(x=read_one)
  read_two <- Biostrings::reverseComplement(Biostrings::DNAStringSet(x=read_two))
  return(list(read1 =read_one, read2 = read_two))
}

