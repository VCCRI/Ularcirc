#' Ularcirc
#'
#' When the function is invoked the Ularcirc shiny app is started. The starting screen has quickstart instructions on how to use the software.
#' Please refer to the Ularcirc vignette for a more detailed workflow.
#' @return
#' Does not return anything
#' @examples
#' # The following commands will load the shiny app either through an RStudio session or
#' # through your internet browser
#'
#' library("Ularcirc")
#' \donttest{ Ularcirc() }
#'
#'
#' @export
Ularcirc <- function()
{
  appDir <- system.file("shiny-app", "circRNA", package = "Ularcirc")

  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `Ularcirc`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode="normal")

}

#' Compatible_Annotation_DBs
#'
#' Interogates bioconductor databases and identifies those that are compatible
#' with Ularcirc. Builds a list of commands that the user can copy to install the
#' required database on their local computer. Once installed the databases are
#' immediately available to Ularcirc upon re-starting the shiny app.
#' This function requires connection to the internet.
#' @return Returns a list of installation commands for all compatible annotation databases
#' @examples
#' # Get all bioconductor annotation databases that are compatible with Ularcirc
#' library('BSgenome')
#' library('httpuv')
#' library('AnnotationHub')
#' all_dbs <- Compatible_Annotation_DBs()
#'
#' # List database IDs
#' names(all_dbs)
#'
#' # Select a database and display the commands needed to install
#' # Use noquote to correctly format output
#' \donttest{ noquote(all_dbs$Hsapiens.UCSC.hg38) }
#' @export
Compatible_Annotation_DBs <- function()
{
  ah <- AnnotationHub::AnnotationHub()
  all_OrgDb <- AnnotationHub::query(ah,"OrgDb")
  Bioconductor_Orgs <- gsub(pattern = ".sqlite", replacement = "",x = all_OrgDb$title[grep(pattern = "db", x = all_OrgDb$title)])

  all_TxDb <- AnnotationHub::query(ah,"TxDb")
  Bioconductor_TxDbs <- gsub(pattern = ".sqlite", replacement = "" ,x = all_TxDb$title[grep(pattern = "UCSC", x = all_TxDb$title)])
  TxDb_Names <- gsub(pattern = "TxDb\\.",replacement = "",x = Bioconductor_TxDbs)
  #all_TxDb$description[grep(pattern = "UCSC", x = all_TxDb$title)]

  Bioconductor_Genomes <- BSgenome::available.genomes(splitNameParts=FALSE, type=getOption("pkgType"))
  BSgenome_Names <- gsub(pattern = "BSgenome.", replacement = "", x = Bioconductor_Genomes)

  GenomeOptions <- {}
  TxDbOptions <- list()
  Org_Annot_Options <- {}
  Available_Organisms <- {}
  Genome_and_TxDb_idx <- 0
  # Parse through each Bioconductor Genome entry and see which also have  matching TxDb
  for (i in 1:length(BSgenome_Names))
  {
    # Determine if there is a genome wide annotation database
    Org_Name <- unlist(strsplit(BSgenome_Names[i],"\\."))[1]
    Species_code <- substr(x = Org_Name,start = 1,stop = 2)
    Org_idx <- grep(pattern= paste("org\\.",Species_code,"\\.",sep=""),x=Bioconductor_Orgs) # USed to test for genome annotation
    TxDb_idx <- grep(pattern = BSgenome_Names[i],x = TxDb_Names)
    if ((length(TxDb_idx) > 0) && (length(Org_idx) > 0))
    { Genome_and_TxDb_idx <- Genome_and_TxDb_idx + 1
      GenomeOptions <-   c(GenomeOptions, BSgenome_Names[i])                # Record entry to display as an option for the user.
      TxDbOptions[[Genome_and_TxDb_idx]]  <- unique(TxDb_Names[TxDb_idx])   #  Record TxDb entries for this genome
      names(TxDbOptions)[Genome_and_TxDb_idx] <-  BSgenome_Names[i]
      Org_Annot_Options <- c(Org_Annot_Options, Bioconductor_Orgs[Org_idx])
      Available_Organisms <- c(Available_Organisms, Org_Name)
    }
  }

  # Construct download commands
  first_instruction <- 'source("http://bioconductor.org/biocLite.R")'
  final_instruction <- "To use this organism with Ularcirc you must download one BSgenome, one TxDb and one Org database"
  download_commands <- list()
  for(i in 1: length(Available_Organisms))
  { download_commands[[i]] <- paste('biocLite("', Org_Annot_Options[[i]],'") # this downloads organism annotations',sep = "")
    download_commands[[i]] <- c(download_commands[[i]], paste('biocLite("BSgenome.',
                                                            GenomeOptions[[i]],'") # this downloads organism genome',sep = ""))
    download_commands[[i]] <- c(download_commands[[i]], paste('biocLite("TxDb.',
                                                            TxDbOptions[[i]],'")  # This downloads a transcript database',sep=""))
    download_commands[[i]] <- c(first_instruction, download_commands[[i]], final_instruction)
  }
  names(download_commands) <- names(TxDbOptions)

  print(names(download_commands))
  invisible(download_commands)
}


