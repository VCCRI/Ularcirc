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
#' Interogates Bioconductor databases and identifies those that are compatible
#' with Ularcirc. Builds a list of commands that the user can copy to install the
#' required database on their local computer. Once installed the databases are
#' immediately available to Ularcirc upon re-starting the shiny app.
#' This function requires connection to the internet.
#' @param list_commands Boolean. By default this is FALSE and will return a dataframe of all
#' compatible databases. When set to TRUE will return a list of commands that can be copied
#' onto console which will download the appropriate databases (see example below).
#' @return Returns a list of compatible annotation databases
#' @examples
#' # Get all Bioconductor annotation databases that are compatible with Ularcirc
#' library('BSgenome')
#' library('httpuv')
#' library('AnnotationHub')
#' # Prepare a dataframe of all compatible annotation databases
#' compatible_DBs <- Compatible_Annotation_DBs()
#'
#' cat(paste(compatible_DBs[1,],collapse=","))
#' @export
Compatible_Annotation_DBs <- function(species ='')
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
  for (i in seq_along(BSgenome_Names))
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
  compatible_databases <- data.frame()
  for(i in seq_along(Available_Organisms))
  {
    download_commands[[i]] <- paste('biocLite("', Org_Annot_Options[[i]],'") # this downloads organism annotations',sep = "")
    download_commands[[i]] <- c(download_commands[[i]], paste('biocLite("BSgenome.',
                                                            GenomeOptions[[i]],'") # this downloads organism genome',sep = ""))
    download_commands[[i]] <- c(download_commands[[i]], paste('biocLite("TxDb.',
                                                            TxDbOptions[[i]],'")  # This downloads a transcript database',sep=""))
    download_commands[[i]] <- c(first_instruction, download_commands[[i]], final_instruction)

    compatible_entry <- data.frame(annotation=Org_Annot_Options[[i]],
                                       genome=paste("BSgenome.",GenomeOptions[[i]],sep=""),
                                       txdb=paste("TxDb.",TxDbOptions[[i]],sep="") )

    if (i == 1)
    {   compatible_databases <- compatible_entry  }
    else
    {   compatible_databases <- rbind(compatible_databases, compatible_entry)  }
  }
  names(download_commands) <- names(TxDbOptions)

 # print(names(download_commands))
  invisible(compatible_databases)
}




## Ularcirc_Annotation_Commands <- function(compatible_databases, species = "None")
##{
##    idx <- grep(pattern = "Hsapiens",x = test[,"genome"])
##    short_list <- compatible_databases[idx,]
##
##
##  download_commands[[i]] <- paste('biocLite("', Org_Annot_Options[[i]],'") # this downloads organism annotations',sep = "")
##  download_commands[[i]] <- c(download_commands[[i]], paste('biocLite("BSgenome.',
##                                                            GenomeOptions[[i]],'") # this downloads organism genome',sep = ""))
##  download_commands[[i]] <- c(download_commands[[i]], paste('biocLite("TxDb.',
##                                                            TxDbOptions[[i]],'")  # This downloads a transcript database',sep=""))
##  download_commands[[i]] <- c(first_instruction, download_commands[[i]], final_instruction)
##}
