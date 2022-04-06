###############################################################################
#' plot_AllJunctions
#'
#' Plots a BSJ, FSJ and transcripts for a nominated gene. Output is combined 
#' onto a single page. This function effectively wraps plotting functions 
#' from plotgardener
#' 
#' @param assembly : Genome assembly
#' @param chrom : chromosome
#' @param chromstart : Starting position of chromosome
#' @param chromend : End position of chromosome
#' @param BSJData : Backsplice junction data table
#' @param BSJ_colors : Backsplice junction assigned colours
#' @param FSJData : Forward junction data table
#' @param FSJ_colors : Forward junction assigned colours
#' @param geneSymbol : Gene symbol
#' @return Returns a list of two DNAstring sets labelled "read1" and "read2" which correspond
#' to forward and reverse read pairs.
#'
#' @examples
#'
#' library('Ularcirc')
#'  # BSJ data.table
#'  BSJ_data <- data.table::data.table(chrom1="chr2", 
#'                 start1=c(40139400, 40160764, 40428472, 40428472), 
#'                 end1=c(40139400, 40160764,40428472, 40428472), 
#'                 chrom2="chr2", start2=c(40178494,40178494,40430302,40430305),
#'                  end2=c(40178494,40178494,40430302,40430305),
#'                score=c(13,20,360,1751))
#'
#'  # FSJ
#'  FSJstarts1 <- c(40115630,40139677,40160865,40164985,40170350,40174721,
#'                               40174843,40175282,40278771,40430302,40430305)
#'  FSJstarts2 <- c(40139400,40160764,40164853,40170280,40174705,40174824,
#'                                40175260,40178386,40428472,40453160,40512348)
#'  FSJ_data <- data.table::data.table(chrom1="chr2", start1=FSJstarts1, end1=FSJstarts1,
#'                      chrom2="chr2", start2=FSJstarts2, end2=FSJstarts2, 
#'                      score=c(225,825,685,666,633,596,517,542,685,101,171))
#'                      
#'  plot_AllJunctions(assembly="hg38", chrom="chr2", 
#'                      chromstart=40096769, chromend=40611554,
#'                      BSJData=BSJ_data, FSJData=FSJ_data, geneSymbol="SLC8A1") 
#'
#' @export
plot_AllJunctions <- function(assembly = "hg38", chrom, chromstart, chromend,
                         BSJData,BSJ_colors='black',
                         FSJData, FSJ_colors='black',
                         geneSymbol)
{
  # Check  parameters
  if(length(BSJ_colors) == 0)
  {
    BSJ_colors <- 'black'
  }
  if(length(FSJ_colors) == 0)
  {
    FSJ_colors <- 'black'
  }
  
  
  ## Create a page
  plotgardener::pageCreate(width = 7.5, height = 6, default.units = "inches", showGuides = FALSE)
 
  
  ## Set the coordinates
  params <- plotgardener::pgParams(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = "hg38",
    width = 7
  )

   
  # Set label coordinates
  yAxis.xpos <- chromstart - ((chromend - chromstart)* 0.05)
  
  ## Plot BSJ
  BSJarchPlot <- plotgardener::plotPairsArches(
    data = BSJData, params = params,
    fill = BSJ_colors,
    #  fill = colorby("length", palette =  colorRampPalette(c("dodgerblue2", "firebrick2"))),
    linecolor = "fill",
    archHeight = BSJData$score,  alpha = 1,
    x = 1, y = 0.5, height = 1.5,
    just = c("left", "top"),
    baseline=TRUE,
    #  bg="white",
    default.units = "inches"
  )
  
  plotgardener::annoYaxis(
    #    plot = archPlot, at = c(0, max(FSJ_junctions$score)),
    plot = BSJarchPlot, at = c(0, max(BSJData$score)), axisLine=TRUE,
    fontsize = 10
  )
  
  plotgardener::annoText(
    label = "Depth", fontsize = 8, plot = BSJarchPlot,
    x = yAxis.xpos, y = max(BSJData$score)/2, rot=90,
    just = "center", default.units = "native"
  )
  plotgardener::plotText(label="Backsplice junctions",x = 3.5, y = 0.4, just = "center", default.units = "inches")
  
  
  ## Plot FSJ
  FSJarchPlot <- plotgardener::plotPairsArches(
    data = FSJData, params = params,
    fill = FSJ_colors,
    #  fill = colorby("length", palette =  colorRampPalette(c("dodgerblue2", "firebrick2"))),
    linecolor = "fill",
    archHeight = FSJData$score,  alpha = 1,
    x = 1, y = 2.25, height = 1.5,
    just = c("left", "top"),
    baseline=TRUE,
    #  bg="white",
    default.units = "inches"
  )
  
  plotgardener::annoYaxis(
    #    plot = archPlot, at = c(0, max(FSJ_junctions$score)),
    plot = FSJarchPlot, at = c(0, max(FSJData$score)), axisLine=TRUE,
    fontsize = 10
  )
  
  plotgardener::annoText(
    label = "Depth", fontsize = 8, plot = FSJarchPlot,
    x = yAxis.xpos, y = max(FSJData$score)/2, rot=90,
    just = "center", default.units = "native"
  )
  
  plotgardener::plotText(label="Linear junctions",x = 3.5, y = 2.1, just = "center", default.units = "inches")
  
  
  transcriptPlot <- plotgardener::plotTranscripts(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = "hg38", labels = "transcript", fontsize = 4,
    x = 1, y = 3.5, width = 7, height = 2,
 #   spaceHeight = 0.1,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Annotate genome label
  plotgardener::annoGenomeLabel(
    plot = transcriptPlot, x = 1, y = 5.5, scale = "Mb",
    just = c("left", "top")
  )  
  
  ## Plot a legend
  plotgardener::plotLegend(
    legend = c("+ strand", "- strand"),
    fill = c("#669fd9", "#abcc8e"), border = FALSE,
    x = 0.5, y = 5.7, width = 0.5, height = 0.25,
    just = c("left", "top")
  )
}