
###########################################################################
## DataPath
##    This simple function returns a directory name to where data is stored.
## There is probably room in the future to allow a more specific directory to
## to be set. Currently all data lives in
DataPath <- function()
{

  extdata_path <- system.file("extdata",package = "Ularcirc")
  if (extdata_path == "") {  # Can not find directory, assume running in outside package mode
    warning("Could not identify data path. May need to reinstall Ularcirc package")
    extdata_path <- "data"
  }
  else
  {
    temp <- read.table(file = paste(extdata_path,"extdata.txt", sep="/"),comment.char = "#",na.strings="NA", stringsAsFactors=FALSE, sep="=", row.names=1, strip.white=TRUE,header=FALSE)
    if (dir.exists(temp["DataPath",]))
    {
      extdata_path <- temp
    }
  }

  return(extdata_path)
}


####################################3
##  List_Genomes
## This function finds all installed libraries and return all genomes that have a linked TxDb
## Need to have a message for when nothing is installed
List_Genomes <- function()
{
  ularcircPath <- as.character(DataPath())
  p <- installed.packages()
  custom_sqlDb <- {}
  if (dir.exists(ularcircPath))
  {
    custom_sqlDb <- list.files(path=ularcircPath, pattern="TxDb.*.sqlite",recursive = TRUE, full.names = TRUE)
    names(custom_sqlDb) <- gsub(pattern="^.*/TxDb.",replacement = "",x = custom_sqlDb) # name of file without directory component
  }

  BSgenome_idx <- grep(pattern = "^BSgenome",x = p[,'Package'],ignore.case = FALSE )   # This identifies BSgenome libraries
  BSgenome_Names<- as.character(gsub(pattern = "BSgenome.", replacement = "", x = p[BSgenome_idx ,'Package'],ignore.case = FALSE))          # Now have list of genome IDs

  TxDB_idx <- grep(pattern = "^TxDb",x = p[,'Package'],ignore.case = FALSE)
  # First list installed databases
  TxDb_Names <- p[TxDB_idx ,'Package']
  names(TxDb_Names) <- as.character(gsub(pattern = "TxDb.", replacement = "", x = p[TxDB_idx ,'Package'],ignore.case = FALSE) )
  # Append on custom installed DB which reside in shiny circRNA directory
  if (length(custom_sqlDb))
  {   TxDb_Names <- c(TxDb_Names, custom_sqlDb)   }


  Org_Annotation_idx <- grep(pattern = "^org.", x = p[,'Package'], ignore.case=FALSE)
  Org_Annotation_Library<- as.character(p[Org_Annotation_idx,'Package'])

  GenomeOptions <- {}
  TxDbOptions <- list()
  Genome_and_TxDb_idx <- 0
  for (i in 1:length(BSgenome_idx))
  {
    if (length(grep(pattern = BSgenome_Names[i],x = TxDb_Names)) > 0)
    {   Genome_and_TxDb_idx <- Genome_and_TxDb_idx + 1
        GenomeOptions <-   c(GenomeOptions, BSgenome_Names[i])                                               # Record entry to display as an option for the user.
        TxDbOptions[[Genome_and_TxDb_idx]]  <- TxDb_Names[grep(pattern = BSgenome_Names[i],x = TxDb_Names)]  #  Record TxDb entried for this genome
        names(TxDbOptions)[Genome_and_TxDb_idx] <-  BSgenome_Names[i]
    }
    else  # no length, therefore no TxDb match
    {
    }
  }
  Org_Annotation_Library <- c(Org_Annotation_Library, "NO_ANNOTATION")
  return (list(GenomeOptions=GenomeOptions, TxDbOptions = TxDbOptions, Org_Annotation_Library = Org_Annotation_Library ))
}




##########################################################################################
## draw.arc
##
##  xc and yc are the X Y center positions of circle
##  r   is radius of arc
##  w1  defines starting position of arc.  270 starts at top of circle
##  w2 defines end position of arc relative to w1.
##      Ularcirc will draw full circle where w1 = 270 and end position is at 630 (i.e. 270 + 360)
##
## This code is taken from https://github.com/Bioconductor-mirror/OmicCircos/blob/release-3.5/R/OmicCircos.R
##
## lend where :
##      0 =  "round" mean rounded line caps [default];
##      1 "butt" mean butt line caps;
##      2 "square" mean square line caps.
##
##
## To use this need to create a plot then call function
## eg:
##
## plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
## draw.arc(400,400,400,270,400,lwd=6)
## draw.arc(400,400,360,270,350,lwd=4, Internal_Labels = c("miRNA"))
draw.arc <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1, draw_tick=FALSE, Internal_Labels = c(""))
{
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 2){
    pix.n <- 2;
  }

  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;


  if ((draw_tick) && (ang.d > 350))   # Make a spiral for ORFs longer than a circRNA length
  { r <-  seq(r * 0.9, r, length.out=pix.n)  }

  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;


  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);

  if (Internal_Labels[1] != "")
  { array_length <- length(fan.i.x)

    text_x <- xc + cos(ang.seq) * (r)
    text_y <- xc - sin(ang.seq) * (r)

    angle =  (360 - w1) - 10 # * 360 / pi * 2 * r
#    text(text_x[array_length], text_y[array_length], paste(Internal_Labels, angle), srt=angle, pos=2, adj = c(0,0)) # 270 down; 90  up; 0 normal ;   180 upside down and left)
    text(text_x[array_length], text_y[array_length], Internal_Labels, srt=angle, pos=2, adj = c(0,0)) # 270 down; 90  up; 0 normal ;   180 upside down and left)
  }

  ## New code added by Dave Humphreys
  if (draw_tick)   ## Draw tick
  {
    tick_length <- r / 10
    tick_x <- xc + cos(ang.seq) * (r-tick_length)
    tick_y <- yc - sin(ang.seq) * (r-tick_length)

    idx <- length(tick_x)  # Originally thought this would be the last pixel position, but actually is first.
    tick_x <- c(fan.i.x[1], tick_x[30])
    tick_y <- c(fan.i.y[1], tick_y[30])
    lines(tick_x, tick_y, col=col, lwd=lwd, type="l", lend=lend);
  }
}


Test_miRNA_circ <-function()
{
  ## source("u:/scripts/VCCRI_git/scripts/RStuff/Shiny/circRNA/Global.R")
  plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
  draw.arc(400,400,400,270,400,lwd=6)
  draw.arc(400,400,360,300,350,lwd=4, Internal_Labels = c("miRNA"))
  draw.arc(400,400,360,360,400,lwd=4, Internal_Labels = c("miRNA"))
  draw.arc(400,400,360,450,480,lwd=4, Internal_Labels = c("miRNA"))
  draw.arc(400,400,360,500,520,lwd=4, Internal_Labels = c("miRNA"))
  draw.arc(400,400,360,590,610,lwd=4, Internal_Labels = c("miRNA"))
}

##########################################################################################
## draw.text.w
##
##  xc and yc are the X Y center positions of circle
##  r   is radius of arc
##  w1  defines starting position of arc.  270 starts at top of circle
##  w2 defines end position of arc relative to w1.
##      Ularcirc will draw full circle where w1 = 270 and end position is at 630 (i.e. 270 + 360)
##
## This code is taken from https://github.com/Bioconductor-mirror/OmicCircos/blob/release-3.5/R/OmicCircos.R
##
## To use this need to create a plot then call function
## eg:
##
## plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
## draw.arc.test(400,400,400,270,400,lwd=4)
## draw.text.w(400,400,400,(270+400)/2, "circRNA")
draw.text.w <- function(xc, yc, r, w, n, col="black", cex=1){
  w <- w%%360;
  w <- w/360*2*pi;
  x <- xc+r*cos(w);
  y <- yc-r*sin(w);
  text(x,y,labels=n, col=col, cex=cex);
}

jsBOXcode <- " shinyjs.collapse = function(boxid)
    { $('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
    }"



List_Species_Files <- List_Genomes()

# The distance used to define boundaries of zoom sliderinput.
# i.e. Gene start/stop +- this distance
zoomOffset <- 50000

FileTypeCounts <- list(CIRI=0,CE2=0, STAR_BSJ=0, STAR_FSJ=0, QORTS=0, REGTOOLS=0)  # Keeps track of loaded data files
# Groupings <- list()


colour_Pallette <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
                       "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
                       "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                       "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3",
                       "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
                       "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
                       "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A",
                       "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                       "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                       "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                       "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
                       "#CCEBC5", "#FFED6F")
