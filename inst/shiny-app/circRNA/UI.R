library(shiny)
library(shinydashboard)
library(shinyjs)
library(Sushi)
library(data.table)
library(shinyFiles)
library(DT)


options(shiny.maxRequestSize=900*1024^2)  # Set upper limit at 900MB

List_Saved_Projects<-function()
{
  # list all RData files
  extdata_path <- as.character(DataPath())  # function from Global.R
  ProjectNames <- list.files(path= extdata_path,pattern = "*.RData",recursive = FALSE)
  return(gsub(pattern = ".RData", x=ProjectNames, replacement = ""))
}

# Define UI for dataset viewer application
shinyUI(
  fluidPage(
    useShinyjs(),
  titlePanel("UlarCirc : Analysing & Visualising Circular and linear RNA"),


  sidebarLayout(
    sidebarPanel(

		conditionalPanel('input.PanelSelect === "Setup"',
		    selectizeInput("Setup_Options",label="Setup configuration",
		                                choices = c('Load transcript database','Load new data', 'CircRNA education')),br(),br(),
		    conditionalPanel('input.Setup_Options == "Load transcript database"',
  				h4('ORGANISM',style="color:red"),
  				selectizeInput("Annotation_lib",label="Annotation library",choices = List_Species_Files$Org_Annotation_Library, selected="NO_ANNOTATION"),  # List_Species_Files$Org_Annotation_Library[1]),
	  			uiOutput("Display_Species_Options"),
		  		uiOutput("Display_Txdb_Options"),br(),
  				textOutput("List_Loaded_TxDB"),  ###  DISPLAY CURRENT LOADED TXDB
  				br()), # 		    conditionalPanel('input.Setup_Options == "Load transcript database"',

		    conditionalPanel('input.Setup_Options == "Load new data"',
  				h4('FILTER OPTIONS:',style="color:red"),
	  			checkboxInput('ChromosomeFilter', 'Same chromosomes:',TRUE),
		  		conditionalPanel(condition ="input.ChromosomeFilter == true",
		  	  sliderInput("GenomicDistance", "Chimeric genomic distance:",min = 10, max = 1000000, value = c(200,100000))),
			  	checkboxInput('StrandFilter', 'Same strand:',TRUE),
				  checkboxInput('CanonicalJuncs', 'DON\'T remove any canonical junctions',TRUE),

	        fileInput('JunctionFile', 'Chimeric junction File(s)',multiple=TRUE),   # accept=c('text/csv', 'text/comma-separated-values,text/plain','.csv'),
		      br()),   # conditionalPanel('input.Setup_Options == "Load new data"',

		    conditionalPanel('input.Setup_Options == "CircRNA education"',
		      h4('Sequencing  Parameters', style='color:orange'),
		      sliderInput("FragmentSize", "Select library fragment size", min=100, max=1000, value=300),
		      sliderInput("ReadLength", "Select sequencing read length", min=50, max=500, value=100),
		      radioButtons("ReadType", "Paired or single end:",choices = c("Paired","Single"), selected=c("Paired")),
		      sliderInput("ReadNumber", "Number of backsplice junction reads recovered", min=1, max=50, value=10),
		      h4('Circular RNA  Parameters', style='color:orange'),
		      sliderInput("CircRNA_Size", "Select length of circular RNA", min=300, max=10000, value=500),
		      br()),
				br()
			),
    conditionalPanel('input.PanelSelect === "Projects"',
        radioButtons("LibraryStrandType", "Library prep (TruSeq = opposing strand):",choices = c("Same Strand","Opposing strand", "Unstranded"), selected=c("Opposing strand")),
        selectizeInput("Load_or_Save",label="Load or save project",choices =  c("Load","Save")),
              conditionalPanel('input.Load_or_Save == "Load"',
                    h4('LOAD',style="color:blue"),
                    selectizeInput("LoadExistingProject",label="Choose an pre-existing project",choices =  List_Saved_Projects()),
                    actionButton("LoadProjectRequest","LOAD"),
                  #  sliderInput("Number_BiologicalSamples", "Number of biological treatments in data set:",min = 1, max = 10, value = 1),
                    uiOutput("DisplayGroupNames"),
                    br()
                      ),
              conditionalPanel('input.Load_or_Save == "Save"',
                    h4('SAVE',style="color:red"),
                    textInput('NewProject_Filename', 'Name of project'),
                    textAreaInput('ProjectNotes','Notes regarding project'),
                    actionButton("SaveProjectRequest","SAVE"),
                    br()
                      )
    ),
		conditionalPanel('input.PanelSelect == "Gene_View" && output.fileUploaded == true',
		    h4('DISPLAY MODE:',style="color:red"),
		    radioButtons("Display_Gene_View_Mode", "",
		    #selectizeInput("Display_Gene_View_Mode", "",
		                 choices = c("Display_Gene_Transcripts","Tabulated_Counts"), selected=c("Tabulated_Counts")), #Display_Gene_Transcripts")),

				conditionalPanel('input.Display_Gene_View_Mode == "Display_Gene_Transcripts"',
				    h4('Gene Display Options',style="color:red"),
				    uiOutput("Display_Gene_Zoom_Coords"),
				   # actionButton("Navigate_Around_Gene","Navigate"),
  				  checkboxInput('ShowTranscriptTable', 'Display transcript table:',FALSE),
				    checkboxInput('ShowBSJunctionCountTable', 'Display Backsplice junction count data',FALSE),
				    checkboxInput('ShowCanonicalCountTable', 'Display Canonical junction count data',FALSE),
				    radioButtons("JunctionType", "Junction Type:",choices = c("Backsplice","Alternative Canonical","All"), selected=c("Backsplice"))
				),
				conditionalPanel('input.Display_Gene_View_Mode == "Tabulated_Counts"',
				    #radioButtons("Annotation_Options",label="Choose how data should be tabulated",
				    h4('Table Display Options',style="color:red"),
				   # selectizeInput("Annotation_Options",label="Data analysis mode",
				    #    choices = c('Selected sample analysis','Grouped')),
				     radioButtons("Annotation_Options",label="Data analysis mode", choices = c('Selected','Grouped'), inline=TRUE),
				    uiOutput("TwoGroupCompareChoices"),             # This displays a selectizeInput menu for the possible group comparison combinations

				    conditionalPanel('input.BSJ_data_source =="STAR"',
				      checkboxInput('Percent_of_Parent', 'Display % parent transcript:',FALSE),
				      checkboxInput('Annotate_with_GeneName', "Annotate with parental gene:",TRUE),
				      checkboxInput('DisplayFilterOptions', "Display filter options:",FALSE)
				    ),
					    conditionalPanel('input.DisplayFilterOptions == true',
  				    checkboxInput('Display_RAD_Score', "Apply RAD filter:",TRUE),
	  			    sliderInput("RAD_filter", "Accepted RAD score", min=0, max=1, value=c(0,1)),
  				    numericInput("RAD_count_threshold", "Minimum count to apply RAD score ", value=9, min = 1, max = 50, step = 1),
  				    sliderInput("Apply_RAD_count", "Minimum count to apply RAD score ", min=0, max=1, value=c(0,1)),
			  	    checkboxInput('Apply_FSJ_Filter', "Apply FSJ filter:",TRUE),
				      sliderInput("FSJ_filter_count_range", "Range to apply FSJ support", min=0, max=50, value=c(0,10))
				      # Perhaps provide options to set min/max count and minimum sample number
				    ),
				    selectizeInput("MAX_BS_juncs_to_annotate", label= "Number of BS junctions to display",choices= as.numeric(c("5","10","20","35","50","75","100","250","500","1000","2000","5000","20000")), selected=10),
				    selectizeInput("Normalisation",label="Raw counts or CPM",choices = c("Raw counts","CPM", "CPM_Gene"), selected=c("Raw counts")),
  				  conditionalPanel('input.Annotation_Options == "Grouped"',
  				          selectizeInput("DisplayMode",label="Display mode",choices = c("Table", "Plots"), selected=c("Table"))
  				          ),
				    conditionalPanel('input.DisplayMode != "Plots"', actionButton("buildTable_Button", "Build table")),

				    conditionalPanel('input.DisplayMode == "Plots"',
				            selectizeInput("Global_Analysis_Plots_Options",label="Please select plot",
				                  choices = c("PCA","Heatmap","Unique Number of circRNAs","Genes producing circRNAs",
				                                  "Cummulative distribution"),
				                                    selected=c("Unique Number of circRNAs"))
				    ),
		        br()
		    ), # conditionalPanel('input.Display_Gene_View_Mode == "Tabulated_Counts"
				br()
			), #conditionalPanel
    conditionalPanel('input.PanelSelect == "Genome_View" && output.fileUploaded == true',
            h4('Display Options for Genome_view'),
            textInput('GenomeChrom_Input', 'Chromosome'),
            textInput('GenomeStart_Input', 'Start'),
            textInput('GenomeEnd_Input', 'End'),
            textInput('GenomeStrand_Input', 'Strand (+/-)'),
            actionButton("Update_Genome_Position","Navigate"),
            br(),
#            p("Slc8a1 in hg38 is captured between chr2:  40108472 - 40530305  "),
#            p("chr8:57303263-57324835"),p("chr5: 77000000 - 77100000"),p("MSTR.5271.1 ch13: 102349720 - 102359752 stand 2"),
#                 radioButtons("JunctionType", "Junction Type:",choices = c("Backsplice","Alternative Canonical","All"), selected=c("Backsplice")),
            checkboxInput('ShowGenomeCanonicalCountTable', 'Display forward canonical junction count data',FALSE),
            checkboxInput('ShowFSJ_Sequence', 'Display splice junction sequence',FALSE),
            br()
      ), #conditionalPanel

		conditionalPanel('input.PanelSelect === "Junction_View" && output.fileUploaded == true',
		    selectizeInput("Junction_View_Mode", "Select junction type to view",choices = c("Backsplice","Canonical"), selected=c("Backsplice")),
		    br(),
    		conditionalPanel('input.Junction_View_Mode == "Backsplice"',

    		    radioButtons('circRNA_Sequence_Analysis', 'Select analysis to perform',
    		                 choices = c("Display backsplice junction sequence","Open reading frame analysis","miRNA binding site analysis"),
    		                 selected=c("Display backsplice junction sequence")),

      		    #			  h5('Analyze flanking intron button'),
      		    #			  checkboxInput('Sequence_motif_analysis', 'Search for sequence motifs',TRUE),
    		    checkboxInput('circRNA_sequence_checkbox','Display circRNA sequence', FALSE),
    		    checkboxInput('Display_STAR_junction_data', 'Display raw data',FALSE),
    		    checkboxInput('Display_FAD', 'Display Distribution of reads across BSJ',FALSE),br(),
    		    sliderInput("FragSize", "Fragment size:",min = 100, max = 500, value = 300),
    		    sliderInput("ReadLength", "Read length:",min = 50, max = 300, value = 100),
    		    actionButton("PE_Fastq_Request", "Generate fastq file (paired end)"),br(),
  				  br()), # conditionalPanel('input.Junction_View_Mode == "Backsplice"',

  	  	conditionalPanel('input.Junction_View_Mode == "Canonical"',
		        h4('Display local gene map'),
		        h4('Display sequence - define squence length'),
		        br()),
		      br()

			)
    ),	# sidebarPanel


    # Show the caption, a summary of the dataset and an HTML
	 # table with the requested number of observations
    mainPanel(

		tabsetPanel(
			id = 'PanelSelect',
			tabPanel('Setup',
				conditionalPanel('input.Setup_Options == "Load transcript database"',
  				p('Welcome to Ularcirc!'),
  				p('If you are new to circular RNAs you may want to select "circRNA Education" from setup options configuration menu on left.',
              'Don\'t forget to come back here when done'),

  				p(h4('Instructions:'),' To get started follow the steps listed below:'),br(),
  				p(strong('STEP 1:'),'Load transcriptional database.',br(),
            'Select appropriate organism, genome and transcript database options and press load in side tab.',br(),
  				  'The database loaded will be listed below LOAD button when ready'),br(),

  				br(),
  				p(strong('STEP 2: '),'Load data.',br(),
              'This can either be an existing project which can be loaded under the Project tab.',br(),
  				    'Alternatively you can uploading new data. Select "Load new data" under setup option configuration on side menu.',
  				    'Select filter options and then click upload file button to load data.
                Once data is loaded you can navidate to project tab to save as a project.
  				      Note that reads that are filtered out are removed permanently'),
  				br(),
  				p(strong('STEP 3: '),'Search for circRNA.',br(),
  				  'Under Gene tab you can either navidate to your favourite gene or build tables of abundant circRNAs.',
  				  'Many of the tables under this tab can be used to select junctions of interest. Select junctions before proceeding to step 4.'),

  				br(),
  				p(strong('STEP 4: '),'Explore junction data.',br(),
  				  'After selecting junction(s) of interest navigate to the Junction tab which will provide a detailed report on the type of junction selected'),

  				br(),br(),
  				p('Keep an eye on the Ularcirc website for future updates and functionality'),
  				br()),   #conditionalPanel('input.Setup_Options == "Load transcript database"',

				conditionalPanel('input.Setup_Options == "Load new data"',

    			h4("Input file details:"),
				  tableOutput("FileNameDataTable"),
				  textOutput("FileNameDataTableDetails"),br(),
				  h4("Filtered junction details:"),
				  br()), # conditionalPanel('input.Setup_Options == "Load new data"',

				conditionalPanel('input.Setup_Options == "CircRNA education"',
				  h4("How are circRNA detected?"),
				  p('By identifying backsplice junctions (BSJ).
				      A BSJ is ultimately  a 2nt sequence which represent a donor and acceptor base from asyncronous exon(s) sequence.
				      Sequencing Reads from circRNA may not always include a BSJ. A type I read is indistinguisable from reads that align to linear RNA.
              Type II/III/IV reads are those that cap capture a BSJ.
				      The graph below demonstrates the longer a circRNA is the less chance of detecting a BSJ.'),
				  plotOutput("circRNA_Read_Distribution"),

				  p('Use the options on the side menu to theoretically estimate how many reads are in your data set for a particular circRNA.
            The table below predicts where the reads were assigned.
            Highlighted cells estimate the coverage of read types that have not been detected.
				    '),
				  DT::dataTableOutput("Predicted_Read_Distribution"),
				  p('* Detection of TypeIV reads may vary depending on the pipeline used and therefore what is displayed above may not be accurate'),
			#    uiOutput("Predicted_Read_Distribution"),
				  ## Table of expected Read type distributions

				  br()), # conditionalPanel('input.Setup_Options == "CircRNA education"',

				br()
				),		# tabPanel 'Setup'
			tabPanel('Projects',
			  h4('Project working directory'),
			  textInput(inputId="Project_WD", label="Directory where projects are saved to", value=as.character(DataPath()), width = '100%'),
			  actionButton("UpdateProjectWorkingDirectory","UPDATE"),
		#	  shinyDirButton("dir","Choose directory/Folder","Select Directory/Folder"),
			  tags$hr(style="border-color:black"),

        h4('"Selected Sample" data sets (sub sample analysis)'),p('Select samples for explatory analysis of a select number of samples'),
			  uiOutput("InputFiles"),		# This lists uploaded files in a checkboxInput format to shiny inputID of SelectedFiles
			  textOutput("Save_Load_Status"),hr(),tags$hr(style="border-color:black"),
	      h4('"Grouped analysis" data sets (whole project analysis)'),
        uiOutput("DisplayGroupings"),	# Allow user to allocate grouping to allow simple differential analysis or trend analysis
        br()
			  ),
			tabPanel('Gene_View',
				conditionalPanel(condition = "output.fileUploaded == false",
					br(),br(),
					h4('No uploaded file detected, please wait or go back to "Setup" tab and load a data set',style="color:red"),
					br()),

				conditionalPanel(condition = "output.fileUploaded == true",
				  # Display option of what datasets to display

				  checkboxInput("DataSourceOptions", "Show data source options", TRUE),
#          checkboxGroupInput("DataSourceOptions", choices=c("Show data source options","Show sample IDs"), choiceValues=TRUE, inLine=TRUE),
				  fluidRow(
				    shinydashboard::box(id="box1", background="red",  solidHeader = F, collapsible = F,
				      uiOutput("DisplayDataSetButtons")
				      ),
				    verbatimTextOutput("ShowDataSets_on_GeneView")
				  ),

				  conditionalPanel('input.Display_Gene_View_Mode == "Display_Gene_Transcripts"',
				      selectInput("GeneListDisplay", NULL, choices = NULL),
				      actionButton(inputId = "Update_Gene_of_Interest",label = "Select Gene"),
					   # uiOutput("DisplayDataSetButtons"),  # To display Ularcirc | circExplorer or other input data sets
					    plotOutput("distPlot"),

					    conditionalPanel(condition = "input.ShowTranscriptTable == true ",
					      h4('Transcript Table'),
					      downloadButton('download_Transcript_Table','Download Transcript table' ),
					      DT::dataTableOutput("TranscriptTable"),
					      h5('Exon Table (populated once a row is selected from transcript table)'),
						    DT::dataTableOutput("ExonTable")
		              ), # conditionalPanel(condition = "output.ShowExonTable == true ",

					    conditionalPanel(condition = "input.ShowBSJunctionCountTable == true ",
						    h5('Backsplice Junction table') ,
						    downloadButton('download_BS_Junc_Count_Table','Download BSJ' ),
						    DT::dataTableOutput("BS_Junction_Count_Table")
					       ), # conditionalPanel(condition = "input.ShowBSJunctionCountTable == true "

				 		  conditionalPanel(condition = "input.ShowCanonicalCountTable == true",
			  		    h5('Canonical Junction table'),
			  		    downloadButton('download_FSJ_Junc_Count_Table','Download FSJ' ),
						    DT::dataTableOutput("CanonicalJunctionCountTable")
		              ), # conditionalPanel(condition = "ShowCanonicalCountTable == true",

					    br() ), #conditionalPanel('input.Display_Gene_View_Mode == "Display_Gene_Transcripts"',


				  # Following conditionalPanel is set up based on previous menu settings shown below
				  # selectizeInput("Annotation_Options",label="Data sets to analyse",
				  #           choices = c('Selected','Grouped analysis'))
				  conditionalPanel('input.Display_Gene_View_Mode == "Tabulated_Counts"',
				      conditionalPanel('input.Annotation_Options == "Selected"',
    					  #h4('Junction table of selected data sets'),
    					  uiOutput("BSJ_count_table_header"),
    					  conditionalPanel('input.BSJ_data_source =="STAR"',
      					  downloadButton('downloadSelectJunctionCountTable','Download' ),
	  				      DT::dataTableOutput("DisplayJunctionCountTable")
      					  ),
    					  conditionalPanel('input.BSJ_data_source !="STAR"',
    					     DT::dataTableOutput("Display_externalBSJ_CountTable")
    					     ),
					      br() ),
				      conditionalPanel('input.Annotation_Options == "Grouped"',
				        # Following conditional Panel is set up based on previous menu setting shown below:
				        # selectizeInput("DisplayMode",label="Display mode",
				        #             choices = c("Table","Heat map", "Volcano plot"), selected=c("Table"))
				        conditionalPanel('input.DisplayMode == "Table"',
				          h4('Junction table of grouped data sets'),
				          conditionalPanel('input.BSJ_data_source =="STAR"',
  				          downloadButton('downloadGroupedJunctionCountTable','Download' ),
	  			          DT::dataTableOutput("DisplayGroupJunctionCountTable")
				          ),
				          conditionalPanel('input.BSJ_data_source !="STAR"',
                  #  downloadButton('downloadGroupedJunctionCountTable','Download' ),
				            DT::dataTableOutput("Display_externalBSJ_GroupCountTable")
				          ),
				          br() ),
				        conditionalPanel('input.DisplayMode == "Heat map"',
				          h4('Heat map of grouped data sets'),
				          plotOutput("DisplayGroupHeatMap"),
				          br()),
				        conditionalPanel('input.DisplayMode == "PCA"',
				          h4('Principle components analysis (PCA)'),
				          plotOutput("DisplayBSJ_PCA"),
				        #  DT::dataTableOutput("DisplayJunctionCountTable"),
				          br() ),
				        conditionalPanel('input.Global_Analysis_Plots_Options == "Heatmap"',
				          selectInput("HeatmapGeneNumber",label="Number of variable genes",choices = seq(from=10,to=50,by=10),selected = 10)
				        ),
				        conditionalPanel('input.DisplayMode == "Plots"',
				          plotOutput("Global_Analysis_Plots")
				          ),
				        br() ), # conditionalPanel('input.Annotation_Options == "Grouped"',
					   br() ), # conditionalPanel('input.Display_Gene_View_Mode == "Tabulated_Counts"',
					br() ) # conditionalPanel(condition = "output.fileUploaded == true",

				),		# tabPanel 'Gene_View'

	    tabPanel('Genome_View',
	         conditionalPanel(condition = "output.fileUploaded == false",
	                          br(),br(),
	                          h4('No uploaded file detected, please wait or go back to "Setup" tab and load a data set',style="color:red"),
	                          br()),

	         conditionalPanel(condition = "output.fileUploaded == true",

	                          verbatimTextOutput("ShowDataSets_on_Genome_View"),

	                          plotOutput("genomePlot"),

	                          conditionalPanel(condition = "input.ShowExonTable == true ",
	                                           h5('Exon Table') #,
	                                           #DT::dataTableOutput("ExonTable")  # Need to make genome version
	                          ),

	                          conditionalPanel(condition = "input.ShowGeneJunctionTable == true ",
	                                           h5('Junction table')
	                                           #DT::dataTableOutput("JunctionTable") # Need to make genome version
	                          ),
	                          conditionalPanel(condition = "input.ShowGenomeCanonicalCountTable == true",
	                                           h5('Canonical Junction table'),
	                                           DT::dataTableOutput("GenomeCanonicalJunctionCountTable")
	                                           ), # conditionalPanel(condition = "ShowGenomeCanonicalCountTable == true",
	                          conditionalPanel(condition = "input.ShowFSJ_Sequence == true",
	                                           h5('Junction sequence'),
	                                           uiOutput("Predicted_Genomic_Junction_Sequence")
	                                            ),

	                          br())



	    ),		# tabPanel 'Genome_View'


			tabPanel('Junction_View',
				conditionalPanel(condition = "output.fileUploaded == false",
					br(),br(),
					h3('No file loaded, go back to "Setup" tab and load a data set',style="color:red"),
					br()),

				conditionalPanel(condition = "output.fileUploaexport(UTR_3p[Multi_UTR_3p_idx], 'c:/temp/multi_3p.gtf')
ded == true",
				  verbatimTextOutput("ShowDataSets_on_JunctionView"),
				  conditionalPanel('input.Junction_View_Mode == "Backsplice"',
  					uiOutput("DisplayBS_sequence_details"),

  					conditionalPanel(condition = 'input.circRNA_Sequence_Analysis == "Display backsplice junction sequence"',
  	  					uiOutput("DisplayBS_sequence"),br()),

  					conditionalPanel(condition = "input.circRNA_sequence_checkbox == true",
  					   uiOutput("Predicted_circRNA_Sequence"), br(),
  					   br()),   # conditionalPanel(condition = "input.circRNA_sequence_checkbox == true",
  					conditionalPanel(condition = "input.Display_STAR_junction_data == true",
  					   h5("Raw junction data from STAR aligner"),
  					   DT::dataTableOutput("JunctionTableOther"),
  					   br()), # conditionalPanel(condition = "input.Display_STAR_junction_data == true",

  					conditionalPanel(condition = 'input.circRNA_Sequence_Analysis == "Open reading frame analysis"',
  					    plotOutput("circRNA_Sequence_Analysis_ORF"),
  					    DT::dataTableOutput("circRNA_Sequence_Analysis_Table"),
  					    br()),

  					conditionalPanel(condition = 'input.circRNA_Sequence_Analysis == "miRNA binding site analysis"',
  					     uiOutput("miRNA_Options"),   # Seed length,  how many duplicates
  					     plotOutput("circRNA_Sequence_Analysis_miRNA"),
  			#		     DT::dataTableOutput("circRNA_Sequence_Analysis_Table"),
  					     br()),

  			conditionalPanel(condition = 'input.Display_FAD == true',
  			        plotOutput("Plot_RAD_Histogram"),
  			        br()),
  			#		uiOutput("circRNA_sequence_analysis"),
  					br(),

					  br()),
	  		  br()),
				  conditionalPanel('input.Junction_View_Mode == "Canonical"',
				  #  uiOutput("DisplayCanonical_sequence"),br(),br(),         #renderUI
				    verbatimTextOutput("DisplayCanonical_sequence"),br(),br(),   #renderText
				    br()),
				  br()
				) 		# tabPanel 'Junction_View'

			)	# tabsetPanel



    ) # main panel
  ) # sidebarLayout
))
