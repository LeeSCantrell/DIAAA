rm(list = ls()) 
#Thank you for using DIAAA (Data-Indpendent Acquisition Analysis App) - If you have any suggestions, e-mail lee.s.cantrell@vanderbilt.edu
#To run app, open in RStudio (https://www.rstudio.com/products/rstudio/download/) and click "Run App" in top right. 
#To open app in browser, select "Open in New Window" in the viewer window of RStudio
#Packages that require GitHub installs include DREAM-AI and DIA-NN R package. Code for install commented below default parameter block

############################################################################################################################################
#Default Parameters
#Set where your data is located, defaults to system root
workingDirectory <- '~/Desktop/'
#Set Min Required Peptides Per Protein (default 2)
minPeptides <- 2
#Set Min Required Quantitative Fragments per Peptide (default 3, EncyclopeDIA only)
minFragments <- 3
#Standard Contaminants List from MaxQuant is annotated below. Edit as appropriate 
############################################################################################################################################

#GitHub Package Installs - remove # before each line, then select hilight code block, select command+enter (Mac) or control+enter (Windows)
#install.packages("devtools")
#library(devtools)
#install_github("https://github.com/vdemichev/diann-rpackage")
#install_github("WangLab-MSSM/DreamAI/Code")

#package load
require(shiny)
require(shinyWidgets)
require(shinyFiles)
require(diann)
require(vroom)
require(shinybusy)
require(NOISeq)
require(hutils)
require(ggplot2)
require(reshape)
require(viridis)
require(varhandle)
require(ggpubr)
require(plotly)
require(GGally)
require(cluster)
require(factoextra)
require(ggrepel)
require(data.table)
require(DreamAI)

#default contaminant list for greedy data filtering
contaminants <- c("P00761","Q32MB2","P19013","Q7RTT2","P15636","P09870","Q9R4J5",
                  "P0C1U8","P00766","P13717","Q9U6Y5","P21578","O76009","O76011",
                  "O76013","O76014","O76015","P08779","Q14525","Q14532","Q15323",
                  "Q92764","Q14533","Q9NSB4","P78385","Q9NSB2","P78386","O43790",
                  "Q6IFU5","Q9UE12","Q8IUT8","Q6NT21","Q6ISB0","Q6NTB9","Q6IFU6",
                  "P04264","P13647","P35908","P13645","P35527","A3EZ79","P02533",
                  "P02538","P48668","P04259","A3EZ82","Q2KIG3","Q0VCM5","Q3SZ57",
                  "Q9N2I2","Q3SZH5","P28800","Q1A7A4","P41361","Q2YDI2","Q3Y5Z3",
                  "P81644","Q2KJ83","Q2KIT0","A2I7N3","Q3SZV7","Q2KJC7","Q3SZR3",
                  "Q28107","P02672","Q1RMN8","Q58D62","P06868","Q2KJF1","P02584",
                  "P02777","Q3SX14","P17697","Q6T181","P34955","P21752","Q32PJ2",
                  "Q28194","P00978","Q5XQN5","Q32PI4","Q9TTE1","Q2KIU3","P67983",
                  "Q28065","Q862S4","Q2KIF2","Q3SX28","Q0V8M9","Q148H6","Q29RQ1",
                  "Q95M17","P07224","Q2HJF0","Q2KIH2","Q04695","A2I7N0","P12763",
                  "P17690","P02769","P02676","P50448","P01030","P01966","P00735",
                  "Q03247","Q3ZBS7","Q2UVX4","Q9TT36","Q28085","Q3SX09","Q3ZBD7",
                  "Q3MHN2","Q9TRI1","P15497","Q95121","Q05443","P02070","Q2KIS7",
                  "Q3MHH8","Q3T052","Q3KUS7","Q1RMK2","Q2TBQ1","Q05B55","A2I7N1",
                  "P04258","Q2KJ62","Q0IIK2","Q3MHN5","P02662","P02663","P02666",
                  "P02668","P31096","P02754","P00711","P62894","Q29443","P19001",
                  "A2AB72","Q8VED5","Q61726","Q3ZAW8","P50446","Q497I4","Q9D312",
                  "Q922U2","Q8BGZ7","A2A4G1","Q9QWL7","Q6IME9","Q6NXH9","A2VCT4",
                  "P07744","Q6IFZ6","Q6IFX2","Q9R0H5","Q3TTY5","Q0VBK2","Q61782",
                  "A2A5Y0","Q99PS0","Q9D646","P05784","Q9DCV7","Q9Z2K1","P07477",
                  "P05787","Q7Z794","Q9BYR9","Q9BYQ5","Q9BYR8","Q9BYQ7","Q3LI72",
                  "Q9BYR4","Q9BYQ8","P60413","P19012","Q2M2I5","O95678","Q01546",
                  "Q99456","Q9H552","P35900","Q3SY84","Q8N1A0","Q5XKE5","P12035",
                  "Q9C075","P08729","Q7Z3Y8","Q7RTS7","Q7Z3Y9","Q7Z3Z0","Q7Z3Y7",
                  "P08727","Q3KNV1","Q86YZ3","P20930","Q5D862")

#Design Matrix Initiate - do not edit
designColumnNames <- c("Run.Id", "Sample.Id", "NumericVar1", "NumericVar2", 
                       "NumericVar3", "NumericVar4", "CategoricVar5",
                       "CategoricVar6", "CategoricVar7", "CategoricVar8",
                       "CategoricVar9", "CategoricVar10")



############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
#User Interface (HTML generator)
ui <- fluidPage(
  setBackgroundColor('GhostWhite'),
  h1(HTML(paste0("DIA Analysis App (DIA",tags$sup("3"), ")"))),
  "Written 9/27/22 - lee.s.cantrell@vanderbilt.edu",
  add_busy_spinner(spin = "orbit", color="steelblue"),
  navlistPanel(
  ############################################################################################################################################
    tabPanel("Data Load",
             "This analysis is protein-centric for multi-sample bioanalytical pipelines.
             For the evaluation of single peptides and post translational modifications,
             alternative tools may be used, or custom R scripts generated. For evaluation
             of EncyclopeDIA data, upload only the peptide level inference set.",
             
             br(),
             
             br(),
             
             'Please select which search algorithm the data was acquired from prior to upload to the server.',
             
             hr(),
             
             fluidRow(
               column(2, radioButtons(inputId = 'algorithm', label = "Search Algorithm Used", c("DIA-NN", "EncyclopeDIA", 'MSFragger'))),
               column(2, shinyFilesButton('rawDataFile', label='DIA Report File', title='Please select a DIA report file', multiple=FALSE)),
               column(2, actionButton(inputId = 'startUpload', 'Start Data Processing')),
               column(6,verbatimTextOutput('rawDataFilePath'),
               )
               
              ), 
             
             fluidRow(
             ),
             
             hr(),
             
             "All data is searched through a design matrix format reference database. 
             A single excel file can be automatically generated in the directory of the report file entered above. All run IDs will be loaded
             in one column. Suggested sample names can be added to the sample ID column, but may be changed. 
             Up to 10 additional variables can be entered in columns specifying traits that can be explored by multiple statistical methods 
             following data loading. Column names may be changed as long as renamed columns remain quantative/continuous and categorical/discrete.", 
             
             hr(),
             
             fluidRow(
               column(2, actionButton("autoFillDesign", "AutoFill Exp. Design",  width = '100%')),
               column(2, shinyFilesButton('designModified', label = 'Filled Design .csv', title = 'Please select the Design report file', multiple = FALSE)), 
               column(2, actionButton(inputId = 'designUpload', 'Start Design Processing')), 
               column(6, verbatimTextOutput('modifiedDesignFilePath'))
             ),
             
  ), 
  ############################################################################################################################################
  tabPanel("Normalization &\nExporatory Data Analysis", 
           "Normalization is performed in customizable steps by combination of Global Sum Normalization, Trimmed Mean of M-Values (TMM), and/or
           Median Normalization of non-log2 transformed data in respective sequence.
           It is recommended to utilize only TMM.",
           br(),
           hr(),

           fluidRow(
             column(2, h5("Global Normalization"), switchInput(inputId = 'globalSum', label = '', value = F)), 
             column(2, h5("TMM"), switchInput(inputId = 'TMMSlider', label = "", value = T)),
             column(2, h5("Median Normalization"), switchInput(inputId = 'medianSlider', label = "", value = F)),
             column(4, h5("Data for EDA"), selectInput(inputId = 'EDAdatachoice', label = NULL, choices = list("Normalized Data" = 2, "Raw Data" = 1))),
             column(2, h5("Display Log2 Treated Data"), switchInput(inputId = 'normLog2', label = "", value = T ))
           ),
           hr(),
           fluidRow(
             uiOutput('uiBoxplot')
           ),
           hr(), 
           fluidRow(
             selectInput(inputId = 'correlationMethod', label = "Correlation Method", choices = list("Pearson" = 1, "Spearman" = 2))
           ),
           fluidRow(
             plotOutput('normCorrelation', height = '600px')
           ),
           hr(),
           fluidRow(
             plotOutput('normCDF')
           ),
           hr(),
           fluidRow(
             uiOutput('uiSum')
           )
           ),
  ############################################################################################################################################
  tabPanel("Imputation",
             "Imputation is not advised for DIA data, and should only be used as a last resort. In case needed, 
             methods from the CPTAC DREAM-AI suite (https://github.com/WangLab-MSSM/DreamAI) are made available below. 
             Imputation may fail depending on data completeness and desired level of imputation. Imputation is computationally
             intensive and may take a few minutes to process. Empirical observation shows that Ensemble method outperforms single
             methods with values missing at random and values missing not at random.",
             hr(),
             column(12, radioButtons(inputId = 'imputeMethod', label = "Selected Imputation Method",
                                                                       choices = list("K-Nearest Neighbors" = 1, "Missing Forest" = 2, 
                                                                                      "Abundance Dependent Missing Imputation (ADMIN)" = 3,
                                                                                      "IRNN-SCAD (Birnn)" = 4, 
                                                                                      "Matrix Factorization (SpectroFM)" =5,
                                                                                      "GLM Ridge Regression (RegImpute)" = 6,
                                                                                      "Ensemble (Average of all)" = 7), selected = 7)),
           fluidRow(sliderInput(inputId = 'imputeFraction', label = "Percent Missing Values to Impute", min = 0, max = 50, value = 25)),
           fluidRow(
             column(2, h5("Use Imputation (Not Recommended)"), switchInput(inputId = 'imputeSlider', label = '', value = F)))
           ),
  ############################################################################################################################################
  tabPanel("Dimensional Reduction & Clustering",
           sidebarLayout(
             sidebarPanel(
               selectInput(inputId = 'pcaCategorical', label = 'Color By', choices = c("K-Means Clustering")), 
               br(),
               sliderInput(inputId <- "kmeansClusters", label = "Number of K-Means Clusters", value = 1, min = 1, max = 2, round = T, step=1)
             ),
             mainPanel(
               tabsetPanel(type = 'tabs',
                           tabPanel("PCA", plotOutput("pcaPlot")),
                           tabPanel("Scree Plot", plotOutput("screePlot")),
                           tabPanel("Optimal Clusters", plotOutput('optimalClustersPlot')),
                           tabPanel("PCA - Interactive", plotlyOutput("pcaPlotInteractive"))
                           )
             )
           )
           ),
  ############################################################################################################################################
  tabPanel("Differential Expression", 
           tabsetPanel(type = "tabs", 
                       tabPanel("Volcano Plot",
                                sidebarPanel(
                                  "Set T-Test Parameters",
                                  selectInput(inputId = "tTestCategorical", label = 'Comparison Variable', choices = c(), selected = NULL), 
                                  selectInput(inputId = 'tTest1', label = "Comparison Group 1", choices = c(), selected = NULL),
                                  selectInput(inputId = 'tTest2', label = "Comparison Group 2", choices = c(), selected = NULL), 
                                  selectInput(inputId = 'tMethod', label = "T-Test Method", 
                                              choices = c("Welch's T-Test" = 0, "T-Test" = 1, "Limma - Bayes Moderated" = 2 )),
                                  sliderInput(inputId <- 'pvalueSlider', label ="Select -Log10 P-Value Threshold", min = 1, max = 10, value = 2, step = .25), 
                                  sliderInput(inputId <- 'log2Slider', label = "Select Log2FC Threshold", min = 0, max = 10, value = 2, step = .25)
                                ), 
                                mainPanel(
                                  plotOutput("volcanoPlot", height = '780px')
                                ),
                                hr(),
                                dataTableOutput("volcanoPlotOutput")
                                ),
                       tabPanel("Interactive Volcano Plot", 
                                plotlyOutput("interactiveVolcanoPlot", height = '800px'))
           )
                       
                       
                           
           ),
  ############################################################################################################################################
  tabPanel("Single Protein Analysis")
  )
)


############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
#Server (Programming)
server <- function(input,output, session){
  volumes = getVolumes()
  masterData <- reactiveValues(rawData = 0)
  
  ############################################################################################################################################
  #Data Load Tab
  #Report Input
  observe({
    setwd(workingDirectory)
    shinyFileChoose(input, "rawDataFile", roots = volumes, defaultPath = getwd(), filetypes=c('', 'txt', 'csv', 'tsv'))
    file_selected <- parseFilePaths(volumes, input$rawDataFile)
    output$rawDataFilePath <- renderText(as.character(file_selected$datapath))
    masterData$filePath <- renderText(as.character(file_selected$datapath))
  }) # load file names into UI
  file_selected <- reactive(parseFilePaths(volumes, input$rawDataFile))
  filePath <- reactive(as.character(file_selected()$datapath)) #save filepath for dataset load
  observeEvent(input$startUpload, {
    rawData <- data.frame(fread(filePath()))
    searchAlgo <- input$algorithm
    if(searchAlgo == "EncyclopeDIA"){
      for(i in 1:length(contaminants)){
        pb <- txtProgressBar(1,length(contaminants), style = 3)
        rawData <- rawData[!grepl(contaminants[i], rawData$Protein, fixed = T),]
        setTxtProgressBar(pb, i)
      }
      masterData$rawData <- rawData
      rownames(rawData) <- rawData[,1]
      rawData <- rawData[rawData[,2] > 1,]
      rawData <- rawData[,-c(1,2,3)]
      masterData$rawProteins <- rawData
    }
    if(searchAlgo == "DIA-NN"){
            rawData <- data.frame(rawData)
      raw.names <- stringr::str_split(rawData$Protein.Ids, ';')
      raw.names <- lapply(raw.names, `length<-`, max(lengths(raw.names)))
      raw.names <- (matrix(data = unlist(raw.names), nrow = nrow(raw.data), ncol = length(raw.names[[1]]), byrow = T))
      raw.names2 <- matrix(raw.names %in% contaminants, nrow <- nrow(raw.names), ncol = ncol(raw.names), byrow = F)
      raw.names <- data.frame(raw.names)
      raw.names$sum <- rowSums(raw.names2)
      raw.names$filter <- raw.names$sum > 0
      dat.filter <- which(raw.names$filter == TRUE) 
      rawData <- rawData[-dat.filter,]
      masterData$rawData <- rawData
      masterData$rawProteins <- diann_maxlfq(rawData[rawData$Q.Value <= 0.01  & rawData$Global.PG.Q.Value <= 0.01,], 
                                                             group.header="Genes", id.header = "Precursor.Id", 
                                                             quantity.header = "Precursor.Quantity", sample.header = 'Run')
    }
  }) # Generate protein matrix
  fileDir <- reactive(dirname(filePath()))
  
  #MatrixUpload
  observeEvent(input$autoFillDesign, {
    rawFile <- colnames(masterData$rawProteins)
    run.ids <- colnames(masterData$rawProteins)
    designGeneric <- matrix(data = NA, nrow = length(run.ids), ncol = length(designColumnNames))
    designGeneric <- data.frame(designGeneric)
    colnames(designGeneric) <- designColumnNames
    designGeneric$Run.Id <- run.ids
    designGeneric$Sample.Id <- trim_common_affixes(run.ids, warn_if_no_prefix = F, warn_if_no_suffix = F)
    masterData$genericDesign <- designGeneric
    setwd(fileDir())
    write.csv(masterData$genericDesign, file = 'GenericDesign.csv', row.names = F)
  }) #create generic design matrix 
  observe({
    shinyFileChoose(input, "designModified", roots = volumes, defaultPath = getwd(), filetypes=c( 'csv'))
    file_selected <- parseFilePaths(volumes, input$designModified)
    output$modifiedDesignFilePath <- renderText(as.character(file_selected$datapath))
  }) # load complete design into UI
  design_selected <- reactive(parseFilePaths(volumes, input$designModified))
  designPath <- reactive(as.character(design_selected()$datapath)) #save filepath for design load
  observeEvent(input$designUpload, {
    design <- (diann_load(designPath()))
    for(i in 3:6){
      design[,i] <- as.numeric(design[,i])
    }
    for(i in 7:length(designColumnNames)){
      design[,i] <- as.factor(design[,i])
    }
    comparable = c()
    for(i in 1:length(designColumnNames)){
      if(length(unique(design[,i])) > 1){a <- TRUE}
      else{a <- FALSE}
      comparable[i] <- a
    }
    design <- design[,comparable]
    masterData$design <- design
    if(sum(comparable) >2 ){
      masterData$comparables <- colnames(design)[3:sum(comparable)]
    }
    namesRawProteins <- colnames(masterData$rawProteins)
    namesDesign <- design$Run.Id
    reorderColumnsIndex <- match(namesDesign, namesRawProteins)
    masterData$rawProteins <- masterData$rawProteins[,reorderColumnsIndex]
    colnames(masterData$rawProteins) <- design[, 2]
  })
  output$data <- renderTable(masterData$rawProteins)
  
  ############################################################################################################################################
  #Normalization, Imputation, and Exploratory Data Analysis Tab
  
  temp <- observe({
    if(!is.null(masterData$rawProteins)){
      temp <- data.frame(masterData$rawProteins)
      if(input$EDAdatachoice == 2 & input$imputeSlider == F){
        if(input$globalSum == T){
          tempSum <- temp
          columnsums <- base::colSums(tempSum, na.rm = T)
          mincolumn <- min(columnsums)
          columnnormfactor <- mincolumn/columnsums
          tempSum <- mapply('*', tempSum, columnnormfactor)
          temp <- data.frame(tempSum)
        }
        if(input$TMMSlider == T){
          tempNoMissingValues <- temp[which(rowSums(is.na(temp))==0), ]
          tempTMM <- tmm(tempNoMissingValues)
          TMMRatios <- (tempTMM[1,]/tempNoMissingValues[1,])
          temp <- mapply('*', temp, TMMRatios)
          temp <- data.frame(temp)
        }
        if(input$medianSlider == T){
          tempNoMissingValues <- temp[which(rowSums(is.na(temp))==0), ]
          medians <- apply(temp, 2, median, na.rm = T)
          medianRatios <- medians[1]/medians
          temp <- mapply('*', temp, medianRatios)
          temp <- data.frame(temp)
        }
        masterData$normalizedProteins <- temp
        masterData$workingProteins <- temp
      }
      if(input$EDAdatachoice == 1 & input$imputeSlider == F){
        masterData$workingProteins <- temp
      }
      if(input$EDAdatachoice == 1 & input$imputeSlider == T){
        if(input$imputeMethod == 1){imputeMethod = 'KNN'}
        if(input$imputeMethod == 2){imputeMethod = "MissForest"}
        if(input$imputeMethod == 3){imputeMethod = "ADMIN"}
        if(input$imputeMethod == 4){imputeMethod = "Birnn"}
        if(input$imputeMethod == 5){imputeMethod = "SpectroFM"}
        if(input$imputeMethod == 6){imputeMethod = "RegImpute"}
        if(input$imputeMethod == 7){imputeMethod = c("KNN","MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute")}
        allowedMissing <- floor(ncol(temp)*.01*input$imputeFraction)
        tempExclude <- temp[rowSums(is.na(temp)) > allowedMissing, ]
        tempImpute <- temp[rowSums(is.na(temp)) <= allowedMissing, ]
        if(sum(is.na(tempImpute)) > 0){
          tempImputed <- DreamAI(tempImpute, k = 10, maxiter_MF = 10, ntree = 100,
                                 maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
                                 gamma_ADMIN = NA, gamma = 50, CV = FALSE,
                                 fillmethod = "row_mean", maxiter_RegImpute = 10,
                                 conv_nrmse = 1e-06, iter_SpectroFM = 40, 
                                 method = imputeMethod,
                                 out = c("Ensemble"))
          temp <- data.frame(rbind(tempImputed$Ensemble, tempExclude))
        }
        masterData$rawProteinsImpute <- temp
        masterData$workingProteins <- temp
      }
      if(input$EDAdatachoice == 2 & input$imputeSlider == T){
        if(input$globalSum == T){
          tempSum <- temp
          columnsums <- base::colSums(tempSum, na.rm = T)
          mincolumn <- min(columnsums)
          columnnormfactor <- mincolumn/columnsums
          tempSum <- mapply('*', tempSum, columnnormfactor)
          temp <- data.frame(tempSum)
        }
        if(input$TMMSlider == T){
          tempNoMissingValues <- temp[which(rowSums(is.na(temp))==0), ]
          tempTMM <- tmm(tempNoMissingValues)
          TMMRatios <- (tempTMM[1,]/tempNoMissingValues[1,])
          temp <- mapply('*', temp, TMMRatios)
          temp <- data.frame(temp)
        }
        if(input$medianSlider == T){
          tempNoMissingValues <- temp[which(rowSums(is.na(temp))==0), ]
          medians <- apply(temp, 2, median, na.rm = T)
          medianRatios <- medians[1]/medians
          temp <- mapply('*', temp, medianRatios)
          temp <- data.frame(temp)
        }
        if(input$imputeMethod == 1){imputeMethod = 'KNN'}
        if(input$imputeMethod == 2){imputeMethod = "MissForest"}
        if(input$imputeMethod == 3){imputeMethod = "ADMIN"}
        if(input$imputeMethod == 4){imputeMethod = "Birnn"}
        if(input$imputeMethod == 5){imputeMethod = "SpectroFM"}
        if(input$imputeMethod == 6){imputeMethod = "RegImpute"}
        if(input$imputeMethod == 7){imputeMethod = c("KNN","MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute")}
        allowedMissing <- floor(ncol(temp)*.01*input$imputeFraction)
        tempExclude <- temp[rowSums(is.na(temp)) > allowedMissing, ]
        tempImpute <- temp[rowSums(is.na(temp)) <= allowedMissing, ]
        if(sum(is.na(tempImpute)) > 0){
          tempImputed <- DreamAI(tempImpute, k = 10, maxiter_MF = 10, ntree = 100,
                                 maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
                                 gamma_ADMIN = NA, gamma = 50, CV = FALSE,
                                 fillmethod = "row_mean", maxiter_RegImpute = 10,
                                 conv_nrmse = 1e-06, iter_SpectroFM = 40, 
                                 method = imputeMethod,
                                 out = c("Ensemble"))
          temp <- data.frame(rbind(tempImputed$Ensemble, tempExclude))
        }
        masterData$normalizedProteinsImpute <- temp
        masterData$workingProteins <- temp
      }
            if(input$normLog2 == T){
        transform <- function(x){log2(x)}
      }
      if(input$normLog2 == F){
        transform <- function(x){x}
      }
      tempViolin <- data.frame(melt(data.frame(temp)))
      tempViolin$variable <- unfactor(tempViolin$variable)
      tempViolin$value <- as.numeric(tempViolin$value)
      colnames(tempViolin) <- c('SampleID', 'Intensity')
      b <- ggplot(tempViolin, aes(SampleID, transform(Intensity), fill = SampleID)) + geom_boxplot(alpha = 0.7) + 
        theme_bw(base_size = 14) + 
        scale_x_discrete("Sample ID") + scale_y_continuous("Log2 Intensity") + 
        ggtitle("Protein Expression Per Sample") + scale_fill_viridis(discrete = T) +
        theme(axis.text.x = element_text(vjust = 0.5, hjust=1)) + 
        theme(legend.position = 'none') + coord_flip()
      if(input$normLog2 == F){b <- b + scale_y_continuous("Intensity")}
      
      runs <- length(temp)
      output$uiBoxplot <- renderUI({
        output$normBoxplot <- renderPlot(b)
        pixels <- 100 + (25*runs)
        plotOutput("normBoxplot", width = '100%', height = pixels)
      })
      
      
      
      tempCorr <- data.frame(temp)
      c <- ggcorr(data = tempCorr)
      if(input$correlationMethod == 2){
        c <- ggcorr(data = tempCorr, method = c("pairwise", "spearman"))
      }
      compares <- c$data["coefficient"]
      compares <- c(compares$coefficient)
      minimum <- min(compares)
      diffCompares <- 1-minimum
      c <- ggcorr(data = tempCorr, limits = c(minimum, 1), midpoint = 1-(diffCompares/2), 
                  label = T, label_round = 2, geom = 'tile', name = "Correlation", method = c('pairwise', 'pearson'))
      if(input$correlationMethod == 2){
        c <- ggcorr(data = tempCorr, limits = c(minimum, 1), midpoint = 1-(diffCompares/2), 
                    label = T, label_round = 2, geom = 'tile', name = "Correlation", method = c('pairwise', 'spearman'))
      }
      output$normCorrelation <- renderPlot(c)
      
      tempCDF <- temp
      tempCDF <- rowSums(is.na(tempCDF))
      colCount <- ncol(temp)
      quantNA <- c()
      for(i in 1:(colCount + 1)){
        quantNA[i] <- sum(grepl(paste0("^", i-1, "$"), tempCDF))
      }
      cdfDataFrame <- data.frame("Missing.Values" = 0:colCount, "QuantProteins" = quantNA)
      cdfQuantNA <- c()
      for(i in 1:length(cdfDataFrame$QuantProteins)){
        cdfQuantNA[i] <- sum(cdfDataFrame$QuantProteins[1:i])
      }
      cdfDataFrame$CDF <- cdfQuantNA
      cdfDataFrame <- cdfDataFrame[-nrow(cdfDataFrame),]
      d  <- ggplot(cdfDataFrame, aes(Missing.Values, CDF)) + geom_col() + theme_classic(base_size = 14) + 
        scale_x_continuous("Missing Values") + scale_y_continuous("Proteins Identified")
      output$normCDF <- renderPlot(d)
      tempSum <- temp
      sums <- colSums(tempSum, na.rm = T)
      sumFrame <- data.frame('file' <- colnames(tempSum), 'sum' <- sums)
      colnames(sumFrame) <- c("file", "sum")
      s <- ggplot(sumFrame, aes(x = file, y = sum)) + geom_bar(stat= 'identity') + theme_classic(base_size = 14) + 
        scale_x_discrete("Sample ID") + scale_y_continuous("Summed Intensity") + coord_flip()

      output$uiSum <- renderUI({
        output$normSummedSignal <- renderPlot(s)
        pixels <- 100 + (25*runs)
        plotOutput("normSummedSignal", width = '100%', height = pixels)
      })
      
      
      
      
    }  
    
    })
  
  ############################################################################################################################################
  #Dimensional Reduction & Clustering
  
  temp <- observe({
    if(!is.null(masterData$rawProteins)){
      design <- masterData$design
      nk <- ncol(data.frame(masterData$rawProteins))
      updateSliderInput(session, inputId = "kmeansClusters", max = nk-1)
      temp <- data.frame(masterData$workingProteins)
      temp <- temp[rowSums(is.na(temp))==0,]
      temp <- data.frame(t(temp))
      pc.out <-prcomp((log2(temp)))
      pc.data <- data.frame(pc.out$x)
      pc.var <- (pc.out$sdev)**2
      pc.var <- c(pc.var/sum(pc.var))
      k2 <- as.factor(kmeans(temp, centers = input$kmeansClusters, nstart = 25)$cluster)
      comparables = c("K-Means Clustering")
      if(length(masterData$comparables) > 0){
        comparables <- c(comparables, masterData$comparables)
      }
      updateSelectInput(session, "pcaCategorical", choices = comparables, selected = input$pcaCategorical)
      category <- input$pcaCategorical
      if(category == "K-Means Clustering"){
        categoryOut <- k2
        pc.data$Categorical <- categoryOut
        pc.data$Sample <- colnames(masterData$rawProteins)
        p <- ggplot(pc.data, aes(PC1, PC2, col = Categorical, sample = Sample)) + geom_point(shape = 10, size = 5) + theme_bw(base_size = 14) + 
          scale_x_continuous(paste0("PC1 ", round(pc.var[1]*100,0), '% PVE')) +
          scale_y_continuous(paste0('PC2 ', round(pc.var[2]*100,0), '% PVE')) + 
          labs(col = "K-Means Cluster")
        
      }
      if(category != "K-Means Clustering"){
        categoryOut <- design[,category]
        pc.data$Categorical <- categoryOut
        pc.data$Sample <- colnames(masterData$rawProteins)
        p <- ggplot(pc.data, aes(PC1, PC2, col = Categorical, sample = Sample)) + geom_point(shape = 10, size = 5) + theme_bw(base_size = 14) + 
          scale_x_continuous(paste0("PC1 ", round(pc.var[1]*100,0), '% PVE')) +
          scale_y_continuous(paste0('PC2 ', round(pc.var[2]*100,0), '% PVE')) + 
          labs(col = category)
      }
      
      output$pcaPlot <- renderPlot(p)
      output$pcaPlotInteractive <- renderPlotly(ggplotly(p))
      
      xlabs <- c()
      labels <- for(i in 1:length(pc.var)){xlabs[i] <- paste0("PC", i)}
      screeData <- data.frame(1:length(pc.var), pc.var)
      colnames(screeData) <- c("PC", "PVE")
      scree <- ggplot(screeData, aes(PC, PVE)) + geom_bar(stat = 'identity') + scale_x_discrete("Principal Component") + scale_y_continuous("Proportion of Variance Explained") + 
        ggtitle("Scree Plot") + theme_bw(base_size = 14)
      output$screePlot <- renderPlot(scree)
      if(nrow(temp)>2){
        kco <- fviz_nbclust(temp, kmeans, k.max = nrow(temp) - 1, method = 'silhouette')
        output$optimalClustersPlot <- renderPlot(kco)
      }
    }
    
  })
  
  
  ############################################################################################################################################
  #Differential Expression/Volcano Plot
  temp <- observe({
    if(!is.null(masterData$rawProteins) & length(masterData$comparables) > 0){
      if(sum(unlist(lapply(masterData$design, is.factor))) > 0){
        design <- masterData$design
        comparableT <- which(unlist(lapply(design, is.factor)))
        updateSelectInput(session, "tTestCategorical", choices = names(comparableT), selected = input$tTestCategorical)
        if(is.null(input$tTestCategorical)){updateSelectInput(session, "tTestCategorical", selected = comparableT[1])}
        if(nchar(input$tTestCategorical) > 0){
          comparableTselect <- input$tTestCategorical
          comparableLevels <- c(levels(design[,comparableTselect]))
          updateSelectInput(session, "tTest1", choices = comparableLevels, selected = input$tTest1)
          updateSelectInput(session, "tTest2", choices = comparableLevels, selected = input$tTest2)
          if(nchar(input$tTest1) > 0){
            tTest1cols <- (which(grepl(input$tTest1, design[,comparableTselect])))
            tTest2cols <- (which(grepl(input$tTest2, design[,comparableTselect])))
            if(length(tTest1cols) > 1 & length(tTest1cols) > 1){
              tTestInput <- masterData$workingProteins
              rownames(tTestInput) <- rownames(masterData$rawProteins)
              rowKeep <- (rowSums(!is.na(tTestInput[,tTest1cols])) > 1) + (rowSums(!is.na(tTestInput[,tTest2cols])) > 1)
              rowKeep <- as.numeric(which(rowKeep == 2))
              tTestInput <- tTestInput[rowKeep, ]
              tempResult <-  tTestInput
              tTestInput <- as.matrix(tTestInput)
              p.value <- c()
              log2FC <- c()
              for(i in 1:nrow(tTestInput)){
                test1dat <- as.numeric(unlist(na.omit(tTestInput[i,tTest1cols])))
                test2dat <- as.numeric(unlist(na.omit(tTestInput[i,tTest2cols])))
                if(input$tMethod == 0){varequal = FALSE}
                if(input$tMethod == 1){varequal = TRUE}

                t.testTemp <- t.test(x = test1dat, y = test2dat, var.equal = varequal)
                p.value[i] <- t.testTemp$p.value
                log2FC[i] <- log2((mean(c(tTestInput[i,tTest1cols]), na.rm = T))/mean(c(tTestInput[i,tTest2cols]), na.rm = T))
              }
              tempResult$p.value <- p.value
              tempResult$log2FC <- log2FC
              tempResult$Significant <- FALSE
              logThreshold <- input$log2Slider
              p.threshold <- input$pvalueSlider
              
              
              tempResult$Significant[intersect(which(log2FC > logThreshold | log2FC < -logThreshold), which(-p.threshold > log10(p.value)))] <- TRUE
              tempResult$Protein.Group <- rownames(tempResult)
              
              
              v <- ggplot(tempResult, aes(log2FC, -log10(p.value), label = Protein.Group)) + 
                geom_point(data = tempResult, aes(col = Significant), shape = 1, size = 2) + 
                theme_bw() + geom_vline(xintercept = logThreshold,  linetype = 'dotted') +
                geom_vline(xintercept = -logThreshold,  linetype = 'dotted') +
                geom_hline(yintercept = p.threshold, linetype = 'dotted') + 
                scale_color_manual(values = c('grey', 'red')) + 
                theme(legend.position = 'none') + 
                scale_x_continuous(paste0("Log2FC\n", nrow(tempResult), ' Protein Groups')) +
                scale_y_continuous('-Log10 P-Value') + ggtitle(paste0(input$tTest1, '/', input$tTest2)) + 
                geom_text_repel(aes(label=ifelse(Significant== TRUE, as.character(Protein.Group),'')))
              
              output$volcanoPlot <- renderPlot(v)
              
              output$interactiveVolcanoPlot <- renderPlotly(ggplotly(v))
              
              tableOutput <- data.frame(tempResult$Protein.Group, tempResult$p.value, tempResult$log2FC, tempResult$Significant, tempResult[,1:(ncol(tempResult)-4)])
              colnames(tableOutput)[1:4] <- c("Protein Group", "p-value", 'Log2 FC', "Significant")
              
              output$volcanoPlotOutput <- renderDataTable(tableOutput)
              
              
              
              
            }
              
               
          }
              
        }
            
      }
    }

  })
  
  
  
  
  
}



#Execute App
shinyApp(ui = ui, server = server)

