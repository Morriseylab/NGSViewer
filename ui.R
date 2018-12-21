library(shinydashboard)
#library(shinyIncubator)
library(shiny)
library(shinyBS)
library(plotly)
library(d3heatmap)
library(shinyjs)
library(rglwidget)
library(SPIA)
options(shiny.sanitize.errors = FALSE)
ui <- dashboardPage(
  dashboardHeader(title = "NGS Data Analysis Web Interface",titleWidth = 350),
  dashboardSidebar(width = 350,
                   div(style="overflow-y: scroll"),
                   tags$style(type="text/css",
                              ".shiny-output-error { visibility: hidden; }",
                              ".shiny-output-error:before { visibility: hidden; }"
                   ),
                   tags$head(tags$style(HTML(".sidebar { height: 250vh; overflow-y: auto; }" ))),
                   
                   sidebarMenu(
                   menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
                   uiOutput("projects"),
                   fluidRow(
                     column(12,uiOutput("contrasts"))
                   ),
                     menuItem('PCA-Plot', tabName = 'pcaplot', icon = icon('hand-o-right'), 
                              menuSubItem("PCA Plot", tabName = "pcaplot"),
                              menuSubItem('Display Variances', tabName = 'var'),
                              menuSubItem('Show 3D plot', tabName = '3dplot')),
                     menuItem('Project Summary and Results', tabName = 'summres', icon = icon('hand-o-right'), 
                              menuSubItem('Results', tabName = 'geneselection'),
                              menuSubItem('View volcano plot', tabName = 'volcanoplot'),
                              menuSubItem('View Limma results of Multiple Contrasts', tabName = 'multilimma'),
                              menuSubItem(icon=NULL,checkboxInput("check", label = "Display Contrast List", value = FALSE)))
                   ),#end of sidebar menu
                   conditionalPanel(
                     condition = "input.check ==true",
                     uiOutput("contrastslimma")
                   ),
                   sidebarMenu(
                     menuItem('View Raw Expression Data', tabName = 'voom', icon = icon('hand-o-right')),
                     menuItem('View Sample Data', tabName = 'phenofile', icon = icon('hand-o-right')),
                     menuItem('Generate Heatmap', tabName = 'heatmap', icon = icon('hand-o-right')),
                     menuItem('GSEA Using Camera', tabName = 'cam', icon = icon('hand-o-right'),
                              menuSubItem('View Camera results', tabName = 'camera'),
                              menuSubItem(icon=NULL,uiOutput("cameradd"))),
                     menuItem('Pathway Analysis using SPIA', tabName = 'spia', icon = icon('hand-o-right')),
                     menuItem('Pathway Analysis using Reactome', tabName = 'reactometab', icon = icon('hand-o-right'),
                              menuSubItem('View ReactomePA results', tabName = 'reactome'),
                              menuSubItem('View CNET plot', tabName = 'cnet'),
                              menuSubItem('View EMAP plot', tabName = 'emap'),
                              menuSubItem('Gene Set Enrichment Analysis', tabName = 'gsea')),
                     menuItem('GO Analysis using GAGE', icon = icon('hand-o-right'),
                              menuSubItem('GAGE Results', tabName = 'gogage')
                     )
                   ),#end of sidebar menu
                   sidebarMenu(
                     menuItem("Help Page", tabName = "help", icon = icon("hand-o-right")))
  ),#end dashboardSidebar
  
  ######################################################################################################################################
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
    tabItems(
      tabItem(tabName = "dashboard",
              box(
                width = 12, status = "primary",solidHeader = TRUE,
                title = "RNA-Seq Data Sets",
                tableOutput("dashdata")
              )
      ),
      ######################################################################################################################################
      tabItem(tabName = "pcaplot",
              box(
                width = 10, status = "primary",solidHeader = TRUE,
                title = "Select Options",
                fluidRow(
                  column(6,uiOutput("pcaxoptions")),
                  column(6,uiOutput("pcayoptions"))
                ),
                br(),textOutput("biplottitle"),br(),
                fluidRow(
                  column(6,uiOutput("pcipslide")),
                  column(6,uiOutput("pcslide"))
                ),br(),
                uiOutput("ellipse")),
              box(
                width = 10, status = "primary",solidHeader = TRUE,
                title = "PCA Plot",
                plotOutput("biplot",height=600)
              ),br(),
              fluidRow(
                column(6,uiOutput("dwldbiplot")))
      ),
      ######################################################################################################################################
      tabItem(tabName = "var",
              box(width = 10, status = "primary",solidHeader = TRUE,title = "Variances of principal Components",
                  textOutput("pcatitle"),plotOutput("pcaplot_ip",width=700,height=400),br()),
              box(width = 12, status = "primary",solidHeader = TRUE,title = "Amount of Variation explained by each Principle Component",DT::dataTableOutput('pcaplot_tab'))),
      
      ######################################################################################################################################
      tabItem(tabName = "3dplot",h4("3D plot"),br(),br(),rglwidgetOutput("pcaplot3d",width = "850px", height = "750px")),
      
      ######################################################################################################################################
      tabItem(tabName = "geneselection",
              box(width = 10, status = "primary",solidHeader = TRUE,title = "Project Description",
                  textOutput("pdesc")),
              
       #Dot plot box
              fluidRow(
                box(width = 8, status = "primary",solidHeader = TRUE,title = "Dot Plot of the gene of interest",
                    fluidRow(
                      column(6,uiOutput("boxplotcol")),
                      column(6,uiOutput("boxplotcol2"))
                    ),
                    uiOutput("minexprline"),
                    plotOutput('dotplot')
                ),
       #Gene selection control box
                box(width = 4, status = "primary",solidHeader = TRUE,title = "Gene Selection",
                    radioButtons("radio", label = h4("Gene Selection"),
                                 choices = c("None" = 'none',"Upregulated" = 'up', "Downregulated" = 'down', "Both" = 'both'),
                                 selected = 'none'),
                    
                    sliderInput("lfc", label = h4("Fold Change"), min = 0.5,max = 6, value = 2),
                    sliderInput("apval", label = h4("P. Value"), min = 0.01,max = 0.2, value =0.05),br(),
                    fluidRow(
                      column(6,downloadButton('dwld','Download results table')),
                      column(6,downloadButton('downloaddotplot', 'Download Dot plot')))
                )),
       #LImma data table       
              box(width = 12, status = "primary",solidHeader = TRUE,title = "Limma data",
                  h5(p(div(span("Note:fc - Fold Change",style="color:red")))),
                  br(),textOutput("contrdesc"),br(),DT::dataTableOutput('table'))),
      
      ######################################################################################################################################
      tabItem(tabName = "volcanoplot",
              box(width = 8, status = "primary",solidHeader = TRUE,title = "Volcano Plot",
                  plotlyOutput("volcanoplot",height=700)),
              box(width = 4, status = "primary",solidHeader = TRUE,title = "Controls",
                  uiOutput("volcdrop"),br(),uiOutput("volcslider"),br(),
                  downloadButton('dwldvolcanoplot', 'Download Volcano plot')),
                br(),DT::dataTableOutput('table_volc')),
      ######################################################################################################################################
      
      tabItem(tabName = "multilimma",DT::dataTableOutput('table_TRUE'),fluidRow(uiOutput("dwldmultitab"))),
      
      ######################################################################################################################################
      
      tabItem(tabName = "voom",DT::dataTableOutput('table3'),
              fluidRow(uiOutput("dwldrawtab"))
      ),
      ######################################################################################################################################
      
      tabItem(tabName = "phenofile",DT::dataTableOutput('phenofile')),
      ######################################################################################################################################
      
      tabItem(tabName = "heatmap",
              box(width = 8, status = "primary",solidHeader = TRUE,title = "Heatmap",
                  textOutput("htitle"),br(),
                  fluidRow(
                    column(6,h3("")),
                    column(width = 3, offset = 2,plotOutput('hmpscale_out',width = 200,height = 65))
                  ),
                  fluidRow(
                    column(6,checkboxInput("hmpsamp", label = "View Heatmap of all samples", value = TRUE)),
                    column(6,h4(""))
                  ),
                  d3heatmapOutput('heatmap',width=550,height=900)
              ),
              
              box(width = 4, status = "primary",solidHeader = TRUE,title = "Controls",
                  selectInput("hmip", "Select Heatmap input type",c('Top number of genes' = "genenum",'Enter Genelist' = "geneli", 'Top Variable genes' = "vargenes")),
                  selectInput("hmpcol", "Select Heatmap Color Palette",c('YlGnBu' = "YlGnBu",'RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'PRGn'="PRGn", 'Blues' = "Blues")),
                  selectInput("clusterby", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none")),
                  checkboxInput("checkbox", label = "Reverse Colors", value = FALSE),
                  
                  conditionalPanel(
                    condition = "input.hmip == 'genenum'",
                    uiOutput("dropdown"),
                    sliderInput("gene", label = h4("Top number of genes"), min = 2,max = 500, value = 50)
                  ),
                  conditionalPanel(
                    condition = "input.hmip == 'geneli'",
                    selectInput(inputId = 'selectidentifier',label='Select Identifier',choices=list('Ensembl ID'='ensembl','Entrez Id'='entrez','Gene Symbol'='genesym')),
                    fileInput('genelistfile', 'Upload Text File',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                    actionButton(inputId = 'ga', label = 'Display Results')
                  ),
                  conditionalPanel(
                    condition = "input.hmip == 'vargenes'",
                    #uiOutput("dropdown"),
                    sliderInput("vgene", label = h4("Top variable genes"), min = 2,max = 500, value = 50)
                  ),br(),
                  downloadButton('downloadheatmap', 'Download')
              )#end box
      ),#end tabItem
      
      ######################################################################################################################################
      tabItem(tabName = "camera",
              fluidRow(
                box(width = 8, status = "primary",solidHeader = TRUE,title = "Camera Heatmap",
                    fluidRow(
                      column(6,uiOutput('hmplimcam')),
                      column(width = 3, offset = 2,plotOutput('hmpscale_out2',width = 200,height = 65))
                    ),
                    fluidRow(
                      column(6,checkboxInput("hmpsamp2", label = "View Heatmap of all samples", value = TRUE)),
                      column(6,h4(""))
                    ),
                    d3heatmapOutput('camheatmap',height=900)),
               
                 box(width = 4, status = "primary",solidHeader = TRUE,title = "Controls",
                    selectInput("hmpcol2", "Select Heatmap Color Palette",c('YlGnBu' = "YlGnBu",'RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'PRGn'="PRGn", 'Blues' = "Blues")),
                    selectInput("clusterby2", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none")),
                    checkboxInput("checkbox2", label = "Reverse Colors", value = FALSE),
                    br(),
                    downloadButton('downloadcamheatmap', 'Download Heatmap'))),
              
              box(width = 12, status = "primary",solidHeader = TRUE,title = "Table",
                  DT::dataTableOutput('tablecam'),textOutput("camdesc"),DT::dataTableOutput('campick3'))
   ),#end tabItem
      ######################################################################################################################################
      tabItem(tabName = "spia", 
              fluidRow(
                box(width = 8, status = "primary",solidHeader = TRUE,title = "SPIA Heatmap",
                    fluidRow(
                      column(6,uiOutput('hmplimspia')),
                      column(width = 3, offset = 2,plotOutput('hmpscale_out2spia',width = 200,height = 65))
                    ),
                    fluidRow(
                      column(6,checkboxInput("hmpsamp2spia", label = "View Heatmap of all samples", value = TRUE)),
                      column(6,h4(""))
                    ),
                    d3heatmapOutput('spiaheatmap',width=550,height=900)),
                box(width = 4, status = "primary",solidHeader = TRUE,title = "Controls",
                    selectInput("hmpcolspia", "Select Heatmap Color Palette",c('YlGnBu' = "YlGnBu",'RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'PRGn'="PRGn", 'Blues' = "Blues")),
                    selectInput("clusterbyspia", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none")),
                    checkboxInput("checkboxspia", label = "Reverse Colors", value = FALSE),
                    br(),
                    downloadButton('downloadspiaheatmap', 'Download Heatmap'),
                    hr(),
                    downloadButton('dwldspia', 'Download SPIA Results'))),
              
              box(width = 12, status = "primary",solidHeader = TRUE,title = "SPIA Result",
                  DT::dataTableOutput('spiaop'),textOutput("spiadesc"),DT::dataTableOutput('spiagenes'))
      ),
     ######################################################################################################################################
   tabItem(tabName = "reactome", 
           box(width = 12, status = "primary",solidHeader = TRUE,title = "Visualize enrichment result",
               fluidRow(
               column(6,radioButtons("enrichradio", label = h4("Plot Type"),choices = c("Barplot" = 'barplot',"Dotplot" = 'dotplot',"Enrichment Map"='enrich'),selected = 'barplot')),
               column(6,sliderInput("ncat", "Number of categories to show:",min = 5, max = 20, value = 10,step=1))
               ),
               plotOutput("enrichplot", height = 700)),
           box(width = 12, status = "primary",solidHeader = TRUE,title = "Enriched Pathways",
           DT::dataTableOutput('enrichpath')),
           box(width = 12, status = "primary",solidHeader = TRUE,title = "Genes in pathway",
               verbatimTextOutput('enrichgenes'))
   ),
   ######################################################################################################################################
   tabItem(tabName = "cnet", 
           box(width = 12, status = "primary",solidHeader = TRUE,title = "Visualize enrichment result using cnetplot",
               plotOutput('cnetplot', height = 800)
               )
   ),
   ######################################################################################################################################
   tabItem(tabName = "emap", 
           box(width = 12, status = "primary",solidHeader = TRUE,title = "Visualize enrichment result using cnetplot",
               plotOutput('plotemap', height = 800)
           )
   ),
   ######################################################################################################################################
   tabItem(tabName = "gsea", 
           box(width = 12, status = "primary",solidHeader = TRUE,title = "Visualize GSEA result",
               plotOutput("plotgsea", height = 800)),
           box(width = 12, status = "primary",solidHeader = TRUE,title = "GSEA Result",
               DT::dataTableOutput('gseares')),
           box(width = 12, status = "primary",solidHeader = TRUE,title = "Visualize GSEA result",
               plotOutput("plotpath", height = 800))
   ),
   ######################################################################################################################################
   
      tabItem(tabName = "gogage",
              fluidRow(
                box(width = 8, status = "primary",solidHeader = TRUE,title = "GO Heatmap",
                    fluidRow(
                      column(6,uiOutput('hmplimgo')),
                      column(width = 3, offset = 2,plotOutput('hmpscale_out3',width = 200,height = 65))
                    ),
                    fluidRow(
                      column(6,checkboxInput("hmpsamp3", label = "View Heatmap of all samples", value = TRUE)),
                      column(6,h4(""))
                    ),
                    d3heatmapOutput('goheatmap',width=550,height=900)),
                box(width = 4, status = "primary",solidHeader = TRUE,title = "Controls",
                    fluidRow(
                      column(6,radioButtons(inputId='gage', label = h5("Select ontology"),
                                            choices = c("Biological Process" = 'BP', "Cellular Component" = 'cc', "Molecular Function" = 'MF'))),
                      column(6,selectInput("go_dd", "GO Selection",c('Upregulated' = "upreg",'Downregulated' = "downreg")))),
                    fluidRow( 
                      column(6,downloadButton('downloadgo', 'Download GO Data')),
                      column(6,downloadButton('downloadgogene', 'Download GO Genelist'))),
                    br(),
                    downloadButton('downloadgoheatmap', 'Download GO Heatmap'),
                    br(),
                    selectInput("hmpcol3", "Select Heatmap Color Palette",c('YlGnBu' = "YlGnBu",'RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'PRGn'="PRGn", 'Blues' = "Blues")),
                    selectInput("clusterby3", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none")),
                    checkboxInput("checkbox3", label = "Reverse Colors", value = FALSE))),
              box(width = 12, status = "primary",solidHeader = TRUE,title = "Table",
                  DT::dataTableOutput('table4'),textOutput("godesc"),DT::dataTableOutput('x4'))
      ),
      
      ######################################################################################################################################
      
      tabItem(tabName = "help",
              h4(p(strong("1. PCA Plot"))),
              h4(p(div("The PCA plot tab displays the biplot by default. You can select the principle component to plot on the x and y axis of the plot from the drop-down menu. You can also specify the number of top genes showing maximum variance to be used as the input for the bioplot as well as the number of genes you want to view in the plot. The ",em("Display variances of PC"),"tab displays the barplot showing the proportion of variance retained by each principle component.","The ",em("3D plot tab"),"displays the 3D plot of the top 3 principle components"))),
              h4("Helpful Links:", a("Click Here for information on PCA biplot", href="http://www.nature.com/nbt/journal/v26/n3/full/nbt0308-303.html")),
              br(),
              h4(p(strong("2. Project Summary and Results"))),
              h4("Select a project and a comparison. (Comparisons are automatically populated in the drop-down menu)"),
              h4(p(div(span("The Project Summary and Results tab",style="color:blue"), "will display the limma (differential expression analysis) output for that comparison in a  table. Clicking on any row will display the dot plot for that gene. Select an Attribute from the drop-dowm menu near the dot-plot to color the plot by that feature."))),
              
              h4("Make a Gene Selection by selecting the radio button to view the list of upregulated and/or downregulated genes in the ", span("Results tab",style="color:blue")),   
              h4(p(div("Type in the Fold Change cutoff and Adjusted PValue cutoff and view the updated table in the same tab. Click on",em("Download Data"),"button to download the table as a csv file"))),
              h4(p(div(span("Note:Make sure the radio button 'None' is not selected when setting FC and P.Value cutoffs",style="color:red")))),
              h4(p(div(span("Click on",em("View volcano plot"),"for the volcano plot.You can adjust the input type and the number of genes to display on the plot")))),
              h4(p(div(span("Click on",em("View Limma results of multiple contrasts")," and select",em("Display Contrast list") ,"checkbox to view the foldchange and adjusted pvalues of multiple comparisons")))),
              br(),
              h4(p(strong("3. Raw Expression Data"))),
              h4("Click on the", em("Click to view the Raw expression data"),"button to view the expression data in the",span("Raw Data tab",style="color:blue")),
              br(),
              h4(p(strong("4. Heatmap"))),
              h4("Select Heatmap type from drop-down menu and give appropriate inputs to view the Heatmap in the",span("Generate Heatmap tab",style="color:blue")),
              h4(p(div("Select color using the",em("Heatmap color palette"),"dropdown and reverse color palette using the",em("Reverse Colors"),"checkbox. Cluster by Rows and/or columns or none by selecting the option from the",em("Cluster By"),"drop-down menu"))),
              h4(p(div("Use the",em("Slider"),"to select number of genes if Heatmap type selected is 'Top number of genes'.Default is 50, minimum is 2 and maximum is 500"))),
              h4(p(div(span("Note:If the limma table has fewer genes, the heatmap will display only those despite the slider value",style="color:red")))),
              h4(p(div("You can also choose to view the heatmap of all the samples and not the samples associated with the selected contrast by unchecking the ",em("View Heatmap of all samples"),"checkbox"))),
              br(),
              h4(p(div("Upload a text file with one gene in each row if the Heatmap type selected is 'Enter Genelist'."))),
              h4(p(div(span("Note:Make sure you select the correct identifier from the drop down menu (ENSEMBL ID, ENTREZ ID, Gene Symbol)",style="color:red")))),
              h4(p(div("Generate Camera and GO data to view the heatmap if the Heatmap type selected is 'Heatmap from Camera' or 'Heatmap from GO'."))),
              br(),
              h4(p(div("You can also view a heatmap of the most variable genes from the PCA"))),
              br(),
              h4(p(strong("5. GSEA using Camera"))),
              #h4("Click on",strong("Click to view Camera results"),"button to view Camera results in the",span("GSEA tab",style="color:blue")),
              h4(p(div("Select a gene set from the ",strong("Select a Gene Set")," dropdown"))),
              h4(p(div("The Camera function in the limma package for testing differential expression, tests whether a set of genes is highly ranked relative to other genes in terms of differential expression. It takes into account the inter-gene correlation.CAMERA, an acronym for Correlation Adjusted MEan RAnk gene set test, is based on the idea of estimating the variance inflation factor associated with inter-gene correlation, and incorporating this into parametric or rank-based test procedures. It returns the number of genes in the set, the inter-gene correlation value, the direction of change (Up or Down), the two-tailed p-value and the Benjamini & Hochberg FDR adjusted P-value"))),
              h4(p(div(span("The GSEA using Camera tab",style="color:blue"), "will display the Camera output in a table. Clicking on any row will display the gene list from the user dataset that belongs to that Gene set category in table below it and the heatmap of those genes above the camera results table."))),
              #h4(p(div("Click on",em("Download Camera Data"),"button to download the table as a csv file"))),
              h4("Helpful Links:", a("Click Here for for information on Camera", href="http://nar.oxfordjournals.org/content/early/2012/05/24/nar.gks461.full")),
              br(),
              h4(p(strong("6. Pathway Analysis using SPIA"))),
              h4("The",span("Pathway Analysis using SPIA tab",style="color:blue"),"displays the SPIA results"),
              h4("Helpful Links: Click", a("here", href="http://www.bioconductor.org/packages/release/bioc/vignettes/SPIA/inst/doc/SPIA.pdf"),"and",a("here", href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1987343/"),"for information on SPIA"),
              br(),
              h4(p(strong("7. Gene Ontology using GAGE"))),
              h4("The default Ontology is set to Biological process. click on the",strong("GAGE Results"),"button in the sidebar to view GAGE results."),
              h4(p(div(span("The GO Analysis Using GAGE tab",style="color:blue"), "will display the GO output for that ontology in a table. Clicking on any row will display the gene list from the user dataset that belong to that GO-term in table below it and the heatmap correspoding to the genelist above the GAGE results table."))),
              
              h4(p(div("Click on",strong("Download GO Data"),"button to download the table as a csv file,",strong("Download GO Genelist")," for the genes associated with each GO category and ",strong("Download GO Heatmap")," for the heatmap associated with each GO term"))),
              h4("Helpful Links:", a("Click Here for information on GAGE", href="http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-161"))
      )
      ######################################################################################################################################
      
      
    )
  )
)


