library(shiny)
library(shinyBS)
#library(pathview)
library("AnnotationDbi")
library("org.Mm.eg.db")
library(gage)
library(gageData)
library(RColorBrewer)
library(NMF)
library(KEGGgraph)
library(Biobase)
library(reshape2)
library(ggplot2)
library(biomaRt)
library(KEGGREST)
library(png)
library(GO.db)
library(d3heatmap)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(FactoMineR)
library(factoextra)
library(shinyRGL)
library(rgl)
library(rglwidget)
library(SPIA)
library(ReactomePA)
library(limma)
library(ggrepel)
library(readxl)
library(biomaRt)

#Specify user-ids and passwords
auth=read.csv("data/authentication.csv")
my_username <- auth$user
my_password <- auth$pwd

#Create a theme for all plots.
plotTheme <-theme_bw() + theme(axis.title.x = element_text(face="bold", size=12),
                               axis.text.x  = element_text(angle=35, vjust=0.5, size=12),
                               axis.title.y = element_text(face="bold", size=12),
                               axis.text.y  = element_text(angle=0, vjust=0.5, size=12))

server <- function(input, output) {

  values <- reactiveValues(authenticated = FALSE)

  # Return the UI for a modal dialog with data selection input. If 'failed'
  # is TRUE, then display a message that the previous value was invalid.
  dataModal <- function(failed = FALSE) {
    modalDialog(
      textInput("username", "Username:"),
      passwordInput("password", "Password:"),
      footer = tagList(
        # modalButton("Cancel"),
        actionButton("ok", "OK")
      )
    )
  }

  # Show modal when button is clicked.
  # This `observe` is suspended only whith right user credential

  obs1 <- observe({
    showModal(dataModal())
  })

  # When OK button is pressed, attempt to authenticate. If successful,
  # remove the modal.

  obs2 <- observe({
    req(input$ok)
    isolate({
      Username <- input$username
      Password <- input$password
    })
    Id.username <- which(my_username == Username)
    Id.password <- which(my_password == Password)
    if (length(Id.username) > 0 & length(Id.password) > 0) {
      if (Id.username == Id.password) {
        Logged <<- TRUE
        values$authenticated <- TRUE
        obs1$suspend()
        removeModal()

      } else {
        values$authenticated <- FALSE
      }
    }
  })


  # output$dataInfo <- renderPrint({
  #   if (values$authenticated) "OK!!!!!"
  #   else "You are NOT authenticated"
  # })
  ##################################################
  ###################################################
  ####### LOAD EXCEL AND POPULATE DROP DOWN #########
  ###################################################
  ###################################################
  
  #Read the parameter file
  readexcel = reactive({
    user=input$username
    file = read.csv(paste("data/allusers_param.csv",sep=""))
    if(user=="allusers"){
      file=file
    }else{
      file=file[file$user==user,]
    }
  })
  
  #Get Project list and populate drop-down
  output$projects = renderUI({
    excel=readexcel()
    prj=excel$projects
    selectInput("projects","Select a project",as.list(sort(as.character(prj))))
  })
  
  ###################################################
  ###################################################
  ####### LOAD RDATA FILE AND GET CONTRASTS##########
  ###################################################
  ###################################################
  
  #Load Rdata
  fileload <- reactive({
    inFile = paste('data/',as.character(input$projects),'.RData',sep = '')
    load(inFile)
    loaddata=results
    return(loaddata)
  })
  
  #Get contrast list and populate drop-down
  output$contrasts = renderUI({
    results=fileload()
    lim=results$limma
    contrasts=as.list(as.character(unlist(lapply((names(lim)),factor))))
    selectInput("contrast","Select a comparison",contrasts,"pick one")
  })

  #Display file in dashboard
  output$dashdata<- renderTable({
    user=input$username
    #read.csv(paste("data/",user,"_param.csv",sep=""))
    file=read.csv('data/allusers_param.csv',stringsAsFactors = F)
    if(user=="allusers"){
      file=file
      colnames(file)=c("Project Name","Project Description","Old?(Y/N)","Type(RNA/Microarray)","Username")
      file=file[order(file$`Project Name`),]
    }else{
      file=file[file$user==user,] %>% dplyr::select(-old:-user)
      colnames(file)=c("Project Name","Project Description")
      file=file[order(file$`Project Name`),]
      
    }
  }, digits = 1)
  ###################################################
  ###################################################
  ##################### PCA PLOT ####################
  ###################################################
  ###################################################
  output$pcaxoptions <- renderUI({
    selectInput("pcaxaxes","Select Principle Component to plot on the X-axis ",c(1:10))
  })
  
  output$pcipslide <- renderUI({
    textInput(inputId = 'pcipslide', label = "Enter top number of input genes that show maximum variance", value = '500')
    #sliderInput("pcipslide", label = h5("Slider: Number of input genes that show maximum variance"), min = 100,max = 1000, value = 200)
  })
  
  output$pcslide <- renderUI({
    textInput(inputId = 'pcslide', label = "Enter number of genes to view in the biplot", value = '0')
    #sliderInput("pcslide", label = h5("Slider: Number of genes to view in the biplot"), min = 2,max = 50, value = 30)
  })
  
  output$pcayoptions <- renderUI({
    selectInput("pcayaxes","Select Principle Component to plot on the Y-axis",c(1:10),selected=2)
  })
  
  output$ellipse <- renderUI({
    checkboxInput("ellipse", label = "Check to view ellipses", value = FALSE)
  })
  
  output$biplottitle <- renderText({
    pcaxaxes=input$pcaxaxes
    pcayaxes=input$pcayaxes
    text=as.character(paste("Dim",pcaxaxes," vs Dim",pcayaxes,sep=""))
    return(text)
  })
  
  plotbiplot = reactive({
    res.pca = res_pca()
    x=as.numeric(input$pcaxaxes)
    y=as.numeric(input$pcayaxes)
    results=fileload()
    v = results$eset
    pData<-phenoData(v)
    validate(
      need(input$pcslide, "Enter number of genes to view in biplot")
    )
    if(input$pcslide==0 & input$ellipse==F){
      fviz_pca_ind(res.pca, repel=T,geom='point',label='var',addEllipses=FALSE, habillage = as.factor(pData$maineffect),pointsize = 3.35,axes=c(x,y))+scale_shape_manual(values = c(rep(19,length(unique(pData$maineffect)))))+theme(axis.title.x = element_text(face="bold", size=14),
                                                                                                                                                                                                                          axis.title.y = element_text(face="bold", size=14),
                                                                                                                                                                                                                          legend.text  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                          legend.title  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                          plot.title  = element_text(angle=0, vjust=0.5, size=16))
    }
    else if(input$pcslide==0 & input$ellipse==T){
        fviz_pca_ind(res.pca, repel=T,geom='point',label='var',addEllipses=T,ellipse.type="confidence",ellipse.alpha=0.2, habillage = as.factor(pData$maineffect),pointsize = 3.35,axes=c(x,y))+scale_shape_manual(values = c(rep(19,length(unique(pData$maineffect)))))+theme(axis.title.x = element_text(face="bold", size=14),
                                                                                                                                                                                                                                       axis.title.y = element_text(face="bold", size=14),
                                                                                                                                                                                                                                       legend.text  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                       legend.title  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                       plot.title  = element_text(angle=0, vjust=0.5, size=16))
        
}
      
    #fviz_pca_ind(res.pca, geom = c("point", "text"))}
    else if(input$pcslide!=0 & input$ellipse==F){fviz_pca_biplot(res.pca,repel=T, label=c("var","ind"),habillage = as.factor(pData$maineffect),pointsize = 3.35,axes=c(x,y),select.var = list(contrib = as.numeric(input$pcslide)))+scale_shape_manual(values = c(rep(19,length(unique(pData$maineffect)))))+theme(axis.title.x = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                                           axis.title.y = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                                           legend.text  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                                           legend.title  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                                           plot.title  = element_text(angle=0, vjust=0.5, size=16))
    }
    
    else{fviz_pca_biplot(res.pca,repel=T, label=c("var","ind"),addEllipses=T,ellipse.type="confidence",ellipse.alpha=0.1,habillage = as.factor(pData$maineffect),pointsize = 3.35,axes=c(x,y),select.var = list(contrib = as.numeric(input$pcslide)))+scale_shape_manual(values = c(rep(19,length(unique(pData$maineffect)))))+theme(axis.title.x = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                                           axis.title.y = element_text(face="bold", size=14),
                                                                                                                                                                                                                                                                           legend.text  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                                           legend.title  = element_text(angle=0, vjust=0.5, size=14),
                                                                                                                                                                                                                                                                           plot.title  = element_text(angle=0, vjust=0.5, size=16))
    }
  })
  
  output$biplot = renderPlot({
    plotbiplot()
  })
  
  output$dwldbiplot = renderUI({
    downloadButton('downloadbiplot', 'Download Biplot')
  }) 
  

  output$downloadbiplot <- downloadHandler(
    filename = function() {
      paste0("biplot.pdf")
    },
    content = function(file){
      pdf(file,width=14,height = 9,useDingbats=FALSE)
      plot(plotbiplot())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ########### VARIANCES OF PCA PLOT #################
  ###################################################
  ###################################################
  
  output$pcatitle <- renderText({
    text="The proportion of variances retained by the principal components can be viewed in the scree plot. The scree plot is a graph of the eigenvalues/variances associated with components"
    return(text)
  })
  
  #get expression data and perform PCA
  res_pca = reactive({
    n=as.numeric(input$pcipslide)
    validate(
      need(as.numeric(input$pcipslide) > 199, "Minimum value of input genes that show maximum variance should at least be 200")
    )
    results=fileload()
    v = results$eset
    keepGenes <- v@featureData@data
    #keepGenes <- v@featureData@data %>% filter(!(seq_name %in% c('X','Y')) & !(is.na(SYMBOL)))
    pData<-phenoData(v)
    v.filter = v[rownames(v@assayData$exprs) %in% rownames(keepGenes),]
    Pvars <- apply(v.filter@assayData$exprs,1,var)
    select <- order(Pvars, decreasing = TRUE)[seq_len(min(n,length(Pvars)))]
    v.var <-v.filter[select,]
    m<-v.var@assayData$exprs
    rownames(m) <- v.var@featureData@data$SYMBOL
    m=as.data.frame(m)
    m=unique(m)
    res.pca = PCA(t(m), graph = FALSE)
  })
  
  #get expression data and perform PCA
  pcaplo_tab = reactive({
    res.pca =res_pca()
    eigenvalues = res.pca$eig
    return(eigenvalues)
  })
  
  output$pcaplot_tab = DT::renderDataTable({
    DT::datatable(pcaplo_tab(),
                  extensions = c('Scroller'),
                  options = list(
                    searchHighlight = TRUE,
                    scrollX = TRUE
                  ))
  })
  
  output$pcaplot_ip = renderPlot({
    res.pca = res_pca()
    fviz_screeplot(res.pca, ncp=10)
  })
  

  ###################################################
  ###################################################
  ##################3D PCA PLOT #####################
  ###################################################
  ###################################################
  output$pcaplot3d = renderRglwidget({
    graphics.off()
    pdf(NULL)
    v=datasetInput3()
    results=fileload()
    pData=pData(results$eset)
    v=t(v)
    v= v[,apply(v, 2, var, na.rm=TRUE) != 0]
    #      pca <- prcomp( v, scale= TRUE )
    #      vars <- apply(pca$x, 2, var)
    #      props <- round((vars / sum(vars))*100,1)
    #      groups=factor(gsub('-','_',pData$maineffect))
    
    pca <- res_pca()
    vars <- apply(pca$var$coord, 2, var)
    props <- round((vars / sum(vars))*100,1)
    groups=factor(gsub('-','_',pData$maineffect))
    
    
    ########
    try(rgl.close())
    open3d()
    # resize window
    par3d(windowRect = c(100, 100, 612, 612))
    palette(c('blue','red','green','orange','cyan','black','brown','pink'))
    plot3d(pca$ind$coord[,1:3], col =as.numeric(groups), type='s',alpha=1.75,axes=F,
           xlab=paste('PC1 (',props[1],'%)',sep=''),
           ylab=paste('PC2 (',props[2],'%)',sep=''),
           zlab=paste('PC3 (',props[3],'%)',sep='')
    )
    axes3d(edges=c("x--", "y--", "z"), lwd=2, expand=10, labels=FALSE,box=T)
    grid3d("x")
    grid3d("y")
    grid3d("z")
    l=length(levels(groups))
    ll=1:l
    y=1+(ll*15)
    #text3d(x=70, y=y, z=0.75,levels(groups) ,col="black")
    #points3d(x=90, y=y, z=0.75, col=as.numeric(as.factor(levels(groups))), size=6)
    legend3d("topright", legend = levels(groups), pch = 16, col=palette(),cex=1, inset=c(0.02))
        #pdf(NULL)
        #rgl.snapshot('./PCA_3d_test.png', fmt = "png", top = TRUE )
        #rgl.postscript('PCA_3D_test.pdf',fmt='pdf')
    
    #movie3d(spin3d(), duration = 5,movie='PCA_movie',dir='./' )
    rglwidget()
  })
  
  ###################################################
  ###################################################
  ######## GET PROJECT DESC AND DISPLAY   ###########
  ###################################################
  ###################################################
  #Read parameter file and get project desc for the project selected
  prjdesc = reactive({
    file = readexcel()
    prj=input$projects
    desc=file$desc[file$projects %in% prj]
    desc=as.character(desc)
  })
  
  #Display text in main panel
  output$pdesc <- renderText({
    desc=prjdesc()
  })
  
  ###################################################
  ###################################################
  ###########LOAD LIMMA FILE AND DISPLAY#############
  ###################################################
  ###################################################
  #Read limma data from eset
  datasetInput0.5 = reactive({
    contrast=input$contrast
    results=fileload()
    k=paste('results$limma$',contrast,sep='')
    limmadata=eval(parse(text = k))
  })
  
  #Update datatable in tab 1 based on gene selection (upregulated, downregulated, both or none)
  datasetInput = reactive({
    contrast=input$contrast #select contrast
    limmadata=datasetInput0.5()
    limmadata=subset(limmadata, select=-c(logFC))
    lfc=as.numeric(input$lfc) #get logFC
    apval=as.numeric(input$apval)#get adjusted P.Vals
    if(is.null(input$radio))
    {
      d = limmadata
    }
    else if(input$radio=='none')
    {
      d=limmadata
    }
    else if(input$radio=='down')
    {
      d=limmadata
      d = d[which(d$fc < (-1*(lfc)) & d$adj.P.Val < apval),]
    }
    else if(input$radio=='up')
    {
      d=limmadata
      d = d[which(d$fc>lfc & d$adj.P.Val < apval),]
    }
    else if(input$radio=='both')
    {
      d=limmadata
      d = d[which(abs(d$fc) > lfc & d$adj.P.Val < apval),]
    }
    #limma = as.data.frame(d) #load limma data
    geneid=d$SYMBOL
    url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
    if(url=="http://www.genecards.org/cgi-bin/carddisp.pl?gene="){
      d$link<-NULL
    }else{
      d$link=paste0("<a href='",url,"'target='_blank'>","Link to GeneCard","</a>")}
    d=as.data.frame(d) 
    #d=as.data.frame(d) %>% select(-t)
    return(d)
  })
  
  #print input file in tab1
  output$table = DT::renderDataTable({
    input$lfc
    input$apval
    input$project
    input$contrast
    DT::datatable(datasetInput(),
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    buttons = list()),
                  rownames=FALSE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  

  
  #download data as excel sheet
  output$dwld <- downloadHandler(
    filename = function() { paste(input$projects, '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file)
    })
  
  ###################################################
  ###################################################
  #######CONDITIONAL PANEL FOR Limma ##################
  ###################################################
  ###################################################
  
  #Create checkboxes with contrasts corresponding to the project (displayed only when multiple contrast checkbox is selected)
  output$contrastslimma <- renderUI({
    results=fileload()
    lim=results$limma
    contrasts=as.list(as.character(unlist(lapply((names(lim)),factor))))
    checkboxGroupInput("multicontrast",label="Pick Contrasts",choices=contrasts)
  })
  
  #create table with p.value and FC value for the contrasts selected
  multilimma = reactive({
    validate(
      need(input$multicontrast, "Please Select at least one comparison ")
    )
    contr=input$multicontrast
    results=fileload()
    full_limma = data.frame(id=as.character())
    for(i in 1:length(contr)){
      k=paste('results$limma$',contr[i],sep='')
      limmadata=eval(parse(text = k))
      limmadata2=data.frame(id=rownames(limmadata),logFC=limmadata$logFC,adj.P.Val=limmadata$adj.P.Val)
      colnames(limmadata2)[-1]=paste(colnames(limmadata2[,c(-1)]),contr[i], sep = "_")
      full_limma=full_join(full_limma,limmadata2,by='id')
    }
    k=data.frame(id=rownames(limmadata),SYMBOL=limmadata$SYMBOL)
    m=full_join(k,full_limma,by='id')
    return(m)
  })
  
  #Display text (contrast name) above limma table
  output$contrdesc <- renderText({
    contrastname=input$contrast
    text=paste('CONTRAST:  ',contrastname,sep="   ")
    
    return(text)
  })
  
  
  
  #update table with the dataframe
  output$table_TRUE = DT::renderDataTable({
    input$project
    input$contrast
    DT::datatable(multilimma(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = list('copy')
                  ),rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  

  output$dwldmultitab = renderUI({
    downloadButton('multidwld','Download Table')
  }) 
  
  output$multidwld <- downloadHandler(
    filename = function() { paste(input$projects, '_multiple_contrasts.csv', sep='') },
    content = function(file) {
      write.csv(multilimma(), file,row.names=FALSE)
    })
  ###################################################
  ###################################################
  ############# DISPLAY VOLCANO PLOT  ###############
  ###################################################
  ###################################################
  datasetInputvol = reactive({
    limmadata=datasetInput()
    #     lfc=as.numeric(input$lfc) #get logFC
    #     apval=as.numeric(input$apval)#get adjusted P.Vals
    #     d = limmadata
    return(limmadata)
  })
  
  output$volcslider <- renderUI({
    conditionalPanel(
      condition = "input.volcdrop == 'signi'",
      fluidRow(
        column(6,sliderInput("volcslider", label = h4("Select top number of genes"), min = 0,max = 25, value = 5))
      ))
    
  })
  
  output$volcdrop <- renderUI({
    selectInput("volcdrop", "Select input type",c('Significant genes' = "signi",'GO genes' = "go"))
    
  })
  
  volcanoplot_out = reactive({
    diff_df=datasetInput0.5()
    
    FDR=input$apval
    lfc=input$lfc
    
    if(input$volcdrop=="signi"){
      diff_df$group <- "NotSignificant"
      # change the grouping for the entries with significance but not a large enough Fold change
      diff_df[which(diff_df['adj.P.Val'] < FDR & abs(diff_df['logFC']) < lfc ),"group"] <- "Filtered by FDR"
      
      # change the grouping for the entries a large enough Fold change but not a low enough p value
      diff_df[which(diff_df['adj.P.Val'] > FDR & abs(diff_df['logFC']) > lfc ),"group"] <- "Filtered by FC"
      
      # change the grouping for the entries with both significance and large enough fold change
      diff_df[which(diff_df['adj.P.Val'] < FDR & abs(diff_df['logFC']) > lfc ),"group"] <- "Significant (Filtered by both FDR and FC)"
      
      # Find and label the top peaks..
      n=input$volcslider
      if(n>0){
      top_peaks <- diff_df[with(diff_df, order(adj.P.Val,logFC)),][1:n,]
      top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(adj.P.Val,-logFC)),][1:n,])
      
      a <- list()
      for (i in seq_len(nrow(top_peaks))) {
        m <- top_peaks[i, ]
        a[[i]] <- list(
          x = m[["logFC"]],
          y = -log10(m[["adj.P.Val"]]),
          text = m[["SYMBOL"]],
          xref = "x",
          yref = "y",
          showarrow = FALSE,
          arrowhead = 0.5,
          ax = 20,
          ay = -40
        )
      }
      p <- plot_ly(data = diff_df, x = diff_df$logFC, y = -log10(diff_df$adj.P.Val),text = diff_df$SYMBOL, mode = "markers", color = diff_df$group) %>% layout(title ="Volcano Plot",xaxis=list(title="Log Fold Change"),yaxis=list(title="-log10(FDR)")) %>%
        layout(annotations = a)
        }
      else{
        p <- plot_ly(data = diff_df, x = diff_df$logFC, y = -log10(diff_df$adj.P.Val),text = diff_df$SYMBOL, mode = "markers", color = diff_df$group) %>% layout(title ="Volcano Plot",xaxis=list(title="Log Fold Change"),yaxis=list(title="-log10(FDR)"))
          }
    }
    else if(input$volcdrop=="go"){
      top_peaks2=GOHeatup()
      diff_df$group <- "All genes"
      diff_df[which(diff_df$SYMBOL %in% top_peaks2$SYMBOL ),"group"] <- "Selected_genes"
      
      # Find and label the top peaks..
      top_peaks <- diff_df[diff_df$SYMBOL %in% top_peaks2$SYMBOL,]
      
      a <- list()
      for (i in seq_len(nrow(top_peaks))) {
        m <- top_peaks[i, ]
        a[[i]] <- list(
          x = m[["logFC"]],
          y = -log10(m[["adj.P.Val"]]),
          text = m[["SYMBOL"]],
          xref = "x",
          yref = "y",
          showarrow = FALSE,
          arrowhead = 0.5,
          ax = 20,
          ay = -40
        )
      }
      p <- plot_ly(data = diff_df, x = diff_df$logFC, y = -log10(diff_df$adj.P.Val),text = diff_df$SYMBOL, mode = "markers", color = diff_df$group) %>% layout(title ="Volcano Plot",xaxis=list(title="Log Fold Change"),yaxis=list(title="-log10(FDR)")) 
      }
    p
  })
  
  volcanoplot_dout = reactive({
    diff_df=datasetInput0.5()
    
    FDR=input$apval
    lfc=input$lfc
    
    if(input$volcdrop=="signi"){
      diff_df$group <- "NotSignificant"
      diff_df[which(diff_df['adj.P.Val'] < FDR & abs(diff_df['logFC']) < lfc ),"group"] <- "Filtered by FDR"
      diff_df[which(diff_df['adj.P.Val'] > FDR & abs(diff_df['logFC']) > lfc ),"group"] <- "Filtered by FC"
      diff_df[which(diff_df['adj.P.Val'] < FDR & abs(diff_df['logFC']) > lfc ),"group"] <- "Significant (Filtered by both FDR and FC)"
      
      # Find and label the top peaks..
      n=input$volcslider
      if(n>0){
        top_peaks <- diff_df[with(diff_df, order(adj.P.Val,logFC)),][1:n,]
        top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(adj.P.Val,-logFC)),][1:n,])
        p <- ggplot(data = diff_df, aes(x = diff_df$logFC, y = -log10(diff_df$adj.P.Val))) + geom_point_rast(aes(color=diff_df$group)) +ggtitle("Volcano Plot") + xlab("Log Fold Change") + ylab("-log10(FDR)") +labs(color="")+ geom_label_repel(data=top_peaks,aes(x = top_peaks$logFC, y = -log10(top_peaks$adj.P.Val),label=top_peaks$SYMBOL)) + theme_bw()
      }
      else{
        p <- ggplot(data = diff_df, aes(x = diff_df$logFC, y = -log10(diff_df$adj.P.Val))) + geom_point_rast(aes(color=diff_df$group)) +ggtitle("Volcano Plot") + xlab("Log Fold Change") + ylab("-log10(FDR)") +labs(color="") + theme_bw()
      }
    }
    else if(input$volcdrop=="go"){
      top_peaks2=GOHeatup()
      diff_df$group <- "All genes"
      diff_df[which(diff_df$SYMBOL %in% top_peaks2$SYMBOL ),"group"] <- "Selected_genes"
      
      # Find and label the top peaks..
      top_peaks <- diff_df[diff_df$SYMBOL %in% top_peaks2$SYMBOL,]
       p <- ggplot(data = diff_df, aes(x = diff_df$logFC, y = -log10(diff_df$adj.P.Val))) + geom_point_rast(aes(color=diff_df$group)) +ggtitle("Volcano Plot") + xlab("Log Fold Change") + ylab("-log10(FDR)") +labs(color="") + theme_bw()
    }
    p
  })
  output$table_volc = DT::renderDataTable({
    DT::datatable(datasetInputvol(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  
  output$volcanoplot = renderPlotly({
   # withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      input$radio
      input$lfc
      input$apval
      input$volcslider
      input$volcdrop
      volcanoplot_out()
   # })
  })
  
 
  output$dwldvolcanoplot <- downloadHandler(
    filename = function() {
      paste0("volcano.pdf")
    },
    content = function(file){
      pdf(file,width=14,height = 9,useDingbats=FALSE)
      plot(volcanoplot_dout())
      dev.off()
    })
  ###################################################
  ###################################################
  ############## DISPLAY VOOM DATA ##################
  ###################################################
  ###################################################
  #load voom data from eset
  datasetInput3 = reactive({
    results=fileload()
    exprsdata=results$eset@assayData$exprs
  })
  
  #annotate voom data using featuresdata 
  datasetInput33 = reactive({
    results=fileload()
    exprsdata=as.data.frame(results$eset@assayData$exprs)
    features=as.data.frame(pData(featureData(results$eset)))
    features$id=rownames(features)
    exprsdata$id=rownames(exprsdata)
    genes <- inner_join(features,exprsdata,by=c('id'='id'))
    return(genes)
  })
  
  #print voom or expression data file in tab2
  output$table3 = DT::renderDataTable({
    DT::datatable(datasetInput33(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,caption= "Voom data")
  })

  output$dwldrawtab = renderUI({
    downloadButton('rawdwld','Download Raw Data')
  }) 
  
  output$rawdwld <- downloadHandler(
    filename = function() { paste(input$projects, '_rawdata.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput33(), file,row.names=FALSE)
    })
  
  ###################################################
  ###################################################
  ############## DISPLAY PHENO DATA ##################
  ###################################################
  ###################################################
  #load pheno from eset
  phenofile = reactive({
    results=fileload()
    pd=pData(results$eset)
    if("minexpr" %in% colnames(pData)){
      pd=pd %>% dplyr::select(-minexpr)
    }
    else{pd=pd}
  })
  
  #print pheno data file 
  output$phenofile = DT::renderDataTable({
    DT::datatable(phenofile(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,caption= "Sample data")
  })
  

  ###################################################
  ###################################################
  ##################### DOT PLOT ####################
  ###################################################
  ###################################################
  #create user interface for dot-plot drop down menu
  output$boxplotcol = renderUI({ 
    results=fileload()
    eset=results$eset
    pData=pData(eset) #get pheno-data
    kc=pData[ , grepl( "var_" , colnames(pData) ) ] #get columns from phenodata that start with "var"
    kt=as.data.frame(t(na.omit(t(kc)))) #omit columns that have only NA's
    kc=data.frame(maineffect=pData$maineffect,sample_name=pData$sample_name,kt) #create a dataframe with maineffect and sample_name and non-null columns starting with var
    bpcols=as.list(as.character(unlist(lapply((colnames(kc)),factor))))
    selectInput("color","Select an Attribute for the X-axis",bpcols) #populate drop down menu with the phenodata columns
  })
  
  output$boxplotcol2 = renderUI({ 
    results=fileload()
    eset=results$eset
    pData=pData(eset) #get pheno-data
    kc=pData[ , grepl( "var_" , colnames(pData) ) ] #get columns from phenodata that start with "var"
    kt=as.data.frame(t(na.omit(t(kc)))) #omit columns that have only NA's
    kc=data.frame(maineffect=pData$maineffect,sample_name=pData$sample_name,kt) #create a dataframe with maineffect and sample_name and non-null columns starting with var
    bpcols=as.list(as.character(unlist(lapply((colnames(kc)),factor))))
    selectInput("color2","Color By",bpcols) #populate drop down menu with the phenodata columns
  })
 
  output$minexprline = renderUI({ 
    tagList(
  checkboxInput("minexprline", label = "Show expression threshold line", value = FALSE),
    bsTooltip("minexprline","Please note that not all projects have this option currently", placement = "bottom", trigger = "hover",options = NULL)
    )
  })
  
  dotplot_out = reactive({
    s = input$table_rows_selected #select rows from table
    dt = datasetInput() #load limma data
    dt$id=rownames(dt)
    dt=data.frame(dt$id,dt[,-ncol(dt)])
    validate(
      need((is.data.frame(dt) && nrow(dt))!=0, "No data in table")
    )
    dt1 = dt[s, , drop=FALSE]#get limma data corresponding to selected row in table
    id = as.character(dt[s,1]) 
    results=fileload()
    eset <- results$eset
    pData=pData(eset) #get pheno-data
    minexpr=pData$minexpr[1]
    e <-data.frame(eset@phenoData@data,signal=eset@assayData$exprs[id,])
    if(is.na(dt1$SYMBOL)) #if gene symbol does not exist,use ENSEMBL id
    {genesymbol=dt1$ENSEMBL}
    else{
      genesymbol=dt1$SYMBOL} #get the gene symbol of the row selected
    if(input$minexprline==T){
    gg=ggplot(e,aes_string(x=input$color,y="signal",col=input$color2))+plotTheme+guides(color=guide_legend(title=as.character(input$color)))+
      labs(title=genesymbol, x="Condition", y="Expression Value") + geom_point(size=5,position=position_jitter(w = 0.1))+ geom_smooth(method=lm,se=FALSE) +
      stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size= 0.3, geom = "crossbar",width=.2) + geom_hline(yintercept=minexpr, linetype="dashed", color = "red")}
    else{
      gg=ggplot(e,aes_string(x=input$color,y="signal",col=input$color2))+plotTheme+guides(color=guide_legend(title=as.character(input$color)))+
        labs(title=genesymbol, x="Condition", y="Expression Value") + geom_point(size=5,position=position_jitter(w = 0.1))+ geom_smooth(method=lm,se=FALSE) +
        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size= 0.3, geom = "crossbar",width=.2)
    }
    gg
  })
  
  # output dotplot
  output$dotplot = renderPlot({
    dotplot_out()
  })
  
  output$downloaddotplot <- downloadHandler(
    filename = function() {
      paste0(input$projects, '_dotplot.jpg', sep='') 
      #paste0("dotplot.jpg")
    },
    content = function(file){
      jpeg(file, quality = 100, width = 800, height = 800)
      plot(dotplot_out())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ############# CAMERA OUTPUT DISPLAY ##############
  ###################################################
  ###################################################
  
  #populate camera dropdown menu with the genesets based on the project RData
  output$cameradd = renderUI({
    results=fileload()
    contrast=input$contrast
    cam=paste("results$camera$",contrast,sep="")
    cam=eval(parse(text=cam))
    cameradd=as.list(as.character(unlist(lapply((names(cam)),factor))))
    selectInput("cameradd","Select a Gene Set",cameradd)
  })
  
  #Get camera data from Rdata file
  geneid = reactive({
#     validate(
#       need(input$camera, "Please Click on the button to view Camera results ")
#     )
    results=fileload()
    cameradd=input$cameradd
    contrast=input$contrast #get user input for contrast/comparison
    c=paste('results$camera$',contrast,'$',cameradd,'$camera_result',sep='') #get camera data corresponding to the contrast chosen
    cam=eval(parse(text = c)) #convert string to variable
    cam=data.frame(name=rownames(cam),cam)
    name=cam$name
    if (cameradd == "GO")
    {
      url= paste("http://amigo.geneontology.org/amigo/term/",name,sep = "") #create link to Gene Ontology Consortium
      cam$link=paste0("<a href='",url,"'target='_blank'>","Link to Gene Ontology Consortium","</a>")
      cam=as.data.frame(cam)
    }else{
      url= paste("http://software.broadinstitute.org/gsea/msigdb/cards/",name,".html",sep = "")
      cam$link=paste0("<a href='",url,"'target='_blank'>","Link to Molecular Dignature Database","</a>")
      cam=as.data.frame(cam)}
    return(cam) # return datatable with camera results
  })
  
  # print out camera data in tab gsea
  output$tablecam = DT::renderDataTable({
    input$camera
    input$cameradd
    input$contrast
    isolate({
      DT::datatable(geneid(),
                    extensions = c('Buttons','Scroller'),
                    options = list(dom = 'Bfrtip',
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                   scrollX = TRUE,
                                   buttons = c('copy', 'print')
                    ),rownames= FALSE,selection = list(mode = 'single', selected =1),escape=FALSE,caption = "Camera Results")
    })
  })
  

  #get the gene-list for every row in camera results table
  campick2 = reactive({
    results=fileload()
    cameradd=input$cameradd
    contrast=input$contrast #get user input for contrast/comparison
    c=paste('results$camera$',contrast,'$',cameradd,'$indices',sep='') #get camera indices corresponding to the contrast chosen
    cameraind=eval(parse(text = c))
    cam=geneid() #get datatable with camera data from reactive
    s=input$tablecam_rows_selected # get  index of selected row from table
    cam=cam[s, ,drop=FALSE]
    res=datasetInput0.5()
    res2=datasetInput33()
    if("ENTREZID" %in% colnames(res2)){
      res2=res2
    }
    else{res2=res}
    
    #get gene list from indices
    if (cameradd == "GO")
    {
      #res2=datasetInput0.5()
      k=paste('res2$ENTREZID[cameraind$`',cam$name,'`]',sep='')}
    else{
      k=paste('res2$ENTREZID[cameraind$',cam$name,']',sep='')
    }
    genes=eval(parse(text = k)) #get entrez id's corresponding to indices
    #genesid=res[match(res$ENTREZID,genes),]
    genesid=res[res$ENTREZID %in% genes,] #get limma data corresponding to entrez id's
    #rownames(genesid)=genesid$id
    return(data.frame(genesid)) #return the genelist
  })
  
  #print data table with gene list corresponding to each row in camera datatable
  output$campick3 = DT::renderDataTable({
    input$cameradd
    input$contrast
    input$projects
    DT::datatable(campick2(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,escape=FALSE,caption="GENE LIST")
  })
  
  
  #download camera datatable
  output$downloadcam <- downloadHandler(
    filename = function() { paste('Camera_',input$projects,'_',input$contrast,'.csv', sep='') },
    content = function(file) {
      write.csv(geneid(), file)
    })
  
  #Generate text title for the gene list table
  output$camdesc <- renderText({
    s = input$tablecam_rows_selected
    dt = geneid() 
    dt = as.character(dt[s, , drop=FALSE]) 
    camname=dt[1]
    text=paste('Gene list for Camera term :',camname,sep="")
    
    return(text)
  })
  
  
  ###################################################
  ###################################################
  ############### ENRICHMENT PLOT ###################
  ###################################################
  ###################################################
  
  output$eplottab = DT::renderDataTable({
    input$eplot
    input$camera
    input$cameradd
    input$contrast
    DT::datatable(geneid(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,escape=FALSE,selection = list(mode = 'single', selected =1))
  })
  
  output$eplotdesc <- renderText({
    res=geneid()
    s=input$eplottab_rows_selected 
    row=res[s, ,drop=FALSE]
    id=row$name
    text=paste('Enrichment plot for :',id,sep="")
    
    return(text)
  })
  
  
  output$en_plot <- renderPlot({
    res=geneid()
    cameradd=input$cameradd
    contrast=input$contrast
    lim=datasetInput0.5()
    stat=data.frame(t=lim$t)
    stat$id=rownames(lim)
    stat=stat[order(stat$id),]
    stat2=stat$t
    names(stat2)=stat$id
    s=input$eplottab_rows_selected 
    row=res[s, ,drop=FALSE]
    id=row$name
    results=fileload()
    contrasts=input$contrast
    y=paste('results$camera$',contrast,'$',cameradd,'$indices$`',id,"`",sep='')
    ind=eval(parse(text =y))
    barcodeplot(stat2,index=ind)
  })
  

  ###################################################
  ###################################################
  ######### CREATE HEATMAP FROM CAMERA ##############
  ###################################################
  ###################################################
  #extract voom expression data of all genes corresponding to selected row in camera datatable
  heatmapcam <- reactive({
    genesid=campick2()  #gene list from camera
    voom=as.data.frame(datasetInput3())#voom data
    genes_cam<-voom[rownames(voom) %in% rownames(genesid),]
    
  })
  
  #create heatmap function
  camheatmap <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    results=fileload()
    pd=pData(results$eset)
    expr <- heatmapcam() #voom expression data of all genes corresponding to selected row in camera datatable
    pval=campick2() #gene list from camera
    file = readexcel()
    prj=input$projects
    old=file$old[file$projects %in% prj]
    old=as.character(old)
    seq=file$seq[file$projects %in% prj]
    seq=as.character(seq)
    if(seq=="R"){
      expr$ENSEMBL=rownames(expr)
      expr=inner_join(expr,pval,by=c('ENSEMBL'='ENSEMBL'))
      
      expr=expr[order(expr$adj.P.Val),]
      rownames(expr)=make.names(expr$SYMBOL,unique=T)
      if(old=="N"){
        expr=expr %>% dplyr::select(-ENSEMBL:-t)}
      else if(old=="Y"){
        expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    else if(seq=="M"){
      expr$id=rownames(expr)
      pval$id=rownames(pval)
      expr=inner_join(expr,pval,by=c('id'='id'))
      expr=expr[order(expr$adj.P.Val),]
      rownames(expr)=make.names(expr$SYMBOL,unique=T)
      if(old=="N"){
        expr=expr %>% dplyr::select(-ENSEMBL:-t)
      }
      else if(old=="Y"){
        expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    hmplim=input$hmplimcam
    top_expr=as.data.frame(expr)
    top_expr=top_expr[1:hmplim,]
    validate(
      need(nrow(top_expr)>1, "No results")
    )
    samples=as.character(rownames(pd))
    top_expr=top_expr %>% dplyr::select(samples)
    if(input$hmpsamp2==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
#       results=fileload()
#       pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    validate(
      need(nrow(top_expr)>1, "No results")
    )
    
    
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkbox2==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby2,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcol2))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby2,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol2)))(30),labRow = sym)}
  }
  ####################################################################################################################################################
  camheatmapalt <- function(){
    results=fileload()
    pd=pData(results$eset)
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    expr <- heatmapcam() #voom expression data of all genes corresponding to selected row in camera datatable
    pval=campick2() #gene list from camera
    #           expr$ENSEMBL=rownames(expr)
    #           #           pval<-expr[rownames(voom) %in% rownames(genesid),]
    #           expr=inner_join(expr,pval,by=c('ENSEMBL'='ENSEMBL'))
    #           #rownames(expr)=expr$SYMBOL
    #           rownames(expr)=make.names(expr$SYMBOL,unique=T)
    #           expr=expr %>% select(-ENSEMBL:-t)
    
    file = readexcel()
    prj=input$projects
    old=file$old[file$projects %in% prj]
    old=as.character(old)
    seq=file$seq[file$projects %in% prj]
    seq=as.character(seq)
    if(seq=="R"){
      expr$ENSEMBL=rownames(expr)
      expr=inner_join(expr,pval,by=c('ENSEMBL'='ENSEMBL'))
      expr=expr[order(expr$adj.P.Val),]
      rownames(expr)=make.names(expr$SYMBOL,unique=T)
      if(old=="N"){
        expr=expr %>% dplyr::select(-ENSEMBL:-t)}
      else if(old=="Y"){
        expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    else if(seq=="M"){
      expr$id=rownames(expr)
      pval$id=rownames(pval)
      expr=inner_join(expr,pval,by=c('id'='id'))
      expr=expr[order(expr$adj.P.Val),]
      rownames(expr)=make.names(expr$SYMBOL,unique=T)
      if(old=="N"){
        expr=expr %>% dplyr::select(-ENSEMBL:-t)
      }
      else if(old=="Y"){
        expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    
    hmplim=input$hmplimcam
    #top_expr=data.frame(expr[,-1])
    #           if(hmplim==0)
    #           {
    #             top_expr=data.frame(expr)}
    #           else{
    top_expr=data.frame(expr)
    top_expr=top_expr[1:hmplim,]
    validate(
      need(nrow(top_expr)>1, "No results")
    )
    # }
    samples=as.character(pd$sample_name)
    top_expr=top_expr[,eval(samples)]
        if(input$hmpsamp2==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
#       results=fileload()
#       pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    validate(
      need(nrow(top_expr)>1, "No results")
    )
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkbox2==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, input$hmpcol2))(30),labRow = sym)}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol2)))(30),labRow = sym)}
  }
  
  ######################################################################################################################################################       

  ##################################################
  ###################################################
  ################ SPIA PATHWAY ANALYSIS#############
  ###################################################
  ###################################################
  
  spia_op <- reactive({
#     validate(
#       need(input$runspia, "Please Click on the button to run SPIA ")
#     )
    results=fileload()
    contrast=input$contrast #get user input for contrast/comparison
    c=paste('results$spia$',contrast,sep='') #get SPIA data corresponding to the contrast chosen
    sp=eval(parse(text = c)) #convert string to variable
    #spia_result=data.frame(name=rownames(sp),sp)
    spia_result=data.frame(sp)
    validate(
      need(nrow(spia_result) > 1, "No Results")
    )
    spia_result$KEGGLINK <- paste0("<a href='",spia_result$KEGGLINK,"' target='_blank'>","Link to KEGG","</a>")
    return(spia_result) 
  })
  
  
  output$spiaop <- DT::renderDataTable({
    input$runspia
    input$contrast
    input$projects
   # withProgress(session = session, message = 'Calculating...',detail = 'This may take a while...',{
      isolate({
        DT::datatable(spia_op(),escape = FALSE,selection = list(mode = 'single', selected =1),
                      extensions = c('Buttons','Scroller'),
                      options = list(
                        dom = 'RMDCT<"clear">lfrtip',
                        searchHighlight = TRUE,
                        pageLength = 10,
                        lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                        scrollX = TRUE,
                        buttons = c('copy', 'print')
                      ),rownames=FALSE)
    #  })
    })
  })
  
  spiagenes = reactive({
    spiaid=spia_op() 
    final_res=datasetInput()
    s=input$spiaop_rows_selected 
    row=spiaid[s, ,drop=FALSE]
    id=paste("mmu",row$ID,sep="")
    #keggid = substr(keggid, start=1, stop=8)
    allgenelist=keggLink("mmu",id) #for each kegg id, get gene list
    p=strsplit(allgenelist,":")
    genes_entrez=sapply(p,"[",2)
    genelist=final_res[final_res$ENTREZID %in% genes_entrez,]
    return(genelist) #return the genelist
  })
  
  output$spiadesc <- renderText({
    s = input$spiaop_rows_selected
    dt = spia_op() 
    dt = dt[s, , drop=FALSE]
    camname=dt$Name
    text=paste('Gene list for SPIA term :',camname,'-',dt[2],sep="")
    
    return(text)
  })
  
  output$spiagenes = DT::renderDataTable({
    DT::datatable(spiagenes(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,escape=FALSE,selection = list(mode = 'single', selected =1,caption="Genelist"))
  })
  
  

  
  output$dwldspia <- downloadHandler(
    filename = function() { paste(input$projects,'_',input$contrast, '_spia.csv', sep='') },
    content = function(file) {
      write.csv(spia_op(), file)
    })
  
  ###################################################
  ###################################################
  ######### CREATE HEATMAP FROM SPIA ##############
  ###################################################
  ###################################################
  #extract voom expression data of all genes corresponding to selected row in camera datatable
  heatmapspia <- reactive({
    genesid=spiagenes()  #gene list from camera
    voom=as.data.frame(datasetInput3())#voom data
    genes_spia<-voom[rownames(voom) %in% rownames(genesid),]
    
  })
  
  #create heatmap function
  spiaheatmap <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    results=fileload()
    pd=pData(results$eset)
    expr <- heatmapspia() #voom expression data of all genes corresponding to selected row in camera datatable
    pval=spiagenes() #gene list from camera
    file = readexcel()
    prj=input$projects
    old=file$old[file$projects %in% prj]
    old=as.character(old)
    seq=file$seq[file$projects %in% prj]
    seq=as.character(seq)
    if(seq=="R"){
      expr$ENSEMBL=rownames(expr)
      expr=inner_join(expr,pval,by=c('ENSEMBL'='ENSEMBL'))
      
      expr=expr[order(expr$adj.P.Val),]
      rownames(expr)=make.names(expr$SYMBOL,unique=T)
      if(old=="N"){
        expr=expr %>% dplyr::select(-ENSEMBL:-t)}
      else if(old=="Y"){
        expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    else if(seq=="M"){
      expr$id=rownames(expr)
      pval$id=rownames(pval)
      expr=inner_join(expr,pval,by=c('id'='id'))
      expr=expr[order(expr$adj.P.Val),]
      rownames(expr)=make.names(expr$SYMBOL,unique=T)
      if(old=="N"){
        expr=expr %>% dplyr::select(-ENSEMBL:-t)
      }
      else if(old=="Y"){
        expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    hmplim=input$hmplimspia
    top_expr=as.data.frame(expr)
    top_expr=top_expr[1:hmplim,]
    validate(
      need(nrow(top_expr)>1, "No results")
    )
    samples=as.character(rownames(pd))
    top_expr=top_expr %>% dplyr::select(samples)
    if(input$hmpsamp2spia==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
      #       results=fileload()
      #       pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    validate(
      need(nrow(top_expr)>1, "No results")
    )
    
    
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkboxspia==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterbyspia,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcolspia))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterbyspia,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcolspia)))(30),labRow = sym)}
  }
  ####################################################################################################################################################
  spiaheatmapalt <- function(){
    results=fileload()
    pd=pData(results$eset)
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    expr <- heatmapspia() #voom expression data of all genes corresponding to selected row in camera datatable
    pval=spiagenes() #gene list from camera
    file = readexcel()
    prj=input$projects
    old=file$old[file$projects %in% prj]
    old=as.character(old)
    seq=file$seq[file$projects %in% prj]
    seq=as.character(seq)
    if(seq=="R"){
      expr$ENSEMBL=rownames(expr)
      expr=inner_join(expr,pval,by=c('ENSEMBL'='ENSEMBL'))
      expr=expr[order(expr$adj.P.Val),]
      rownames(expr)=make.names(expr$SYMBOL,unique=T)
      if(old=="N"){
        expr=expr %>% dplyr::select(-ENSEMBL:-t)}
      else if(old=="Y"){
        expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    else if(seq=="M"){
      expr$id=rownames(expr)
      pval$id=rownames(pval)
      expr=inner_join(expr,pval,by=c('id'='id'))
      expr=expr[order(expr$adj.P.Val),]
      rownames(expr)=make.names(expr$SYMBOL,unique=T)
      if(old=="N"){
        expr=expr %>% dplyr::select(-ENSEMBL:-t)
      }
      else if(old=="Y"){
        expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    
    hmplim=input$hmplimspia
    top_expr=data.frame(expr)
    top_expr=top_expr[1:hmplim,]
    validate(
      need(nrow(top_expr)>1, "No results")
    )
    samples=as.character(pd$sample_name)
    top_expr=top_expr[,eval(samples)]
    if(input$hmpsamp2spia==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
      #       results=fileload()
      #       pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    validate(
      need(nrow(top_expr)>1, "No results")
    )
    sym=rownames(top_expr)
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkboxspia==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, input$hmpcolspia))(30),labRow = sym)}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcolspia)))(30),labRow = sym)}
  }
  
  ##################################################
  ##################################################
  ############## GAGE GENE ONTOLOGY ################
  ##################################################
  ##################################################
  
  datasetInput7 = reactive({
#     validate(
#       need(input$gage, "Please Select Ontology")
#     )
    final_res=datasetInput0.5() #get limma data
    logfc=final_res$fc #get FC values from limma data
    names(logfc)=final_res$ENTREZID # get entrez ids for each row
    results=fileload()
    pd=pData(results$eset)
    organism=pd$organism
    prjs=c("DS_FalcorFoxA2","YT_mir302","RJ_ESC_Laminin","RJ_CardiacHdac7_updated","DS_FalcorKO")
    prj2=c("DK_IPSC_lungepi","ZA_Boa_PKM2")
    if(!input$projects %in% prjs){
      if(!input$projects %in% prj2){
        validate(
          need(length(unique(organism))==1,"Please check pData file for errors in organism column. Does it have more than one organism or is it empty?")
        )
        organism=unique(pd$organism)[1]
      }}
    if(input$projects %in% prjs){
      organism="mouse"
    }
    else if(input$projects %in% prj2){
      organism="human"
    }
    
    #organism=pd$organism[1]
    
    if(organism=="human")
    {
      data(go.sets.hs) #load GO data from gage
      data(go.subs.hs)
      
      if(input$gage=='BP')
      {
        gobpsets = go.sets.hs[go.subs.hs$BP]
        go_res = gage(logfc, gsets=gobpsets)
      }
      else if(input$gage=='cc')
      {
        goccsets = go.sets.hs[go.subs.hs$CC]
        go_res = gage(logfc, gsets=goccsets, same.dir=TRUE)
      }
      else if(input$gage=='MF')
      {
        gomfsets = go.sets.hs[go.subs.hs$MF]
        go_res = gage(logfc, gsets=gomfsets, same.dir=TRUE)
      }}
    else if(organism=="Rat")
    {
      data(go.sets.rn) #load GO data from gage
      data(go.subs.rn)
      
      if(input$gage=='BP')
      {
        gobpsets = go.sets.rn[go.subs.rn$BP]
        go_res = gage(logfc, gsets=gobpsets)
      }
      else if(input$gage=='cc')
      {
        goccsets = go.sets.rn[go.subs.rn$CC]
        go_res = gage(logfc, gsets=goccsets, same.dir=TRUE)
      }
      else if(input$gage=='MF')
      {
        gomfsets = go.sets.rn[go.subs.rn$MF]
        go_res = gage(logfc, gsets=gomfsets, same.dir=TRUE)
      }
    }
    else 
    {
      data(go.sets.mm) #load GO data from gage
      data(go.subs.mm)
      
      if(input$gage=='BP')
      {
        gobpsets = go.sets.mm[go.subs.mm$BP]
        go_res = gage(logfc, gsets=gobpsets)
      }
      else if(input$gage=='cc')
      {
        goccsets = go.sets.mm[go.subs.mm$CC]
        go_res = gage(logfc, gsets=goccsets, same.dir=TRUE)
      }
      else if(input$gage=='MF')
      {
        gomfsets = go.sets.mm[go.subs.mm$MF]
        go_res = gage(logfc, gsets=gomfsets, same.dir=TRUE)
      }
    }
    
    return(go_res)
  })
  
  #Get all GO terms based on user-selection (upregulated/downregulated)
  datasetInput8 = reactive({
    go_res=datasetInput7()
    go_dd=input$go_dd
    if(go_dd=="upreg"){
      res=data.frame(go_res$greater)} #load limma data
    else if(go_dd=="downreg"){
      res=data.frame(go_res$less)
    }
    res = data.frame(GOterm=rownames(res),res)
    
    #Get GO id from GO terms
    row=data.frame(lapply(res,as.character),stringsAsFactors = FALSE)
    p=strsplit(row[,1], " ")
    m=sapply(p,"[",1)
    go_up=data.frame(GO_id=m,res)
    go_term=go_up$GO_id
    url= paste("http://amigo.geneontology.org/amigo/term/",go_term,sep = "") #create link to Gene Ontology Consortium
    go_up$link=paste0("<a href='",url,"'target='_blank'>","Link to Gene Ontology Consortium","</a>")
    go_up=as.data.frame(go_up)
    return(go_up)
  })
  
  #Print GO data in datatable
  output$table4 = DT::renderDataTable({
    #input$ga
    input$go_dd
    input$gage
    input$radio
    input$project
    input$contrast
   # withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      isolate({
        DT::datatable(datasetInput8(),
                      extensions = c('Buttons','Scroller'),
                      options = list(dom = 'Bfrtip',
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                     scrollX = TRUE,
                                     buttons = c('copy','print')
                      ),rownames=FALSE,escape=FALSE,selection = list(mode = 'single', selected =1))
    #  })
    })
  })
  

  
  output$downloadgo <- downloadHandler(
    filename = function() { paste('GO_',input$projects,'_',input$contrast,'_',input$gage,'_',input$go_dd,'.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput8(), file)
    })
  
  ###################################################
  ###################################################
  ############## GET GENES FROM  GO #################
  ###################################################
  ###################################################
  # get GO associated genes
  GOHeatup = reactive({
    s = input$table4_rows_selected
    dt = datasetInput8() #load GO data
    dt = dt[s, , drop=FALSE] #get GO data corresponding to selected row in table
    results=fileload()
    pd=pData(results$eset)
    organism=pd$organism[1]
    prjs=c("DS_FalcorFoxA2","YT_mir302","RJ_ESC_Laminin","RJ_CardiacHdac7_updated","DS_FalcorKO")
    prj2=c("DK_IPSC_lungepi","ZA_Boa_PKM2")
    if(input$projects %in% prjs){
      organism="mouse"
    }
    else if(input$projects %in% prj2){
      organism="human"
    }
    goid=dt$GO_id
    if(organism=="human"){
      enterezid=paste("go.sets.hs$`",goid,"`",sep="")
    }
    else if(organism=="Rat"){
      enterezid=paste("go.sets.rn$`",goid,"`",sep="")
    }
    else{
      enterezid=paste("go.sets.mm$`",goid,"`",sep="")
    }
    entrezid=eval(parse(text=enterezid))
    limma=datasetInput0.5()
    lim_vals=limma[limma$ENTREZID %in% entrezid,]
  })
  
  #Print datatable with gene list
  output$x4 = DT::renderDataTable({
    input$gage
    input$go_dd
    #input$ga
    input$radio
    input$project
    input$contrast
    goheatup=GOHeatup()
  },caption="Gene List",escape=FALSE)
  
  #Text title for gene list table
  output$godesc <- renderText({
    s = input$table4_rows_selected
    dt = datasetInput8() #load GO data
    dt = dt[s, , drop=FALSE] #get GO data corresponding to selected row in table
    goid=dt$GO_id
    text=paste('Gene list for GO term :',goid,sep="")
    
    return(text)
  })
  
  output$downloadgogene <- downloadHandler(
    filename = function() { paste('GO_',input$projects,'_',input$contrast,'_',input$gage,'_',input$go_dd,'.csv', sep='') },
    content = function(file) {
      write.csv(GOHeatup(), file)
    })
  ###################################################
  ###################################################
  ########## MAKE HEATMAP WITH GO ###################
  ###################################################
  ###################################################
  #plot heatmap
  goheatmapup <- function(){
    results=fileload()
    pd=pData(results$eset)
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol3
    pval=GOHeatup()
    hmplim=input$hmplimgo
    top_expr=datasetInput3() 
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]#voom expression data of all genes corresponding to selected row in GO datatable
    top_expr=as.data.frame(top_expr)
    file = readexcel()
    prj=input$projects
    old=file$old[file$projects %in% prj]
    old=as.character(old)
    seq=file$seq[file$projects %in% prj]
    seq=as.character(seq)
    if(seq=="R"){
      top_expr$ENSEMBL=rownames(top_expr)
      top_expr=inner_join(top_expr,pval,by=c('ENSEMBL'='ENSEMBL'))
      top_expr=top_expr[order(top_expr$adj.P.Val),]
      #rownames(top_expr)=top_expr$SYMBOL
      rownames(top_expr)=make.names(top_expr$SYMBOL,unique=T)
      if(old=="N"){
        top_expr=top_expr %>% dplyr::select(-ENSEMBL:-t)}
      else if(old=="Y"){
        top_expr=top_expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    else if(seq=="M"){
      top_expr$id=rownames(top_expr)
      pval$id=rownames(pval)
      top_expr=inner_join(top_expr,pval,by=c('id'='id'))
      top_expr=top_expr[order(top_expr$adj.P.Val),]
      rownames(top_expr)=make.names(top_expr$SYMBOL,unique=T)
      if(old=="N"){
        top_expr=top_expr %>% dplyr::select(-ENSEMBL:-t)
      }
      else if(old=="Y"){
        top_expr=top_expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    top_expr=top_expr[1:hmplim,]
    
    validate(
      need(nrow(top_expr) >1 , "No results")
    )
    samples=as.character(rownames(pd))
    top_expr=top_expr %>% dplyr::select(samples)
    if(input$hmpsamp3==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
#       results=fileload()
#       pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
    sample=as.character(rownames(sample))
    top_expr=top_expr %>% dplyr::select(sample)}
    
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    sym=rownames(top_expr)
    if(input$checkbox3==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby3,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby3,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }
  
  goheatmapupalt <- function(){
    results=fileload()
    pd=pData(results$eset)
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol3
    pval=GOHeatup() #genelist from GO
    hmplim=input$hmplimgo
    top_expr=datasetInput3() 
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]#voom expression data of all genes corresponding to selected row in GO datatable
    top_expr=as.data.frame(top_expr)
    #sym=pval$SYMBOL
    #     hmplim=input$hmplim
    #     top_expr=datasetInput3() #voom expression data of all genes corresponding to selected row in GO datatable
    #     top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    #     top_expr=as.data.frame(top_expr)
    #     top_expr$ENSEMBL=rownames(top_expr)
    #     top_expr=inner_join(top_expr,pval,by=c('ENSEMBL'='ENSEMBL'))
    #     rownames(top_expr)=top_expr$SYMBOL
    #     rownames(top_expr)=make.names(top_expr$SYMBOL,unique=T)
    #     top_expr=top_expr %>% select(-ENSEMBL:-link)
    file = readexcel()
    prj=input$projects
    old=file$old[file$projects %in% prj]
    old=as.character(old)
    seq=file$seq[file$projects %in% prj]
    seq=as.character(seq)
    if(seq=="R"){
      top_expr$ENSEMBL=rownames(top_expr)
      top_expr=inner_join(top_expr,pval,by=c('ENSEMBL'='ENSEMBL'))
      top_expr=top_expr[order(top_expr$adj.P.Val),]
      #rownames(top_expr)=top_expr$SYMBOL
      rownames(top_expr)=make.names(top_expr$SYMBOL,unique=T)
      if(old=="N"){
        top_expr=top_expr %>% dplyr::select(-ENSEMBL:-t)}
      else if(old=="Y"){
        top_expr=top_expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    else if(seq=="M"){
      top_expr$id=rownames(top_expr)
      pval$id=rownames(pval)
      top_expr=inner_join(top_expr,pval,by=c('id'='id'))
      top_expr=top_expr[order(top_expr$adj.P.Val),]
      rownames(top_expr)=make.names(top_expr$SYMBOL,unique=T)
      if(old=="N"){
        top_expr=top_expr %>% dplyr::select(-ENSEMBL:-t)
      }
      else if(old=="Y"){
        top_expr=top_expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
      }
    }
    top_expr=top_expr[1:hmplim,]
    sym=rownames(top_expr)
    validate(
      need(nrow(top_expr) >1 , "No results")
    )
    samples=as.character(pd$sample_name)
    top_expr=top_expr[,eval(samples)]
    if(input$hmpsamp3==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
#       results=fileload()
#       pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    #Remove rows that have variance 0 (This will avoid the Na/Nan/Inf error in heatmap)
    ind = apply(top_expr, 1, var) == 0
    top_expr <- top_expr[!ind,]
    if(input$checkbox3==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv =TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv=TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }

  ###################################################
  ###################################################
  ####### CREATE HEATMAP FOR LIMMA DATA##############
  ###################################################
  ###################################################
  
  # Get gene list from user
  datasetInput41 = reactive({
    #     validate(
    #       need(input$genelist, "Please input Genelist")
    #     )
    file=input$genelistfile
    genes=read.table(file$datapath) #get complete gene list as string
    df=as.vector(genes$V1)
    df=tolower(df)
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    df=firstup(df)
    results=fileload()
    pd=pData(results$eset)
    org=unique(pd$organism)
    
    if(org=="human"){
      dataset="hsapiens_gene_ensembl"
    }
    else if(org=="Rat"){
      dataset="rnorvegicus_gene_ensembl"
    }
    else{
      dataset="mmusculus_gene_ensembl"
    }
    ensembl = useEnsembl(biomart="ensembl", dataset=dataset)
    #     p=strsplit(genes,",") #split string by ','
    #     genelist=sapply(p,"[") #get array of strings with gene id's
    
    #load limma and voom data
    limma=datasetInput()
    voom=datasetInput3()
    #get expression values of the genes in the gene list
    # user-defined identifier for the gene list
    if(input$selectidentifier=='ensembl')
    {
      #ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
      genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters ='ensembl_gene_id', values =df, mart = ensembl)
      genelist=genes$ensembl_gene_id
      
    }
    else if(input$selectidentifier=='entrez')
    {
      #ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
      genes <- getBM(attributes=c('ensembl_gene_id','entrezgene'), filters ='entrezgene', values =df, mart = ensembl)
      genelist=genes$ensembl_gene_id
      
      
      #       limma$id<-rownames(limma)
      #       ensembleid=limma$id[limma$ENTREZID %in% genelist]
      #       sym=limma[limma$ENTREZID %in% genelist,]
      #       sym=sym[,c("ENSEMBL","SYMBOL")]
      #       expr_vals=voom[rownames(voom) %in% ensembleid,]
      #       expr_vals=merge(expr_vals,sym,by="row.names")
      #       rownames(expr_vals)=expr_vals$Row.names
      #       expr_vals=data.frame(expr_vals[,-c(1,(ncol(expr_vals)-1))])
    }
    else if(input$selectidentifier=='genesym')
    {
      #ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
      genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters ='external_gene_name', values =df, mart = ensembl)
      genelist=genes$ensembl_gene_id
      
      #       limma$id<-rownames(limma)
      #       ensembleid=limma$id[limma$SYMBOL %in% genelist]
      #       sym=limma[limma$SYMBOL %in% genelist,]
      #       sym=sym[,c("ENSEMBL","SYMBOL")]
      #       expr_vals=voom[rownames(voom) %in% ensembleid,]
      #       expr_vals=merge(expr_vals,sym,by="row.names")
      #       rownames(expr_vals)=expr_vals$Row.names
      #       expr_vals=data.frame(expr_vals[,-c(1,(ncol(expr_vals)-1))])
    }
    
    sym=limma[limma$ENSEMBL %in% genelist,] %>% dplyr::select(ENSEMBL,SYMBOL)
    #sym=sym[,c("ENSEMBL","SYMBOL")]
    expr_vals=merge(voom,sym,by="row.names")
    rownames(expr_vals)=expr_vals$Row.names
    expr_vals=data.frame(expr_vals[,-c(1,(ncol(expr_vals)-1))])
    
    validate(
      need(nrow(expr_vals) > 1, "Please Check Identifier chosen or Select genelist from Raw Expression Data tab")
    )
    #     validate(
    #       need(nrow(expr_vals)== nrow(genelist), "One or more of the genes entered does not exist. Check list for typos and verify using Raw Expression Data tab")
    #     )
    return(expr_vals)
  })
  
  var.genes = reactive({
    n=as.numeric(input$vgene)
    results=fileload()
    v = results$eset
    keepGenes <- v@featureData@data
    #keepGenes <- v@featureData@data %>% filter(!(seq_name %in% c('X','Y')) & !(is.na(SYMBOL)))
    pData<-phenoData(v)
    v.filter = v[rownames(v@assayData$exprs) %in% rownames(keepGenes),]
    Pvars <- apply(v.filter@assayData$exprs,1,var)
    select <- order(Pvars, decreasing = TRUE)[seq_len(min(n,length(Pvars)))]
    v.var <-v.filter[select,]
    m<-v.var@assayData$exprs
    rownames(m) <- v.var@featureData@data$SYMBOL
    m=as.data.frame(m)
    m=unique(m)
    return(m)
  })
  ###################################################
  ###################################################
  ####### TOP VARIABLE GENES  #######################
  ###################################################
  ###################################################
  varheatmap <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol #user input-color palette
    #hmplim=input$hmplim
    top_expr <- var.genes()
    # results=fileload()
    # v = results$eset
    # keepGenes <- v@featureData@data
    # pval <- datasetInput4()
    # #get expression values of genes with highest pvals
    # top_expr=expr[match(rownames(pval),rownames(expr)),]
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$hmpsamp==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
      if(input$projects =="MoGene2_Dave"){
        ct1="Hopx"
        ct2="SPC"
      }
      results=fileload()
      pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    #top_expr=top_expr[1:hmplim,]
    #sym=pval$SYMBOL
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30))}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30))}
  }
  
  varheatmapalt <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol #user input-color palette
    #hmplim=input$hmplim
    top_expr <- var.genes()
    #pval <- datasetInput4()
    #get expression values of genes with highest pvals
    #top_expr=expr[match(rownames(pval),rownames(expr)),]
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$hmpsamp==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
      results=fileload()
      pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    #top_expr=top_expr[1:hmplim,]
    #sym=pval$SYMBOL
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$checkbox==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv = TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, hmpcol))(30))}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv = TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30))}
  }
  
  ###################################################
  ###################################################
  ####### ENTER GENELIST ############################
  ###################################################
  ###################################################
  
  #create heatmap function for gene-list given by user
  heatmap2 = function(){
    dist2 = function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol#user input-color palette
    expr = datasetInput41()
    sym=expr$SYMBOL
    expr2=data.frame(expr[,-ncol(expr)])
    validate(
      need(nrow(expr2)>1, "No results")
    )
    if(input$hmpsamp==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
      results=fileload()
      pd=pData(results$eset)
      if(input$projects =="MoGene2_Dave"){
        ct1="Hopx"
        ct2="SPC"
      }
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      expr2=expr2[,eval(sample)]}
    
    #rownames(expr2)=expr[,1]
    validate(
      need(nrow(expr2)>1, "No results")
    )
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(expr2),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{d3heatmap(as.matrix(expr2),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }
  
  heatmap2alt = function(){
    dist2 = function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol#user input-color palette
    expr = datasetInput41()
    sym=expr$SYMBOL
    expr2=data.frame(expr[,-ncol(expr)])
    validate(
      need(nrow(expr2)>1, "No results")
    )
    if(input$hmpsamp==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
      results=fileload()
      pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      expr2=expr2[,eval(sample)]}
    
    #rownames(expr2)=expr[,1]
    validate(
      need(nrow(expr2)>1, "No results")
    )
    if(input$checkbox==TRUE){
      aheatmap(as.matrix(expr2),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{aheatmap(as.matrix(expr2),distfun=dist2,scale="row",Rowv=TRUE,Colv=TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }
  
  ###################################################
  ###################################################
  #################### TOP GENES ####################
  ###################################################
  ###################################################
  output$dropdown <- renderUI({
    radio=input$radio
    if(radio=="none"){
      selectInput("sortby", "Sort By",c('FDR'="sortnone",'Absolute Fold Change' = "sortab",'Positive Fold Change' = "sortpos",'Negative Fold Change' = "sortneg"))
    }
    else if(radio=="up"){
      selectInput("sortby", "Sort By",c('FDR'="sortnone",'Fold Change' = "sortab"))
    }
    else if(radio=="down"){
      selectInput("sortby", "Sort By",c('FDR'="sortnone",'Fold Change' = "sortab"))
    }
    else if(radio=="both"){
      selectInput("sortby", "Sort By",c('FDR'="sortnone",'Absolute Fold Change' = "sortab",'Positive Fold Change' = "sortpos",'Negative Fold Change' = "sortneg"))
    }
    
  })
  
  #create heatmap function for top number of genes as chosen from the slider
  datasetInput4 <- reactive({
    validate(
      need(input$gene, "Please Enter number of genes to plot heatmap ")
    )
    #sort by pval
    n<-input$gene #number of genes selected by user (input from slider)
    d<-datasetInput()
    sortby=input$sortby
    
    if(sortby=='sortnone'){
      res<-d[order(d$adj.P.Val),]
    }
    else if(sortby=='sortab'){
      res<-d[order(-abs(d$fc)),]
    }
    else if(sortby=='sortpos'){
      res<-d[order(-d$fc),]
    }
    else if(sortby=='sortneg'){
      res<-d[order(d$fc),]
    }
    
    #res<-d[order(d$adj.P.Val),]
    if(n>nrow(d)){
      reqd_res=res[1:nrow(d),]} #get top n number of genes
    else{
      reqd_res=res[1:n,]
    }
    return(reqd_res)
  })
  
  
  #create heatmap function for top n genes
  heatmap <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol #user input-color palette
    #hmplim=input$hmplim
    expr <- datasetInput3()
    pval <- datasetInput4()
    #get expression values of genes with highest pvals
    top_expr=expr[match(rownames(pval),rownames(expr)),]
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$hmpsamp==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
      if(input$projects =="MoGene2_Dave"){
        ct1="Hopx"
        ct2="SPC"
      }
      results=fileload()
      pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    #top_expr=top_expr[1:hmplim,]
    sym=pval$SYMBOL
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }
  
  heatmapalt <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol #user input-color palette
    #hmplim=input$hmplim
    expr <- datasetInput3()
    pval <- datasetInput4()
    #get expression values of genes with highest pvals
    top_expr=expr[match(rownames(pval),rownames(expr)),]
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$hmpsamp==F){
      contrast=input$contrast
      contr=strsplit(contrast,"_vs_")
      ct1=sapply(contr,"[",1)
      ct2=sapply(contr,"[",2)
      results=fileload()
      pd=pData(results$eset)
      sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
      sample=as.character(sample)
      top_expr=top_expr[,eval(sample)]}
    #top_expr=top_expr[1:hmplim,]
    sym=pval$SYMBOL
    validate(
      need(nrow(top_expr) > 1, "No results")
    )
    if(input$checkbox==TRUE){
      aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv = TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{aheatmap(as.matrix(top_expr),distfun=dist2,scale="row",Rowv = TRUE,Colv = TRUE,fontsize = 10,color = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }
  
  #manually create scale (colorkey) for heatmap
  
  hmpscale <- reactive({
    hmpcol=input$hmpcol #user input-color palette
    voom=datasetInput3()
    min=min(voom)
    max=max(voom)
    val=sort(c(seq(min,max, length.out = 5)),decreasing=TRUE)
    df <- data.frame(x = rep(1, 5),y = val,z = factor(1:5))
    if(input$checkbox==FALSE){
      ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    }
    else{
      ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol,direction = -1)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    }
  })
  
  output$hmpscale_out = renderPlot({
    hmpscale()
  })
  
  hmpscale2 <- reactive({
    hmpcol=input$hmpcol2 #user input-color palette
    voom=datasetInput3()
    min=min(voom)
    max=max(voom)
    val=sort(c(seq(min,max, length.out = 5)),decreasing=TRUE)
    df <- data.frame(x = rep(1, 5),y = val,z = factor(1:5))
    if(input$checkbox2==FALSE){
      ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    }
    else{
      ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol,direction = -1)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    }
  })
  
  hmpscale3 <- reactive({
    hmpcol=input$hmpcol3 #user input-color palette
    voom=datasetInput3()
    min=min(voom)
    max=max(voom)
    val=sort(c(seq(min,max, length.out = 5)),decreasing=TRUE)
    df <- data.frame(x = rep(1, 5),y = val,z = factor(1:5))
    if(input$checkbox3==FALSE){
      ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    }
    else{
      ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol,direction = -1)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    }
  })
  
  hmpscalespia <- reactive({
    hmpcol=input$hmpcolspia #user input-color palette
    voom=datasetInput3()
    min=min(voom)
    max=max(voom)
    val=sort(c(seq(min,max, length.out = 5)),decreasing=TRUE)
    df <- data.frame(x = rep(1, 5),y = val,z = factor(1:5))
    if(input$checkboxspia==FALSE){
      ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    }
    else{
      ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol,direction = -1)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    }
    
    
  })
  
  output$hmpscale_out2 = renderPlot({
    hmpscale2()
  })
  
  output$hmpscale_out2spia = renderPlot({
    hmpscalespia()
  })
  
  output$hmpscale_out3 = renderPlot({
    hmpscale3()
  })
  #Set limit for number of genes that can be viewed in the heatmap
  output$hmplimcam <- renderUI({
    #expr <- heatmapcam() 
    pval=campick2() 
    #       rownames(expr)=pval$SYMBOL
    #       top_expr=as.data.frame(expr)
    top_expr=datasetInput3()
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    mx=nrow(top_expr)
    sliderInput("hmplimcam", label = h5("Select number of genes to view in the heatmap"), min = 2,max =mx, value = mx)
  })
  
  #Set limit for number of genes that can be viewed in the heatmap
  output$hmplimspia <- renderUI({
    pval=spiagenes() 
    top_expr=datasetInput3()
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    mx=nrow(top_expr)
    sliderInput("hmplimspia", label = h5("Select number of genes to view in the heatmap"), min = 2,max =mx, value = mx)
  })
  
  #Set limit for number of genes that can be viewed in the heatmap
  output$hmplimgo <- renderUI({
    pval=GOHeatup()
    top_expr=datasetInput3()
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    mx=nrow(top_expr)
    sliderInput("hmplimgo", label = h5("Select number of genes to view in the heatmap"), min = 2,max =mx, value = mx)
  })
  
  output$hmpsamp <- renderUI({
    checkboxInput("hmpsamp", label = "View Heatmap of all samples", value = TRUE)
  })
  output$hmpsamp2 <- renderUI({
    checkboxInput("hmpsamp2", label = "View Heatmap of all samples", value = TRUE)
  })
  output$hmpsamp3 <- renderUI({
    checkboxInput("hmpsamp3", label = "View Heatmap of all samples", value = TRUE)
   })
  
  output$hmpsamp2spia <- renderUI({
    checkboxInput("hmpsamp2spia", label = "View Heatmap of all samples", value = TRUE)
  })
  
  #Text title for type of heatmap being displayed in the heatmap tab
  output$htitle <- renderText({
    hmip=input$hmip
    if(input$hmip=="genenum"){text="Heatmap of Top Genes "}
    else if(input$hmip=="geneli"){text="Heatmap of Genelist "}
  })
  
  # make heatmap for genes
  output$heatmap <- renderD3heatmap({
    input$hmpcol #user input-color palette
    input$clusterby #user input-cluster by
    input$checkbox #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$hmip
    input$makeheat
    input$gage
    input$go_dd
    input$ga
    input$table4_rows_selected
    input$tablecam_rows_selected
    input$radio
    input$projects
    input$contrast
    input$cameradd
    input$hmpsamp
    input$hmplim
    input$lfc
    input$apval
    input$sortby
    input$vgene
    #if user selected enter n num of genes, call heatmap() and if user entered genelist, call heatmap2()
    isolate({
      if(input$hmip == 'genenum'){heatmap()}
      else if(input$hmip == 'geneli'){heatmap2()}
      else if(input$hmip == 'vargenes' ){varheatmap()}
      #else if(input$hmip == 'hmpgo'){goheatmapup()}
    })
  })
  
  # make heatmap for genes
  output$camheatmap <- renderD3heatmap({
    input$hmpcol #user input-color palette
    input$clusterby #user input-cluster by
    input$checkbox #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$makeheat
    input$gage
    input$go_dd
    #input$ga
    input$table4_rows_selected
    input$tablecam_rows_selected
    #input$radio
    input$projects
    input$contrast
    input$cameradd
    input$hmpsamp2
    input$hmplimcam
    #input$lfc
    #input$apval
    #input$sortby
    #if user selected enter n num of genes, call heatmap() and if user entered genelist, call heatmap2()
    #isolate({
    camheatmap()
    #})
  })
  
  # make heatmap for genes
  output$spiaheatmap <- renderD3heatmap({
    input$hmpcolspia #user input-color palette
    input$clusterbyspia #user input-cluster by
    input$checkboxspia #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$spiaop_rows_selected
    input$projects
    input$contrast
    input$hmpsamp2spia
    input$hmplimspia
    spiaheatmap()
    
  })
  
  # make heatmap for genes
  output$goheatmap <- renderD3heatmap({
    input$hmpcol #user input-color palette
    input$clusterby #user input-cluster by
    input$checkbox #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$makeheat
    input$gage
    input$go_dd
    #input$ga
    input$table4_rows_selected
    input$tablecam_rows_selected
    #input$radio
    input$projects
    input$contrast
    input$cameradd
    input$hmpsamp3
    input$hmplimgo
    #input$lfc
    #input$apval
    #input$sortby
    #if user selected enter n num of genes, call heatmap() and if user entered genelist, call heatmap2()
    #isolate({
    goheatmapup()
    #})
  })
  
  #Download heatmaps 
  output$downloadheatmap <- downloadHandler(
    filename = function(){
      paste0('heatmap','.pdf',sep='')
    },
    content = function(file){
      #png(file)
      #jpeg(file, quality = 100, width = 800, height = 1300)
      pdf(file,width=9,height =14,useDingbats=FALSE)
      if(input$hmip == 'genenum'){heatmapalt()}
      else if(input$hmip == 'geneli'){heatmap2alt()}
             else if(input$hmip == 'vargenes' ){varheatmapalt()}
      #       else if(input$hmip == 'hmpgo'){goheatmapupalt()}
      dev.off()
    })
  
  #Download heatmaps 
  output$downloadcamheatmap <- downloadHandler(
    filename = function(){
      paste0('camera_heatmap','.pdf',sep='')
    },
    content = function(file){
      # png(file)
      # jpeg(file, quality = 100, width = 800, height = 1300)
      pdf(file,width=9,height = 14,useDingbats=FALSE)
      camheatmapalt()
      dev.off()
    })
  
  #Download heatmaps 
  output$downloadgoheatmap <- downloadHandler(
    filename = function(){
      paste0('GO_heatmap','.pdf',sep='')
    },
    content = function(file){
      # png(file)
      # jpeg(file, quality = 100, width = 800, height = 1300)
      pdf(file,width=9,height = 14,useDingbats=FALSE)
      goheatmapupalt()
      dev.off()
    })

  #Download heatmaps 
  output$downloadspiaheatmap <- downloadHandler(
    filename = function(){
      paste0('SPIA_heatmap','.pdf',sep='')
    },
    content = function(file){
      # png(file)
      # jpeg(file, quality = 100, width = 800, height = 1300)
      pdf(file,width=9,height = 14,useDingbats=FALSE)
      spiaheatmapalt()
      dev.off()
    })
  
  

  
}#end of server
