# NGSViewer
R Shiny website for viewing RNA-Seq and microarray data analysed using our [pipeline.](https://github.com/Morriseylab/scripts) 

## Introduction
NGSViewer reads in the expression data, sample data, feature annotation and GSEA results for both RNA-Seq and Microarray datasets as an RData object and enables users to view and interact with their data

## Requirements
- R (version > 3.4)
- RStudio Server
- Shiny Server (if you need to host it online)

If you need help installing the above or getting started, refer to [this](https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#install-r)

## Installation
Run the following commands in R to install all required packages
```
install.packages(c("devtools","shiny","shinydashboard","shinyjs","shinyBS","RColorBrewer","reshape2","ggplot2","ggrepel",
                   "dplyr","tidyr","plotly","htmlwidgets","DT","shinyRGL","rgl","rglwidget","readxl","png","FactoMineR","factoextra",
                    "data.table","NMF"))

## try http:// if https:// URLs are not supported
BiocManager::install(c("biomaRt","Biobase","SPIA","AnnotationDbi","org.Mm.eg.db","gage","gageData","KEGGgraph","KEGGREST",
                  "GO.db","limma","ReactomePA"))
                  
devtools::install_github("rstudio/d3heatmap")
```
An alternative is to use the Dockerfile using [Shinyproxy](https://github.com/openanalytics/shinyproxy)

## Input Data format
Analyse your data using the pipeline and save the results as an RData file. Please note that the object should always be saved as **results** and filename should match the project name specified in the param.csv file. Also, the eset created by limma has a slot for phenodata extracted using the PData function. The phenodata must absolutely have the **sample_name** and **maineffect** column. Sample names must be same and in the same order as expression data. maineffect column refers to the variable or the effect that is being tested.

### Adding your dataset

Add your data to the param.csv file and move it to the data directory. You can find an example dataset [here.](http://165.123.69.6/NGSViewer/example_data.RData) Please note that the data directory must be in the same location as your server.R, ui.R and function.R files (rename the Example data folder into data). The param.csv file should also be saved in the data directory as the RData files.

### NOTE
Please note that this script requires a username and a password. Before running it, either comment out the Authentication section in server.R or add the username and password in authentication.csv file in the data folder. The username has to be entered in the param.csv file as well so that the user can view only specific datasets.
