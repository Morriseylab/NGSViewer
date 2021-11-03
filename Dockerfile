# Install R version 4.0
FROM r-base:4.0.0

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y keyboard-configuration
RUN apt-get install -y libcurl4-openssl-dev libssl-dev
RUN apt-get install -y xorg libx11-dev mesa-common-dev libglu1-mesa-dev
RUN apt-get install -y libxml2-dev
RUN apt-get install -y libftgl2 freetype2-demos libfreetype6-dev
RUN apt-get install -y libhdf5-dev
RUN apt-get install -y r-cran-rcppeigen
RUN apt-get install -y libgit2-dev
RUN apt-get install -y libcairo2-dev
RUN apt-get install -y libxt-dev


RUN R -e 'install.packages(c("devtools","shiny","shinydashboard","shinyjs","shinyBS","RColorBrewer","reshape2","ggplot2","ggrepel","dplyr","tidyr","plotly","htmlwidgets","DT","shinyRGL","rgl","rglwidget","readxl","png","FactoMineR","factoextra","data.table","NMF"))'

#Install packages from bioconductor
RUN R -e 'BiocManager::install(c("biomaRt","Biobase","SPIA","AnnotationDbi","org.Mm.eg.db","gage","gageData","KEGGgraph","KEGGREST","GO.db","limma","ReactomePA"))'


RUN R -e 'devtools::install_github("rstudio/d3heatmap")'


RUN mkdir /srv/shiny-server
ADD .  /srv/shiny-server
# RUN mv /srv/shiny-server/Example\ data /srv/shiny-server/data
WORKDIR /srv/shiny-server

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp(port=3838,host='0.0.0.0')"]
