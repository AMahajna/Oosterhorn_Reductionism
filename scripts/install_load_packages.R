#readr::local_edition(1)
###############################################################################
# install packages

if( !require("corrplot") ){
  install.packages("corrplot")
  library("corrplot")
}

if( !require("plotly") ){
  install.packages("plotly", type = "source")
  library("plotly")
}

if( !require("ggcorrplot") ){
  install.packages("ggcorrplot")
  library("ggcorrplot")
}

if( !require("BiocManager") ){
  install.packages("BiocManager")
  library("BiocManager")
}

if( !require("Biostrings") ){
  BiocManager::install("Biostrings")
  library("Biostrings")
}

if( !require("GGally") ){
  install.packages("GGally")
  library("GGally")
}

if( !require("miaViz") ) {
  BiocManager::install("miaViz")
  library("miaViz")
}

if( !require("mia") ) {
  BiocManager::install("mia")
  library("mia")
}

if( !require("scater") ) {
  BiocManager::install("scater")
  library("scater")
}

if( !require("devtools") ) {
  install.packages("devtools")
  library("devtools")
}

if( !require("tidyverse") ) {
  install.packages("tidyverse")
  library("tidyverse")
}

if( !require("readxl") ) {
  install.packages("readxl")
  library("readxl")
}

if( !require("astsa") ) {
  install.packages("astsa")
  library("astsa")
}

if( !require("readr") ) {
  install.packages("readr")
  library("readr")
}

if( !require("dplyr") ) {
  install.packages("dplyr")
  library("dplyr")
}

if( !require("ggplot2") ) {
  install.packages("ggplot2")
  library("ggplot2")
}

if( !require("forecast") ) {
  install.packages("forecast")
  library("forecast")
}

if( !require("GGally") ) {
  install.packages("GGally")
  library("GGally")
}
