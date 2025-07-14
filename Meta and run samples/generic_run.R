

devtools::load_all("c:/xxxxx/eoiroar")  #change to pathname
library(readxl)
library(scrutiny)
library(dplyr)
library(openxlsx)


pathname <- "xxxx"  #change to path name of excel files

filename <- "all 2025"
exname <- ".xlsx"
addon <-"_analysis"
filename2 <-paste0(filename,addon,exname)

pathload = paste0(pathname,filename,exname)


data <- read_excel(pathload)
Study <- data$STUDY_ID

av <- data$EVENTS_1
bv <- data$TOTAL_1 - av
cv <- data$EVENTS_2
dv <- data$TOTAL_2 - cv

eoioutput <- eoimeta(av,bv,cv,dv,atal = 1)




