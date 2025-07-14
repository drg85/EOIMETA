
setwd("z:/Dropbox/RStudioStuff/VitaminDmodel/")


devtools::load_all("z:/dropbox/RStudiostuff/eoiroar")
library(readxl)
library(scrutiny)
library(dplyr)
library(openxlsx)


pathname <- "z:/Dropbox/RStudioStuff/VitaminDmodel/"

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
rrout <- roarfunc(AV,BV,CV,DV, vis = 1)



#This is on unadjusted! Be wary



#atal vector stuff
jt <- eoioutput$atalvec
aa <- jt$EVENTS_1
ba <- jt$TOTAL_1 - aa
ca <- jt$EVENTS_2
da <- jt$TOTAL_2 - ca

v1c <- aa - av
v2c <- ba - bv
v3c <- ca - cv
v4c <- da - dv











#additional analyses
# Load both packages
library(epitools)
library(epiR)

AV <- sum(av)
BV <- sum(bv)
CV <- sum(cv)
DV <- sum(dv)




# 2x2 table
tab <- matrix(c(AV, BV, CV, DV), nrow = 2, byrow = TRUE,
              dimnames = list(Exposure = c("Placebo", "VitD"),
                              Outcome = c("Dead", "Alive")))

cat("=== 2x2 Table ===\n")
print(tab)

# Wald RR from epitools
cat("\n=== Relative Risk (Wald CI, epitools) ===\n")
print(riskratio(tab, method = "wald"))

# boot RR from epitools
cat("\n=== Relative Risk (Wald CI, epitools) ===\n")
print(riskratio(tab, method = "boot"))


# Chi2 and Fisher
cat("\n=== Chi-squared Test ===\n")
print(chisq.test(tab, correct = FALSE))
cat("\n=== Fisher's Exact Test ===\n")
print(fisher.test(tab))







