library(plyr)
library(tidyverse)
library(invctr)
library(casnet)


# dummy data --------------------------------------------------------------
df <- sample(c(1:9), 100, replace = TRUE)


df_s<-ts_symbolic(df$value)

df <- rio::import("~/Documents/Projects/RNGproject/wouteR/Data_totaal.csv")

df <- rio::import("~/Documents/Projects/RNGproject/wouteR/random sequence.tsv")
y <- df$value
y <- c(2,1)
randseqR::PhiIndex(y = y, minScale = 1, maxScale = 9, responseAlternatives = c(1,2,3,4,5,6,7,8,9))


# general -----------------------------------------------------------------
# response alternatives (a) max = 9
# RNG <- function(df, a) {}

#dff <- data.frame(df)
# if(!is.null(dim(dff))){
#   if(dim(dff)[2]!=1){stop("Must have 1 column based timeseries")}
# }

minScale <- 1
maxScale <- 9













# symHelper <- function(sym_num){
#   out <- sym_num
#   i <- 1
#   same <- 0
#   while(i<=length(sym_num)){
#     if(sym_num[i]%in%c("increase","decrease")){
#       samesame <- TRUE
#       r <- i+1
#       same <- 0
#       while(samesame){
#         if(sym_num[r]%in%"same"){
#           same <- same+1
#           r <- r+1
#         } else {
#           samesame <- FALSE
#         }
#       }
#       if(same>0){
#
#         if(all(!sym_num[i]%in%c("increase"),sym_num[i+same+1]%in%c("peak","increase"))){
#           out[i:(i+same)] <- "trough"
#         } else {
#           if(all(!sym_num[i]%in%c("decrease")&sym_num[i+same+1]%in%c("trough","decrease"))){
#             out[i:(i+same)] <- "peak"
#           }
#         }
#       }
#     }
#     if(same>0){
#       i <- (i+same)
#     } else {
#       i <- i+1
#     }
#   }
#   return(out)
# }
