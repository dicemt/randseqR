
# dummy data --------------------------------------------------------------
df <- sample(c(1:9), 100, replace = TRUE)

df_s<-ts_symbolic(df$value)

df <- rio::import("~/Documents/Projects/RNGproject/wouteR/random sequence.tsv")
y <- df$value
# general -----------------------------------------------------------------
# response alternatives (a) max = 9
# RNG <- function(df, a) {}

#dff <- data.frame(df)
# if(!is.null(dim(dff))){
#   if(dim(dff)[2]!=1){stop("Must have 1 column based timeseries")}
# }

minScale <- 1
maxScale <- 9