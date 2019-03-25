
# HELPERS ----

#' Title
#'
#' @param sym_num
#'
#'
#' @return
#' @export
#'
#' @examples
symHelper <- function(sym_num){
  out <- sym_num
  i <- 1
  same <- 0
  while(i<=length(sym_num)){
    if(sym_num[i]%in%c("increase","decrease")){
      samesame <- TRUE
      r <- i+1
      same <- 0
      while(samesame){
        if(sym_num[r]%in%"same"){
          same <- same+1
          r <- r+1
        } else {
          samesame <- FALSE
        }
      }
      if(same>0){

        if(all(!sym_num[i]%in%c("increase"),sym_num[i+same+1]%in%c("peak","increase"))){
          out[i:(i+same)] <- "trough"
        } else {
          if(all(!sym_num[i]%in%c("decrease")&sym_num[i+same+1]%in%c("trough","decrease"))){
            out[i:(i+same)] <- "peak"
          }
        }
      }
    }
    if(same>0){
      i <- (i+same)
    } else {
      i <- i+1
    }
  }
  return(out)
}


check_y <- function(y, minScale=NULL, maxScale=NULL, noZero = TRUE){

  if(any(is.na(minScale),is.na(maxScale))){
    stop("Must provide minimum and maximum values!")
  }

  if(noZero){
    y[y==0] <- y[y==0]+.Machine$double.eps
  }

  invisible(return(y))
}


first_occurence <- function(y, minScale, maxScale) {
  return(max(plyr::laply(seq(minScale,maxScale), function(x) min(which(y == x))), na.rm = TRUE))
}


#'  get all RNG
#'
#' @param y
#' @param minScale
#' @param maxScale
#' @export
#' @examples
#'
#' RNGmeasures(y = round(runif(100,1,9)), minScale = 1, maxScale = 9)
getallRNG <- function(y, minScale=NA, maxScale=NA){



}


#' RNG measures
#'
#' @param y A number sequence
#' @param minScale Minimum value possible
#' @param maxScale Maximum value paossible
#'
#' @name RNGmeasures
#' @return Output
#'
#' @export
#'
#' @examples
#'
#'
NULL



# Redundancy --------------------------------------------------------------


#' Redundancy
#' @export
#' @rdname RNGmeasures
Redundancy<- function(y, maxScale=NA){
  # Hsingle <- log2(n) - 1/n * sum(table(y)*log2(table(y)))
  # Hmax    <- log2(maxScale)
  # R <- 100 * (1 - (Hsingle/Hmax))

  y <- check_y(y, maxScale = maxScale)

  R <- 100 * (1 - ((log2(NROW(y)) - 1/NROW(y) * sum(table(y)*log2(table(y)))) / log2(maxScale)))

  attr(R,'Name') <- 'Redundancy'
  attr(R,'y') <- y
  attr(R,'maxScale') <- maxScale

  return(R)
}

# RNG ---------------------------------------------------------------------

#' RNG
#' @export
#' @rdname RNGmeasures
RNG <- function(y,minScale=NA, maxScale=NA){

  check_y(y, minScale=NA, maxScale=NA)

  combi  <- DescTools::CombSet(x=minScale:maxScale,m=2,ord=TRUE, repl=TRUE)
  combis <- paste0(combi[,1],".",combi[,2])
  observed <- paste(interaction(y, dplyr::lead(y, 1, default = y[1])))
  Afreq <- table(y)
  freqs <- plyr::ldply(combis, function(c){data.frame(pair=c,freq=sum(observed%in%c))})
  RNG   <- sum(freqs$freq * log(freqs$freq),na.rm = TRUE) / sum(Afreq*log(Afreq), na.rm = TRUE)

  attr(RNG,'Name') <- 'RNG'
  attr(RNG,'y') <- y
  attr(RNG,'minScale') <- minScale
  attr(RNG,'maxScale') <- maxScale

  return(RNG)
}

# RNG2 --------------------------------------------------------------------

#' RNG2
#' @export
#' @rdname RNGmeasures
RNG2 <- function(y, minScale=NA, maxScale=NA){

  check_y(y, minScale=NA, maxScale=NA)

  combi  <- DescTools::CombSet(x=minScale:maxScale,m=2,ord=TRUE, repl=TRUE)
  combis <- paste0(combi[,1],".",combi[,2])
  observed <- paste(interaction(y[1:(NROW(y)-2)], y[3:NROW(y)])) # lead(y, 1,default = y[1])))
  Afreq <- table(y)
  freqs <- plyr::ldply(combis, function(c){data.frame(pair=c,freq=sum(observed%in%c))})
  RNG2 <- sum(freqs$freq * log(freqs$freq),na.rm = TRUE) / sum(Afreq*log(Afreq), na.rm = TRUE)

  attr(RNG2,'Name') <- 'RNG2'
  attr(RNG2,'y') <- y
  attr(RNG2,'minScale') <- minScale
  attr(RNG2,'maxScale') <- maxScale

  return(RNG2)

}


# Response frequencies ----------------------------------------------------


#' RF
#' @export
#' @rdname RNGmeasures
RF <- function(y){
  check_y(y)

  RF  <- dplyr::as_tibble((table(y))
  return(RF)
}



# Coupon ------------------------------------------------------------------


#' Coupon
#' @export
#' @rdname RNGmeasures
Coupon <- function(y, minScale = NA, maxScale = NA){

  check_y(y, minScale = minScale, maxScale = maxScale)

  cnt <- 0
  stp <- FALSE
  cpn <- dplyr::tibble(x = vector(mode = "integer"))

  while (stp == FALSE) {
    cnt <- cnt + 1
    cpn[[cnt,1]] <- first_occurence(y,minScale,maxScale)
    yc <- yc[(cpn[[cnt,1]]+1):length(y)]
    chk <- map_lgl(seq(minScale,maxScale), function(x) any(which(yc == x)))
    stp <- any(chk == FALSE)
  }

  Coupon <- mean(cpn$x, na.rm = TRUE)
  rm(yc, cnt, stp, cpn, chk)

}

# NSQ ---------------------------------------------------------------------
# uses objects from RNG to calculate
NS  <- length(combis)-length(unique(observed))
NSQ <- 100 * (NS / (maxScale^2 - 1))

# First-Order difference --------------------------------------------------
# Uses objects from RNG to calculate
FOD <- as_tibble(table(diff(y),dnn=c("FOD")))


# Adjacency ---------------------------------------------------------------
# uses objects from RNG/FOD to calculate
Ad_a <- 100 * FOD$n[FOD$FOD == 1] / length(y)
Ad_d <- 100 * FOD$n[FOD$FOD == -1] / length(y)
Ad_c <- 100 * sum(FOD$n[FOD$FOD == 1],
                  FOD$n[FOD$FOD == -1], na.rm = TRUE) / length(y)

# TPI ---------------------------------------------------------------------
# TPI depends starting direction
# may straddle in case of response repetitions (using <= and >=)

yc <- ts_symbolic(y,usePlateaus = TRUE)
#yc <- symHelper(as.character(ys))

#yc <- ts_symbolic(as.numeric(tmp),usePlateaus = TRUE)

startWith <- which.min(c(which(yc %in% c("peak"))[1],which(yc %in% c("trough"))[1]))

yPeak   <- yc %in% c("peak")
yTrough <- yc %in% c("trough")
sIndices <- list(which(yPeak[1:99] - yPeak[2:100]>0),which(yTrough[1:99] - yTrough[2:100]>0))[[startWith]]
eIndices <- list(which(yPeak[1:99] - yPeak[2:100]>0),which(yTrough[1:99] - yTrough[2:100]>0))[[-startWith]]

TPIs <- length(sIndices) + length(eIndices)

TPI <- 100 * (TPIs) / (2/3 * (n - 2))

phases <- c(eIndices-sIndices,sIndices[2:length(sIndices)]-eIndices[1:(length(eIndices)-1)])

#table(phases)
runs <- var(phases)

# tmp<-import("/Users/fredhasselman/Documents/Projects/RNGproject/wouteR/temp.txt")

# std <- if_else(y[2] < y[1], ">", "<")
# TPf <- function(y, std) {
#   switch(std,
#          "<" = as_tibble(y[-1L] <= y[-length(y)]),
#          ">" = as_tibble(y[-1L] >= y[-length(y)])
#          ) %>%
#   mutate(y = lead(.$value, n = 1L)) %>%
#   mutate(z = abs(.$value - .$y))
# }
#
# tpo <- TPf(y, std)
# TPI = 100 * sum(tpo$z, na.rm = T) / (2/3 * (n - 2))
#

# TPs <- as_tibble(y[-1L] < y[-length(y)]) %>%
#   mutate(y = lead(.$value, n = 1L)) %>%
#   mutate(z = abs(.$value - .$y))

#
# # Phase length ------------------------------------------------------------
# # uses objects from TPI to calculate
# # excluding first and last run
# # TP at end of response repetitions instead of halfway!
# PhL <- as_tibble(
#   table(
#     rle(
#       tpo$value[head(which(tpo$value == T), 1):tail(which(tpo$value == T), 1)])
#     )
#   )
#
# PLi <- PhL %>%
#   group_by(lengths) %>%
#   summarise(PL_n = sum(n))
#
# }
#
# # Runs --------------------------------------------------------------------
# # uses objects from Phase length to calculate
# # repetitions are treated as breaks in sequences
# Runs <- if_else(std == ">", var(PhL$n[PhL$values == T]),
#                 var(PhL$n[PhL$values == F]))
#

# Repetition Distance -----------------------------------------------------
tmp <- llply(1:maxScale,function(n) diff(which(y==n)))
RD  <- table(unlist(tmp))

#
# Ry <- function(yr, rd_1) {
#   x <- if_else(yr$values[1] == T & yr$values[length(yr$values)] == T, "both",
#        if_else(yr$values[1] == T & yr$values[length(yr$values)] == F, "first",
#        if_else(yr$values[1] == F & yr$values[length(yr$values)] == T, "last",
#          "none")))
#   switch(x,
#          "both"  = rd_1,
#          "first" = rd_1[1:(length(rd_1)-1)],
#          "last"  = rd_1[2:(length(rd_1))],
#          "none"  = rd_1[2:(length(rd_1)-1)]
#   )
# }
#
# rd_t <- list()
# cnt <- 0
#
# for(i in seq_len(maxScale)) {
#   cnt  <- cnt + 1
#   yr  <- rle(y == i)
#   rd_1 <- yr$lengths[yr$values == F] + 1
#   rd_2 <- yr$lengths[yr$values == T & yr$lengths >= 2] - 1
#   rd_t[[cnt]] <- as_tibble(c(Ry(yr, rd_1), rep.int(1, times = sum(rd_2))))
#   rm(yr, rd_1, rd_2)
# }
#
# RD  <- as_tibble(unlist(rd_t))
# RDT <- as_tibble(table(RD$value))
# rm(cnt, i, rd_t)


# Repetition Gap ----------------------------------------------------------
# uses Repetition Distance information to calculate
RD_median <- median(unlist(tmp))
RD_mean <- mean(unlist(tmp))
RD_mode <- as.numeric(names(which.max(RD)))

# getmode <- function(x) {
#  uniqx <- unique(x)
#  uniqx[which.max(tabulate(match(x, uniqx)))]
# }
#
# RDmean <- mean(RD$value, na.rm = T)
# RDmed  <- median(RD$value, na.rm = T)
# RDmode <- as.double(RDT$Var1[RDT$Var1%in%which(RDT$n == max(RDT$n))])

# Phi index ---------------------------------------------------------------
f_exp <- table(y[1:3-1])*table(y[2:3]) / table(y[2:3-1])


# p <- 6
#
# for(i in seq_len(p)) {
#   tmp <- as_tibble(y) %>% mutate(x = lead(y, n = 1)) %>%
#     mutate(y = if_else(.$value == .$x, 0, 1))
#   for(j in seq_len(maxScale)) {
#     tmp <- tmp %>% mutate(z = if_else(.$value == j, 0, 1))
#   }
# }

return(data.frame(RNG = RNG))
}


