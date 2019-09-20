
# HELPERS ----

#' Find turning points
#'
#' @param y Numeric vector.
#'
#' @return A list object with `start` and `end` indices of turning points.
#'
#' @export
#'
findTPIs <- function(y){

  y <- casnet::ts_symbolic(y, usePlateaus = TRUE)
  startWith <- which.min(c(which(y %in% c("peak"))[1],which(y %in% c("trough"))[1]))

  yPeak   <- y %in% c("peak")
  yTrough <- y %in% c("trough")

  sIndices <- list(which(yPeak[1:NROW(y)] - yPeak[2:(NROW(y)-1)]>0), which(yTrough[1:(NROW(y)-1)] - yTrough[2:NROW(y)]>0))[[startWith]]
  eIndices <- list(which(yPeak[1:(NROW(y)-1)] - yPeak[2:NROW(y)]>0), which(yTrough[1:(NROW(y)-1)] - yTrough[2:NROW(y)]>0))[[-startWith]]

  return(list(sIndices = sIndices,
              eIndices = eIndices))

}


#' Check vector
#'
#' @inheritParams allRNG
#' @param y A vector. If not a numeric vector and `toNumeric = TRUE`, unique elements will be labelled by an integer value.
#' @param noZero Can the series have value `0`?
#' @param toNumeric If `TRUE`, the unique elements in a non-numeric vector will be labelled by an integer value.
#'
#' @return Trasformed `y` (if necessary), or, an error message.
#'
#' @export
#'
#' @examples
#'
#' y <- letters
#' check_y(y)
#'
check_y <- function(y, minScale=NA, maxScale=NA, noZero = FALSE, toNumeric = TRUE, responseAlternatives = NA){

  transformation <- character(0)

  if(any(is.na(minScale),is.na(maxScale))&any(is.na(responseAlternatives))){
    stop("Must provide either minScale & maxScale, or, responseAlternatives (or both)")
  }

  if(!is.numeric(y)){
    if(toNumeric){
      y <- casnet::as.numeric_discrete(y, sort.unique = TRUE)
      transformation <- c(transformation,'toNumeric')
    } else {
      stop("y is a non-numeric vector! Assign numbers to unique elements of y.")
    }
  }

  if(noZero){
    y[y==0] <- y[y==0]+.Machine$double.eps
    transformation <- c(transformation,'noZero')
  }

  if(any(is.na(responseAlternatives))){
    responseAlternatives <- minScale:maxScale
  } else {
    responseAlternatives <- casnet::as.numeric_discrete(responseAlternatives)
  }

  if(any(diff(responseAlternatives)!=1)){
    if(any(diff(responseAlternatives)==0)){
      stop("Elements of responseAlternatives must be unqiue!")
    } else {
      warning("Gaps in responseAlternatives. Will create a sequence based on unique elements")
      responseAlternatives <- casnet::as.numeric_discrete(as.character(responseAlternatives), sort.unique = TRUE)
      transformation <- c(transformation,"toCharacter")
    }
  }

   if(!all(y%in%names(responseAlternatives))){
      stop("Elements in y do not correspond to elements of responseAlternatives!")
   }

  if(is.na(minScale)&is.na(maxScale)){
    minScale <- responseAlternatives[1]
    maxScale <- max(responseAlternatives, na.rm = TRUE)
  }

  if(length(minScale:maxScale)!=length(responseAlternatives)){
    stop("length of minScale:maxScale must be equal to lenght of responseAltenatives")
  }

  if(length(transformation)==0){transformation <- "none"}

    attr(y,'transformation') <- transformation
    attr(y,'minScale')       <- minScale
    attr(y,'maxScale')       <- maxScale
    attr(y,"responseAlternatives") <- responseAlternatives

  invisible(return(y))
}


#' Find first occurrence
#'
#' Find first occurrence in vector `minScale:maxScale`
#'
#' @param y Numeric vector
#' @param minScale Minimum expected value.
#' @param maxScale Maximum expected value.
#'
#' @return vector with ocurrences
#' @export
#' @keywords internal
#'
first_occurence <- function(y, minScale, maxScale) {
  return(max(plyr::laply(seq(minScale,maxScale), function(x) which.min(y == x)), na.rm = TRUE))
}


#' Observed and expected combination pairs
#'
#' @param y Numeric vector
#' @param minScale Minimum expected value.
#' @param maxScale Maximum expected value.
#' @param lag Lag the series
#'
#' @return list with observed and expected paired combinations
#' @export
#' @keywords internal
#'
observed_expected <- function(y, minScale, maxScale, lag = 1){

  combi    <- DescTools::CombSet(x=minScale:maxScale,m=2,ord=TRUE, repl=TRUE)
  combis   <- paste0(combi[,1],".",combi[,2])
  observed <- paste(interaction(y[1:(NROW(y)-lag)], y[(lag+1):NROW(y)]))

  #observed <- paste(interaction(y, dplyr::lead(y, 1, default = y[1])))

  return(list(combis   = combis,
              observed = observed))
}

#' Get all classical RNG measures
#'
#' @param y A sequence of symbols. If `y` is non-numeric, unique elements will be labelled by an integer value.
#' @param minScale Minimum expected value. If `y` is a character vector this should refer to the lowest numeric code used.
#' @param maxScale Maximum expected value. If `y` is a character vector this should refer to the highest numeric code used.
#' @param responseAlternatives An optional vector of possible response alternatives. If `NA`, `responseAlternatives` will be set to `seq(minScale,maxScale)`
#' @param R Redundancy (default = `TRUE`)
#' @param RNG RNG (default = `TRUE`)
#' @param RNG2 RNG2 (default = `TRUE`)
#' @param RF Response frequencies (default = `TRUE`)
#' @param Coupon Coupon (default = `TRUE`)
#' @param NSQ NSQ (default = `TRUE`)
#' @param Adjacency Adjacency (default = `TRUE`)
#' @param TPI TPI (default = `TRUE`)
#' @param PhaseLength PhaseLength (default = `TRUE`)
#' @param Runs Runs (default = `TRUE`)
#' @param RepDistance Repetition Distance (default = `TRUE`)
#' @param RepGap Repetition Gap (default = `TRUE`)
#' @param Phi Phi index (default = `TRUE`)
#'
#' @return A list object of classical RNG task measures.
#'
#' @export
#'
#' @examples
#'
#' y <- round(runif(100,1,9))
#'
#' # Omit the Frequecny tables
#' allRNG(y, minScale = 1, maxScale = 9, RF = FALSE)
#'
allRNG <- function(y, minScale=NA, maxScale=NA, responseAlternatives = NA,
                      R      = TRUE,
                      RNG    = TRUE,
                      RNG2   = TRUE,
                      RF     = TRUE,
                      Coupon = TRUE,
                      NSQ    = TRUE,
                      Adjacency   = TRUE,
                      TPI         = TRUE,
                      PhaseLength = TRUE,
                      Runs        = TRUE,
                      RepDistance = TRUE,
                      RepGap      = TRUE,
                      Phi         = TRUE
                      ){

#disp('not implemented yet')

  y <- check_y(y, minScale = minScale, maxScale = maxScale, responseAlternatives = responseAlternatives)
  minScale <- attr(y,"minScale")
  maxScale <- attr(y,"maxScale")
  responseAlternatives <- attr(y,"responseAlternatives")

  Nresp <- NROW(y)

  if(R){
    Redundancy <- Redundancy(y = y, maxScale = maxScale)
  }

  if(RNG){
    RNG <- RNG(y = y,minScale = minScale, maxScale = maxScale)
  }


  if(Phi){
    PhiIndex <- PhiIndex(y = y,minScale = minScale, maxScale = maxScale)
  }

}


#' Measures of the classical RNG task
#'
#' Individual funcations to get the measures of the classical RNG task. Use function `allRNG()` to get a list with all (or selected) measures.
#'
#' @param y A sequence of symbols. If `y` is non-numeric, unique elements will be labelled by an integer value.
#' @param minScale Minimum expected value. If `y` is a character vector this should refer to the lowest numeric code used.
#' @param maxScale Maximum expected value. If `y` is a character vector this should refer to the highest numeric code used.
#'
#' @name classicalRNG
#'
#' @return Output
#'
#' @seealso [allRNG()] to get *all*, or, a selected list of measures.
#'
#' @examples
#'
#' y <- round(runif(100,1,9))
#'
#' R(y, maxScale = 9)
#'
NULL



# Redundancy --------------------------------------------------------------


#' Redundancy
#' @rdname classicalRNG
Redundancy<- function(y, maxScale = NA){
  # Hsingle <- log2(n) - 1/n * sum(table(y)*log2(table(y)))
  # Hmax    <- log2(maxScale)
  # R <- 100 * (1 - (Hsingle/Hmax))

  y <- check_y(y, maxScale = maxScale, noZero = TRUE)
  maxScale <- attr(y,'maxScale')

  R <- 100 * (1 - ((log2(NROW(y)) - 1/NROW(y) * sum(table(y)*log2(table(y)))) / log2(maxScale)))

  attr(R,'Name') <- 'Redundancy'
  attr(R,'y') <- y
  attr(R,'maxScale') <- maxScale

  return(R)
}

# RNG ---------------------------------------------------------------------

#' RNG
#' @rdname classicalRNG
RNG <- function(y, minScale=NA, maxScale=NA){

  y <- check_y(y, minScale=minScale, maxScale=maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')

  o_e <- observed_expected(y, minScale = minScale, maxScale = maxScale, lag = 1)

  combis   <- o_e$combis
  observed <- o_e$observed

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
#' @rdname classicalRNG
RNG2 <- function(y, minScale=NA, maxScale=NA){

  y <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')

  o_e <- observed_expected(y, minScale = minScale, maxScale = maxScale, lag = 2)

  combis   <- o_e$combis
  observed <- o_e$observed

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
#' @rdname classicalRNG
RF <- function(y){

  y <- check_y(y)

  RF  <- dplyr::as_tibble(table(y))
  return(RF)
}

# Coupon ------------------------------------------------------------------

#' Coupon
#' @rdname classicalRNG
Coupon <- function(y, minScale = NA, maxScale = NA){

  y <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')

  cnt <- 0
  stp <- FALSE
  cpn <- dplyr::tibble(x = vector(mode = "integer"))

  while (stp == FALSE) {
    cnt <- cnt + 1
    cpn[[cnt,1]] <- first_occurence(y,minScale,maxScale)
    y <- y[(cpn[[cnt,1]]+1):length(y)]
    chk <- map_lgl(seq(minScale,maxScale), function(x) any(which(y == x)))
    stp <- any(chk == FALSE)
  }

  Coupon <- mean(cpn$x, na.rm = TRUE)
  rm(y, cnt, stp, cpn, chk)

  return(Coupon)
}


# NSQ ---------------------------------------------------------------------

#' NSQ
#' @rdname classicalRNG
NSQ <- function(y, minScale=NA, maxScale=NA){

  y <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')

  o_e <- observed_expected(y, minScale = minScale, maxScale = maxScale, lag = 1)

  combis   <- o_e$combis
  observed <- o_e$observed

  NS  <- length(combis)-length(unique(observed))
  NSQ <- 100 * (NS / (maxScale^2 - 1))

  return(NSQ)
}

# First-Order difference --------------------------------------------------

# Adjacency ---------------------------------------------------------------
# uses objects from RNG/FOD to calculate

#' Adjacency
#' @rdname classicalRNG
Adjacency <- function(y){

  y <- check_y(y)

  FOD <- as_tibble(table(diff(y),dnn=c("FOD")))
  Ad_a <- 100 * FOD$n[FOD$FOD == 1] / length(y)
  Ad_d <- 100 * FOD$n[FOD$FOD == -1] / length(y)
  Ad_c <- 100 * sum(FOD$n[FOD$FOD == 1],
                    FOD$n[FOD$FOD == -1], na.rm = TRUE) / length(y)

  return(list(Ad_a = Ad_a,
              Ad_b = Ad_b,
              Ad_c = Ad_c))
}

# TPI ---------------------------------------------------------------------
# TPI depends starting direction
# may straddle in case of response repetitions (using <= and >=)

#' TPI
#' @rdname classicalRNG
TPI <- function(y){

  y <- check_y(y)

  ind      <- findTPIs(y)
  sIndices <- ind$sIndices
  eIndices <- ind$eIndices

  TPIs <- length(sIndices) + length(eIndices)

  TPI <- 100 * (TPIs) / (2/3 * (NROW(y) - 2))
  return(TPI)
}


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

#' Phase Lengths
#' @rdname classicalRNG
PhL <- function(y){
  y        <- check_y(y)
  ind      <- findTPIs(y)
  sIndices <- ind$sIndices
  eIndices <- ind$eIndices
  PhL      <- c(eIndices-sIndices, sIndices[2:length(sIndices)]-eIndices[1:(length(eIndices)-1)])
  return(PhL)
}

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


#' Runs
#' @rdname classicalRNG
Runs <- function(y){
  phases <- PhL(y)
  return(var(phases))
}

# Runs <- if_else(std == ">", var(PhL$n[PhL$values == T]),
#                 var(PhL$n[PhL$values == F]))
#



# Repetition Distance -----------------------------------------------------

#' Repetition Distance
#' @rdname classicalRNG
repDistance <- function(y, minScale = NA, maxscale = NA){

  y   <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')


  tmp <- llply(minScale:maxScale,function(v) diff(which(y==v)))
  RD  <- table(unlist(tmp))

  return(RD)
}


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

#' Repetition Gap
#' @rdname classicalRNG
repGap <- function(y, minScale = NA, maxScale = NA){

  y   <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')

  tmp <- llply(minScale:maxScale,function(v) diff(which(y==v)))
  RD  <- table(unlist(tmp))

  RD_median <- stats::median(unlist(tmp))
  RD_mean   <- mean(unlist(tmp))
  RD_mode   <- as.numeric(names(which.max(RD)))

  rm(tmp,RD)

  return(list(RD_median = RD_median,
              RD_mean   = RD_mean,
              RD_mode   = RD_mode))
}

# getmode <- function(x) {
#  uniqx <- unique(x)
#  uniqx[which.max(tabulate(match(x, uniqx)))]
# }
#
# RDmean <- mean(RD$value, na.rm = T)
# RDmed  <- median(RD$value, na.rm = T)
# RDmode <- as.double(RDT$Var1[RDT$Var1%in%which(RDT$n == max(RDT$n))])

# Phi index ---------------------------------------------------------------

#' Phi Index
#'
#' @inheritParams allRNG
#' @param maxOrder
#'
#' @rdname classicalRNG
#' @export
#'
PhiIndex <- function(y, minScale = NA, maxScale = NA, responseAlternatives = NA, maxOrder = 7){

  y <- check_y(y, minScale = minScale, maxScale = maxScale, responseAlternatives = responseAlternatives)
  minScale <- attr(y,"minScale")
  maxScale <- attr(y,"maxScale")
  responseAlternatives <- attr(y,"responseAlternatives")

  Nresp   <- NROW(y)
  setSize <- length(responseAlternatives)

  if(maxOrder>Nresp-1){
    stop("maxOrder cannot be larger than N-1")
  }

  PhiArray <- cbind(expand.grid(Phi = 2:maxOrder, RespAlt = minIndex:maxIndex),matrix(0,ncol=length(1:maxIndex)))

  lookup <- data.frame(DescTools::CombSet(c(0,1),7,repl=TRUE, ord = TRUE),0)
  colnames(lookup) <- c(paste0("pos",1:7),"freq")

  phiArrayList <- phiObsPredList <- phiObsFreqList <- phiList <- list()

  PhiFreq <- data.frame(respIns = minIndex:maxIndex, freq = as.numeric(table(y)))
  PhiFreq$invFreq <-  Nresp-PhiFreq$freq

  phiObsPred <- matrix(0,2,2, dimnames = list(c("repeat.pred","alternate.pred"),c("repeat.obs","alternate.obs")))

  for(A in 2:maxOrder){
    phiObsFreq <- list()
    phiArray <- data.frame(DescTools::CombSet(c(0,1),A,repl=TRUE, ord = TRUE),0)
    colnames(phiArray) <- c(paste0("pos",1:A),"freq")
    sameInd <- phiArray[1]==phiArray[A]
    diffInd <- phiArray[1]!=phiArray[A]
    for(B in 1:setSize){
      lookup  <- data.frame(DescTools::CombSet(c(0,1),A,repl=TRUE, ord = TRUE),0)
      colnames(lookup) <- c(paste0("pos",1:A),"freq")
      for(C in 1:((Nresp-A)+1)){
        #O <- numeric(maxOrder)
        O <- numeric(A)
        Pos <- A
        for(D in 0:(A-1)){
          if(y[C+D]==B){O[Pos] <- 1}
          Pos <- Pos-1
        } #D
        lookup$freq[rowSums(t(t(lookup[1:A]) == O))==A] <- lookup$freq[rowSums(t(t(lookup[1:A]) == O))==A] + 1
      } # C series length -A
      phiObsFreq[[B]] <- lookup
      phiArray$freq <- (phiArray$freq + lookup$freq )

      for(r in 1:NROW(lookup)){
        BinString <- lookup[r,1:A]
        if(A==2){
          freq1   <- PhiFreq[B,ifelse(O[1]==1,2,3)]
          freq2   <- PhiFreq[B,ifelse(O[2]==1,2,3)]
          freqMid <- Nresp
        } else {
          tmp     <- phiObsFreqList[[(A-1)]][[B]]
          string1   <- as.numeric(BinString[1:(A-1)])
          string2   <- as.numeric(BinString[2:A])
          freq1 <- tmp$freq[rowSums(t(t(tmp[1:(A-1)]) == string1))==(A-1)]
          freq2 <- tmp$freq[rowSums(t(t(tmp[1:(A-1)]) == string2))==(A-1)]
          rm(tmp)
          if(A==3){
            freqMid <- PhiFreq[B,ifelse(BinString[2]==1,2,3)]
          } else {
            tmp       <- phiObsFreqList[[(A-2)]][[B]]
            stringMid <- as.numeric(BinString[2:(A-1)])
            freqMid   <- tmp$freq[rowSums(t(t(tmp[1:(A-2)]) == stringMid))==(A-2)]
            rm(tmp)
          }
        } # if A==2 else
        if(freqMid>0){
          result <- freq1*freq2/freqMid
        } else {
          result <- 0
        }
        if(BinString[1]==BinString[A]){
          phiObsPred[2,1] <- phiObsPred[2,1] + result
        } else {
          phiObsPred[2,2] <- phiObsPred[2,2] + result
        }
        # rm(BinString,string1,string2,stringMid,freqMid,freq1,freq2,result)
      } #r
      rm(lookup)
    } # B response options
    phiArrayList[[A]]   <- phiArray
    phiObsFreqList[[A]] <- phiObsFreq

    phiObsPred[1,1] <- sum(phiArray$freq[sameInd])
    phiObsPred[1,2] <- sum(phiArray$freq[diffInd])
    phiObsPredList[[A]] <- phiObsPred

    a1 <- phiObsPred[1,1]
    a3 <- phiObsPred[1,2]
    a2 <- phiObsPred[2,1]
    a4 <- phiObsPred[2,2]
    Den <- sum(phiObsPred, na.rm = TRUE)

    if(Den > 0){
      AA1 <- ((a1 + a2) * (a1 + a3)) / Den
      AA2 <- ((a2 + a4) * (a1 + a2)) / Den
      AA3 <- ((a3 + a4) * (a1 + a3)) / Den
      AA4 <- ((a2 + a4) * (a3 + a4)) / Den
    }

    B1 <- ((a1-AA1)**2 / AA1)%00%0
    B2 <- ((a2-AA2)**2 / AA2)%00%0
    B3 <- ((a3-AA3)**2 / AA3)%00%0
    B4 <- ((a4-AA4)**2 / AA4)%00%0

    s <- 1
    if(a1>a4){s <- -1}

    phiList[[A]] <- s*(sqrt((B1 + B2 + B3 + B4)/(setSize*Nresp)))

    rm(a1,a2,a3,a4,AA1,AA2,AA3,AA4,B1,B2,B3,B4,s, Den)

  } # A

  names(phiList) <- paste0("phi",2:maxOrder)

  return(phiList)
}

# p <- 6
#
# for(i in seq_len(p)) {
#   tmp <- as_tibble(y) %>% mutate(x = lead(y, n = 1)) %>%
#     mutate(y = if_else(.$value == .$x, 0, 1))
#   for(j in seq_len(maxScale)) {
#     tmp <- tmp %>% mutate(z = if_else(.$value == j, 0, 1))
#   }
# }




