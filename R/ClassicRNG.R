
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
#  startWith <- which.min(c(which(y %in% c("peak"))[1],which(y %in% c("trough"))[1]))
  startWith <- which.min(c(which(y == 5)[1],which(y == 1)[1]))

  yPeak   <- y == 5 # %in% c("peak")
  yTrough <- y == 1 # %in% c("trough")
  ind     <- list(which(yPeak[1:(NROW(y)-1)] - yPeak[2:NROW(y)]>0), which(yTrough[1:(NROW(y)-1)] - yTrough[2:NROW(y)]>0))

  return(list(sIndices = ind[[startWith]],
              eIndices = ind[[-startWith]]))
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

  if(any(is.na(minScale),is.na(maxScale))&&any(is.na(responseAlternatives))){
    stop("Must provide either minScale & maxScale, and/or, responseAlternatives")
  }


  if(any(is.na(minScale), is.na(maxScale))){
    minScale <- responseAlternatives[1]
    maxScale <- max(responseAlternatives, na.rm = TRUE)
  }


  if(any(is.na(responseAlternatives))){
    responseAlternatives <- minScale:maxScale
  } else {
    responseAlternatives <- casnet::as.numeric_discrete(responseAlternatives)
  }


#  if(is.na(responseAlternatives)&&all(!is.na(minScale),!is.na(maxScale))){
#    responseAlternatives <- minScale:maxScale
#  } else {
#    stop("Must provide minScale and maxScale to generate responseAlternatives")
#  }


  if(!is.numeric(y)){
    if(toNumeric){
      y <- casnet::as.numeric_discrete(y, sortUnique = TRUE)
      transformation <- c(transformation,'toNumeric')
    } else {
      stop("y is a non-numeric vector! Assign numbers to unique elements of y.")
    }
  }


  if(noZero){
    y[y==0] <- y[y==0]+.Machine$double.eps
    transformation <- c(transformation,'noZero')
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


   if(!all(y%in%responseAlternatives)){
      stop("Elements in y do not correspond to elements of responseAlternatives!")
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
  return(max(purrr::map_int(seq(minScale,maxScale), function(x) min(which(y == x), na.rm = TRUE)), na.rm = TRUE))
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
#' # Omit the Frequency tables
#' allRNG(y, minScale = 1, maxScale = 9, RF = FALSE)
#'
allRNG <- function(y,
                   minScale=NA,
                   maxScale=NA,
                   responseAlternatives = NA,
                   R           = TRUE,
                   RNG         = TRUE,
                   RNG2        = TRUE,
                   RF          = TRUE,
                   Coupon      = TRUE,
                   NSQ         = TRUE,
                   Adjacency   = TRUE,
                   TPI         = TRUE,
                   PhL         = TRUE,
                   Runs        = TRUE,
                   repDistance = TRUE,
                   repGap      = TRUE,
                   Phi         = TRUE
                   ){

#disp('not implemented yet')
  y        <- check_y(y, minScale = minScale, maxScale = maxScale, responseAlternatives = NA)
  minScale <- attr(y,"minScale")
  maxScale <- attr(y,"maxScale")
  responseAlternatives <- attr(y,"responseAlternatives")

  if(R){
    Redundancy <- Redundancy(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(RNG){
    RNG <- RNG(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(RNG2){
    RNG2 <- RNG2(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(RF){
    RF <- RF(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(Coupon){
    Coupon <- Coupon(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(NSQ){
    NSQ <- NSQ(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(Adjacency){
    Adjacency <- Adjacency(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(TPI){
    TPI <- TPI(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(PhL){
    PhL <- PhL(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(Runs){
    Runs <- Runs(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(repDistance){
    repDistance <- repDistance(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(repGap){
    repGap <- repGap(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(Phi){
    phiIndex <- phiIndex(y = y, minScale = minScale, maxScale = maxScale)
  }

  out <- list("Reduncany"      = Redundancy,
              "RNG"            = RNG,
              "RNG2"           = RNG2,
              "RF"             = RF,
#              "RF_1"           = RF$`1`,
#              "RF_2"           = RF$`2`,
#              "RF_3"           = RF$`3`,
#              "RF_4"           = RF$`4`,
#              "RF_5"           = RF$`5`,
#              "RF_6"           = RF$`6`,
#              "RF_7"           = RF$`7`,
#              "RF_8"           = RF$`8`,
#              "RF_9"           = RF$`9`,
              "Coupon"         = Coupon,
              "NSQ"            = NSQ,
              "Adjacency"      = Adjacency,
#              "Adj_ascending"  = Adjacency$Adj_ascending,
#              "Adj_descending" = Adjacency$Adj_descending,
#              "Adj_combined"   = Adjacency$Adj_combined,
              "TPI"            = TPI,
              "phaselength"    = PhL,
              "Runs"           = Runs,
              "repDistance"    = repDistance,
              "RD_median"      = repGap$RD_median,
              "RD_mean"        = repGap$RD_mean,
              "RD_mode"        = repGap$RD_mode,
              "Phi_index"      = phiIndex
#              "phi2"           = phiIndex$phi2,
#              "phi3"           = phiIndex$phi3,
#              "phi4"           = phiIndex$phi4,
#              "phi5"           = phiIndex$phi5,
#              "phi6"           = phiIndex$phi6,
#              "phi7"           = phiIndex$phi7
              )

  return(out)

}

#' Measures of the classical RNG task
#'
#' Individual functions to get the measures of the classical RNG task. Use function `allRNG()` to get a list with all (or selected) measures.
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
Redundancy <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  R <- 100 * (1 - ((log2(NROW(y)) - 1/NROW(y) * sum(table(y)*log2(table(y)))) / log2(maxScale)))

  attr(R,'Name')     <- 'Redundancy'
  attr(R,'y')        <- y
#  attr(R,'maxScale') <- maxScale

  return(R)
}

# RNG ---------------------------------------------------------------------

#' RNG
#' @rdname classicalRNG
RNG <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  o_e      <- observed_expected(y, minScale = minScale, maxScale = maxScale, lag = 1)

  combis   <- o_e$combis
  observed <- o_e$observed

  Afreq <- table(y)
  freqs <- plyr::ldply(combis, function(c){data.frame(pair=c,freq=sum(observed%in%c))})
  RNG   <- sum(freqs$freq * log(freqs$freq),na.rm = TRUE) / sum(Afreq*log(Afreq), na.rm = TRUE)

  attr(RNG,'Name')     <- 'RNG'
  attr(RNG,'y')        <- y
#  attr(RNG,'minScale') <- minScale
#  attr(RNG,'maxScale') <- maxScale

  return(RNG)
}

# RNG2 --------------------------------------------------------------------

#' RNG2
#' @rdname classicalRNG
RNG2 <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  o_e      <- observed_expected(y, minScale = minScale, maxScale = maxScale, lag = 2)

  combis   <- o_e$combis
  observed <- o_e$observed

  Afreq <- table(y)
  freqs <- plyr::ldply(combis, function(c){data.frame(pair=c,freq=sum(observed%in%c))})
  RNG2  <- sum(freqs$freq * log(freqs$freq),na.rm = TRUE) / sum(Afreq*log(Afreq), na.rm = TRUE)

  attr(RNG2,'Name')     <- 'RNG2'
  attr(RNG2,'y')        <- y
#  attr(RNG2,'minScale') <- minScale
#  attr(RNG2,'maxScale') <- maxScale

  return(RNG2)
}

# Response frequencies ----------------------------------------------------

#' RF
#' @rdname classicalRNG
RF <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  RF <- as.list(table(y))

  return(RF)
}

# Coupon ------------------------------------------------------------------

#' Coupon
#' @rdname classicalRNG
Coupon <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  stp <- FALSE
  cnt <- 0
  cpn <- vector(mode = "integer")

  if(all(purrr::map_lgl(responseAlternatives, function(x) any(which(y == x))))) {
    while (stp == FALSE) {
      cnt      <- cnt + 1
      cpn[cnt] <- first_occurence(y, minScale, maxScale)
      y        <- y[(cpn[cnt]+1):length(y)]
      chk      <- purrr::map_lgl(responseAlternatives, function(x) any(which(y == x)))
      stp      <- any(chk == FALSE)
    }
    Coupon     <- mean(cpn, na.rm = TRUE)
    rm(y, cnt, stp, cpn, chk)
  } else {
      Coupon   <- NULL
    }

  attr(Coupon,'Name')     <- 'Coupon'
#  attr(Coupon,'y')        <- y
  attr(Coupon,'minScale') <- minScale
  attr(Coupon,'maxScale') <- maxScale

  return(Coupon)
}

# NSQ ---------------------------------------------------------------------

#' NSQ
#' @rdname classicalRNG
NSQ <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  o_e      <- observed_expected(y, minScale=minScale, maxScale=maxScale, lag = 1)

  combis   <- o_e$combis
  observed <- o_e$observed

  NS  <- length(combis)-length(unique(observed))
  NSQ <- 100 * (NS / (maxScale^2 - 1))

  return(NSQ)
}

# Adjacency ---------------------------------------------------------------

#' Adjacency
#' @rdname classicalRNG
Adjacency <- function(y, minScale=NA, maxScale=NA){

  y    <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  FOD  <- table(diff(y))
  if(any(diff(y) ==  1)) {ascending  <- 100 * FOD[["1"]]  / length(y)} else {ascending <- 0}
  if(any(diff(y) == -1)) {descending <- 100 * FOD[["-1"]] / length(y)} else {descending <- 0}
#  combined   <- 100 * sum(FOD[["1"]], FOD[["-1"]], na.rm = TRUE) / length(y)
  combined <- sum(ascending, descending)

  return(list(Adj_ascending  = ascending,
              Adj_descending = descending,
              Adj_combined   = combined))
}

# TPI ---------------------------------------------------------------------

#' TPI
#' @rdname classicalRNG
TPI <- function(y, minScale = NA, maxScale = NA){

  y <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  ind      <- findTPIs(y)
  sIndices <- ind$sIndices
  eIndices <- ind$eIndices

  TPIs     <- length(sIndices) + length(eIndices)
  TPI      <- 100 * (TPIs) / (2/3 * (NROW(y) - 2))

  return(TPI)
}

# Phase length ------------------------------------------------------------

#' Phase Lengths
#' @rdname classicalRNG
PhL <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  ind      <- findTPIs(y)
  if(length(ind[[1]]) != length(ind[[2]])) {
    ind <- purrr::map(ind, function(x) {length(x) <- max(lengths(ind)); x})
  }

  sIndices <- ind$sIndices
  eIndices <- ind$eIndices
  PhL      <- table(c(eIndices-sIndices, sIndices[2:length(sIndices)]-eIndices[1:(length(eIndices)-1)]))

  return(PhL)
}

# Runs --------------------------------------------------------------------

#' Runs
#' @rdname classicalRNG
Runs <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  ind      <- findTPIs(y)
  if(length(ind[[1]]) != length(ind[[2]])) {
    ind <- purrr::map(ind, function(x) {length(x) <- max(lengths(ind)); x})
  }

  sIndices <- ind$sIndices
  eIndices <- ind$eIndices
  phases   <- c(eIndices-sIndices, sIndices[2:length(sIndices)]-eIndices[1:(length(eIndices)-1)])

  return(var(phases, na.rm = TRUE))
}

# Repetition Distance -----------------------------------------------------

#' Repetition Distance
#' @rdname classicalRNG
repDistance <- function(y, minScale = NA, maxScale = NA){

  y   <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  tmp <- purrr::map(responseAlternatives, function(v) diff(which(y==v)))
  RD  <- table(unlist(tmp))

  return(RD)
}

# Repetition Gap ----------------------------------------------------------

#' Repetition Gap
#' @rdname classicalRNG
repGap <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  tmp <- purrr::map(responseAlternatives, function(v) diff(which(y==v)))
  RD  <- table(unlist(tmp))

  RD_median <- stats::median(unlist(tmp))
  RD_mean   <- mean(unlist(tmp))
  RD_mode   <- as.numeric(names(which.max(RD)))

  rm(tmp,RD)

  return(list(RD_median = RD_median,
              RD_mean   = RD_mean,
              RD_mode   = RD_mode))
}

# Phi index ---------------------------------------------------------------

#' Phi Index
#'
#' @inheritParams allRNG
#' @param maxOrder

#'
#' @rdname classicalRNG
#' @export
#'
phiIndex <- function(y, minScale = NA, maxScale = NA, responseAlternatives = NA, maxOrder = 7){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale, responseAlternatives = responseAlternatives)
  minIndex <- attr(y,"minScale")
  maxIndex <- attr(y,"maxScale")
  responseAlternatives <- attr(y,"responseAlternatives")

  Nresp   <- NROW(y)
  setSize <- length(responseAlternatives)
  `%00%`  <- invctr::`%00%`

  if(maxOrder>Nresp){
    stop("maxOrder cannot be larger than N")
  }

  PhiArray <- cbind(expand.grid(Phi = 2:maxOrder, RespAlt = minIndex:maxIndex),matrix(0,ncol=maxOrder))

  lookup <- data.frame(DescTools::CombSet(c(0,1),maxOrder,repl=TRUE, ord = TRUE),0)
  colnames(lookup) <- c(paste0("pos",1:maxOrder),"freq")

  phiArrayList    <- phiObsPredList <- phiObsFreqList <- phiList <- list()
  PhiFreq         <- data.frame(respIns = minIndex:maxIndex, freq = 0)
  PhiFreq$freq    <- plyr::laply(PhiFreq$respIns,function(n) sum(y==n, na.rm = TRUE))
  PhiFreq$invFreq <- Nresp-PhiFreq$freq

  for(A in 2:maxOrder){
    phiObsPred <- matrix(0,2,2, dimnames = list(c("obs","pred"),c("repeat","alternate")))
    phiObsFreq <- list()
    phiArray   <- data.frame(DescTools::CombSet(c(0,1),A,repl=TRUE, ord = TRUE),0)
    colnames(phiArray) <- c(paste0("pos",1:A),"freq")
    sameInd <- phiArray[1]==phiArray[A] #phiArray[1]==1&phiArray[A]==1&(phiArray[1]==phiArray[A])
    diffInd <- phiArray[1]!=phiArray[A] #&(phiArray[1]!=0|phiArray[A]!=0)
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
        O <- rev(O)
        lookup$freq[rowSums(t(t(lookup[1:A]) == O))==A] <- lookup$freq[rowSums(t(t(lookup[1:A]) == O))==A] + 1
      } # C series length -A
      phiObsFreq[[B]] <- lookup
      phiArray$freq <- (phiArray$freq + lookup$freq )

      for(r in 1:NROW(lookup)){
        BinString <- lookup[r,1:A]
        if(A==2){
          freq1   <- PhiFreq[B,ifelse(BinString[1]==1,2,3)]
          freq2   <- PhiFreq[B,ifelse(BinString[2]==1,2,3)]
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
        if(BinString[1]==BinString[A]){ #all(BinString[1]==1,BinString[A]==1,BinString[1]==BinString[A])
          phiObsPred[2,1] <- phiObsPred[2,1] + result
        } else {
         # if(all(BinString[1]!=BinString[A])){
            phiObsPred[2,2] <- phiObsPred[2,2] + result
        #  }
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
    a2 <- phiObsPred[1,2]
    a3 <- phiObsPred[2,1]
    a4 <- phiObsPred[2,2]
    Den <- sum(phiObsPred, na.rm = TRUE)

    if(Den > 0){
      AA1 <- ((a1 + a2) * (a1 + a3)) / Den
      AA2 <- ((a2 + a4) * (a1 + a2)) / Den
      AA3 <- ((a3 + a4) * (a1 + a3)) / Den
      AA4 <- ((a2 + a4) * (a3 + a4)) / Den
    }

    B1 <- ((a1-AA1)^2 / AA1)%00%0
    B2 <- ((a2-AA2)^2 / AA2)%00%0
    B3 <- ((a3-AA3)^2 / AA3)%00%0
    B4 <- ((a4-AA4)^2 / AA4)%00%0

    s <- 1
    if(a2>a4){s <- -1}

    phiList[[A]] <- s*(sqrt((B1 + B2 + B3 + B4)/(setSize*Nresp))*100)

    rm(a1,a2,a3,a4,AA1,AA2,AA3,AA4,B1,B2,B3,B4,s, Den)

  } # A

  names(phiList) <- c(NA,paste0("phi",2:maxOrder))
  phiList[[1]] <- NULL

  return(phiList)
}
