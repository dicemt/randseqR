
# HELPERS ----

#' Find peaks/troughs
#'
#' @param y Numeric vector
#'
#' @return converted y with values 1 (trough); 2 (decrease); 3 (same); 5 (increase); 5 (peak)
#'
#' @export
#'
TPIhelper <- function(y){

  TT <- length(y)
  su <- vector(mode = "numeric", length = length(y))
  su[1] <- NA

  # sym_num conversion ------------------------------------------------------
  for (t in 2:(TT)) {
    if(t != TT) {
      if (y[t] == y[t - 1] | y[t] == y[t + 1]) {
        su[t] <- 3
      }
      if (y[t] > y[t - 1]) {
        if (y[t] > y[t + 1]) {
          su[t] <- 5
        }
        else {
          su[t] <- 4
        }
      }
      if (y[t] < y[t - 1]) {
        if (y[t] < y[t + 1]) {
          su[t] <- 1
        }
        else {
          su[t] <- 2
        }
      }
    }
    else {
      if (y[t - 1] > y[t]) {
        su[t] <- 2
      }
      if (y[t - 1] < y[t]) {
        su[t] <- 4
      }
      if (y[t - 1] == y[t]) {
        su[t] <- 3
      }
    }
  }

  # sym_label conversion with plateaus -------------------------------------
  sa <- dplyr::case_when(su == 1 ~ "trough", su == 2 ~ "decrease", su == 3 ~ "same",
                         su == 4 ~ "increase", su == 5 ~ "peak", is.na(su) ~ NA_character_)

  ou   <- su
  oa   <- sa
  i    <- 1
  same <- 0

  while (i <= length(sa)) {
    if (sa[i] %in% c("increase", "decrease")) {
      samesame <- TRUE
      r <- i + 1
      same <- 0
      while (samesame) {
        if (sa[r] %in% "same") {
          same <- same + 1
          r <- r + 1
        }
        else {
          samesame <- FALSE
        }
      }
      if (same > 0) {
        if (all(!sa[i] %in% c("increase"), sa[i + same + 1] %in% c("peak", "increase"))) {
          oa[i:(i + same)] <- "trough"
          ou[i:(i + same)] <- 1
        }
        else {
          if (all(!sa[i] %in% c("decrease"), sa[i + same + 1] %in% c("trough", "decrease"))) {
            oa[i:(i + same)] <- "peak"
            ou[i:(i + same)] <- 5
          }
        }
      }
    }
    if (same > 0) {
      i <- (i + same)
      same <- 0
    }
    else {
      i <- i + 1
    }
  }

  return(ou)

}

#' Find turning points
#'
#' @param y Numeric vector.
#'
#' @return A list object with `start` and `end` indices of turning points.
#'
#' @export
#'
findTPIs <- function(y){


  # y <- casnet::ts_symbolic(y, usePlateaus = TRUE)
  # y <- attr(y, "sym_numeric")

  z <- TPIhelper(y)

  startWith <- which.min(c(which(z == 5)[1], which(z == 1)[1]))
  # startWith <- which.min(c(which(y %in% c("peak")[1]), which(y %in% c("trough")[1])))

  yPeak          <- z == 5 # %in% c("peak")
  yTrough        <- z == 1 # %in% c("trough")
  # names(yPeak)   <- paste0(1:length(yPeak), "peak")
  # names(yTrough) <- paste0(1:length(yTrough), "trough")
  ind            <- list(which(yPeak[1:(NROW(y)-1)] - yPeak[2:NROW(y)]>0),
                          which(yTrough[1:(NROW(y)-1)] - yTrough[2:NROW(y)]>0))


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
#' @return Transformed `y` (if necessary), or, an error message.
#'
#' @examples
#'
#' y <- sample(letters,10)
#' check_y(y, responseAlternatives=letters)
#'
check_y <- function(y, minScale=NA, maxScale=NA, responseAlternatives = NA,
                    noZero = FALSE, toNumeric = TRUE){

  transformation <- character(0)

  if(any(is.na(minScale),is.na(maxScale))&&any(is.na(responseAlternatives))){
    stop("Must provide either minScale & maxScale, and/or, responseAlternatives")
  }


  #  if(is.na(responseAlternatives)&&all(!is.na(minScale),!is.na(maxScale))){
  #    responseAlternatives <- minScale:maxScale
  #  } else {
  #    stop("Must provide minScale and maxScale to generate responseAlternatives")
  #  }


  if(any(is.na(responseAlternatives))){
    responseAlternatives <- minScale:maxScale
  } else {
    responseAlternatives <- casnet::as.numeric_discrete(responseAlternatives)
  }


  if(any(is.na(minScale), is.na(maxScale))){
    minScale <- as.numeric(responseAlternatives[1])
    maxScale <- max(responseAlternatives, na.rm = TRUE)
  }


  if(!is.numeric(y)){
    if(toNumeric){
      y <- casnet::as.numeric_discrete(y)
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
      responseAlternatives <- casnet::as.numeric_discrete(as.character(responseAlternatives))
      transformation <- c(transformation,"toCharacter")
    }
  }


   #if(!all(y%in%responseAlternatives)){
   #  stop("Elements in y do not correspond to elements of responseAlternatives!")
   #}


   #if(length(minScale:maxScale)!=length(responseAlternatives)){
   #  stop("length of minScale:maxScale must be equal to length of responseAltenatives")
   #}


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
observed_expected <- function(y, minScale, maxScale, lag){

  combi    <- DescTools::CombSet(x=minScale:maxScale, m=2, ord=TRUE, repl=TRUE)
  combis   <- paste0(combi[,1],".",combi[,2])
  observed <- paste(interaction(y[1:(NROW(y)-lag)], y[(lag+1):NROW(y)]))

  #observed <- paste(interaction(y, dplyr::lead(y, 1, default = y[1])))

  return(list(combis   = as.numeric(combis),
              observed = as.numeric(observed)))
}


#' Get all classical RNG measures
#'
#' @param y A sequence of symbols. If `y` is non-numeric, unique elements will be labelled by an integer value.
#' @param minScale Minimum expected value. If `y` is a character vector this should refer to the lowest numeric code used.
#' @param maxScale Maximum expected value. If `y` is a character vector this should refer to the highest numeric code used.
#' @param responseAlternatives An optional vector of possible response alternatives. If `NA`, `responseAlternatives` will be set to `seq(minScale,maxScale)`
#' @param results Either 'randseqR' (default) or 'classical'.
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
                   minScale = NA,
                   maxScale = NA,
                   responseAlternatives = NA,
                   results = c("classical", "randseqR")[2],
                   Redundancy  = TRUE,
                   RNG         = TRUE,
                   RNG2        = TRUE,
                   RF          = TRUE,
                   Coupon      = TRUE,
                   NSQ         = TRUE,
                   FOD         = TRUE,
                   Adjacency   = TRUE,
                   TPI         = TRUE,
                   PhL         = TRUE,
                   Runs        = TRUE,
                   repDistance = TRUE,
                   repGap      = TRUE,
                   PhiIndex    = TRUE
                   ){

#disp('not implemented yet')
  y        <- check_y(y, minScale = minScale, maxScale = maxScale, responseAlternatives = responseAlternatives)
  minScale <- attr(y,"minScale")
  maxScale <- attr(y,"maxScale")
  responseAlternatives <- attr(y,"responseAlternatives")

  if(Redundancy){
    Redundancy <- Redundancy(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(RNG){
    RNG <- RNG(y = y, minScale = minScale, maxScale = maxScale, results = results)
  }

  if(RNG2){
    RNG2 <- RNG2(y = y, minScale = minScale, maxScale = maxScale, results = results)
  }

  if(RF){
    RF <- RF(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(Coupon){
    Coupon <- Coupon(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(NSQ){
    NSQ <- NSQ(y = y, minScale = minScale, maxScale = maxScale, results = results)
  }

  if(FOD){
    FOD <- FOD(y = y, minScale = minScale, maxScale = maxScale, results = results)
  }

  if(Adjacency){
    Adjacency <- Adjacency(y = y, minScale = minScale, maxScale = maxScale,
                           results = results)
  }

  if(TPI){
    TPI <- TPI(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(PhL){
    PhL <- PhL(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(Runs){
    Runs <- Runs(y = y, minScale = minScale, maxScale = maxScale, results = results)
  }

  if(repDistance){
    repDistance <- repDistance(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(repGap){
    repGap <- repGap(y = y, minScale = minScale, maxScale = maxScale)
  }

  if(PhiIndex){
    phiIndex <- phiIndex(y = y, minScale = minScale, maxScale = maxScale)
  }

  out <- list("Results"     = results,
              "N"           = length(y),
              "Redundancy"  = Redundancy,
              "RNG"         = RNG,
              "RNG2"        = RNG2,
              "RF"          = RF,
              "Coupon"      = Coupon,
              "NSQ"         = NSQ,
              "FOD"         = FOD,
              "Adjacency"   = Adjacency,
              "TPI"         = TPI,
              "phaselength" = PhL,
              "Runs"        = Runs,
              "repDistance" = repDistance,
              "repGap"      = repGap,
              "PhiIndex"    = phiIndex
              )

  return(out)

}


#' Measures of the classical RNG task
#'
#' Individual functions to get the measures of the classical RNG task (see Details below). Use function `allRNG()` to get a list with all (or selected) measures.
#'
#' @param y A sequence of symbols. If `y` is non-numeric, unique elements will be labelled by an integer value.
#' @param minScale Minimum expected value. If `y` is a character vector this should refer to the lowest numeric code used.
#' @param maxScale Maximum expected value. If `y` is a character vector this should refer to the highest numeric code used.
#' @param results either 'randseqR' (default) or 'classical. randseqR gives a better consistency among RNG measures, while classical gives output comparible to RGcalc by Towse and Neil (1998)
#'
#' @details Avialable classical measures:
#' - **Redundancy**: Measures how,..
#' - **RNG**: Measures of
#' - **RNG2**: blah
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
# > NULL

# Redundancy --------------------------------------------------------------


#' Redundancy
#'
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
Redundancy <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  R <- 100 * (1 - ((log2(NROW(y)) - 1/NROW(y) * sum(table(y)*log2(table(y)))) / log2(maxScale)))

  attr(R,'Name')     <- 'Redundancy'
  attr(R,'y')        <- y

  return(R)
}

# RNG ---------------------------------------------------------------------

#' RNG

#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
RNG <- function(y, minScale = NA, maxScale = NA,
                results = c("classical", "randseqR")[2]){


  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  o_e      <- observed_expected(y, minScale = minScale, maxScale = maxScale, lag = 1)
  combis   <- o_e$combis


  if(results == "randseqR") {
    observed <- c(o_e$observed)
    Afreq <- table(y)
    freqs <- purrr::map_dfr(combis, function(c) {data.frame(pair=c, freq=sum(observed%in%c))})
    RNG   <- 100 * sum(freqs$freq * log(freqs$freq),na.rm = TRUE) / sum(Afreq*log(Afreq), na.rm = TRUE)
  } else {
    observed <- c(o_e$observed, paste0(dplyr::last(y),".",dplyr::first(y)))
    Afreq <- table(y)
    freqs <- purrr::map_dfr(combis, function(c) {data.frame(pair=c, freq=sum(observed%in%c))})
    RNG   <- sum(freqs$freq * log(freqs$freq),na.rm = TRUE) / sum(Afreq*log(Afreq), na.rm = TRUE)
  }

  attr(RNG,'Name')     <- 'RNG'
  attr(RNG,'y')        <- y


  attr(RNG,'results')  <- results

#  attr(RNG,'minScale') <- minScale
#  attr(RNG,'maxScale') <- maxScale


  return(RNG)
}

# RNG2 --------------------------------------------------------------------

#' RNG2
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
RNG2 <- function(y, minScale = NA, maxScale = NA,
                 results = c("classical", "randseqR")[2]){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  o_e      <- observed_expected(y, minScale = minScale, maxScale = maxScale, lag = 2)

  combis   <- o_e$combis
  observed <- o_e$observed

  Afreq <- table(y)
  freqs <- purrr::map_dfr(combis, function(c) {data.frame(pair=c,freq=sum(observed%in%c))})

  if(results == "randseqR") {
    RNG2  <- 100 * sum(freqs$freq * log(freqs$freq),na.rm = TRUE) / sum(Afreq*log(Afreq), na.rm = TRUE)
  } else {
    RNG2  <- sum(freqs$freq * log(freqs$freq),na.rm = TRUE) / sum(Afreq*log(Afreq), na.rm = TRUE)
  }

  attr(RNG2,'Name')     <- 'RNG2'
  attr(RNG2,'y')        <- y
  attr(RNG2, 'results') <- results

  return(RNG2)
}

# Response frequencies ----------------------------------------------------

#' RF
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
RF <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  RFList        <- as.list(table(y))
  names(RFList) <- c(paste0("RF_", minScale:maxScale))

  attr(RFList, 'Name') <- 'RF'
  attr(RFList, 'y')    <- y

  return(RFList)
}

# Coupon ------------------------------------------------------------------

#' Coupon
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
Coupon <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  stp <- FALSE
  cnt <- 0
  cpn <- vector(mode = "integer")
  tmp <- y

  if(all(purrr::map_lgl(responseAlternatives, function(x) any(which(y == x))))) {
    while (stp == FALSE) {
      cnt      <- cnt + 1
      cpn[cnt] <- first_occurence(tmp, minScale, maxScale)
      tmp      <- tmp[(cpn[cnt]+1):length(tmp)]
      chk      <- purrr::map_lgl(responseAlternatives, function(x) any(which(tmp == x)))
      stp      <- any(chk == FALSE)
    }
    Coupon     <- mean(cpn, na.rm = TRUE)
    rm(tmp, cnt, stp, cpn, chk)
  } else {
      Coupon   <- NaN
    }

  attr(Coupon,'Name')     <- 'Coupon'
  attr(Coupon,'y')        <- y

  return(Coupon)
}

# NSQ ---------------------------------------------------------------------

#' NSQ
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
NSQ <- function(y, minScale = NA, maxScale = NA,
                results = c("classical", "randseqR")[2]){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  o_e      <- observed_expected(y, minScale=minScale, maxScale=maxScale, lag = 1)

  combis   <- o_e$combis

  if(results == "randseqR"){
    observed <- o_e$observed
  } else {
    observed <- c(o_e$observed, paste0(dplyr::last(y),".",dplyr::first(y)))
  }

  NS  <- length(combis)-length(unique(observed))
  NSQ <- 100 * (NS / (maxScale^2 - 1))

  attr(NSQ,'Name')     <- 'NSQ'
  attr(NSQ,'y')        <- y
  attr(NSQ, 'results') <- results

  return(NSQ)
}


# First-order difference --------------------------------------------------

#' First-order difference
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
FOD <- function(y, minScale = NA, maxScale = NA,
                results = c("classical", "randseqR")[2]) {

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')

  if(results == "randseqR") {
    FOD            <- table(diff(y))
    FODList        <- as.list(FOD)
    names(FODList) <- c(paste0("Diff_", names(FOD)))
  } else {
    FOD            <- table(c(diff(y), dplyr::first(y) - dplyr::last(y)))
    FODList        <- as.list(FOD)
    names(FODList) <- c(paste0("Diff_", names(FOD)))
  }


  attr(FODList, 'Name')     <- 'FOD'
  attr(FODList, 'y')        <- y
  attr(FODList, 'results')  <- results

  return(FODList)
}

# Adjacency ---------------------------------------------------------------

#' Adjacency
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
Adjacency <- function(y, minScale = NA, maxScale = NA,
                      results = c("classical", "randseqR")[2]){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  if(results == "randseqR") {
    FOD <- table(diff(y))
  } else {
    FOD <- table(c(diff(y), dplyr::first(y) - dplyr::last(y)))
  }

  if(any(diff(y) ==  1)) {ascending  <- 100 * FOD[["1"]]  / length(y)} else {ascending <- 0}
  if(any(diff(y) == -1)) {descending <- 100 * FOD[["-1"]] / length(y)} else {descending <- 0}
#  combined   <- 100 * sum(FOD[["1"]], FOD[["-1"]], na.rm = TRUE) / length(y)
  combined <- sum(ascending, descending)

  attr(Adjacency,'Name')     <- 'Adjacency'
  attr(Adjacency,'y')        <- y
  attr(Adjacency, 'results') <- results

  return(list(Ascending  = ascending,
              Descending = descending,
              Combined   = combined))
}

# TPI ---------------------------------------------------------------------

#' TPI
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
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

  attr(TPI,'Name') <- 'TPI'
  attr(TPI,'y')    <- y

  return(TPI)
}

# Phase length ------------------------------------------------------------

#' Phase Lengths
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
PhL <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  ind     <- findTPIs(y)
  if(diff(lengths(ind)) == 0) {
    PhL <- c(ind$eIndices-ind$sIndices, ind$sIndices[2:length(ind$sIndices)]-ind$eIndices[1:(length(ind$eIndices)-1)])
    PhL <- table(PhL)
  } else {
    idMax        <- which.max(c(length(ind$sIndices), length(ind$eIndices)))
    idMin        <- which.min(c(length(ind$sIndices), length(ind$eIndices)))
    ind[[idMin]] <- c(ind[[idMin]], ind[[idMax]][length(ind[[idMax]])])

    PhL <- c(ind$eIndices-ind$sIndices, ind$sIndices[2:length(ind$sIndices)]-ind$eIndices[1:(length(ind$eIndices)-1)])
    PhL <- table(PhL[! PhL%in%0])
  }

  PhLList        <- as.list(PhL)
  names(PhLList) <- c(paste0("PhL_", names(PhL)))

  attr(PhLList,'Name') <- 'PhL'
  attr(PhLList,'y')    <- y

  return(PhLList)
}

# Runs --------------------------------------------------------------------

#' Runs
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
Runs <- function(y, minScale = NA, maxScale = NA,
                 results = c("classical", "randseqR")[2]){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  ind     <- findTPIs(y)
  if(results == "randseqR") {
    if(diff(lengths(ind)) == 0) {
      R    <- c(ind$eIndices-ind$sIndices, ind$sIndices[2:length(ind$sIndices)]-ind$eIndices[1:(length(ind$eIndices)-1)])
      Runs <- var(R, na.rm = TRUE)
      } else {
        idMax        <- which.max(c(length(ind$sIndices), length(ind$eIndices)))
        idMin        <- which.min(c(length(ind$sIndices), length(ind$eIndices)))
        ind[[idMin]] <- c(ind[[idMin]], ind[[idMax]][length(ind[[idMax]])])

        R    <- c(ind$eIndices-ind$sIndices, ind$sIndices[2:length(ind$sIndices)]-ind$eIndices[1:(length(ind$eIndices)-1)])
        R    <- R[! R%in%0]
        Runs <- var(R, na.rm = TRUE)
        }
    } else {
      warning("Unable to replicate classical Runs output, returns NA")
      Runs <- NA
    }

  attr(Runs,'Name')     <- 'Runs'
  attr(Runs,'y')        <- y
  attr(Runs, 'results') <- results

  return(Runs)

}

# Repetition Distance -----------------------------------------------------

#' Repetition Distance
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
repDistance <- function(y, minScale = NA, maxScale = NA){

  y   <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  tmp <- purrr::map(responseAlternatives, function(v) diff(which(y==v)))
  RD  <- table(unlist(tmp))

  RDList        <- as.list(RD)
  names(RDList) <- c(paste0("RD_", names(RD)))

  attr(RDList,'Name') <- 'repetition Distance'
  attr(RDList,'y')    <- y

  return(RDList)
}

# Repetition Gap ----------------------------------------------------------

#' Repetition Gap
#'
#' @inheritParams allRNG
#'
#' @rdname classicalRNG
#' @export
#'
repGap <- function(y, minScale = NA, maxScale = NA){

  y        <- check_y(y, minScale = minScale, maxScale = maxScale)
  minScale <- attr(y,'minScale')
  maxScale <- attr(y,'maxScale')
  responseAlternatives <- attr(y,"responseAlternatives")

  tmp <- purrr::map(responseAlternatives, function(v) diff(which(y==v)))
  RD  <- table(unlist(tmp))

  median <- stats::median(unlist(tmp))
  mean   <- mean(unlist(tmp))
  mode   <- as.numeric(names(which.max(RD)))

  rm(tmp,RD)

  attr(repGap,'Name') <- 'repGap'
  attr(repGap,'y')    <- y

  return(list(RG_median = median,
              RG_mean   = mean,
              RG_mode   = mode))
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
phiIndex <- function(y, minScale = NA, maxScale = NA, responseAlternatives = NA,
                     maxOrder = 7){

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

  names(phiList) <- c(NA,paste0("Phi_",2:maxOrder))
  phiList[[1]] <- NULL

  attr(phiList,'Name') <- 'phiIndex'
  attr(phiList,'y')  <- y

  return(phiList)
}

