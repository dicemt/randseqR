# RNG_functions

get.sessions <- function(listID = 1:10, flist, angularTS = TRUE, latencyTS = TRUE, sessionN, subjectN){
  out <- list()
  ang <- list()
  rt  <- list()

  for(s in listID){

    sessionNR  <- sessionN[s]
    subjectNR  <- subjectN[s]

    df         <- import(file = flist[[s]])
    id         <- diff(c(0,which(!(is.na(df$ButtonPressed)))))
    trial      <- sapply(seq_along(id), function(l) rep(l,length=id[[l]]))
    df$subject <- subjectNR
    df$session <- sessionNR
    df$trial   <- unlist(trial)

    id  <- which(!(is.na(df$ButtonPressed)))
    if(latencyTS){
      rt[[s]] <- cbind.data.frame(
        subjectNR,sessionNR,
        trial   = df$trial[id],
        number  = as.numeric(df$ButtonPressed[id]),
        rt      = diff(df$TimeStamp[c(1,id)])
      )
    }

    df  <- df[-id,]

    if(angularTS){
      ang[[s]] <- ldply(1:550, function(t) cbind.data.frame(
        df[df$trial==t&df$session==sessionNR,],
        angularV(x=df$MouseX[df$trial==t&df$session==sessionNR],
                 y=df$MouseY[df$trial==t&df$session==sessionNR],
                 t=df$TimeStamp[df$trial==t&df$session==sessionNR])))
    }

    out[[s]] <- df

  }

  return(list(xy   = ldply(out),
              rt   = ldply(rt),
              angV = ldply(ang)
  )
  )
}

angularV <- function(x,y,t=NULL){

  if(is.null(t)){t <- rep(1,nrow(as.matrix(x)))}

  x <- diff(x) #out$MouseX[out$trial==1])
  y <- diff(y) #out$MouseY[out$trial==1])
  t <- diff(t) #out$TimeStamp[out$trial==1])

  # To polar
  radius <- sqrt(x^2 + y^2)
  theta  <- atan2(y, x)

  omega <- theta
  omega[theta>0] <-theta[theta>0]/t[theta>0]

  return(cbind.data.frame(d.x = c(0,x), d.y = c(0,y), d.t = c(0,t),
                          r= c(0,radius), angle = c(0,theta), d.angle = c(0,omega))
  )
}

get.entropy <- function(data){

  out <- list()

  for(s in unique(data$session)){

    out[[s]] <- ldply(1:550, function(t) sample_entropy(data$angle[data$trial==t&data$session==s]))
  }
  return(ldply(out))
}
