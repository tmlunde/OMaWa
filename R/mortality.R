calcArabMort <- function (tMax, bovine, human, HBI){
 aCof <- c(1, 3, 6, 10, 15, 21, 28, 36, 45) - 1
 bCof <- c(2, 5, 9, 14, 20, 27, 35, 44, 200) - 1 
 cCof <- bCof-aCof
 d <- 1:200

 m1 <- m2 <- m3 <- m4 <- m5 <- m6 <- m7 <- m8 <- m9 <- tMax*0
 mrt <- d*0
 nR <- length(tMax)
 mR <- length(d)
 #dyn.load("src/arabMort.so")

 outMrt <- .C("arabMort", as.double(tMax), as.integer(nR), as.double(cCof), 
              as.double(d), as.integer(mR), m1=as.double(m1), m2=as.double(m2), 
              m3=as.double(m3), 
              m4=as.double(m4), m5=as.double(m5), m6=as.double(m6), 
              m7 = as.double(m7), m8 = as.double(m8), m9=as.double(m9), 
              as.double(mrt), as.double(bovine), as.double(human), as.double(HBI))

 mort.arabTmp <- data.frame(m1=outMrt$m1, m2=outMrt$m2, m3=outMrt$m3, m4=outMrt$m4, m5=outMrt$m5,
           m6=outMrt$m6, m7=outMrt$m7, m8=outMrt$m8, m9=outMrt$m9)
 for (i in 1:9) {
  mort.arabTmp[,i][tMax>40] <- 1
  mort.arabTmp[,i][tMax>39] <- .5

 }

 #dyn.unload("src/arabMort.so")
 
 return(mort.arabTmp)
}


################################################################################

calcGambMort2 <- function (t2, RH, bovine, human, HBI, d = 0:199, mSiz = NULL, mType = "a"){
 aCof <- c(1, 3, 6, 10, 15, 21, 28, 36, 45) - 1
 bCof <- c(2, 5, 9, 14, 20, 27, 35, 44, 200) - 1 
 cCof <- bCof-aCof

 m1 <- m2 <- m3 <- m4 <- m5 <- m6 <- m7 <- m8 <- m9 <- t2*0
 mrt <- d*0
 nR <- length(t2)
 mR <- length(d)
 f <- factorial(0:5)
 kR <- length(f)

 if (is.null(mSiz)) {
  mSiz <- rep(3.050182, length(t2))
 }

# dyn.load("src/gambMort.so")

outMrt <- .C("gambMort", as.double(t2), as.double(RH), as.integer(nR), as.double(cCof), 
              as.double(d), as.integer(mR), m1=as.double(m1), m2=as.double(m2), 
              m3=as.double(m3), 
              m4=as.double(m4), m5=as.double(m5), m6=as.double(m6), 
              m7 = as.double(m7), m8 = as.double(m8), m9=as.double(m9), 
              mrt = as.double(mrt), as.double(bovine), as.double(human), as.double(HBI), 
              f = as.double(f), kR = as.integer(kR), myK = as.double(0:5), mSiz = mSiz)
#outMrt[t2>40] <- 1
#outMrt[outMrt>1] <- 1

mort.arabTmp <- data.frame(m1=outMrt$m1, m2=outMrt$m2, m3=outMrt$m3, m4=outMrt$m4, m5=outMrt$m5,
          m6=outMrt$m6, m7=outMrt$m7, m8=outMrt$m8, m9=outMrt$m9)

if (mType == "a") {
 for (i in 1:9) {
  mort.arabTmp[,i][t2>40] <- 1
  mort.arabTmp[,i][t2>39] <- .5
 }
} else if (mType == "g") {
 for (i in 1:9) {
  mort.arabTmp[,i][t2>40] <- 3
  mort.arabTmp[,i][t2>39] <- 1
}


}

#dyn.unload("src/gambMort.so")

#return(mort.arabTmp)
return(outMrt)
}

################################################################################

mortGambAq <- function (tslb, arab, gamb) {
    myFrac <- arab/(gamb)
  myFrac[myFrac > 1] <- 1
  out <- (0.002404075*(tslb)^2  -0.1127944*(tslb) + 1.337783)
  out[tslb>=25 & tslb <= 35] <- out[tslb>=25 & tslb <= 35]*
                               (.4+.6*(1+sin(-10.9956+0.3142*tslb[tslb>=25 & tslb <= 35])))^myFrac[tslb>=25 & tslb <= 35]
  out[is.na(out)] <- 3
  out
}


################################################################################

mortArabAq <- function (tslb, arab, gamb) {
  myFrac <- gamb/(arab)
  myFrac[myFrac > 1] <- 1
  out <- (0.0006556736*(tslb)^2  -0.02980226*(tslb) + 0.3587285)
  out[tslb < 21.91209] <- (0.002404075*(tslb[tslb < 21.91209])^2  -0.1127944*(tslb[tslb < 21.91209]) + 1.337783)
  out[tslb>=25 & tslb <= 35] <- out[tslb>=25 & tslb <= 35]*
                               ((2+cos(-18.8496 + 0.6283*tslb[tslb>=25 & tslb <= 35]))^(log10(2.842150)/log10(3)))^myFrac[tslb>=25 & tslb <= 35]
  out[is.na(out)] <- 3

  out
}

################################################################################


aqDev <- function (Temperature, stage) {
  t.devel <- matrix(c(1.018064, 4.381969, 6.387254, 8.722533, 10.200527, 10.637424,
                    20.148304, 28.382355, 162.568299, 1136.370873, 13.585902, 20.387559,
                    12.096000, 14.134000, 9.402000, 7.038000, 20.742000, 19.759000,
                     4.839000, 5.425000, 4.658000, 5.174000, 8.946000, 6.827000), 6, 4)

  rownames(t.devel) <- c("Egg", "First instar", "Second instar", "Third instar",
                      "Fourth instar", "Pupa")
  colnames(t.devel) <- c("a.t", "b.t", "c.t", "d.t")

  if (stage == 1) { dev.time <- t.devel[stage,1] + t.devel[stage,2]/( 1 + ((Temperature) / t.devel[stage,3] ) ^ t.devel[stage,4])
  } else {dev.time <- (t.devel[stage,1] + t.devel[stage,2]/( 1 + ((Temperature) / t.devel[stage,3] ) ^ t.devel[stage,4])) -
            (t.devel[stage-1,1] + t.devel[stage-1,2]/( 1 + ((Temperature) / t.devel[stage-1,3] ) ^ t.devel[stage-1,4]))}
  return(dev.time)
}


################################################################################

gambExtraDev <- function (arab, gamb) {
  myFrac <- 100*arab/(gamb+arab)
  myFrac[myFrac > 75] <- 75
  gambDelay <- c(9.5, 8.9)/9.5
  1/(1+myFrac*(-0.0008421))
}

################################################################################

arabExtraDev <- function (arab, gamb) {
  myFrac <- 100*gamb/(gamb+arab)
  myFrac[myFrac > 75] <- 75
  arabDelay <- c(10.6, 12.3)/10.6
  1/(1+myFrac*0.002138)
}

################################################################################


################################################################################

gonHosh <- function (T){
  Dd <- 37
  Tc <- 7.7
  out <- 1/(1+Dd/(T-Tc))
if (any(T<=7.7)){
  out[T<=7.7] <- 1/(1+Dd/(7.701-Tc))
}
  out
}


gonot3 <- function (T){
a1 = 1.713413
a2 = 544347.6
a3 = 3.929327
  frac <- function (TT) sapply(TT, function (x) min(max(-0.66667+ 0.03333 * x, 0), .5))
  out <- a1 + a2*((T^(-a3)))
  if (any(T<=0)) out[T <= 0] <- 1/gonHosh(7.701)
  if (any(T>=20)){
    out[T >= 20] <- 1/gonHosh(T[T >= 20]) * frac(T[T >= 20]) + (a1 + a2*((T[T >= 20]^(-a3))))*(1-frac(T[T >= 20]))
  }
  1/out
}


################################################################################



#  pfF1 <- function (T) {
#    a <- 19.2088799
#    b <- 39.6831875
#    c <- 0.1459128
#    log((T - a)/b)/(-c)
#  }

#pfF2 <- function (T) ( 149.1201/(T-13.7435))

#pfF2 <- function (T) (111 / (T - 16))

pfF <- function (T) 9.5907 + 0.0051029/(exp(T^0.7349 - 17.0325))


################################################################################

































































