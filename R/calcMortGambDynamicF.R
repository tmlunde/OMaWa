calcGambMortF <- function (t2, RH, bovine, human, HBI, d = 0:199, mSiz = NULL, gamma=6){
aCof <- c(1, 3, 6, 10, 15, 21, 28, 36, 45) - 1
bCof <- c(2, 5, 9, 14, 20, 27, 35, 44, 200) - 1 
cCof <- bCof-aCof
#t2 <- runif(10000, 14, 40)
m1 <- m2 <- m3 <- m4 <- m5 <- m6 <- m7 <- m8 <- m9 <- t2*0
mrt <- d*0
nR <- length(t2)
mR <- length(d)
f <- factorial(0:(gamma-1))
kR <- length(f)

if (is.null(mSiz)) {
 mSiz <- rep(3.050182, length(t2))
}

#dyn.load("src/gambMort.so")

outMrt <- .C("gambMort", as.double(t2), as.double(RH), as.integer(nR), as.double(cCof), 
              as.double(d), as.integer(mR), m1=as.double(m1), m2=as.double(m2), m3=as.double(m3), 
              m4=as.double(m4), m5=as.double(m5), m6=as.double(m6), 
              m7 = as.double(m7), m8 = as.double(m8), m9=as.double(m9), mrt = as.double(mrt), as.double(bovine), as.double(human), as.double(HBI), 
              f = as.double(f), kR = as.integer(kR), myK = as.double(0:(gamma-1)), mSiz = mSiz)
#outMrt[t2>40] <- 1
#outMrt[outMrt>1] <- 1

mort.arabTmp <- data.frame(m1=outMrt$m1, m2=outMrt$m2, m3=outMrt$m3, m4=outMrt$m4, m5=outMrt$m5,
          m6=outMrt$m6, m7=outMrt$m7, m8=outMrt$m8, m9=outMrt$m9)
for (i in 1:9) {
mort.arabTmp[,i][t2>40] <- 1
mort.arabTmp[,i][t2>39] <- .5

}

#dyn.unload("src/gambMort.so")

#return(mort.arabTmp)
return(outMrt)
}

