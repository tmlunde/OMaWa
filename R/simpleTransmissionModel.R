
simpleMalaria <- function (Tair = 20, RH = 80, siz = 3, fracHumInf = .01, times = seq(0, 70, .05)){

 if (length(Tair) > 1 | length(RH)>1 | length(siz)>1 ) {
  stop("Only one temperature, RH and size at the time")
 }

 if (RH > 100){
  stop("RH > 100%. Supersatiration?")
 }

 if (Tair > 70){
  stop("Air temperature should be given in C")
 }

 if (Tair < 10){
  warning("Air temperature is low. Might have problems solving this")
 }

 if (fracHumInf > 1 | fracHumInf < 0){
  stop("Fraction of infectious humans should be between 0 and 1")
 }

if (siz < 2 | siz > 5) {
 stop("Size should be given in mm measured as wing length. Are you providing reasonable values?")
}



prms <- calcGambMort(Tair, RH, 1, 1, 1, d = 0:199, mSiz = siz)

aCof <- c(1, 3, 6, 10, 15, 21, 28, 36, 45) - 1
bCof <- c(2, 5, 9, 14, 20, 27, 35, 44, 200) - 1 
cCof <- bCof-aCof
pras <- 1/cCof

     ## Parameters 
     parms  <- c(m1 = prms$m1, m2 = prms$m2, m3 = prms$m3,
                 m4 = prms$m4, m5 = prms$m5, m6 = prms$m6,
                 m7 = prms$m7, m8 = prms$m8, m9 = prms$m9, 
                 a1 = pras[1], a2 = pras[2], a3 = pras[3],
                 a4 = pras[4], a5 = pras[5], a6 = pras[6],
                 a7 = pras[7], a8 = pras[8], a9 = pras[9], gonRate = gonot3(Tair), pf = 1/pfF(Tair), 
                 fracHumInf = fracHumInf)
     
     ## vector of timesteps
     y <- xstart <- c(M1 = 1000, M2 = 0, M3 = 0, M4 = 0, M5 = 0, M6 = 0, M7 = 0, M8 = 0, M9 = 0, 
                      Md2 = 0, Md3 = 0, Md4 = 0, Md5 = 0, Md6 = 0, Md7 = 0, Md8 = 0, Md9 = 0, 
                      Mi2 = 0, Mi3 = 0, Mi4 = 0, Mi5 = 0, Mi6 = 0, Mi7 = 0, Mi8 = 0, Mi9 = 0)

     
     ## Solving
     out3 <-  lsoda(xstart, times, simpleOMaWa, parms) 
out3
}



simpleOMaWa <- function(t, x, parms) {
       with(as.list(c(parms, x)), {
         dM1 <- -(m1+a1)*M1     
         dM2 <- a1*M1 -  (m2+a2 +gonRate*fracHumInf)*M2   
         dM3 <- a2*M2 -  (m3+a3+gonRate*fracHumInf)*M3   
         dM4 <- a3*M3 -  (m4+a4+gonRate*fracHumInf)*M4   
         dM5 <- a4*M4 -  (m5+a5+gonRate*fracHumInf)*M5   
         dM6 <- a5*M5 -  (m6+a6+gonRate*fracHumInf)*M6   
         dM7 <- a6*M6 -  (m7+a7+gonRate*fracHumInf)*M7   
         dM8 <- a7*M7 -  (m8+a8+gonRate*fracHumInf)*M8   
         dM9 <- a8*M8 -  (m9+a9+gonRate*fracHumInf)*M9   
         
         
         dMd2 <- (gonRate*fracHumInf)*M2 - (m2+a2 +pf)*Md2   
         dMd3 <- (gonRate*fracHumInf)*M3 + a2*Md2 -  (m3+a3 +pf)*Md3   
         dMd4 <- (gonRate*fracHumInf)*M4 + a3*Md3 -  (m4+a4 +pf)*Md4   
         dMd5 <- (gonRate*fracHumInf)*M5 + a4*Md4 -  (m5+a5 +pf)*Md5   
         dMd6 <- (gonRate*fracHumInf)*M6 + a5*Md5 -  (m6+a6 +pf)*Md6   
         dMd7 <- (gonRate*fracHumInf)*M7 + a6*Md6 -  (m7+a7 +pf)*Md7   
         dMd8 <- (gonRate*fracHumInf)*M8 + a7*Md7 -  (m8+a8 +pf)*Md8   
         dMd9 <- (gonRate*fracHumInf)*M9 + a8*Md8 -  (m9+a9 +pf)*Md9   

         dMi2 <- (pf)*Md2 - (m2+a2)*Mi2   
         dMi3 <- (pf)*Md3 + a2*Mi2 -  (m3+a3)*Mi3   
         dMi4 <- (pf)*Md4 + a3*Mi3 -  (m4+a4)*Mi4   
         dMi5 <- (pf)*Md5 + a4*Mi4 -  (m5+a5)*Mi5   
         dMi6 <- (pf)*Md6 + a5*Mi5 -  (m6+a6)*Mi6   
         dMi7 <- (pf)*Md7 + a6*Mi6 -  (m7+a7)*Mi7   
         dMi8 <- (pf)*Md8 + a7*Mi7 -  (m8+a8)*Mi8   
         dMi9 <- (pf)*Md9 + a8*Mi8 -  (m9+a9)*Mi9   

         
         res <- c(dM1 ,  dM2, dM3, dM4, dM5, dM6, dM7, dM8, dM9, 
                 dMd2, dMd3, dMd4, dMd5, dMd6, dMd7, dMd8, dMd9,
                 dMi2, dMi3, dMi4, dMi5, dMi6, dMi7, dMi8, dMi9)

         list(res)
       })
     }






