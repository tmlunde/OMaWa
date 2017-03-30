morA <- function (T){
c = 0.000203
Tm = 42.3
T0 = 11.7
T[T < T0] <- NA
T[T > Tm] <- NA
res <- c*T*(T - T0)*(Tm - T)^(1/2)
res[res<0] <- 0
res
}



morBc <- function (T){
q = -0.54
r = 25.2
s = -206
res <- q*T^2 + r*T + s
res[res<0] <- 0
res
}


morP <- function (T){
q = -0.000828
r = 0.0367
s = 0.522
q*T^2 + r*T + s

}


morPdr <- function (T){

c = 0.000111
Tm = 34.4
T0 = 14.7
res <- c*T*(T - T0)*(Tm - T)^(1/2)
res[res<0] <- 0
res
}


morPea <- function (T){
q = -0.00924
r = 0.453
s = -4.77 
res <- q*T^2 + r*T + s
res[res<0] <- 0
res
}


morMdr <- function (T){
c = 0.000111
Tm = 34
T0 = 14.7
res <- c*T*(T - T0)*(Tm - T)^(1/2)
res[res<0] <- 0
res

}

morEfd <- function (T){
q = -0.153
r = 8.61
s = -97.7 
res <- q*T^2 + r*T + s
res[res<0] <- 0
res

}

morR0 <- function (T, N=1, r=1){
res <- (morA(T)^2 * morBc(T) * exp(log(morP(T))/morPdr(T)) * morEfd(T) * morPea(T)*morMdr(T)/(N*r*(-log(morP(T)))^3))^.5
res[is.na(res)] <- 0
res
}


