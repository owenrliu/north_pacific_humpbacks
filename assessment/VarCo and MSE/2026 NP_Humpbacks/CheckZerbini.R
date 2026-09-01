SA <- 0.984
SC <- 0.669
Apar <- 5.9
Amat <- Apar-1.0
birth <- 0.41
q <- 0.5
r <- 0.073
r <- 0.09
T1 <- exp(Apar*r)
T2 <- SA*exp((Apar-1)*r)
T3 <- birth*q*SC*SA^(Apar-1)
print(T1)
print(T2+T3)


birth <- (exp(r*(Amat+1.0))-SA*exp(r*Amat))/(SC*SA^(Amat))
print(birth)

T1 <- exp(Apar*r)
T2 <- SA*exp((Apar-1)*r)
T3 <- birth*SC*SA^(Apar-1)
print(T1)
print(T2+T3)
###
print("part2")
SA <- 0.96
SC <- 0.66
Apar <- 9
Amat <- Apar-1.0
birth <- 0.41
q <- 0.5
r <- 0.0905
T1 <- exp(Apar*r)
T2 <- SA*exp((Apar-1)*r)
T3 <- birth*q*SC*SA^(Apar-1)
print(T1)
print(T2+T3)


birth <- (exp(r*(Amat+1.0))-SA*exp(r*Amat))/(SC*SA^(Amat))
print(birth)

T1 <- exp(Apar*r)
T2 <- SA*exp((Apar-1)*r)
T3 <- birth*SC*SA^(Apar-1)
print(T1)
print(T2+T3)
