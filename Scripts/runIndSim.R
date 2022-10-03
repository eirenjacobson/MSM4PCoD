
source("./Scripts/initPop.R")
source("./Scripts/projPop.R")
source("./Scripts/redistributePopInd.R")

set.seed(20221003)
Zinit <- initPop()

Ka <- rep(1000, 100)
Kb <- c(rep(1000, 50), rep(500, 50))

plot(x = NA, y = NA, xlim = c(0, 100), ylim = c(-10, 10))

Na_t <- rep(NA, 100)
Nb_t <- rep(NA, 100)

for (t in 1:100){
  if (t == 1){
    Za_t <- projPop(Zinit = Zinit, nyears = 1)
    Zb_t <- projPop(Zinit = Zinit, nyears = 1)
  } else {
    Za_t <- projPop(Zinit = Za_tplus1, nyears = 1, K = Ka[t])
    Zb_t <- projPop(Zinit = Zb_tplus1, nyears = 1, K = Kb[t])
  }
  Z_new <- redistributePop(Za = Za_t, Zb = Zb_t, Ka = Ka[t], Kb = Kb[t], c = 1)
  Za_tplus1 <- Z_new$Za_new
  Zb_tplus1 <- Z_new$Zb_new
  
  Na_t[t] <- sum(Za_tplus1[,ncol(Za_tplus1)])
  Nb_t[t] <- sum(Zb_tplus1[,ncol(Zb_tplus1)])
  
#  points(x= t, y = sum(Za_t[,ncol(Za_t)]), col = "blue", pch = 1)
  points(x = t, y = sum(Za_tplus1[,ncol(Za_tplus1)]) - sum(Za_t[,ncol(Za_t)]), col = "blue", pch = 16)
  
 # points(x = t, sum(Zb_t[,ncol(Zb_t)]), col = "green", pch = 2)
  points(x = t, sum(Zb_tplus1[,ncol(Zb_tplus1)]) - sum(Zb_t[,ncol(Zb_t)]), col = "green", pch = 17)
  
} # end for t

plot(Na_t, type = "l", ylim = c(0, 2000))
lines(Nb_t)
