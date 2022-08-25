
source("./Scripts/initPop.R")
source("./Scripts/projPop.R")
source("./Scripts/redistributePopInd.R")

Zinit <- initPop()

Ka <- rep(1000, 100)
Kb <- c(rep(1000, 50), rep(500, 50))

for (t in 1:100){
  if (t == 1){
    Za_t <- projPop(Zinit = Zinit, nyears = 1)
    Zb_t <- projPop(Zinit = Zinit, nyears = 1)
  } else {
    Za_t <- projPop(Zinit = Za_tplus1, nyears = 1, K = Ka[t])
    Zb_t <- projPop(Zinit = Zb_tplus1, nyears = 1, K = Kb[t])
  }
  Z_new <- redistributePop(Za = Za_t, Zb = Zb_t, Ka = Ka[t], Kb = Kb[t])
  Za_tplus1 <- Z_new$Za_new
  Zb_tplus1 <- Z_new$Zb_new
  
} # end for t

plot(colSums(Za_tplus1)[50:150], type = "l", ylim = c(0, 2000))
lines(colSums(Zb_tplus1)[50:150])
