# ---------------------------------------------
# Collect data for timing plots (ALBI v score)
# --------------------------------------------
# setwd("TimingComparison")
nseq <- c(seq(200, 2000, length=9), c(2250, 2500, 2750, 3000))       # size of K
ours.times01 <- matrix(0, 100, length(nseq))
fiesta.times01 <- matrix(0, 100, length(nseq))
albi.times01 <- matrix(0, 100, length(nseq))
rlrt.times01 <- matrix(0, 100, length(nseq))

for(uu in 1:100){
  t0 <- readRDS(paste("Results_R1/h2_01_Replicate", uu, ".RDS", sep=""))
  ours.times01[uu,] <- t0$time$proposed
  fiesta.times01[uu,] <- t0$time$fiesta
  albi.times01[uu,] <- t0$time$albi
  rlrt.times01[uu,] <- t0$time$rlrt
  cat(uu, "\n")
}


ours.times1 <- matrix(0, 100, length(nseq))
fiesta.times1 <- matrix(0, 100, length(nseq))
albi.times1 <- matrix(0, 100, length(nseq))
rlrt.times1 <- matrix(0, 100, length(nseq))

for(uu in 1:100){
  t0 <- readRDS(paste("Results_R1/h2_1_Replicate", uu, ".RDS", sep=""))
  ours.times1[uu,] <- t0$time$proposed
  fiesta.times1[uu,] <- t0$time$fiesta
  albi.times1[uu,] <- t0$time$albi
  rlrt.times1[uu,] <- t0$time$rlrt
  cat(uu, "\n")
}




ours.times5 <- matrix(0, 100, length(nseq))
fiesta.times5 <- matrix(0, 100, length(nseq))
albi.times5 <- matrix(0, 100, length(nseq))
rlrt.times5 <- matrix(0, 100, length(nseq))

for(uu in 1:100){
  t0 <- readRDS(paste("Results_R1/h2_5_Replicate", uu, ".RDS", sep=""))
  ours.times5[uu,] <- t0$time$proposed
  fiesta.times5[uu,] <- t0$time$fiesta
  albi.times5[uu,] <- t0$time$albi
  rlrt.times5[uu,] <- t0$time$rlrt
  cat(uu, "\n")
}


result <- list(
  "ours.times01" = ours.times01,
  "ours.times1" = ours.times1,
  "ours.times5" = ours.times5,
  "fiesta.times01" = fiesta.times01,
  "fiesta.times1" = fiesta.times1,
  "fiesta.times5" = fiesta.times5,
  "albi.times01" = albi.times01,
  "albi.times1" = albi.times1,
  "albi.times5" = albi.times5,
  "rlrt.times01" = rlrt.times01,
  "rlrt.times1" = rlrt.times1,
  "rlrt.times5" = rlrt.times5
)

saveRDS(result, "TimingResults.RDS")

