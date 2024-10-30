#library(MASS)

# find boundaries in a spatial dataset given sequences of row and column
findboundary <- function(rw, cl) {
  ur <- sort(unique(rw))
  uc <- sort(unique(cl))
  index_bd <- NULL
  # for each row, find which index has column values in boundary
  for (i in 1:length(ur)) {
    colsinrowi <- cl[which(rw == ur[i])]     # column values in row i
    index_bd <- c(index_bd, which(cl == min(colsinrowi) & rw == ur[i]), which(cl == max(colsinrowi) & rw == ur[i]))
  }
  # for each col, find which index has row values in boundary
  for (i in 1:length(uc)) {
    rowsincoli <- rw[which(cl == uc[i])]     # row values in row i
    index_bd <- c(index_bd, which(rw == min(rowsincoli) & cl == uc[i]), which(rw == max(rowsincoli) & cl == uc[i]))
  }
  return(unique(index_bd))               # remove duplicate
}

# compute spatial coordinate distance matrix K
spatialCor <- function(coords, l = 2, scale = FALSE) {
  if (scale == TRUE) {
    maxrange <- max(apply(coords, 2, function(x) diff(range(x))))
    coords <- apply(coords, 2, function(x) (x - min(x)) / maxrange)
  }
  n <- nrow(coords)
  S <- matrix(0, nrow = n, ncol = n)
  for (kk in 1:n) {
    for (jj in 1:n) {
      dist <- sqrt((coords[kk, 1] - coords[jj, 1]) ^ 2 + (coords[kk, 2] - coords[jj, 2]) ^ 2)
      S[jj,kk] <- exp(-dist / l)
    }
  }
  return(S)
}