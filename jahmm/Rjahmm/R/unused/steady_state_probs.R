steady_state_probs <- function(Q) {
   if (any(is.na(Q)) || any(is.infinite(Q))) {
      return(eig$vectors[,1] / sum(eig$vectors[,1]))
   }
   eig <- eigen(t(Q))
   if (is.complex(eig$values[1]) || eig$values[1] < 0.99) {
      return(rep(1/nrow(Q),nrow(Q)))
   }
   return(eig$vectors[,1] / sum(eig$vectors[,1]))
}
