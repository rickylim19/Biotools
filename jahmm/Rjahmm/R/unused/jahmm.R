BaumWelch <- function(m, y, pi_, alpha, C, blocks, index,
      verbose=TRUE, maxiter=1000, tol=1e-6) {

   nas <- is.na(colSums(yz, na.rm=FALSE))

   r <- nrow(y)
   n <- ncol(y)

   n_na.rm <- sum(!nas)
   y_na.rm <- y[,!nas]

   index <- rep(-1L, n)
   index_na.rm <- rep(-1L, n_na.rm)

   Q <- matrix(.05/(m-1), nrow=m, ncol=m)
   diag(Q) <- ifelse(m > 1, .95, 1.0)

   p <- matrix(NA, ncol=m, nrow=r+1)
   ybar <- mean(y[1,], na.rm=TRUE)
   # Fill in p row-wise.
   p[1,] <- ybar * C
   if (m == 1) mult <- 1
   if (m == 2) mult <- c(1,2)
   if (m == 3) mult <- c(.5,1,2)
   for (i in 2:r) {
      zbar <- mean(y[i,], na.rm=TRUE)
      p[i+1,] <- zbar * mult
   }
   # Assert that values are well-defined.
   if (any(p < .Machine$double.eps)) {
      stop('p undefined: check input to jahmm')
   }
   p <- scale(p, center=FALSE, scale=colSums(p))

   new.p <- p

   for (iter in 1:maxiter) {

      if (verbose) cat(paste("iteration:", iter, "\r"), file=stderr())
      oldparams <- p

      # E-step.

      # Compute (unnormalized) emission probabilities.
      C_call_1 <- .C(
         zinm_prob,
         # input #
         as.integer(m),
         as.integer(n),
         as.integer(r),
         as.integer(y),
         # params #
         as.double(pi_),
         as.double(alpha),
         as.double(p),
         # index #
         as.integer(index),
         # control #
         as.integer(0+4),  # No warning, linear space unless underflow.
         # output #
         double(m*n),      # Emission probabilities.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE
      )

      initialProb <- steady_state_probs(Q)

      # Compute smoothing phis by the forward-backward algorithm.
      C_call_2 <- .C(
         block_fwdb,
         # input #
         as.integer(m),
         as.integer(length(blocks)),
         as.integer(blocks),
         # params #
         as.double(Q),
         as.double(initialProb),
         # output #
         C_call_1[[11]],  # Forward alphas.
         double(m*n),     # Smoothing phis.
         double(m*m),     # Sum of transitions.
         double(1),       # Log-likelihood.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE
      )

      # M-step.
      Q <- matrix(C_call_2[[8]], nrow=m)
      Q <- Q / rowSums(Q)

      phi <- matrix(C_call_2[[7]], ncol=m, byrow=TRUE)[!nas,]

      # TODO: Update this part. # 
      for (j in 1:10) {

         old.p <- new.p
         old.q <- new.q

#         C_call_1 <- .C(
#            zinm_prob,
#            # input #
#            as.integer(m),
#            as.integer(n_na.rm),
#            as.integer(r),
#            as.integer(yz_na.rm),
#            # params #
#            as.double(p_),
#            as.double(alpha),
#            as.double(old.p),
#            # index #
#            as.integer(index_na.rm),
#            # control #
#            as.integer(2+4),     # No warning, ratio of mixture states.
#            # output #
#            double(m*n_na.rm),   # Probability ratios.
#            # extra '.C()' arguments #
#            NAOK = FALSE,
#            DUP = FALSE
#         )

         theta_i <- matrix(C_call_1[[11]], byrow=TRUE, ncol=m)
         zi1 <- scale(yz_na.rm %*% (phi * theta_i), center=FALSE,
                  scale=colSums(phi * theta_i))
         zi0 <- scale(yz_na.rm %*% (phi * (1-theta_i)), center=FALSE,
                  scale=colSums(phi * (1-theta_i)))

         # Add 'alpha' to first row.
         zi1[1,] <- alpha + zi1[1,]
         zi0[1,] <- alpha + zi0[1,]
         new.p <- scale(zi1, center=FALSE, scale=colSums(zi1))
         new.q <- scale(zi0, center=FALSE, scale=colSums(zi0))
         new.p[1,] <- new.p[1,] / (C1+1)
         new.q[1,] <- new.q[1,] / (C2+1)
         new.p <- rbind(C1*new.p[1,], new.p)
         new.q <- rbind(C2*new.q[1,], new.q)

         if (any(is.na(c(new.p,new.q,old.p,old.q)))) {
            if (verbose) {
               cat("Baum-Welch algorithm failed\n", file=stderr())
            }
            return (list(Q=Q, p=old.p, q=old.q, index=index,
                         iter=iter, failed=TRUE))
         }
         if (all(abs(c(new.p,new.q)-c(old.p,old.q))< 1e-3)) break

      }

      p <- new.p
      q <- new.q

      if (all(abs(oldparams - c(p,q)) < tol)) break

   }

   if (verbose) cat("\n", file=stderr())

   return (list(Q=Q, p=new.p, q=new.q, index=index,
                iter=iter, failed=FALSE))

}


jahmm <- function (data, PSO=FALSE, verbose=TRUE, threshold=0.08, ...) {

###############################################
#              OPTION PROCESSING              #
###############################################

   n <- nrow(data)

   # EM for mixture distribution of the baseline.
   baseline_params <- zinb(data[,2])

   pi_ <- baseline_params[1]
   alpha <- baseline_params[2]
   C <- baseline_params[3] / (1-baseline_params[3])

   # Coerce input data to matrix for C calls.
   y <- t(as.matrix(data[,-1]))


###############################################
#                 MAIN LOOP                   #
###############################################


   blocks <- tapply(X=rep(1,n), INDEX=as.character(data[,1]), FUN=sum)
   sorted <- match(unique(data[,1]), names(blocks))
   blocks <- blocks[sorted]

   loglik <- rep(NA, 3)
#   for (m in c(1,3)) {
   for (m in 3) {
      BW <- BaumWelch(m, y, pi_, alpha, C, blocks, verbose=verbose, ...)

      Q <- BW$Q
      p <- BW$p

      index <- BW$index

      if ((PSO || BW$failed) && m > 1) {
         if (verbose) cat("running PSO\n", file=stderr())
         # Run PSO 5 times in case the Baum-Welch algorithm failed.
         ncycles <- ifelse(BW$failed, 5, 1)
         for (cycle in 1:ncycles) {
            pso_call <- .C(
               pso,
               # input #
               as.integer(m),
               as.integer(n),
               as.integer(nrow(y)),
               as.integer(y),
               # fixed params #
               as.double(pi_),
               as.double(alpha),
               # params #
               as.double(Q),
               as.double(p),
               # index #
               as.integer(index),
               # output #
               double(1),
               # extra '.C()' arguments #
               NAOK = TRUE,
               DUP = FALSE
            )
            Q <- array(pso_call[[7]], dim=dim(Q))
            p <- array(pso_call[[8]], dim=dim(p))
         }
      }

      C_call_1 <- .C(
         zinm_prob,
         # input #
         as.integer(m),
         as.integer(n),
         as.integer(nrow(y)),
         as.integer(y),
         # params #
         as.double(pi_),
         as.double(alpha),
         as.double(p),
         # index #
         as.integer(index),
         # control #
         as.integer(1+4),     # No warning, log space.
         # output #
         double(m*n),         # Emission probabilities.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE
      )

      initialProb <- steady_state_probs(Q)

      if (m > 1) {
         vitC <- .C(
            block_viterbi,
            # input #
            as.integer(m),
            as.integer(length(blocks)),
            as.integer(blocks),
            as.double(log(Q)),
            as.double(log(initialProb)),
            as.double(C_call_1[[11]]),
            # control #
            as.integer(1),
            # output #
            integer(n),
            # extra '.C()' arguments #
            NAOK = TRUE,
            DUP = FALSE
         )

         vPath <- vitC[[8]]

      }

      C_call_2 <- .C(
         block_fwdb,
         # input #
         as.integer(m),
         as.integer(length(blocks)),
         as.integer(blocks),
         # params #
         as.double(Q),
         as.double(initialProb),
         # output #
         C_call_1[[11]],  # Forward alphas.
         double(m*n),     # Smoothing phis.
         double(m*m),     # Sum of transitions.
         double(1),       # Log-likelihood.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE
      )

      loglik[m] <- C_call_2[[9]]

   }

   phi <- matrix(C_call_2[[7]], ncol=m, byrow=TRUE)

   return(list(
      vPath = vPath,
      alpha = alpha,
      pi_ = pi_,
      Q = Q,
      p = p,
      score = score,
      loglik = loglik[3],
      phi = phi,
      iter = BW$iter
   ))

}
