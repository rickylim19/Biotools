EM.mnb <- function(x, a, theta, p, q, tol=1e-6, maxiter=100) {
   # Does maximum likelihood estimation of negative binomial
   # parameters for sample 'x'.

   # Ignore NAs.
   n <- sum(!is.na(x))
   m <- mean(x, na.rm=TRUE)
   tab <- tabulate(x+1L)
   u <- 0:(length(tab)-1L)
   u <- u[tab > 0]
   tab <- tab[tab > 0]

   if (missing(a)) {
      a <- 1
   }
   if (missing(theta)) {
      # Empty bins are usually over-represented. Get an initial
      # value of theta as the probability that a bin is not empty.
      theta <- 1 - tab[1] / n
   }
   usplit <- 1
   if (missing(p)) {
      R <- sum(u[u > usplit]*tab[u > usplit]) / sum(tab[u > usplit])
      p <- a/(a+R)
   }
   if (missing(q)) {
      L <- sum(u[u <= usplit]*tab[u <= usplit]) / sum(tab[u <= usplit])
      q <- a/(a+L)
   }
      
   f <- function(a) {
      return(-n*digamma(a) + 
         sum(tab*(digamma(a+u) - w*log(1+x1/a)-(1-w)*log(1+x2/a))))
   }   
   df <- function(a) {
      return(-n*trigamma(a) +
         sum(tab*(trigamma(a+u) + w*x1/(a*(x1+a))+(1-w)*x2/(a*(x2+a)))))
   }   

   old_params <- c(a, theta, p, q)

   # Start the EM cycles.
   for (iter in 1:maxiter) {
      # Compute the weights in log space for numeric stability.
      # This is equivalent to w = L1 / (L1 + L2), where
      #    L1 = theta * p^a * (1-p)^u and
      #    L2 = (1-theta) * q^a * (1-q)^u,
      # avoiding the undefined case 0/0.
      logterms <- log(1-theta)-log(theta) + a*(log(q)-log(p)) +
         u*(log(1-q)-log(1-p))
      w <- 1 / (1 + exp(logterms))

      # If 'q' or 'p' approaches 1 the weights will be
      # undetermined for u = 0.
      if (is.na(w[1])) w[1] = 0

      x1 <- sum(w*u*tab) / sum(w*tab)
      x2 <- sum((1-w)*u*tab) / sum((1-w)*tab)

      if (f(a) < 0) {
         a <- a / 2
         while(f(a) < 0) { a <- a / 2 }
         a. <- a
         .a <- a * 2
      }
      else {
         a <- a * 2
         while(f(a) >= 0) { a <- a * 2 }
         a. <- a / 2
         .a <- a
      }

      # Failed...
      if (a. > 128) { return(c(loglik=-Inf, a=Inf, theta=0, p=0, q=0)) }

      new.a <- a + 2*tol

      # Solve equation in alpha (a) numerically with
      # Newton-Raphson method.
      for (whileiter in 1:maxiter) {
         a <- ifelse(new.a < a. || new.a > .a, (a. + .a) / 2, new.a)
         fa <- f(a)
         if (fa < 0) { .a <- a }
         else { a. <- a }
         new.a <- a - fa / df(a)
         if (abs(new.a - a) < tol) break
      }

      a <- new.a
      theta <- sum(w*tab) / n
      p <- a / (x1+a)
      q <- a / (x2+a)

      if (all(abs(c(a, theta, p, q) - old_params) < tol)) break
      old_params <- c(a, theta, p, q)

   }

   # Swap parameters so that 'p' is always smaller than 'q'.
   if (p > q) {
      tmp <- p; p <- q; q <- tmp;
      theta <- 1-theta
   }

   # The general log-likelihood equation is given by
   # 'sum(lgamma(x+a) - lgamma(a) - lgamma(x+1) + log(L1+L2))', where
   #    L1 = theta * p^a * (1-p)^x, and
   #    L2 = (1-theta) * q^a * (1-q)^x
   # We rewrite the last term of the equation as
   #    log(L1+L2) = log(1+L2/L1) + log(L1)
   # Because p < q, the term 'L2/L1' is bounded, which prevents
   # overflow. Because 'L1' is a product, the term 'log(L1)' will
   # not underflow.
   logL1L2 <- log(1 + (1-theta)/theta * (q/p)^a * ((1-q)/(1-p))^u) + 
      log(theta) + a*log(p) + u*log(1-p)
   loglik <- sum(tab*lgamma(a+u))/n - lgamma(a) -
      sum(tab*lgamma(u+1))/n + sum(tab*(logL1L2))/n

   return(c(loglik=loglik, a=a, theta=theta, p=p, q=q))

}
