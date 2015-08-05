zinb <- function(x, tol=1e-6) {
   stopifnot(is.integer(x))
   MAXITER <- 50
   # Luckily, 'tabulate()' does not count zeros, neither 'NA'.
   # This allows us to compute everything conditionally on
   # x being non zero.
   tab <- tabulate(x)
   u <- 1:length(tab)
   u <- u[tab > 0]
   tab <- tab[tab > 0]
   
   n0 <- sum(x == 0)
   n <- sum(tab)
   N <- n + n0
   S <- sum(tab*u)

   # Define convenience functions.
   ll <- function(a, p, theta) {
      n0*log(theta*p^a + 1-theta) - n*lgamma(a) + n*log(theta) +
         a*n*log(p) + S*log(1-p) + sum(tab*lgamma(a+u))
   }

   f <- function(a) {
      -n*digamma(a) + n_*log(a/(m_+a)) + sum(tab*(digamma(a+u)))
   }

   df <- function(a) {
      -n*trigamma(a) + n_*m_/(a*(m_+a)) + sum(tab*(trigamma(a+u)))
   }

   start.a <- c(rep(NA, 11), 1)
   start.p <- c(rep(NA, 11), .5)

   for (deficit in seq(from=0, to=1, by=.1)) {
      
      n_ <- as.integer(N - n0 * deficit)
      m_ <- S / n_

      a <- 1
      # Find a lower and an upper bound for 'a'.
      if (f(a) < 0) {
         a <- a / 2
         while(f(a) < 0) { a <- a / 2 }
         a_lo <- a
         a_hi <- a * 2
      }
      else {
         a <- a * 2
         while(f(a) >= 0) { a <- a * 2 }
         a_lo <- a / 2
         a_hi <- a
      }

      # Failed...
      if (a_lo > 128) next

      new.a <- (a_lo + a_hi) / 2 + 2*tol

      # Solve equation numerically with Newton-Raphson method.
      for (j in 1:MAXITER) {
         a <- ifelse(new.a < a_lo || new.a > a_hi,
               (a_lo + a_hi) / 2, new.a)
         fa <- f(a)
         if   (fa < 0) { a_hi <- a }
         else          { a_lo <- a }
         new.a <- a - fa / df(a)
         if ((a_hi - a_lo) < tol) break
      }

      start.a[10*deficit+1] <- new.a
      start.p[10*deficit+1] <- new.a / (m_+new.a)

   }

   # For convenience, the functions to compute are defined below.
   f <- function(a,p) {
      n*a/(p*(1-p^a)) - S/(1-p)
   }
   g <- function(a,p) {
      n*log(p)/(1-p^a) - n*digamma(a) + sum(tab*(digamma(a+u)))
   }
   dfdp <- function(a,p) {
      -(n*a*(1-(a+1)*p^a)/(p*(1-p^a))^2 + S/(1-p)^2)
   }
   dgda <- function(a,p) {
      n*(log(p))^2*p^a/(1-p^a)^2 - n*trigamma(a) + sum(tab*(trigamma(a+u)))
   }
   dfda <- function(a,p) {
      n * (1-p^a + a*p^a*log(p)) / (p*(1-p^a)^2)
   }


   old.l <- -Inf
   for (idx in 1:12) {
      a <- start.a[idx]
      p <- start.p[idx]

      if (any(is.na(c(a,p)))) next

      compf <- f(a,p)
      compg <- g(a,p)
      oldgrad <- compf^2 + compg^2

      # Newton-Raphson.
      for (i in 1:MAXITER) {

         if (!is.na(oldgrad) && oldgrad < tol) break

         # Compute Hessian.
         A <- suppressWarnings(dfdp(a,p))
         B <- suppressWarnings(dfda(a,p))
         C <- B
         D <- suppressWarnings(dgda(a,p))
         H <- matrix(c(A,B,C,D), nrow=2)

         dpa <- NULL
         dpa <- try(- solve(H) %*% c(compf, compg), silent=TRUE)
         if (!is.null(dpa)) {
            dp <- dpa[1]
            da <- dpa[2]
         }
         else {
            dp <- 1
            da <- 1
         }

         compf <- suppressWarnings(f(a+da,p+dp))
         compg <- suppressWarnings(g(a+da,p+dp))
         newgrad <- compf^2 + compg^2

         for (j in 1:10) {
            if (!is.na(newgrad) && newgrad < oldgrad) break
            da <- da/2
            dp <- dp/2
            compf <- suppressWarnings(f(a+da,p+dp))
            compg <- suppressWarnings(g(a+da,p+dp))
            newgrad <- compf^2 + compg^2
         }

         p <- p + dp
         a <- a + da
         oldgrad <- compf^2 + compg^2

      }

      theta <- min(n / N / (1-p^a), 1.00)
      l <- suppressWarnings(ll(a,p,theta))

      if (!is.na(l) && l > old.l) {
         old.l <- l
         best.a <- a
         best.p <- p
         best.theta <- theta
      }

   }

   return(c(theta=best.theta, a=best.a, p=best.p))

}
