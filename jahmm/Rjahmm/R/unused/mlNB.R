mlNB <- function(x, tol=1e-6) {
# Maximum likelihood estimation of negative binomial parameters.

   n <- length(x)
   m <- mean(x, na.rm=TRUE)
   tab <- tabulate(x+1L)
   u <- 0:(length(tab)-1L)
   a <- mean(x)
   new.a <- a + 2*tol

   while (abs(new.a - a) > tol) {
      a <- new.a
      num <- sum(tab*digamma(a+u)) - n*digamma(a) + n*log(a/(m+a))
      denom <- sum(tab*trigamma(a+u)) - n*trigamma(a) + n*m/(a*(m+a))
      new.a <- a - num/denom
   }

   return(c(a=new.a, p=m/(m+new.a)))

}
