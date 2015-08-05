\name{wig2jahmm}
\alias{wig2jahmm}
\title{wig2jahmm}
\usage{
wig2jahmm(directory)
}
\arguments{
  \item{directory}{directory where the wig files are. ChIP-seq Input profiles must contain the strings \code{input} or \code{control} (case-insensitive) in their file name in order to be recognized as Input.}
}
\value{
  A fitted HMM using the \code{jahmm} discretizer. (See \code{jahmm})
}
\description{
  Wrapper for the jahmm discretizer. Run jahmm directly from wig files.
}
\examples{
fit <- wig2jahmm("./exampleDirectory")
}
\author{
  Pol Cuscó. date: May 23, 2014.
}
\seealso{
  \code{\link{jahmm}}
}
\details{
  The ChIP-seq profiles to discretize must be in the same directory. They consist in one or more Input profiles that will be used as baseline, and one or more Target profiles (these must have been obtained from the same DNA-binding protein and, ideally, they should be experimental replicates). Only \code{fixedStep} formatted wig files are supported.
}
