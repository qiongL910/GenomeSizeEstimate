\name{b}
\alias{b}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("b")}
\format{
  Real dateset with 10001 observations on the following 2 variables. This is data is the output from function \code{count} in \code{jellyfish} with k-mer length of 30.
  \describe{
    \item{\code{frequency}}{The frequency of kmer}
    \item{\code{counts}}{The counts of kmer at different frequencies}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Real data calculated from \emph{Glycine latifolia} sequencing reads.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(b)
## maybe str(b) ; plot(b) ...
}
\keyword{datasets}
