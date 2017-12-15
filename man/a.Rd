\name{a}
\alias{a}
\docType{data}
\title{
Example data for package \code{GenomeSizeEstmate}
}
\description{
Object a contains real counted kmer frequency data from Glycine latifolia whole genome sequencing project. This data is the output from function \code{count} in \code{jellyfish} with k-mer length of 19.
}
\usage{data("a")}
\format{
  Real dateset with 10001 observations on the following 2 variables. 
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
data(a)
## maybe str(a) ; plot(a) ...

}
\keyword{datasets}
