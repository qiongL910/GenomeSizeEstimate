\name{b}
\alias{b}
\docType{data}
\title{
Example data for package \code{GenomeSizeEstmate}
}
\description{
Object \code{b} contains real counted kmer frequency data from \emph{Glycine latifolia} whole genome sequencing project.
This data is the output from function \code{count} in \code{jellyfish} with k-mer length of 30.
}
\usage{data("b")}
\format{
  Real dateset with 10001 observations on the following 2 variables.
  \describe{
    \item{\code{frequency}}{The frequency of kmer}
    \item{\code{counts}}{The counts of kmer at different frequencies}
  }
}
\details{
The data set comes with 10001 observations with 2 variable. The First column is named with "frequency" whch records the frequency of kmer, and the second column is named with "counts" which records the number of kmers at certain frequency in the data set.
}
\source{
Real data calculated from \emph{Glycine latifolia} sequencing reads sequenced on two lanes of Hiseq 2500.
}
\references{
More details about \code{jellyfish} software and its function could be found at:

http://www.cbcb.umd.edu/software/jellyfish/
}
\examples{
data(b)
## maybe str(b) ; plot(b) ...
}
\keyword{datasets}
