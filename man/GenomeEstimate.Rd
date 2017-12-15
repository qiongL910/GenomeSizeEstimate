\name{GenomeEstimate}
\alias{GenomeEstimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate genome size and single copy region in genome
}
\description{
Estimated genome size is an important statistic to collect in genomic studies, especially for the evaluation purpose for de novo whole genome assembly. Besides flow cytometric analysis, genome size can also be simply estimated with comparable accuracy based on k-mer analysis using NGS short-read sequencing data.
}
\usage{
GenomeEstimate(file, kmer_len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
Counted kmer frequency file from \code{jellyfish count} function. The first and second columns have to be names as "frequency" and "counts" respectively.
}
  \item{kmer_len}{
The length of kmers
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns the mean coverage of kmers, estimated genome size based on kmer analysis, the end point of single copy region, estimated single copy region and the percentage of single copy region in genome. Among which, the mean coverage of kmers, single copy region end point are needed in \code{PlotKmerFrequency} function.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Qiong Liu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load an example data called a
data(glat_2lane_kmer19)
GenomeEstimate(a,19)

# load an example data called b
data(glat_2lane_kmer30)
GenomeEstimate(b,30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
