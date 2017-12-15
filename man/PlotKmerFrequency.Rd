\name{PlotKmerFrequency}
\alias{PlotKmerFrequency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the kmer frequency distribution
}
\description{

}
\usage{
PlotKmerFrequency(file, kmer_len, start_point, peak, end_point)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
Counted kmer frequency file from \code{jellyfish count} function. The first and second columns have to be names as "frequency" and "counts" respectively.
}
  \item{kmer_len}{
  The length of kmers.
}
  \item{start_point}{
  The starting point of trusted kmers. Get this value trhough function \code{FindTrustKmer(file,kmer_len)}.
}
  \item{peak}{
  The mean coverage of kmers. Get this value through function \code{GenomeEstimate(file,kmer_len)}.
}
  \item{end_point}{
  The end point of single copy region. Get this value through function \code{GenomeEstimate(file,kmer_len)}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function will return a plot of counted kmer frequency distribution. Function \code{PlotKmerFrequency} also give a theoratical poisson distribution with a mean of mean coverage of kmers.
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
# load an example data called b
data(b)

# Get startingpoint of trusted kmer
FindTrustKmer(b,30)

# Get the mean covergae of kmer, end point of single copy region
GenomeEstimate(b,30)

# Plot the figure
PlotKmerFrequency(file=b,kmer_len=30,start_point = 13,peak=36,end_point = 80)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
