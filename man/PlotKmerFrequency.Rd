\name{PlotKmerFrequency}
\alias{PlotKmerFrequency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the kmer frequency distribution
}
\description{
This function will utilize part of the output from both function \code{GenomeEstimate} and \code{FindTrustKmer} and plot the distribution of counted kmers at various frequencies.

Additionally, it will indicate the mean coverage of kmers, as well as the theoretical poisson distribution of kmer frequency with the mean of kmer coverage.
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
This function takes the output from \code{jellyfish count} function and the column of the data have to be renamed into "frequency" and "counts" respectively.

The function will utilize the output (The starting point of trusted kmer, mean coverage of kmer, and ending point of single copy region) from function \code{FindTrustKmer} and \code{GenomeEstimate}. The theoretical poission distribution of counted kmer frequency is based on the calcualted single copy region and the mean coverage of kmer.
}
\value{
This function will return a plot of counted kmer frequency distribution. Function \code{PlotKmerFrequency} also give a theoratical poisson distribution with a mean of mean coverage of kmers.
}
\references{
More details about calculation can be referred at:

http://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html

https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
}
\author{
Qiong Liu
}
\note{
N/A
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Function \code{FindTrustKmer} and \code{GenomeEstimate}
}
\examples{
# load an example data called b
data(b)

# Get startingpoint of trusted kmer
FindTrustKmer(b,30)

# Get the mean covergae of kmer, ending point of single copy region
GenomeEstimate(b,30)

# Plot the figure
PlotKmerFrequency(file=b,kmer_len=30,start_point = 13,peak=36,end_point = 80)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
