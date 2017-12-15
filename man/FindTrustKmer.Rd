\name{FindTrustKmer}
\alias{FindTrustKmer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate trusted and untrusted kmer from counted kmer frequency file
}
\description{
The counted kmer frequency file calculated from whole genome sequencing data using \code{jellyfish count} function contains large amount of kmers with low frequency. These kmers are considered untrusted kmers from sequencing errors. 
\code{FrindTrustKmer} detect the starting point of trusted kmer and calculate the total number of untrusted and trusted kmers in the dataset. 
}
\usage{
FindTrustKmer(file, kmer_len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
Counted kmer frequency file from \code{jellyfish count} function. The first and second columns have to be named as "frequency" and "counts" respectively.
}
  \item{kmer_len}{
The length of kmers.
}
}
\details{
This function takes the output from \code{jellyfish count} function and the column of the data have to be renamed into "frequency" and "counts" respectively.

The function will first detect the trust kmers starting point, and subsequently calculate the total number of trusted and untrusted kmers at various frequency, thus generating the percentage of trusted and untrusted kmers in the data set.
}
\value{
This function will return two values.
The first value is the starting point of trusted kmers. This value will be useful in the \code{PlotKmerFrequency} function.
The second value is a vector which contains the total number of all kmers, trusted kmers, and untrusted kmers, as well as the percentage of trusted and untrusted kmers.
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
Function \code{GenomeEstimate} and \code{PlotKmerFrequency}.
}
\examples{
# load the example data. This gives you an example data called a with kmer length 19bp
data(a)

FindTrustKmer(a,19)

# load the example data. This gives you an example data called b with kmer length 30bp
data(b)

FindTrustKmer(b,30)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
