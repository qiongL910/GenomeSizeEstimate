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
This function takes the output from \code{jellyfish count} function and the column of the data have to be renamed into "frequency" and "counts" respectively.

The function will first detect the trust kmer starting point, and subsequently identify the mean coverage of kmer after that. It then estimate the gneome size based on the counted trusted kmer and its mean coverage.

The end point of single copy region is detected at the proximal end of bell shape distribution of kmer frequency.
}
\value{
This function returns the mean coverage of kmers, estimated genome size based on kmer analysis, the end point of single copy region, estimated single copy region and the percentage of single copy region in genome. Among which, the mean coverage of kmers, single copy region end point are needed in \code{PlotKmerFrequency} function.
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
Function \code{FindTrustKmer} and \code{PlotKmerFrequency}
}
\examples{
# load an example data called a with kmer length 19bp
data(a)
GenomeEstimate(a,19)

# load an example data called b with kmer length 30bp
data(b)
GenomeEstimate(b,30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
