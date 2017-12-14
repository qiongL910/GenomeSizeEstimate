\name{GenomeEstimate}
\alias{GenomeEstimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GenomeEstimate(file, kmer_len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{kmer_len}{
%%     ~~Describe \code{kmer_len} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file, kmer_len) 
{
    for (q in 1:length(file$counts)) {
        if (file[q, 2] < file[q + 1, 2]) {
            start_point <- q
            break
        }
    }
    for (h in start_point:length(file$counts)) {
        if (file[h, 2] > file[h + 1, 2]) {
            mean_coverage <- h
            cat(paste("The mean coverage of ", kmer_len, "-mer is ", 
                mean_coverage, "\n", sep = ""))
            break
        }
    }
    trusted_kmer <- sum(as.numeric(file$frequency[start_point:length(file$frequency)]) * 
        as.numeric(file$counts[start_point:length(file$counts)]))
    genome_size <- round(trusted_kmer/(mean_coverage * 10^9), 
        3)
    cat(paste("Estimated genome size based on ", kmer_len, "-mer analysis: ", 
        genome_size, " GB\n", sep = ""))
    i <- mean_coverage + start_point
    subset_file <- file$counts[i:length(file$counts)]
    decrease <- numeric(length(subset_file))
    for (h in 1:length(subset_file)) {
        decrease[h] <- (subset_file[h] - subset_file[h + 1])/subset_file[h]
    }
    for (k in 1:length(decrease)) {
        if (decrease[k] < 0.03) {
            stop_point <- k
            break
        }
    }
    single_copy_stop <- i + stop_point
    cat(paste("Single copy region end point: ", single_copy_stop, 
        "\n", sep = ""))
    single_copy_region <- round(sum(as.numeric(file$frequency[start_point:single_copy_stop]) * 
        as.numeric(file$counts[start_point:single_copy_stop]))/(mean_coverage * 
        1000000000), 3)
    single_copy_percentage <- round(single_copy_region * 100/genome_size, 
        2)
    cat("Estimated single copy region:", single_copy_region, 
        " Gb\n", sep = "")
    cat("Percentage of single copy region in genome: ", single_copy_percentage, 
        "\%\n", sep = "")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
