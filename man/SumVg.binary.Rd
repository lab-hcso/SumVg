\name{SumVg.binary}
\alias{SumVg.binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sum of variance explained by all SNPs in GWAS (for a binary trait), with standard error estimate
}
\description{
This function estimates the sum of variance explained (i.e. heritbaility) by all SNPs in GWAS of a binary trait,
with standard error estimate
}
\usage{
SumVg.binary(zall, method = "paraboot", d = 1, repl = 200, out = "unconditional", SE, K, MAF )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zall}{
A vector of z-statistics from a GWAS (or meta-analyses of GWAS) for a binary outcome; z-statistics are assumed to be LD-pruned
}
  \item{method}{
Describes the method for standard error estimation.
Options include:
\cr "jack"(jackknife)
\cr "paraboot"(parametric bootstrap using the corrected z value as mean for simulation)
\cr "fdrboot1"(parametric bootstrap using the observed z value as mean with consideration of local fdr)
\cr "fdrboot2"(parametric bootstrap using the corrected z value as mean with consideration of local fdr)
}
  \item{d}{
The number of observations to be left out in each jackknife replicate
}
  \item{repl}{
Number of replications for resampling
}
  \item{out}{
Describes the type of sum of Vg. Options include "conditional" (conditional on H1) and "unconditional".
}
  \item{SE}{
Standard error of the logistic regression coefficients (a vector)
}
  \item{K}{
Prevalence of disease
}
  \item{MAF}{
Vector of minor allele frequencies for the corresponding SNPs (does not need to match to the effect allele)
}
}
\details{
This function estimates the sum of variance explained (i.e. heritbaility) by all SNPs in GWAS of a binary trait.
The standard error is estimated using various resampling approaches.
}
\value{

\item{Est.SumVg}{The estimated sum of variance explained by all SNPs in GWAS}
\item{SE.SumVg}{The standard error of the estimated sum of variance}

}
\references{
So, H.C., Li, M. and Sham, P.C. (2011) Uncovering the total heritability explained by all true susceptibility variants in a genome-wide association study
Genet Epidemiol, 35, 447-456.
\cr So, H.C., Xue, X. and Sham, P.C. (2023) SumVg: Total heritability explained by all variants in genome-wide association studies based on summary statistics with standard error estimates.
}
\author{
Hon-Cheong So
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##simulate z-statistics under the complete null
## SE and MAF are made up for reference only
zall = rnorm(n=1000)
SumVg.binary(zall=zall, method="paraboot", d=1, repl=5, out="unconditional", SE= rep(0.1,1000), K=0.01, MAF=rep(0.2, 1000) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
