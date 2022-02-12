\name{SumVg}
\alias{SumVg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sum of variance explained by all SNPs in GWAS (for a continuous trait), with standard error estimate
}
\description{
This function estimates the sum of variance explained (i.e. heritbaility) by all SNPs in GWAS of a continuous trait, 
with standard error estimate
}
\usage{
SumVg(zall, totalN, method = "paraboot", d = 1, repl = 50, 
out = "unconditional")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zall}{
A vector of z-statistics from a GWAS (or meta-analyses of GWAS)
}
  \item{totalN}{
Total sample size
}
  \item{method}{
Describes the method for standard error estimation. Options include: 
\cr "jack"(jackknife); 
\cr "paraboot"(parametric bootstrap using the corrected z value as mean for simulation);  
\cr "fdrboot1"(parametric bootstrap using the observed z value as mean with consideration of local fdr) ; 
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
}
\details{
This function estimates the sum of variance explained (i.e. heritbaility) by all SNPs in GWAS of a continuous trait. 
The standard error is estimated using various resampling approaches. 
}
\value{
\item{Est.SumVg}{The estimated sum of variance explained by all SNPs in GWAS}
\item{SE.SumVg}{The standard error of the estimated sum of variance}
}
\references{
So, H.C., Li, M. and Sham, P.C. (2011) Uncovering the total heritability explained by all true susceptibility variants in a genome-wide association study
Genet Epidemiol, 35, 447-456.
\cr So, H.C. and Sham, P.C. (2014) Standard error estimation of the total heritability explained by all variants in genome-wide association studies based 
on summary statistics. Submitted. 

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
zall = rnorm(10000)
SumVg(zall=zall, totalN=5000, method="jack", d=2000, repl=5,out="unconditional" )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
