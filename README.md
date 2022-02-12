### Estimating the standard error of the total heritability explained by all variants in GWAS

The R package SumVg provides estimates of the sum of heritability explained by all true susceptibility variants in GWAS. We have also recently derived methods to estimate the standard error (SE) based on re-sampling approaches. These methods have been implemented in the R package SumVg. 

Example:

\# install SumVg first

\# installation of sfsmisc and locfdr 

install.packages("sfsmisc")

install.packages("locfdr")

library(sfsmisc)

library(locfdr)

library(SumVg)

\#\# simulate z-statistics under the complete null for testing

zall = rnorm(n=10000, mean=0, sd = 1)

\#\#examples using delete-d-jackknife 

\#\#assume outcome is continuous

SumVg(zall=zall, totalN=10000, method="jack", d=2000, repl=5,out="unconditional") 

\#\# assume outcome is binary

SumVg.binary(zall=zall, method="jack", d=2000, repl=5, out="unconditional", caseNo=10000, ctrlNo=10000, K=0.01) 

#### Citation

To cite this project in publications use:

- [H.-C. So and P.C. Sham, “Improving polygenic risk prediction from summary statistics by an empirical Bayes approach,” Scientific Reports, vol. 7, no. 1, 2017, pp. 41262; DOI 10.1038/srep41262.](https://doi.org/10.1038/srep41262)

