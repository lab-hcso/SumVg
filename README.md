## SumVg: Estimating the standard error of the total heritability explained by all variants in GWAS
Hon-Cheong SO

### Introduction
The R package SumVg provides estimates of the sum of heritability explained by all true susceptibility variants in GWAS. We have also recently derived methods to estimate the standard error (SE) based on re-sampling approaches. These methods have been implemented in the R package SumVg. 

### Installation
To install the SumVg R package from GitHub using `devtools`, run

```R
devtools::install_github("lab-hcso/SumVg")
```

To install the SumVg R package by file,

downloading the tar file *SumVg_1.0.tar.gz* from [GitHub](https://github.com/lab-hcso/SumVg), run

```R
install.packages("C:/SumVg_1.0.tar.gz", repos=NULL)
```

### Usage

```R
library(SumVg)

## Simulate z-statistics under the complete null for testing
zall = rnorm(n=10000, mean=0, sd = 1)

## Examples using delete-d-jackknife 
## Assume the outcome is continuous
SumVg(zall=zall, totalN=10000, method="jack", d=2000, repl=5,out="unconditional") 

## Assume the outcome is binary
## SE and MAF are made up for reference only
SumVg.binary(zall=zall, method="paraboot", d=1, repl=5, out="unconditional", SE= rep(0.1,10000), K=0.01, MAF=rep(0.2, 10000) )
```

### Citation

To cite this project in publications use:

- [H.-C. So, et al., “Uncovering the total heritability explained by all true susceptibility variants in a genome-wide association study,” Genetic Epidemiology, vol. 35, no. 6, 2011, pp. 447-456; DOI https://doi.org/10.1002/gepi.20593.](https://doi.org/10.1002/gepi.20593)

- [So, H. C., Xue, X. & Sham, P. C. (2023). SumVg: Total heritability explained by all variants in genome-wide association studies based on summary statistics with standard error estimates](https://arxiv.org/abs/2306.14200)
