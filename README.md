## COMBAT [![CRAN](http://www.r-pkg.org/badges/version/COMBAT)](https://cran.r-project.org/package=COMBAT) [![Downloads](http://cranlogs.r-pkg.org/badges/COMBAT?color=brightgreen)](http://www.r-pkg.org/pkg/COMBAT) [![Total downloads]( https://cranlogs.r-pkg.org/badges/grand-total/COMBAT?color=brightgreen)](http://www.r-pkg.org/pkg/COMBAT)

### Description
Genome-wide association studies (GWAS) have been widely used for identifying common variants associated with complex diseases. Due to the small effect sizes of common variants, the power to detect individual risk variants is generally low. Complementary to SNP-level analysis, a variety of gene-based association tests have been proposed. However, the power of existing gene-based tests is often dependent on the underlying genetic models, and it is not known a priori which test is optimal.  Here we proposed COMBined Association Test (`COMBAT`) to incorporate strengths from multiple existing gene-based tests, including VEGAS, GATES and simpleM. Compared to individual tests, `COMBAT` shows higher overall performance and robustness across a wide range of genetic models.

### Installation
`COMBAT` is available from `CRAN` so the simplest way to install in `R` is by running `install.packages("COMBAT")`.

To install the latest update from here in `github`, run `devtools::install_github("mw201608/COMBAT")` in `R`.

### Usage
Run `?COMBAT` to see example usage after loading the package in R.

### Reference
[Minghui Wang, Jianfei Huang, Yiyuan Liu, Li Ma, James B Potash, Shizhong Han (2017) OMBAT: A Combined Association Test for Genes Using Summary Statistics. *Genetics* 207(3): 883-891.](https://doi.org/10.1534/genetics.117.300257)
