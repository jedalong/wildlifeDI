
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wildlifeDI

[![](https://cranlogs.r-pkg.org/badges/wildlifeDI)](https://cran.r-project.org/package=wildlifeDI)

The wildlifeDI package facilitates the calculation of indices of dynamic
interaction for wildlife telemetry data. There are also functions for
more advanced contact analysis. For more information on the methods used
within see the documentation but also see:

Long, J.A., Nelson, T.A., Webb, S.L., Gee, K.L. (2014) A critical
examination of indices of dynamic interaction for wildlife telemetry
studies. Journal of Animal Ecology. 83(5):1216-1233.

## Installation

You can install the latest (under development version) of wildlifeDI
from github with:

``` r
devtools::install_github("jedalong/wildlifeDI")
```

The version that can currently be downloaded from CRAN is wildlifeDI
v0.4.1:

``` r
install.packages('wildlifeDI')
```

## Example

A small dataset is provided to demonstrate the functionality of the
package.

``` r
library(wildlifeDI)
data(deer)
deer
#> 
#> *********** List of class ltraj ***********
#> 
#> Type of the traject: Type II (time recorded)
#> * Time zone unspecified: dates printed in user time zone *
#> Irregular traject. Variable time lag between two locs
#> 
#> Characteristics of the bursts:
#>   id burst nb.reloc NAs          date.begin            date.end
#> 1 37    37      551   0 2005-03-07 19:03:00 2005-03-13 18:47:00
#> 2 38    38      567   0 2005-03-07 19:02:00 2005-03-13 18:47:00
#> 
#> 
#>  infolocs provided. The following variables are available:
#> [1] "pkey"
```

Different tests of dynamic interaction can be straightforwardly
computed. Most require a distance threshold (dc) to define when two
animals are proximal. Similarly, most methods require a time threshold
(tc) to define when two fixes are simultaneous.

Here we use dc = 50 and tc = 8 mins (8\*60).

A few examples are given below:

``` r
dc <- 50
tc <- 8*60
#Proximity Analysis
Prox(deer[1],deer[2])
#> [1] 0.4293948
#Doncaster's test
Don(deer[1],deer[2],tc=tc,dc=dc,plot=F)
#> $conTable
#>            below.crit above.crit Totals
#> Paired            226        320    546
#> Non-Paired       3925     293637 297570
#> Totals           4159     293957 298116
#> 
#> $p.value
#> [1] 0
#Kenward's Coefficient of Sociality
Cs(deer[1],deer[2],tc=tc)
#> $Do
#> [1] 422.3163
#> 
#> $De
#> [1] 873.9819
#> 
#> $Cs
#> [1] 0.3484272
#> 
#> $p.Attract
#> [1] 7.403068e-53
#> 
#> $p.Avoid
#> [1] 1
#Shirabe's correlation coefficient
Cr(deer[1],deer[2],tc=tc)
#> [1] 0.3706059
#Dynamic Interaction index
DI(deer[1],deer[2],tc=tc)
#> $DI
#> [1] 0.1510688
#> 
#> $DI.theta
#> [1] 0.1735282
#> 
#> $DI.d
#> [1] 0.5910381
#> 
#> $P.positive
#> [1] 0.001834862
#> 
#> $P.negative
#> [1] 1
#Benhamou's interaction statistic
IAB(deer[1],deer[2],tc=tc,dc=dc)
#> $IAB.obs
#> [1] 0.3986007
#> 
#> $IAB.exp
#> [1] 0.01694352
#> 
#> $P.attract
#> [1] 0.001831502
#> 
#> $P.avoid
#> [1] 1
```

For much more detailed information on the package please see the
documentation and the vignette.

— END —
