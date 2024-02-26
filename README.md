
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wildlifeDI

<!-- badges: start -->

![wildlifeDI CRAN
version](http://www.r-pkg.org/badges/version/wildlifeDI) ![wildlifeDI
Age](https://rpkg.net/pub-age/wildlifeDI/blue) [![wildlifeDI
downloads](https://cranlogs.r-pkg.org/badges/last-month/wildlifeDI)](https://cran.r-project.org/package=wildlifeDI)
<!-- badges: end -->

The wildlifeDI package facilitates the calculation of indices of dynamic
interaction for wildlife telemetry data. There are also functions for
more advanced contact analysis. There are two main streams of analysis
that it facilitates: 1) the calculation of a number of ‘dynamic
interaction’ or ‘spatial-temporal association’ indices betwewn dyads
(i.e., pairs of individuals), and 2) the identification and analysis of
contacts within large tracking datasets. There are currently two
vignettes, one for each of these streams of analysis:

1)  [Dynamic interaction analysis amoung
    dyads](https://cran.r-project.org/package=wildlifeDI/vignettes/wildlifeDI-vignette.html)
2)  [Contact analysis in large
    datasets](https://cran.r-project.org/package=wildlifeDI/vignettes/wildlifeDI-vignette-contact_analysis.html)

For more information on the methods used within see the documentation
which references the methods cited within. The two papers that can be
cited when using the package are:

Long, J.A., Nelson, T.A., Webb, S.L., Gee, K.L. (2014) A critical
examination of indices of dynamic interaction for wildlife telemetry
studies. Journal of Animal Ecology. 83(5):1216-1233.
[Link](https://doi.org/10.1111/1365-2656.12198)

Long, J.A., Webb, S.L., Harju, S.M., Gee, K.L. (2022) Analyzing contacts
and behavior from high frequency tracking data using the wildlifeDI R
package. Geographical Analysis. 54(3):648-663.
[Link](https://doi.org/10.1111/gean.12303)

## Installation

You can install the latest (under development version) of wildlifeDI
from github:

``` r
devtools::install_github("jedalong/wildlifeDI")
```

To download the latest version from CRAN:

``` r
install.packages('wildlifeDI')
```

# Important Update - wildlifeDI version \>= 1.0

In wildlifeDI version 1.0, I have moved from the
[adehabitat](https://cran.r-project.org/package=adehabitatLT)
[ltraj](https://rdrr.io/cran/adehabitatLT/man/as.ltraj.html) objects to
the newer [move2](https://cran.r-project.org/package=move2) class of
objects. The rationale here is that the *move2* objects extend the *sf*
class of spatial objects and therefore easily integrate common data
processing tools such as the *tidyverse* and are easily mapped and
integrated with other spatial analysis workflows. This means all
wildlifeDI functions now expect a *move2* object as input. There is a
new helper function *ltraj_move2* that can help you to convert an
*ltraj* to a *move2* tracking object.

All the analytical functions now accept tracking data in an identical
manner. This is a change to promote consistency across the methods in
wildlifeDI. The user can specify a single tracking dataset with multiple
individuals (2 or more) and the functions will calculate the chosen
metric for all dyad pairs in the dataset with a non-zero temporal
overlap (Carefully read the documentation for the function *checkTO* to
understand how this works). Alternatively, analytical functions can
accept two tracking datasets which allows one to study interactions from
the first tracking dataset to the second. This can be very useful in
multi-species study areas (e.g., comparing interactions of predators and
prey) or in scenarios where there are two groupings within a dataset
that you wish to compare between only (e.g., male and female).

The format of *move2* objects provides numerous efficiencies for quickly
mapping and summarizing tracking data over the *ltraj* object type. For
that reason a number of the helper functions in wildlifeDI have been
removed. These are summarized below:

| Deleted Function | Alternative Functions or Examples         |
|------------------|-------------------------------------------|
| sf2ltraj         | move2_ltraj                               |
| ltraj2sf         | ltraj_move2                               |
| filterTraj       | move2::mt_filter_per_interval             |
| conSpatial       | move2::mt_segments, move2::mt_track_lines |
| conTemporal      | see contact analysis vignette             |
| conSummary       | see contact analysis vignette             |
| conPairs         | conProcess(return=‘contact’)              |
| conMatrix        | see contact analysis vignette             |

— END —
