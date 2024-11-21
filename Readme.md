<img src="inst/EGAnet_hex_2.png" width = 250 />

<div id="badges"><!-- pkgdown markup -->
<a href="https://CRAN.R-project.org/package=EGAnet"><img border="0" src="https://www.r-pkg.org/badges/version/EGAnet?color=blue" alt="CRAN version"/></a>
<a href="https://github.com/hfgolino/EGAnet/releases"><img src="https://img.shields.io/github/v/release/hfgolino/EGAnet" alt="GitHub version"/></a>
<a href="https://github.com/hfgolino/EGAnet/actions/workflows/r.yml"><img border="0" src="https://github.com/hfgolino/EGAnet/actions/workflows/r.yml/badge.svg" alt="R CMD check status"/></a> </br>
<a href="https://www.repostatus.org/#active"><img border="0" src="https://www.repostatus.org/badges/latest/active.svg" alt="Project Status: Active – The project has reached a stable, usable state and is being actively developed."/></a>
<a href="https://r-ega.net"><img border="0" src="https://cranlogs.r-pkg.org/badges/grand-total/EGAnet?color=blue" alt="Downloads Total"/></a>
<a href="https://r-ega.net"><img border="0" src="http://cranlogs.r-pkg.org/badges/EGAnet?color=blue" alt="Downloads per month"/></a>
<a href="https://r-ega.net"><img border="0" src="http://cranlogs.r-pkg.org/badges/last-day/EGAnet" alt="Downloads Yesterday"/></a>
</div>

### To get started, check out the website: [r-ega.net](https://r-ega.net) ###

How to Install
=============
```
if(!"devtools" %in% row.names(installed.packages())){
  install.packages("devtools")
}

devtools::install_github("hfgolino/EGAnet")
```

Exploratory Graph Analysis
=============
### Hudson F. Golino ###
### Associate Professor of Quantitative Methods, Department of Psychology, University of Virginia. ###
### Contact: <hfg9s@virginia.edu> ###

### Alexander P. Christensen ###
### Assistant Professor of Quantitative Methods, Department of Psychology and Human Development, Vanderbilt University
###
### Contact: <alexander.christensen@vanderbilt.edu> ###


News
============

The EGAnet package version 2.0.0 (Blue Metallic) is now available on [CRAN](https://cran.r-project.org/web/packages/EGAnet/index.html).

You can read about some of the changes in this dedicated [Wiki](https://github.com/hfgolino/EGAnet/wiki/What's-Changed%3F). There's also a [Wiki for the new plotting system](https://github.com/hfgolino/EGAnet/wiki/Plots-in-%7BEGAnet%7D), and a series of posts about the package on Twitter: [some new functions](https://twitter.com/GolinoHudson/status/1691800126866829739?s=20), [new stability analysis](https://twitter.com/GolinoHudson/status/1684912389194436610?s=20), and [more](https://twitter.com/GolinoHudson/).

The EGAnet package is currently supported by two University of Virginia grants, one from the STAR - [Support Transformative Autism Research](https://curry.virginia.edu/faculty-research/centers-labs-projects/supporting-transformative-autism-research-star) initiative and one from the [Democracy Initiative](http://democracyinitiative.virginia.edu).

References
============

Christensen, A. P. (2023). Unidimensional community detection: A Monte Carlo simulation, grid search, and comparison. *PsyArXiv*. doi:[10.31234/osf.io/ep3vx](https://doi.org/10.31234/osf.io/ep3vx)
+ Related functions: `community.unidimensional`

Christensen, A. P., Garrido, L. E., & Golino, H. (2023). Unique variable analysis: A network psychometrics method to detect local dependence. *Multivariate Behavioral Research*. doi:[10.1080/00273171.2023.2194606](https://doi.org/10.1080/00273171.2023.2194606)
+ Related functions: `UVA`

Christensen, A. P., Garrido, L. E., Guerra-Pena, K., & Golino, H. (2023). Comparing community detection algorithms in psychometric networks: A Monte Carlo simulation. *Behavior Research Methods*. doi:[10.31234/osf.io/hz89e](https://doi.org/10.31234/osf.io/hz89e)
+ Related functions: `EGA`

Christensen, A. P., & Golino, H. (2021).
Estimating the stability of psychological dimensions via Bootstrap Exploratory Graph Analysis: A Monte Carlo simulation and tutorial. *Psych*, *3*(3), 479-500.
doi:[10.3390/psych3030032](https://doi.org/10.3390/psych3030032)
+ Related functions: `bootEGA`, `dimensionStability`, and `itemStability`

Christensen, A. P., & Golino, H. (2021). Factor or network model? Predictions from neural networks. *Journal of Behavioral Data Science*, *1*(1), 85-126. doi:[10.35566/jbds/v1n1/p5](https://doi.org/10.35566/jbds/v1n1/p5)
+ Related functions: `LCT`

Christensen, A. P., & Golino, H. (2021). On the equivalency of factor and network loadings. *Behavior Research Methods*, *53*, 1563–1580. doi:[10.3758/s13428-020-01500-6](https://doi.org/10.3758/s13428-020-01500-6)
+ Related functions: `LCT` and `net.loads`

Christensen, A. P., Golino, H., & Silvia, P. J. (2020). A psychometric network perspective on the validity and validation of personality trait questionnaires. *European Journal of Personality*, *34*, 1095-1108. doi:[10.1002/per.2265](https://doi.org/10.1002/per.2265)
+ Related functions: `bootEGA`, `dimensionStability`, `EGA`, `itemStability`, and `UVA`

Christensen, A. P., Gross, G. M., Golino, H., Silvia, P. J., & Kwapil, T. R. (2019). Exploratory graph analysis of the Multidimensional Schizotypy Scale. *Schizophrenia Research*, *206*, 43-51. doi:[10.1016/j.schres.2018.12.018](https://doi.org/10.1016/j.schres.2018.12.018)
+ Related functions: `CFA` and `EGA`

Garcia-Pardina, A., Abad, F. J., Christensen, A. P., Golino, H., & Garrido, L. E. (2022). Dimensionality assessment in the presence of wording effects: A network psychometric and factorial approach. *PsyArXiv*. doi:[10.31234/osf.io/7yqau](https://doi.org/10.31234/osf.io/7yqau)
+ Related functions: `riEGA`

Golino, H., Christensen, A. P., Moulder, R. G., Kim, S., & Boker, S. M. (2022). Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections. *Psychometrika*, *87*(1), 156-187. doi:[10.1007/s11336-021-09820-y](https://doi.org/10.1007/s11336-021-09820-y)
+ Related functions: `dynEGA`, `net.loads`, and`simDFM`

Golino, H., & Demetriou, A. (2017). Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis. *Intelligence*, *62*, 54-70. doi:[10.1016/j.intell.2017.02.007](https://doi.org/10.1016/j.intell.2017.02.007)
+ Related functions: `EGA`

Golino, H., & Epskamp, S. (2017). Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research. *PLoS ONE*, *12*, e0174035. doi:[10.1371/journal.pone.0174035](https://doi.org/10.1371/journal.pone.0174035)
+ Related functions: `CFA`, `bootEGA`, and `EGA`

Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Neito, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020). Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables. *Multivariate Behavioral Research*. doi:[10.31234/osf.io/mtka2](https://doi.org/10.31234/osf.io/mtka2)
+ Related functions: `entropyFit`, `tefi`, and `vn.entropy`

Golino, H., Nesselroade, J. R., & Christensen, A. P. (2022). Towards a psychology of individuals: The ergodicity information index and a bottom-up approach for finding generalizations. *PsyArXiv*. doi:[10.31234/osf.io/th6rm](https://doi.org/10.31234/osf.io/th6rm)
+ Related functions: `boot.ergoInfo`, `ergoInfo`, `jsd`, and `infoCluster`

Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., Thiyagarajan, J. A., & Martinez-Molina, A. (2020). Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial. *Psychological Methods*, *25*, 292-320. doi:[10.1037/met0000255](https://doi.org/10.1037/met0000255)
+ Related functions: `EGA`

Golino, H., Thiyagarajan, J. A., Sadana, M., Teles, M., Christensen, A. P., & Boker, S. M. (under review). Investigating the broad domains of intrinsic capacity, functional ability, and environment: An exploratory graph analysis approach for improving analytical methodologies for measuring healthy aging. *PsyArXiv*. doi:[10.31234/osf.io/hj5mc](https://doi.org/10.31234/osf.io/hj5mc)
+ Related functions `EGA.fit` and `tefi`

Jamison, L., Christensen, A. P., & Golino, H. (2021). Optimizing Walktrap's community detection in networks using the Total Entropy Fit Index. *PsyArXiv*. doi:[10.31234/osf.io/9pj2m](https://doi.org/10.31234/osf.io/9pj2m)
+ Related functions: `EGA.fit` and `tefi`
  
Jamison, L., Golino, H., & Christensen, A. P. (2023). Metric invariance in exploratory graph analysis via permutation testing. *PsyArXiv*. doi:[10.31234/osf.io/j4rx9](https://doi.org/10.31234/osf.io/j4rx9)
+ Related functions: `invariance`

Jiménez, M., Abad, F. J., Garcia-Garzon, E., Golino, H., Christensen, A. P., & Garrido, L. E. (2023). Dimensionality assessment in bifactor structures with multiple general factors: A network psychometrics approach. *Psychological Methods*. doi:[10.1037/met0000590](https://doi.org/10.1037/met0000590)
+ Related functions: `hierEGA` and `net.scores`
  
Shi, D., Christensen, A. P., Day, E., Golino, H., & Garrido, L. E. (2023). A Bayesian approach for dimensionality assessment in psychological networks. *PsyArXiv*. doi:[10.31234/osf.io/9rcev](https://doi.org/10.31234/osf.io/9rcev)
+ Related functions: `EGA`
