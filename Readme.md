### CRAN 1.2.3 | GitHub 1.2.4

<img src="inst/EGAnet_hex.png" width = 250 />

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![Downloads Total](https://cranlogs.r-pkg.org/badges/grand-total/EGAnet?color=brightgreen)](https://cran.r-project.org/package=EGAnet) [![Downloads per month](http://cranlogs.r-pkg.org/badges/EGAnet?color=brightgreen)](https://cran.r-project.org/package=EGAnet) 

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
### Assistant Professor of Quantitative Methods, Department of Psychology and Human Development, Vanderbilt University ###
### Contact: <alexander.christensen@vanderbilt.edu> ###


News
============

The EGAnet package is currently supported by two University of Virginia grants, one from the STAR - [Support Transformative Autism Research](https://curry.virginia.edu/faculty-research/centers-labs-projects/supporting-transformative-autism-research-star) initiative and one from the [Democracy Initiative](http://democracyinitiative.virginia.edu).


The old EGA package is now EGAnet, and is now available in The Comprehensive R Archive Network (CRAN): https://cran.r-project.org/web/packages/EGAnet/index.html


References
============

Chrisensen, A. P., Garrido, L. E., & Golino, H. (in press). Unique variable analysis: A network psychometrics method to detect local dependence. *Multivariate Behavioral Research*. doi:[10.31234/osf.io/4kra2](https://doi.org/10.31234/osf.io/4kra2)
+ Related functions: `UVA`

Christensen, A. P., Garrido, L. E., Guerra-Pena, K., & Golino, H. (in press). Comparing community detection algorithms in psychometric networks: A Monte Carlo simulation. *Behavior Research Methods*. doi:[10.31234/osf.io/hz89e](https://doi.org/10.31234/osf.io/hz89e)
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

Golino, H., Christensen, A. P., Moulder, R. G., Kim, S., & Boker, S. M. (2022). Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections. *Psychometrika*, *87*(1), 156-187. doi:[10.1007/s11336-021-09820-y](https://doi.org/10.1007/s11336-021-09820-y)
+ Related functions: `dynEGA` and `simDFM`

Golino, H., & Demetriou, A. (2017). Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis. *Intelligence*, *62*, 54-70. doi:[j.intell.2017.02.007](https://www.sciencedirect.com/science/article/pii/S0160289616302240)
+ Related functions: `EGA`

Golino, H., & Epskamp, S. (2017). Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research. *PLoS ONE*, *12*, e0174035. doi:[10.1371/journal.pone.0174035](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0174035)
+ Related functions: `EGA`

Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Neito, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020). Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables. *Multivariate Behavioral Research*. doi:[10.31234/osf.io/mtka2](https://doi.org/10.31234/osf.io/mtka2)
+ Related functions: `entropyFit`, `tefi`, and `vn.entropy`

Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., Thiyagarajan, J. A., & Martinez-Molina, A. (2020). Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial. *Psychological Methods*, *25*, 292-320. doi:[10.1037/met0000255](https://psycnet.apa.org/doiLanding?doi=10.1037/met0000255)
+ Related functions: `EGA`

Golino, H., Thiyagarajan, J. A., Sadana, M., Teles, M., Christensen, A. P., & Boker, S. M. (under review). Investigating the broad domains of intrinsic capacity, functional ability, and environment: An exploratory graph analysis approach for improving analytical methodologies for measuring healthy aging. *PsyArXiv*. doi:[10.31234/osf.io/hj5mc](https://doi.org/10.31234/osf.io/hj5mc)
+ Related functions `EGA.fit` and `tefi`

