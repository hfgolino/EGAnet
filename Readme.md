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

### Exploratory Graph Analysis: a framework for estimating the number of dimensions in multivariate data using network psychometrics ###

### To get started, check out the website: [r-ega.net](https://r-ega.net) ###

What is EGAnet?
=============

{EGAnet} implements the Exploratory Graph Analysis (EGA) framework for dimensionality and psychometric assessment. Instead of fitting a factor model, EGA represents measured variables as a network — items are nodes, their partial correlations are edges — and applies community detection algorithms to estimate how many dimensions organize the data and which items belong to each one. On top of this core method, the package provides bootstrap procedures for assessing the stability of dimensions and items, Unique Variable Analysis for detecting local dependence (redundancy) between items, network loadings that parallel factor loadings and can be used to compute network scores, configural and metric invariance testing, Hierarchical EGA for nested structures, and Dynamic EGA for time series and intensive longitudinal data at the individual, group, or population level.

How to Install
=============

From CRAN (recommended):
```r
install.packages("EGAnet")
```

Development version from GitHub:
```r
if(!"devtools" %in% row.names(installed.packages())){
  install.packages("devtools")
}

devtools::install_github("hfgolino/EGAnet")
```

Quick Start
=============

```r
library(EGAnet)

# Wiener Matrizen Test 2 (WMT-2): 18 fluid-intelligence items
colnames(wmt2)
ega_wmt <- EGA(wmt2[, 7:24])
plot(ega_wmt)
```

`EGA()` estimates the network, detects communities (dimensions), and returns a plot showing which items group together. For a full walkthrough and next steps (bootstrapping, plotting options, invariance testing, dynamic data, and more), see the [Quick Start guide](https://r-ega.net/articles/quick-start.html) and [Workflows](https://r-ega.net) on the website.

Key Features
=============

+	**EGA core functions** — `EGA`, `bootEGA`, `dynEGA`, `hierEGA`, `riEGA`: estimate and validate the number of dimensions in cross-sectional, hierarchical, and longitudinal data.
+	**Psychometric tools** — `itemStability`, `invariance`, `net.loads`, `net.scores`, `UVA`, `CFA`: assess item/dimension stability, measurement invariance, network loadings and scores, local dependence, and confirmatory comparisons.
+	**Exploratory Graph Model** — `simEGM`: simulate data from an exploratory graph model for methods research and power analysis.
+	**Information theory** — `tefi`, `entropyFit`, `ergoInfo`, `totalCor`: fit indices and information-theoretic measures for evaluating dimensional structure.
+	**Visualization** — `plot_clusters`, `compare.EGA.plots`, `color_palette_EGA`: publication-ready plots for networks, dimension comparisons, and dynamic clusters.
+	**Foundational network methods** — `TMFG`, `polychoric.matrix`, `community.detection`, `community.consensus`: the underlying network estimation and community detection building blocks used throughout the package.

See the full [function reference](https://r-ega.net/reference/index.html) on the website for every exported function.

Learn More
=============

+	[Quick Start guide](https://r-ega.net/articles/quick-start.html)
+	[Workflows](https://r-ega.net/workflows/index.html) — step-by-step guides through the analysis process
+	[Plotting](https://r-ega.net/plotting.html) — customizing and combining network plots
+	[Articles](https://r-ega.net/articles/index.html) — tutorials and methods write-ups
+	[Publications](https://r-ega.net/publications.html) — papers behind the package
+	[News](https://r-ega.net/news/index.html) — release notes and announcements

Authors & Contact
=============

### Hudson F. Golino ###
Associate Professor of Quantitative Methods, Department of Psychology, University of Virginia

Contact: <hfg9s@virginia.edu>

### Alexander P. Christensen ###
Assistant Professor of Quantitative Methods, Department of Psychology and Human Development, Vanderbilt University

Contact: <alexander.christensen@vanderbilt.edu>

Funding
=============

The EGAnet package is currently supported by two University of Virginia grants, one from the STAR - [Support Transformative Autism Research](https://curry.virginia.edu/faculty-research/centers-labs-projects/supporting-transformative-autism-research-star) initiative and one from the [Democracy Initiative](http://democracyinitiative.virginia.edu).

For release notes and updates, see the [NEWS](NEWS) file or the website's [News](https://r-ega.net/news/index.html) page.

Citing EGAnet
=============

If you use {EGAnet} in your research, please cite:

> Golino, H., & Christensen, A. P. (2026). EGAnet: Exploratory Graph Analysis – A framework for estimating the number of dimensions in multivariate data using network psychometrics. doi:[10.32614/CRAN.package.EGAnet](https://doi.org/10.32614/CRAN.package.EGAnet)

Or, in R, run `citation("EGAnet")` for the citation matching your installed version.

References
=============

<details>
<summary>Full reference list</summary>

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

</details>

Contributing / Issues
=============

Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/hfgolino/EGAnet/issues) — please use the provided [bug report](.github/ISSUE_TEMPLATE/bug-report.md) or [feature request](.github/ISSUE_TEMPLATE/feature-request.md) templates.

License
=============

{EGAnet} is licensed under [AGPL (>= 3.0)](LICENSE).
