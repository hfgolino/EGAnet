Exploratory Graph Analysis: An introduction using the EGA package.
=============
### Hudson F. Golino ###
### Assistant Professor of Quantitative Methods, Department of Psychology, University of Virginia. ###
### Contact: <hfg9s@virginia.edu> ###


Introduction
============

Estimating the correct number of dimensions in psychological and educational instruments is a long-standing problem in psychometrics (Preacher & MacCallum, 2003; Ruscio & Roche, 2012; Velicer & Jackson, 1990). There are many techniques and methods to accomplish this task, such as Horn’s (PA; 1965) parallel analysis, the Kaiser-Guttman eigenvalue greater than-one-rule (Guttman, 1954; Kaiser, 1960), Velicer’s minimum average partial procedure (MAP; 1976), the very simple structure (VSS; Revelle & Rocklin, 1979), the use of the Bayesian information criterion (BIC; Schwarz, 1978) or the extended Bayesian information criterion (EBIC; Chen & Chen, 2008) to investigate the best fit over a series of structural models with different number of factors, among other techniques. 
There are several studies investigating the performance of parallel analysis (Buja & Eyubuglu, 1992; Crawford et al., 2010; Green, Redell, Thompson, & Levy, 2016; Keith, Caemmerer, & Reynolds, 2016; Timmerman & Lorenzo-Seva, 2011; Velicer, Eaton, & Fava, 2000, Velicer, 1986; Zwick & Velicer, 1986), multiple average partial procedure (Garrido, Abad, & Ponsoda, 2011; Keith, Caemmerer, & Reynolds, 2016; Velicer, Eaton, & Fava, 2000), BIC (Dziak, Coffman, Lanza, & Li, 2012; Lopes & West, 2004; Preacher, Zhang, Kim, & Mels, 2013; Song & Belin, 2008) and the Kaiser-Guttman eigenvalue rule (Hakstian, Rogers, & Catell, 1982; Keith, Caemmerer, & Reynolds, 2016; Ruscio & Roche, 2012; Velicer et al., 2000; Zwick & Velicer, 1982, 1986) in estimating the correct number of factors. In sum, PA and MAP work quite well when there is a low or moderate correlation between factors, when the sample size is equal to or greater than 500 and when the factor loadings are from moderate to high (Crawford et al., 2010; Green, et al., 2011; Keith, Caemmerer, & Reynolds, 2016; Zwick & Velicer, 1986). However, they tend to underestimate the number of factors when the correlation between factors are high, when the sample size is small and when there is small number of indicators per factor (Crawford et al., 2010; Green, et al., 2011; Keith, Caemmerer, & Reynolds, 2016; Ruscio and Roche, 2012). Regarding BIC, Preacher, Zhang, Kim and Mels (2013) showed that it performs well when the sample size is small, but it tends to overestimate the number of factors in large datasets. Dziak, Coffman, Lanza and Li (2012), by the other side, showed that BIC increases its accuracy in estimating the number of factors when the sample size is greater than 200 cases. The Kaiser-Guttman rule is the default method for choosing the number of factors in many commercial softwares (Bandalos & Boehm-Kaufman, 2009), but simulation studies show that it overestimates the number of factors, especially with a large number of items and a large sample size (Hakstian, Rogers, & Catell, 1982; Keith, Caemmerer, & Reynolds, 2016;  Ruscio & Roche, 2012; Velicer et al., 2000; Zwick & Velicer, 1982, 1986). Ruscio and Roche (2012) provided a startling evidence in this direction, since the Kaiser-Gutman rule overestimated the number of factors in 89.87% of the 10,000 simulated datasets, generated with different number of factors, sample size, number of items, number of response categories per item and correlation between factors. In face of the evidences from the simulation studies, some researchers strongly recommend not to use this method (Bandalos, & Boehm-Kaufman, 2009; Velicer et al., 2000). 

These simulation studies highlight a very complicated problem within psychology, since it is very common to find areas in which the correlation between factors is high, the sample size is up to 500 cases, especially in the intelligence field (Keith, Caemmerer, & Reynolds, 2016). Thus, in such situations, neither PA, MAP or comparing different number of factors via BIC can lead to the correct number of dimensions, proving that estimating the number of factors is still a non-trivial task, in spite of the past decades’ developments. Recently a new technique termed *Exploratory Graph Analysis (EGA)* was proposed by Golino and Epskamp (2017), that seems to overcome the limitations presented by the tradicional methods pointed above. 

What is *Exploratory Graph Analysis*
------------------------------------

Exploratory Graph Analysis is part of a new area called network psychometrics (see Epskamp, Maris, Waldorp & Borsboom, in press), that focuses on the estimation of undirected network models (i.e. Markov Random Fields; Lauritsen, 1996) to psychological datasets. This area has been applied in different areas of psychology, from psychopathology (e.g., Borsboom et al, 2011; Borsboom & Cramer, 2013; Fried et al., 2015) to developmental psychology (van der Maas, 2006), passing through quality of life research (Kossakowski et al., 2015). In network psychometrics, variables are represented as nodes and their relations as edges. When analyzing data from psychological instruments, one may be interest in observing if many nodes are connected with each other, forming clusters, as it is argued clusters of nodes may emerge due to underlying latent variables. If a latent variable model is the true underlying causal model, we would expect indicators in a network model to form strongly connected clusters for each latent variable. Network models can be shown to be mathematically equivalent under certain conditions to latent variable models in both binary (Epskamp, Maris, Waldorp & Borsboom, in press) and Gaussian datasets (Chandrasekaran, Parrilo & Willsky, 2010).  

By defining a cluster as a group of connected nodes regardless of edge weight, Golino and Epskamp (2017) pointed to a fundamental rule of network psychometrics: clusters in network equals latent variables. This is not only a philosophical interpretation of networks, nor solely an empirical finding (Cramer et al., 2010; Cramer et al., 2012; Costantini et al., 2014; Borsboom et al., 2011; Epskamp et al., 2012; van der Maas et al., 2006), but a mathematical characteristic of networks (Golino & Epskamp, 2017). Golino & Epskamp (2017) showed that:

1. If the latent factors are orthogonal, the resulting network model consists of unconnected clusters.
2. Assuming factor loadings and residual variances are reasonably on the same scale for every item, the off-diagonal blocks of the variance-covariance matrix will be scaled closer to zero than the diagonal blocks of the variance-covariance matrix. Hence, the resulting network model will contain weighted clusters for each factor.

Estimating network models can be done by using penalized maximum likelihood estimation, such as the least absolute shrinkage and selection operator (LASSO; Tisbhirani, 1996), one of the most used methods for network estimation on psychological datasets (van Borkulo et al., 2014; Kossakowski et al., 2015; Fried et al., 2015). The LASSO technique avoid overfitting and may result in many parameters to be estimated to exactly equal zero. This indicates conditional independence and facilitates the interpretability of the network structure.

The *Exploratory Graph Analysis* works as follows. Firstly it estimates the correlation matrix of the observable variables, then proceeds to use the *graphical LASSO* estimation to obtain the sparse inverse covariance matrix, with the regularization parameter defined via EBIC over 100 different values. In the last step the *walktrap* algorithm (Pons & Latapy, 2004) is used to find the number of clusters of the partial correlation matrix. In *EGA*, the number of clusters identified equals the number of latent factors in a given dataset. For more details see Golino and Epskamp (2017) and Golino and Demetriou (2017).

The EGA package
===============

The `EGA` was developed as a simple and easy way to implement the *Exploratory Graph Analysis* technique. The package has three main functions:

- `EGA`: Estimates the number of dimensions of a given dataset/instrument using graphical lasso and a random walk algorithm. The glasso regularization parameter is set via EBIC. 
- `bootEGA`: Estimates the number of dimensions of *n* **bootstraps** from the empirical correlation matrix, and returns a typical network (i.e. the network formed by the median pairwise partial correlations over the *n* **bootstraps**) and its dimensionality. 
- `CFA`: Verifies the fit of the structure suggested by `EGA` using confirmatory factor analysis.


References
==========


Borsboom, D., & Cramer, A. O. J. (2013). Network analysis: An integrative approach to the structure of psychopathology. Annual Review of Clinical Psychology, 9, 91-121.

Borsboom, D., Cramer, A. O. J., Schmittmann, V. D., Epskamp, S., & Waldorp, L. J. (2011). The small world of psychopathology. PLoS ONE, 11, e27407

Buja, A., & Eyuboglu, N. (1992). Remarks on parallel analysis. Multivariate Behavioral Research, 27, 509-540.
Chandrasekaran, V., Parrilo, P. A., & Willsky, A. S. (2010, September). Latent variable graphical model selection via convex optimization. In Communication, Control, and Computing (Allerton), 2010 48th Annual Allerton Conference on (pp. 1610-1613). IEEE.

Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.

Condon, D. M. & Revelle, W. (2014). The International Cognitive Ability Resource: Development and initial validation of a public-domain measure. Intelligence, 43, 52–64, DOI: http://dx.doi.org/10.1016/j.intell.2014.01.004

Condon, D.M. & Revelle, W. (2016). Selected ICAR Data from the SAPA-Project: Development and Initial Validation of a Public-Domain Measure. Journal of Open Psychology Data, 4(1), DOI: http://doi.org/10.5334/jopd.25

Costantini, G., Epskamp, S., Borsboom, D., Perugini, M., Mõttus, R., Waldorp, L. J., & Cramer, A. O. (2014). State of the aRt personality research: A tutorial on network analysis of personality data in R. Journal of Research in Personality.

Crawford, A. V., Green, S. B., Levy, R., Lo, W. J., Scott, L., Svetina, D., & Thompson, M. S. (2010). Evaluation of parallel analysis methods for determining the number of factors. Educational and Psychological Measurement, 70(6), 885–901. 

Dinno, A. (2009). Exploring the sensitivity of Horn's parallel analysis to the distributional form of random data. Multivariate Behavioral Research, 44(3), 362-388.

Dziak, J. J., Coffman, D. L., Lanza, S. T., & Li, R. (2012). Sensitivity and specificity of information criteria. PeerJ PrePrints, 3, e1350.

Epskamp, S. (2014). semPlot: Path diagrams and visual analysis of various SEM packages' output. R package version 1.0.1. http://CRAN.R-project.org/package=semPlot

Epskamp, S., Cramer, A. O. J., Waldorp, L.J., Schmittmann, V.D., & Borsboom, D. (2012). qgraph: Network Visualizations of Relationships in Psychometric Data. Journal of Statistical Software, 48(4), 1-18. URL  http://www.jstatsoft.org/v48/i04/.

Epskamp, S., Maris, G., Waldorp, L. J., & Borsboom, D. (in press). Network Psychometrics. To appear in: Irwing, P., Hughes, D., & Booth, T. (Eds.), Handbook of Psychometrics. New York: Wiley.

Formann, A. K., Waldherr, K., & Piswanger, K. (2011). Wiener Matrizen-Test 2 (WMT-2): Ein Rasch-skalierter sprachfreier Kurztest zur Erfassung der Intelligenz. Beltz Test.

Foygel, R., & Drton, M. (2010). Extended bayesian information criteria for gaussian graphical models. In Advances in Neural Information Processing Systems (pp. 604-612).

Fried, E. I., Bockting, C., Arjadi, R., Borsboom, D., Tuerlinckx, F., Cramer, A., Epskamp, S., Amshoff, M., Carr, D., & Stroebe, M. (2015). From Loss to Loneliness: The Relationship Between Bereavement and Depressive Symptoms. Journal of Abnormal Psychology, 124, 256-265. 

Garrido, L. E., Abad, F. J. & Ponsoda, V. (2011). Performance of Velicer’s Minimum Average Partial factor retention method with categorical variables. Educational and Psychological Measurement, 71(3), 551-570.

Garrido, L. E., Abad, F. J., & Ponsoda, V. (2013). A new look at Horn’s parallel analysis with ordinal variables. Psychological Methods, 18(4), 454.

Golino, H. F., & Demetriou, A. (2017). Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis. Intelligence.

Golino, H. F., & Epskamp, S. (2017). Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research. PloS one, 12(6), e0174035.

Green, S. B., Redell, N., Thompson, M. S., & Levy, R. (2016). Accuracy of Revised and Traditional Parallel Analyses for Assessing Dimensionality with Binary Data. Educational And Psychological Measurement, 76(1), 5-21.

Haggerty, M., Terman, L., Thorndike, E., Whipple, G., & Yerkes, R (1920). National Intelligence Tests. Manual of Directions. For Use with Scale A, Form 1 and Scale B, Form 1. New York: World Book Company.  

Horn, J. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30, 179-185.

Horn, J. L. and Engstrom, R. (1979). Cattell’s scree test in relation to bartlett’s chi-square test and other observations on the number of factors problem. Multivariate Behavioral Research, 14(3), 283–300.

Keith, T. Z., Caemmerer, J. M., & Reynolds, M. R. (2016). Comparison of methods for factor extraction for cognitive test-like data: Which overfactor, which underfactor?. Intelligence, 54, 37-54.

Kossakowski, J. J., Epskamp, S., Kieffer, J. M., Borkulo, C. D. van, Rhemtulla, M., & Borsboom, D. (2015). The Application of a Network Approach to Health-Related Quality of Life: Introducing a New Method for Assessing HRQoL in Healthy Adults and Cancer Patients. Quality of Life Research.

Latapy, M., & Pons, P. (2004). Computing communities in large networks using random walks. arXiv preprint cond-mat/0412368.

Lauritzen, S. L. (1996). Graphical models. Clarendon Press.

Must, O. & Must, A., (2014). Data from “Changes in test-taking patterns over time” concerning the Flynn Effect in Estonia. Journal of Open Psychology Data. 2(1), p.e2. DOI: http://doi.org/10.5334/jopd.ab

Preacher, K. J., Zhang, G., Kim, C., & Mels, G. (2013). Choosing the optimal number of factors in exploratory factor analysis: A model selection perspective. Multivariate Behavioral Research, 48(1), 28-56.

R Core Team (2016). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.

Revelle, W., & Rocklin, T. (1979). Very simple structure: An alternative procedure for estimating the optimal number of interpretable factors. Multivariate Behavioral Research, 14(4), 403-414.

Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. URL  http://www.jstatsoft.org/v48/i02/.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society. Series B (Methodological), 267-288.

Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment of ordered polytomous items with parallel analysis. Psychological Methods, 16, 209-220.

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for constructing networks from binary data. Scientific Reports, 4, 5918.

Van Der Maas, H. L., Dolan, C. V., Grasman, R. P., Wicherts, J. M., Huizenga, H. M., & Raijmakers, M. E. (2006). A dynamical model of general intelligence: the positive manifold of intelligence by mutualism. Psychological review, 113(4), 842.

Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations. Psychometrika, 41(3), 321-327.

Velicer, W. F., & Jackson, D. N. (1990). Component analysis versus common factor analysis: Some issues in selecting an appropriate procedure. Multivariate Behavioral Research, 25(1), 1-28.
