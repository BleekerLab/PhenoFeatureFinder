---
title: 'PhenoFeatureFinder: a python package for linking developmental phenotypes to omics features'
tags:
  - insect development
  - phenotyping
  - metabolomics
  - omics
  - feature selection
  - preprocessing
authors:
  - name:
      given-names: Lissy-Anne M.
      surname: Denkers
    orcid: 0009-0004-8222-8689
    affiliation: 1
  - name:
      given-names: Marc D.
      surname: Galland
    orcid: 0000-0003-2161-8689
    affiliation: 2
  - name:
      given-names: Annabel
      surname: Dekker
    affiliation: 3
  - name:
      given-names: Valerio
      surname: Bianchi
    orcid: 0000-0001-9025-7923
    affiliation: 3
  - name:
      given-names: Petra M.
      surname: Bleeker
    orcid: 0000-0003-2222-9263
    affiliation: 1
affiliations:
  - name: University of Amsterdam, Department of Plant Physiology, Green Life Science Research Theme, Swammerdam Institute for Life Sciences, Amsterdam, The Netherlands
    index: 1
  - name: INRAE, Institute of Genetics, Environment and Plant Protection (IGEPP—Joint Research Unit 1349), Le Rheu, France
    index: 2
  - name: Enza Zaden R&D B.V., BTR-BM Bioinformatics, Enkhuizen, The Netherlands
    index: 3
date: 5 June 2024
bibliography: paper.bib
---

# Summary

`PhenoFeatureFinder` is designed to facilitate the analyses required to analyse quantitative and/or progressive phenotypic- and omics data, and link those using Machine Learning with the aim to identify causal features, in one package. It can be used for 1) evaluation and visualisation of phenotype progression over multiple stages and between groups (e.g. treatments, genotypes), 2) pre-processing of omics data, and 3) prediction of features that explain the phenotypic classification. To facilitate usability, each step in the pipeline can also be performed independently, hence has been assigned a class in the package (\autoref{fig1}). We provide an example of implementation below that focuses on insect development through time and the selection of metabolic features causal to the observed phenotype, but different input data could be used, provided it has a similar structure. This could be any phenotype that is scored in progressive stages over time. Also, `PhenoFeatureFinder` was developed initially with metabolomics data, but users can evaluate its fit applying other types of omics data.

![Overview of the package, consisting of three classes that can be used separately or as a workflow. Class 1: analysing and visualising the phenotype, Class 2: preprocessing and visualising omics datasets, and Class 3: feature selection through a Machine Learning approach.  \label{fig1}](./package_figure.png)

# Statement of need

The analysis of developmental phenotypes can be challenging, due to the many variables involved (e.g. time, developmental stages, replicates, treatments), especially for researchers whose strength or interest does not lie in data analysis. The same goes for the pre-processing of omics data and linking the omics data and developmental phenotypes. With `PhenoFeatureFinder`, we aim to support such research by combining the nessecary functionalities in one package with easy to follow manuals and examples. 

In R, the package `drc` is available for fitting dose-response curves [@Ritz2015], offering an extensive and versatile set of functionalities. However, for the purposes described here `drc` poses some limitations, such as the options for custom pre-processing and analyses of multiple experimental groups simultaneously. Here we implemented pre-processing steps and aimed to decrease the amount of coding needed to obtain a fitted development curve.

# Use case example

Plants interact with their (a)biotic environment through a range of specialised metabolites and deal with pathogens and pest attack through constitutive or inducible production of those defence molecules [@Erb2020; @García-Olmedo1998]. High-throughput “-omics” tools including (untargeted) metabolomics have been successfully implemented in plant biology [@Dalio2021], but the accompanying resistance phenotyping often lacks in robustness [@Song2021]. 

Proliferation of an insect population is affected by various factors, including the chemical composition of the host, and/or the environment [@Ma2022]. In particular, host resistance via hampered larval development is noteworthy, because reducing the speed at which larvae reach the adult stage and produce offspring negatively affects pest-population development [@Maharijaya2019; @Muema2016; @Vengateswari2022]. However, evaluating larval development results in a complex dataset that is challenging to process. Developmental success is based on the number of larvae throughout various larval stages, as well as on the speed of development. 

To identify underlying mechanisms of resistance, the chemical or molecular composition of a plant can be investigated. Proteins and metabolites are commonly analysed through untargeted Mass-Spectrometry, yielding exhaustive profiles generally consisting of many thousands of unannotated features. Often such data displays sparsity, i.e. missing values between datasets, and a low sample-to-feature ratio, adding to the complexity of the analysis [@Kortbeek2021; @Liebal2020]. Tree-based Machine-Learning algorithms (e.g., random forest) are suitable for the analysis of, and feature selection from, untargeted data [@Liebal2020] computing the contribution of each feature in the phenotypic classification. 

## Class I: `PhenotypeAnalysis`

A binary classification of plants into "resistant" or "susceptible" helps to extract relevant features especially when threshold effects or sparsity (presence/absence) effects are at play. Here we firstly assess performance over different developmental stages of larvae on different host plants. The number of individuals in each stage at a given time is recorded. When plotted, the cumulative data of these bioassays resemble a growth- or dose-response curve that can be used to manually assign a binary phenotype (e.g., resistant/non-resistant), a resistance classification used as input for `FeatureSelection` (Class 3). 

To account for missing data when individuals that reached the final developmental stage are removed from the experiment, we implemented an automated correction step. The count data can be transformed to cumulative data to analyse the maximum of individuals that reach each of the developmental stages. Next, the time to reach a specific stage can be compared between treatments by fitting a 3-parameter log-logistic curve [@Muse2021; @Seefeldt1995; @Vliet2013] to the cumulative data for each treatment, with the function:   

$$ f(x) = \frac{m}{1 + \exp(s \times (\log(x) - \log(e_{50})))} $$

where $x$ is time, $m$ is the upper limit (or maximum of individuals that developed to the stage of interest), $s$ is the slope of the linear part of the curve and $e_{50}$ is the EmT50 (the timepoint at which 50% of the individuals have developed to the stage of interest). We added the possibility to compare performance between treatments by fitting a curve with the function:

$$ f(x) = \frac{a \times \frac{s}{m} \times (\frac{x}{m})^{s-1}}{1 + (\frac{x}{m})^{s}} $$

Here, $x$ is time, $a$ the area under the curve, $s$ is the shape of the curve and $m$ the median time point. Both functions output a table with the model parameters, confidence intervals and the model fit, together with a plot displaying the observed data and the fitted model. For both functions it is possible to predict the potential maximum beyond the final experimental measurements.

## Class II: `OmicsAnalysis`

Untargeted omics results in large datasets that tend to contain background noise and unreliable features. To clean the data, multiple filtering methods are implemented in the `OmicsAnalysis` class, including the removal of contaminants present in blank samples, filtering to decrease sparsity and other quality control steps. The structure of the data can subsequently be visualised with a PCA and an UpSet plot.  

## Class III: `FeatureSelection`

Combining the output of Classes 1 and 2, i.e. the binary phenotype classification and the tidied untargeted metabolomics, `FeatureSelection` is set up to predict features that can explain the phenotypic observation under study. This part of the pipeline was built as a wrapper around the Python libraries `scikit-learn` and `TPOT` [@Olson2016; @Pedregosa2011]. The `FeatureSelection` wrapper is designed to select optimal pipelines for data preprocessing and identification of the most suitable Machine Learning model. One characteristic of metabolomics data is strongly correlated features (linear dependencies between variables) that make it difficult to extract individual feature importance. Therefore, this method implements a PCA as dimensionality reduction method before searching for the best fitting pipeline. Finally, the importance of the Principal Components and their most related features (high loadings) can be retrieved to select features with predicted importance to the phenotypic classification.  

# Acknowledgements

The Amsterdam Data Science Centre is acknowledged for their input and Frans van der Kloet for his statistical support. Gea van der Lee and Piet Verdonschot are kindly thanked for providing the data from [@Lee2020] used in the example in GitHub. 

# Author contributions

The software was written by Lissy-Anne Denkers (LD) and Marc Galland (MG), with input imput from Petra Bleeker (PB) and tested by Annabel Dekker (AD) and Valerio Bianchi (VB). The manuals and examples were written by LD with major input from AD and VB. The manuscript was designed, written and revised by LD, MG, AD, VB and PB.

# References