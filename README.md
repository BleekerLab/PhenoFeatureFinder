# PhenoFeatureFinder

**Linking developmental phenotypes to metabolic features**

`PhenoFeatureFinder` is divided into three classes:
* `PhenotypeAnalysis`
* `OmicsAnalysis`
* `FeatureSelection`

![Overview of the package](documentation/paper/package_figure.png)

`PhenotypeAnalysis` was designed to analyse the development of insects through their larval stages over time in multiple treatments (e.g. hosts). These types of phenotyping analyses can be challenging, due to the many variables involved (e.g. time, developmental stages, replicates, treatments), especially for researchers whose strength does not lie in data analysis. `PhenotypeAnalysis` offers a set of functions to visualise the development while taking into account those different variables, and to perform the necessary data preprocessing steps. From the output, it is easy to manually assign binary phenotypes to your groups (treatments, genotypes, etc), if you want to use it as input for `FeatureSelection`. 

With `OmicsAnalysis`, you can filter large untargeted metabolomics datasets and visualise the structure of the data. The filtered data and a corresponding set of binary phenotypes can then be used as input for `FeatureSelection`. With only a few lines of code, the best fitting pipeline to link the phenotypes to metabolic features is created using Automated Machine Learning with `TPOT` and `scikit-learn`.

Although `OmicsAnalysis` and `FeatureSelection` are designed for metabolomics data, they might also be used for other types of omics data. The user would have to keep in mind that the functions were written for the specifics of metabolomics data (high sparsity, strongly correlated features) and first assess the fit for other types of data. 

## Installation

**Needs update**

```bash
$ pip install phloemfinder
```

## Usage

For each of the classes, you can find a manual with an explanation for all of their functions in the [manuals folder](documentation/manuals/).

If you want to see an example of how `PhenoFeatureFinder` can be used for real-world data, you can take a look at one of the two [examples](documentation/examples/). The first example showcases the use of the [`PhenotypeAnalysis` class](documentation/examples/caddisfly/) for the analysis of the development of caddisfly larvae in four freshwater streams. In the second example, the [`OmicsAnalysis` and `FeatureSelection` classes](documentation/examples/MicroMass/) are used to analyse and select interesting features from a mass spectrometry dataset of a panel of bacterial species.

### Dependencies

Required for all classes:
- NumPy
- pandas
- Matplotlib
- seaborn

Additionally required for PhenotypeAnalysis:
- SciPy

Additionally required for OmicsAnalysis:
- scikit-learn
- UpSetPlot

Additionally required for FeatureSelection:
- scikit-learn
- TPOT
- auto-sklearn (auto-sklearn is made for Linux operating systems. On macOS it needs to be installed manually with brew and pip. You can do this by following [these instructions](https://gist.github.com/simonprovost/051952533680026b67fa58c3552b8a7b).)

## Citation

Insert citation option when ready

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`PhenoFeatureFinder` was created by Lissy-Anne Denkers and Marc Galland, with input from Annabel Dekker, Valerio Bianchi and Petra Bleeker. It is licensed under the terms of the Apache License 2.0 license.

## Credits

`PhenoFeatureFinder` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).

## Useful reading

- [Auto-sklearn talks](https://github.com/automl/auto-sklearn-talks)
- [NumPy docstring examples](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html#example-numpy)
