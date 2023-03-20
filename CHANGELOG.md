# Changelog

<!--next-version-placeholder-->

## v0.4.0 (17-03-2023)

- Since features are highly correlated, it is very difficult to estimate the individual feature importance using permutation_importance directly on the feature data. Rather a dimensionality reduction method is being used and importance is calculated directly on the Principal Components. This still gives an idea of the model performance. 
- Add a function to extract the features most important per Principal Component. function is called `get_names_of_top_n_features_from_selected_pc()`
- Add a unit test for the `get_names_of_top_n_features_from_selected_pc()` function. 

## v0.3.0 (09-01-2023)

- Add a function to compute metabolome matrix sparsity.
- Add a function to create a histogram/density plot of feature peak areas for each group e.g. genotype. 
- Add a function to filter the features per group based on a selected peak area percentile. 

## v0.2.0 (06-01-2023)

- Add search for best model and computes feature importances.
- Fix small bug with PCA scree plot.

## v0.1.0 (31-10-2022)

- First release of `phloemfinder`!

