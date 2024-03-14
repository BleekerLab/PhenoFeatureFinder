# README

This example dataset comes from the MicroMass dataset published at UCI Machine Learning Repositiory. 

The data used here is the reference panel of 20 Gram positive and negative bacterial species. The samples are clustered per species, with each instance of a species as a replicate. 

The Gram types are used as phenotypes.

![taxonomy](./taxonomy.png)


| species                       | replicates | Gram type |
|-------------------------------|------------|-----------|
| *Bacillus cereus*             | 26         | positive  |
| *Bacillus thuringiensis*      | 11         | positive  |
| *Citrobacter braakii*         | 26         | negative  |
| *Citrobacter freundii*        | 28         | negative  |
| *Clostridium difficile*       | 14         | positive  |
| *Clostridium glycolicum*      | 16         | positive  |
| *Enterobacter asburiae*       | 29         | negative  |
| *Enterobacter cloacae*        | 52         | negative  |
| *Escherichia coli*            | 60         | negative  |
| *Haemophilus influenzae*      | 50         | negative  |
| *Haemophilus parainfluenzae*  | 21         | negative  |
| *Listeria ivanovii*           | 29         | positive  |
| *Listeria monocytogenes*      | 31         | positive  |
| *Shigella boydii*             | 18         | negative  |
| *Shigella flexneri*           | 32         | negative  |
| *Shigella sonnei*             | 31         | negative  |
| *Streptococcus mitis*         | 26         | positive  |
| *Streptococcus oralis*        | 24         | positive  |
| *Yersinia enterocolitica*     | 27         | negative  |
| *Yersinia frederiksenii*      | 20         | negative  |


## Workflow

This example starts with two datasets: 'feature_data_unfiltered.csv' with the unfiltered mass spectrometry data and 'phenotype.csv' with Gram types of each sample as explained above.

In 'feature_filtering.ipynb', the mass spectrometry dataset is filtered based on the presence of each feature in at least 11 replicates of at least one species. This threshold is chosen, because this is the smallest number of replicates in the data (the sample size of *Bacillus thuringiensis*). There are no blanks in the dataset, so the blank filtering step is omitted. The filtered dataset is exported as 'filtered_features.csv'.

In 'feature_selection.ipynb' the filtered dataset and 'phenotype.csv' are used as input to build a Machine Learning pipeline for the classification of Gram positive and negative bacteria. This results in the best fitting pipeline, outputted as 'best_fitting_pipeline.py', and the option to extract the most important features for the classification.

## Reference 

Mah,Pierre and Veyrieras,Jean-Baptiste. (2013). MicroMass. UCI Machine Learning Repository. https://doi.org/10.24432/C5T61S.


**Link:** [https://doi.org/10.24432/C5T61S](https://doi.org/10.24432/C5T61S)

