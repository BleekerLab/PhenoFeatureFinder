# README

This example dataset comes from a publication on Arabidopsis metabolomics from a set of 309 accessions subjected to control or stress conditions. 

Since the control and stress conditions are easily separated (see PCA below) based on the metabolomic profiles, it is possible to use this dataset to benchmark the Machine Learning module of the phloemfinder package for the classification method. The method should be able to identify metabolites able to separate the two conditions. 

![PCA](./Wu_2017_pca.png)


## Workflow

Since no biological replicate is present, the accessions were clustered in groups following their metabolomic profile in order to make artificial biological replicates. 

__Control condition__ the 309 accessions were grouped in 34 clusters using the `create_bio_replicates_from_ecotypes_control_condition.R` (see at the bottom of this README file for exact number of replicates per cluster). Only clusters with at least 4 "replicates" were kept. 

__Stress condition__ the 309 accessions were grouped in XX clusters using the `create_bio_replicates_from_ecotypes_stress_condition.R` (see at the bottom of this README file for exact number of replicates per cluster). Only clusters with at least 4 "replicates" were kept. 




## Reference 

Wu S., Tohge T., Cuadros-Inostroza A´ ., Tong H., Tenenboim H., Kooke R., Me´ ret M., Keurentjes J.B.,
Nikoloski Z., Fernie A.R., Willmitzer L., and Brotman Y. (2018). Mapping the Arabidopsis Metabolic Landscape
by Untargeted Metabolomics at Different Environmental Conditions. Mol. Plant. 11, 118–134.

**Link:** [https://pubmed.ncbi.nlm.nih.gov/28866081/](https://pubmed.ncbi.nlm.nih.gov/28866081/)

## Number of accessions per group and per condition

### Control condition

|         cluster_id    	|                 n_accessions              	|
|-----------------------	|-------------------------------------------	|
|     1                 	|     8                                     	|
|     2                 	|     15                                    	|
|     3                 	|     11                                    	|
|     4                 	|     28                                    	|
|     5                 	|     4                                     	|
|     6                 	|     10                                    	|
|     7                 	|     11                                    	|
|     8                 	|     17                                    	|
|     9                 	|     12                                    	|
|     10                	|     5                                     	|
|     11                	|     9                                     	|
|     12                	|     16                                    	|
|     13                	|     8                                     	|
|     14                	|     8                                     	|
|     15                	|     9                                     	|
|     16                	|     4                                     	|
|     17                	|     1                                     	|
|     18                	|     9                                     	|
|     19                	|     3                                     	|
|     20                	|     3                                     	|
|     21                	|     8                                     	|
|     22                	|     2                                     	|
|     23                	|     3                                     	|
|     24                	|     9                                     	|
|     25                	|     11                                    	|
|     26                	|     17                                    	|
|     27                	|     7                                     	|
|     28                	|     4                                     	|
|     29                	|     14                                    	|
|     30                	|     14                                    	|
|     31                	|     9                                     	|
|     32                	|     8                                     	|
|     33                	|     1                                     	|
|     34                	|     3                                     	|

### Stress condition


