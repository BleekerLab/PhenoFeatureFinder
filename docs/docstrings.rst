*************
Documentation
*************


PhenotypeAnalysis
#################

A class to analyse data from developmental bioassays and group the samples in distinct phenotypic classes. 


Parameters
----------

bioassay_csv: string
    A path to a .csv file with the bioassay count data.


Attributes
----------

bioassay_csv
    The path to the csv file with the input data.
bioassay
    A pandas dataframe with the development data.
sample_id
    The name of the column that contains the sample identifiers.
cumulative_data
    A pandas dataframe with cumulative development data.
max_counts
    A pandas dataframe with the maximum number of nymphs developed to or past each developmental stage per group.
survival_data
    A pandas dataframe with the number of living nymphs per stage at each time point per group.
order_of_groups
    A string with the order in which the groups should be plotted.
absolute_counts
    A pandas dataframe with values from max_counts in long format.
max_relative
    A pandas dataframe based on absolute_counts, with values relative to values of specified developmental stage.


Methods
-------

reshape_to_wide
    Reshapes the dataframe from a long to a wide format to make the data accessible for pre-processing. With the counts of each developmental stage in a seperate columns.
combine_seperately_counted_versions_of_last_recorded_stage
    Calculates the total number of nymphs developed to the final developmental stage per sample on each timepoint.
convert_counts_to_cumulative
    Calculates the total number of nymphs developed to or past each stage on each timepoint.
correct_cumulative_counts
    Inner function for convert_counts_to_cumulative().
create_df_with_max_counts_per_stage
    Inner function for convert_counts_to_cumulative(). 
prepare_for_plotting
    Prepare the order in which the groups should be plotted.
plot_counts_per_stage
    Plots the counts per nymphal stage in boxplots.
plot_development_over_time_in_fitted_model
    Fits a 3 parameter log-logistic curve to the development over time to a specified stage.
ll3
    A three parameter log-logistic function.
plot_survival_over_time_in_fitted_model
    Fits a 3 parameter log-normal curve to the number of living nymphs over time.
hazard
    A three parameter log-normal function.


reshape_to_wide
***************

Reshapes the dataframe from a long to a wide format to make the data accessible for pre-processing.
with the counts of each developmental stage in a seperate columns.


Usage
-----

reshape_to_wide(
    self,
    sample_id='sample_id',
    grouping_variable='genotype',
    developmental_stages='stage',
    count_values='number',
    time='day')


Parameters
----------

sample_id: string, default='sample_id'
    The name of the column that contains the sample identifiers.
grouping_variable: string, default='genotype'
    The name of the column that contains the names of the grouping variables.
    Examples are genotypes or treatments
developmental_stages: string, default='stage'
    The name of the column that contains the developmental stages that were scored during the bioassay.
count_values: string, default='numbers'
    The name of the column that contains the counts.
time: string, default='day'
    The name of the column that contains the time at which bioassay scoring was performed.
    Examples are the date or the number of days after infection.

Examples
--------

Example of an input dataframe

+-----------+-----------+-------+---------------+-------+
| sample_id | genotype  | day   | stage         | count |
+===========+===========+=======+===============+=======+
| mm_1      |   mm      | 5     | eggs          | 45    |
+-----------+-----------+-------+---------------+-------+
| mm_1      |   mm      | 5     | first_instar  | 0     |
+-----------+-----------+-------+---------------+-------+
| mm_1      |   mm      | 5     | second_instar | 0     |
+-----------+-----------+-------+---------------+-------+


Example of a reshaped output dataframe

+-----------+-----------+-------+-------+---------------+---------------+---------------+
| sample_id | genotype  | day   | eggs  | first_instar  | second instar | third_instar  |
+===========+===========+=======+=======+===============+===============+===============+
| mm_1      |   mm      | 5     | 45    | 0             | 0             | 0             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+
| mm_1      |   mm      | 9     | NA    | 10            | 5             | 0             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+
| mm_1      |   mm      | 11    | NA    | 15            | 17            | 4             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+


combine_seperately_counted_versions_of_last_recorded_stage
**********************************************************

Calculates the total number of nymphs developed to the final developmental stage per sample on each timepoint.
This is used when nymphs in the (late) final nymph stage were removed after each counting moment and/or
when exuviea and last instar stage nymphs were counted seperately.
Removal of late last stage nymphs could for example be used to prevent adults from emerging and escaping.


Usage
-----

combine_seperately_counted_versions_of_last_recorded_stage(
    self,
    exuviea='exuviea',
    late_last_stage='late_fourth_instar',
    early_last_stage='early_fourth_instar',
    new_last_stage='fourth_instar',
    seperate_exuviea=True,
    late_last_stage_removed=True,
    early_last_stage_kept=True,
    remove_individual_stage_columns=True)


Parameters
----------

exuviea: string, default='exuviea'
    The name of the column that contains the exuviea counts. 
late_last_stage: string, default='late_fourth_instar'
    The name of the column that contains the counts of the last developmental stage recorded in the bioassay.
early_last_stage: string, default='early_fourth_instar'
    The name of the column that contains the counts of the nymphs in early last developmental stage.
    Is used when nymphs counted in late_last_stage were removed after each counting moment during the bioassay.
new_last_stage: string, default='fourth_instar'
    Name for new column with the returned total final stage data
seperate_exuviea: boolean, default=True
    If True, sums exuviea and late_last_stage per sample per timepoint.
    If exuviea were counted seperately from late_last_stage, set to True.
    If exuviea count was included in late_last_stage, set to False
late_last_stage_removed: boolean, default=True
    If True, returns the cumulative number of late_last_stage(+exuviea) per sample over time.
    If nymphs counted in late_last_stage (and exuviea if counted seperately) were removed after each counting 
    moment, set to True.
    If nymphs counted in late_last_stage (and exuviea if counted seperately) were left on the sample until
    ending the bioassay, set to False.
early_last_stage_kept: boolean, default=True
    If True, sums the early and late last stage counts per sample per timepoint
    If late last stage nymphs were removed after each counting moment, but early last stage nymphs were left on
    sample, set to True.
    If early and late last stage nymphs were not counted seperately, set to False
remove_individual_stage_columns: boolean, default=True
    If True, removes exuviea, late_last_stage, early_last_stage columns from dataframe after returning 
    new_last_stage column.


Examples
--------

Example of an input dataframe

+-----------+-----------+-------+-------+-----+---------------+-----------+---------------------+--------------------+
| sample_id | genotype  | day   | eggs  | ... | third_instar  | exuviea   | early_fourth_instar | late_fourth_instar |
+===========+===========+=======+=======+=====+===============+===========+=====================+====================+
| mm_1      |   mm      | 5     | 45    | ... | 0             | 0         | 0                   | 0                  |
+-----------+-----------+-------+-------+-----+---------------+-----------+---------------------+--------------------+
| mm_1      |   mm      | 9     | NA    | ... | 0             | 1         | 5                   | 0                  |
+-----------+-----------+-------+-------+-----+---------------+-----------+---------------------+--------------------+
| mm_1      |   mm      | 11    | NA    | ... | 4             | 0         | 7                   | 4                  |
+-----------+-----------+-------+-------+-----+---------------+-----------+---------------------+--------------------+


Example of an output dataframe

+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| sample_id | genotype  | day   | eggs  | first_instar  | second instar | third_instar  | fourth_instar |
+===========+===========+=======+=======+===============+===============+===============+===============+
| mm_1      |   mm      | 5     | 45    | 0             | 0             | 0             | 0             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| mm_1      |   mm      | 9     | NA    | 10            | 5             | 0             | 6             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| mm_1      |   mm      | 11    | NA    | 15            | 17            | 4             | 12            |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+


convert_counts_to_cumulative
****************************

Calculates the total number of nymphs developed to or past each stage on each timepoint.
Cumulative counts make the analysis of development over time and the comparison of number of nymphs past a stage easier.
If nymphs in the (late) final nymph stage were removed after each counting moment and/or
when exuviea and/or early and late last instar stage nymphs were counted seperately, 
total_last_stage() should be used first.

Usage
-----

def convert_counts_to_cumulative(
    self,
    n_developmental_stages=4,
    sample_id='sample_id',
    eggs='eggs',
    first_stage='first_instar',
    second_stage='second_instar',
    third_stage='third_instar',
    fourth_stage='fourth_instar',
    fifth_stage='fifth_instar',
    sixth_stage='sixth_instar')


Parameters
----------

n_developmental_stages: integer, default=4
    The number of developmental stages which were recorded seperately. 
    Can range from 2 to 6.
sample_id: string, default='sample_id'
    The name of the column that contains the sample identifiers.
eggs: string, default='eggs'
    The name of the column that contains the counts of the eggs.
first_stage: string, default='first_instar'
    The name of the column that contains the counts of the first developmental stage recorded in the bioassay.
second_stage: string, default='second_instar'
    The name of the column that contains the counts of the second developmental stage recorded in the bioassay.
third_stage: string, default='third_instar'
    The name of the column that contains the counts of the third developmental stage recorded in the bioassay.
fourth_stage: string, default='fourth_instar'
    The name of the column that contains the counts of the fourth developmental stage recorded in the bioassay.
fifth_stage: string, default='fifth_instar'
    The name of the column that contains the counts of the fifth developmental stage recorded in the bioassay.
sixth_stage: string, default='sixth_instar'
    The name of the column that contains the counts of the sixth developmental stage recorded in the bioassay. 


correct_cumulative_counts
*************************

Inner function for convert_counts_to_cumulative(). If nymphs die during the bioassay, 
they should be included in the cumulative count for the stages it had passed. 
Otherwise, the cumulative count could go down over time. This function corrects the cumulative
count if it is lower than the previous count.

Usage
-----

correct_cumulative_counts(
    self, 
    current_stage,
    grouping_variable)


create_df_with_max_counts_per_stage
***********************************

Inner function for convert_counts_to_cumulative(). 
With the maximum number of nymphs developed to or past each developmental stage per plant, 
making graphs becomes easier.

Usage
-----

create_df_with_max_counts_per_stage(
    self, 
    egg_column,
    last_stage,
    grouping_variable):


prepare_for_plotting
********************

Prepare the order in which the groups should be plotted.


Usage
-----

prepare_for_plotting(
self,
order_of_groups)


Parameters
----------

order_of_groups: string
    List of the group names in the prefered order for plotting
    For example: ['MM', 'LA', 'PI']


plot_counts_per_stage
*********************

Plots the counts per nymphal stage in boxplots. The nymph counts are given as the absolute number of nymphs that 
developed to or past each stage at the last timepoint and as a fraction of nymphs that developed to or past each 
stage at the last timepoint relative to another developmental stage. The other developmental stage to which the 
data is made relative defaults to the first instar stage, because this represents the number of hatched eggs. This
means that in this case only the succes of the development is compared between groups (e.g. genotypes or 
treatments) and the hatching rate of the eggs is not taken into acount.

The imput dataframe 'max_counts' is created with convert_counts_to_cumulative.


Usage
-----

plot_counts_per_stage(
    self,
    grouping_variable='genotype',
    sample_id='sample_id',
    eggs='eggs',
    first_stage='first_instar',
    second_stage='second_instar',
    third_stage='third_instar',
    fourth_stage='fourth_instar',
    absolute_x_axis_label='genotype',
    absolute_y_axis_label='counts (absolute)',
    relative_x_axis_label='genotype',
    relative_y_axis_label='relative number of nymphs',
    make_nymphs_relative_to='first_instar')


Parameters
----------

grouping_variable: string, default='genotype'
    The name of the column that contains the names of the grouping variables.
    Examples are genotypes or treatments
sample_id: string, default='sample_id'
    The name of the column that contains the sample identifiers.
eggs: string, default='eggs'
    The name of the column that contains the counts of the eggs.
first_stage: string, default='first_instar'
    The name of the column that contains the counts of the first developmental stage recorded in the bioassay.
second_stage: string, default='second_instar'
    The name of the column that contains the counts of the second developmental stage recorded in the bioassay.
third_stage: string, default='third_instar'
    The name of the column that contains the counts of the third developmental stage recorded in the bioassay.
fourth_stage: string, default='fourth_instar'
    The name of the column that contains the counts of the fourth developmental stage recorded in the bioassay.
absolute_x_axis_label: string, default='genotype'
    Label for the x-axis of the boxplots with count data.
absolute_y_axis_label: string, default='counts (absolute)'
    Label for the y-axis of the boxplots with count data.
relative_x_axis_label: string, default='genotype'
    Label for the x-axis of the boxplots with relative development.
relative_y_axis_label: string, default='relative number of nymphs'
    Label for the y-axis of the boxplots with relative development.
make_nymphs_relative_to: string, default='first_instar'
    The name of the column that contains the counts of the developmental stage which should be used to calculate 
    the relative development to all developmental stages.


Examples
--------

Example of an input dataframe

+-----------+-----------+-------+-------+---------------+---------------+--------------+---------------+
| sample_id | genotype  | day   | eggs  | first_instar  | second_instar | third_instar | fourth_instar |
+===========+===========+=======+=======+===============+===============+==============+===============+
| mm_1      |   mm      | 28    | 45    | 34            | 30            | 30           | 29            |
+-----------+-----------+-------+-------+---------------+---------------+--------------+---------------+
| mm_2      |   mm      | 28    | 50    | 39            | 33            | 28           | 26            |
+-----------+-----------+-------+-------+---------------+---------------+--------------+---------------+
| LA_1      |   LA      | 28    | 42    | 30            | 25            | 17           | 4             |
+-----------+-----------+-------+-------+---------------+---------------+--------------+---------------+


plot_development_over_time_in_fitted_model
******************************************

Fits a 3 parameter log-logistic curve to the development over time to a specified stage. The fitted curve and the
observed datapoints are plotted and returned with the model parameters. 
The reduced Chi-squared is provided to asses the goodness of fit for the fitted models for each group (genotype, 
treatment, etc.). Optimaly, the reduced Chi-squared should approach the number of observation points per sample. A
much larger reduced Chi-squared indicates a bad fit. A much smaller reduced Chi-squared indicates overfitting of 
the model.


Usage
-----

plot_development_over_time_in_fitted_model(
    self, 
    grouping_variable='genotype',
    sample_id='sample_id',
    time='day',
    x_axis_label='days after infection',
    y_axis_label='development to 4th instar stage (relative to 1st instars)',
    stage_of_interest='fourth_instar',
    use_relative_data=True,
    make_nymphs_relative_to='first_instar',
    predict_for_n_days=0)


Parameters
----------

grouping_variable: string, default='genotype'
    The name of the column that contains the names of the grouping variables.
    Examples are genotypes or treatments
sample_id: string, default='sample_id'
    The name of the column that contains the sample identifiers.
time: string, default='day'
    The name of the column that contains the time at which bioassay scoring was performed.
    Examples are the date or the number of days after infection.
x_axis_label: string, default='days after infection'
    Label for the x-axis
y_axis_label: string, default='development to 4th instar stage (relative to 1st instars)'
    Label for the y-axis
stage_of_interest: string, default='fourth_instar'
    The name of the column that contains the data of the developmental stage of interest.
use_relative_data: boolean, default=True
    If True, the counts for the stage of interest are devided by the stage indicated at 'make_nymphs_relative_to'.
    The returned relative rate is used for plotting and curve fitting.
make_nymphs_relative_to: string, default='first_instar'
    The name of the column that contains the counts of the developmental stage which should be used to calculate 
    therelative development to all developmental stages.
predict_for_n_days: default=o
    Continue model for n days after final count.


Examples
--------

Example of an input dataframe

+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| sample_id | genotype  | day   | eggs  | first_instar  | second instar | third_instar  | fourth_instar |
+===========+===========+=======+=======+===============+===============+===============+===============+
| mm_1      |   mm      | 5     | 45    | 15            | 7             | 0             | 0             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| mm_1      |   mm      | 9     | NA    | 24            | 14            | 6             | 3             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| mm_1      |   mm      | 11    | NA    | 38            | 27            | 16            | 12            |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+


ll3
***

A three parameter log-logistic function.

Usage
-----

ll3(x,slope,maximum,emt50):


Parameters
----------

slope: 
    the slope of the curve
maximum: 
    the maximum value of the curve
emt50: 
    the EmT50, the timepoint at which 50% of nymphs has developed to the stage of interest


Model
-----

y(x) = maximum/(1+np.exp(slope*(np.log(x)-np.log(emt50))))


plot_survival_over_time_in_fitted_model
***************************************

Fits a 3 parameter log-normal curve to the number of living nymphs over time. The fitted curve and the
observed datapoints are plotted and returned with the model parameters. 
The reduced Chi-squared is provided to asses the goodness of fit for the fitted models for each group (genotype, 
treatment, etc.). Optimaly, the reduced Chi-squared should approach the number of observation points per sample. A
much larger reduced Chi-squared indicates a bad fit. A much smaller reduced Chi-squared indicates overfitting of 
the model.


Usage
-----

plot_survival_over_time_in_fitted_model(
    self,
    grouping_variable='genotype',
    sample_id='sample_id',
    time='day',
    x_axis_label='days after infection',
    y_axis_label='number of nymphs per plant',
    stage_of_interest='first_instar',
    use_relative_data=False,
    make_nymphs_relative_to='eggs',
    predict_for_n_days=0)


Parameters
----------

grouping_variable: string, default='genotype'
    The name of the column that contains the names of the grouping variables.
    Examples are genotypes or treatments
sample_id: string, default='sample_id'
    The name of the column that contains the sample identifiers.
time: string, default='day'
    The name of the column that contains the time at which bioassay scoring was performed.
    Examples are the date or the number of days after infection.
x_axis_label: string, default='days after infection'
    Label for the x-axis
y_axis_label: string, default='development to 4th instar stage (relative to 1st instars)'
    Label for the y-axis
stage_of_interest: string, default='first_instar'
    The name of the column that contains the data of the developmental stage of interest.
use_relative_data: boolean, default=False
    If True, the counts for the stage of interest are devided by the stage indicated at 'make_nymphs_relative_to'.
    The returned relative rate is used for plotting and curve fitting.
make_nymphs_relative_to: string, default='eggs'
    The name of the column that contains the counts of the developmental stage which should be used to calculate 
    the relative development to all developmental stages.
predict_for_n_days: default=o
    Continue model for n days after final count.


Examples
--------
Example of an input dataframe

+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| sample_id | genotype  | day   | eggs  | first_instar  | second instar | third_instar  | fourth_instar |
+===========+===========+=======+=======+===============+===============+===============+===============+
| mm_1      |   mm      | 5     | 45    | 15            | 7             | 0             | 0             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| mm_1      |   mm      | 9     | NA    | 24            | 14            | 6             | 3             |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+
| mm_1      |   mm      | 11    | NA    | 38            | 27            | 16            | 12            |
+-----------+-----------+-------+-------+---------------+---------------+---------------+---------------+


hazard
******

A three parameter log-normal function.


Usage
-----

hazard(x,auc,median,shape)


Parameters
----------

auc: 
    area under the curve
median: 
    median time point
shape: 
    shape of the curve


model
-----

y(x) = (auc*(shape/median)*pow(x/median,shape-1))/(1+pow(x/median,shape))
   

OmicsAnalysis
#############

A class to streamline the filtering and exploration of a metabolome dataset.   


Parameters
----------

metabolome_csv: str
    A path to a .csv file with the metabolome data (scaled or unscaled).
    Shape of the dataframe is usually (n_samples, n_features) with n_features >> n_samples

metabolome_feature_id_col: str, optional
    The name of the column that contains the feature identifiers (default is 'feature_id').
    Feature identifiers should be unique (=not duplicated).


Attributes
----------

metabolome: `pandas.core.frame.DataFrame`, (n_samples, n_features)
    The metabolome Pandas dataframe imported from the .csv file. 
    metabolome_validated: `bool`
    Is the metabolome dataset validated?
    Default is False.
blank_features_filtered: `bool`
    Are the features present in blank samples filtered out from the metabolome data?
    Default by False.
filtered_by_percentile_value: bool
    Are the features filtered by percentile value?
unreliable_features_filtered: `bool`
    Are the features not reliably present within one group filtered out from the metabolome data?
pca_performed: `bool`
    Has PCA been performed on the metabolome data?
    Default is False. 
exp_variance: `pandas.core.frame.DataFrame`, (n_pc, 1)
    A Pandas dataframe with explained variance per Principal Component.
    The index of the df contains the PC index (PC1, PC2, etc.).
    The second column contains the percentage of the explained variance per PC.
metabolome_pca_reduced: `numpy.ndarray`, (n_samples, n_pc)
    Numpy array with sample coordinates in reduced dimensions.
    The dimension of the numpy array is the minimum of the number of samples and features. 
sparsity: float
    Metabolome matrix sparsity.


Methods
-------

validate_input_metabolome_df
    Check if the provided metabolome file is suitable. Turns attribute metabolome_validated to True. 
discard_features_detected_in_blanks
    Removes features only detected in blank samples. 
impute_missing_values_with_median
    Impute missing values with the median value of the feature.
filter_out_unreliable_features
    Filter out features not reliably detectable in replicates of the same grouping factor. 
    For instance, if a feature is detected less than 4 times within 4 biological replicates, it is discarded with argument nb_times_detected=4.  
filter_features_per_group_by_percentile
    Filter out features whose abundance within the same grouping factor is lower than a certain percentile value.
    For instance, features lower than the 90th percentile within a single group are discarded with argument percentile=90. 
compute_metabolome_sparsity
    Computes the sparsity percentage of the metabolome matrix (percentage of 0 values e.g. 100% for an matrix full of 0 values)
write_clean_metabolome_to_csv
    Write the filtered and analysis-ready metabolome data to a .csv file.  


Examples
-----

Example of an input metabolome input format (from a csv file)

+----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
| feature_id           | blank_1 | blank_2 | blank_3 | blank_4 | MM_1  | MM_2  | MM_3  | MM_4  | LA1330_1 | LA1330_2 | LA1330_3 | LA1330_4 |
+======================+=========+=========+=========+=========+=======+=======+=======+=======+==========+==========+==========+==========+
| rt-0.04_mz-241.88396 | 280     | 694     | 502     | 604     | 554   | 678   | 674   | 936   | 824      | 940      | 794      | 828      |
+----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
| rt-0.05_mz-143.95911 | 1036    | 1566    | 1326    | 1490    | 1364  | 1340  | 1692  | 1948  | 1928     | 1956     | 1730     | 1568     |
+----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
| rt-0.06_mz-124.96631 | 1308    | 992     | 1060    | 1010    | 742   | 990   | 0     | 888   | 786      | 668      | 762      | 974      |
+----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
| rt-0.08_mz-553.45905 | 11340   | 12260   | 10962   | 11864   | 10972 | 11190 | 12172 | 11820 | 12026    | 11604    | 11122    | 11260    |
+----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
| rt-0.08_mz-413.26631 | 984     | 1162    | 1292    | 1104    | 1090  | 1106  | 1290  | 1170  | 1282     | 924      | 1172     | 1062     |
+----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+


validate_input_metabolome_df
****************************

Validates the dataframe containing the feature identifiers, metabolite values and sample names.
Will place the 'feature_id_col' column as the index of the validated dataframe. 
The validated metabolome dataframe is stored as the 'validated_metabolome' attribute. 

Usage
-----

validate_input_metabolome_df(
    self, 
    metabolome_feature_id_col='feature_id')


Parameters
----------

metabolome_feature_id: str, optional 
    The name of the column that contains the feature identifiers (default is 'feature_id').
    Feature identifiers should be unique (=not duplicated).
    

Returns
-------

self: object
    Object with attribute metabolome_validated set to True if tests are passed. 


Examples
-----

Example of a valid input metabolome dataframe

+-------------+----------------+----------------+----------------+----------------+
| feature_id  | genotypeA_rep1 | genotypeA_rep2 | genotypeA_rep3 | genotypeA_rep4 |
+=============+================+================+================+================+
| metabolite1 |   1246         | 1245           | 12345          | 12458          |
+-------------+----------------+----------------+----------------+----------------+
| metabolite2 |   0            | 0              | 0              | 0              |
+-------------+----------------+----------------+----------------+----------------+
| metabolite3 |   10           | 0              | 0              | 154            |
+-------------+----------------+----------------+----------------+----------------+


discard_features_detected_in_blanks
***********************************

Removes features present in blanks.
Steps:

#. Sum the abundance of each feature in the blank samples.
#. Makes a list of features to be discarded (features with a positive summed abundance).
#. Returns a filtered Pandas dataframe with only features not detected in blank samples


Usage
-----

discard_features_detected_in_blanks(
    self, 
    blank_sample_contains='blank')


Parameters
----------

blank_sample_contains: str, optional.
    Column names with this name will be considered blank samples.
    Default is='blank'

Returns
-------

metabolome: pandas.core.frame.DataFrame
    A filtered Pandas dataframe without features detected in blank samples and with the blank samples removed. 


create_density_plot
*******************

For each grouping variable (e.g. genotype), creates a histogram and density plot of all feature peak areas.
This plot helps to see whether some groups have a value distribution different from the rest. 
The percentage is indicated on the y-axis (bar heights sum to 100).

Usage
-----

create_density_plot(
    self, 
    name_grouping_var="genotype", 
    n_cols=3, 
    nbins=1000)


Parameters
----------

name_grouping_var: str, optional
    The name used when splitting between replicate and main factor.
    For example "genotype" when splitting MM_rep1 into 'MM' and 'rep1'.
    Default is 'genotype'. 
n_cols: int, optional
    The number of columns for the final plot.
nbins: int, optional
    The number of bins to create. 

Returns
-------

matplotlib Axes
    Returns the Axes object with the density plots drawn onto it.


filter_features_per_group_by_percentile
***************************************

Filter metabolome dataframe based on a selected percentile threshold.
Features with a peak area values lower than the selected percentile will be discarded. 
The percentile value is calculated per grouping variable. 

For instance, selecting the 50th percentile (median) will discard 50% of the features with a peak area
lower than the median/50th percentile in each group. 


Usage
-----

filter_features_per_group_by_percentile(
self, 
name_grouping_var="genotype",
separator_replicates="_",
percentile=50)


Parameters
----------

name_grouping_var: str, optional
    The name of the grouping variable (default is "genotype")
separator_replicates: str, optional
    The character used to separate the main grouping variable from biological replicates. 
    Default is "_: (underscore)
percentile: float, optional
    The percentile threshold. Has to be comprised 0 and 100.


Returns
-------

self: object
    The object with the .metabolome attribute filtered and the filtered_by_percentile_value set to True. 


.. seealso:: Use create_density_plot() method to decide on a suitable percentile value. 


filter_out_unreliable_features
******************************

Removes features not reliably detectable in multiple biological replicates from the same grouping factor. 

Takes a dataframe with feature identifiers in index and samples as columns.
Steps:

#. First melt and split the sample names to generate the grouping variable
#. Count number of times a metabolite is detected in the groups. If number of times detected in a group = number of biological replicates then it is considered as reliable. Each feature receives a tag  'reliable' or 'not_reliable'
#. Discard the 'not_reliable' features and keep the filtered dataframe. 


Usage
-----

filter_out_unreliable_features(
    self,
    name_grouping_var="genotype", 
    nb_times_detected=4,
    separator_replicates='_')


Parameters
------

name_grouping_var: str, optional
    The name used when splitting between replicate and main factor.
    For example "genotype" when splitting MM_rep1 into 'MM' and 'rep1'.
    Default is 'genotype'. 
nb_times_detected: int, optionaldefault=4
    Number of times a metabolite should be detected to be considered 'reliable'. 
    Should be equal to the number of biological replicates for a given group of interest (e.g. genotype)
separator_replicates: string, default="_"
    The separator to split sample names into a grouping variable (e.g. genotype) and the biological replicate number (e.g. 1)


Returns
-------

metabolome: ndarray
    A Pandas dataframe with only features considered as reliable, sample names and their values. 

Examples
-----
Example of an input dataframe

+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| feature_id            | MM_1  	| MM_2  	| MM_3  	| MM_4  	| LA1330_1 	| LA1330_2 	|
+=======================+===========+===========+===========+===========+===========+===========+
| rt-0.04_mz-241.88396 	| 554   	| 678   	| 674   	| 936   	| 824      	| 940      	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| rt-0.05_mz-143.95911 	| 1364  	| 1340  	| 1692  	| 1948  	| 1928     	| 1956     	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| rt-0.06_mz-124.96631 	| 0      	| 0     	| 0     	| 888   	| 786      	| 668      	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| rt-0.08_mz-553.45905 	| 10972 	| 11190 	| 12172 	| 11820 	| 12026    	| 11604    	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+


Example of an output df (rt-0.06_mz-124.96631 is kicked out because 3x0 and 1x888 in MM groups)

+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| feature_id            | MM_1  	| MM_2  	| MM_3  	| MM_4  	| LA1330_1 	| LA1330_2 	|
+=======================+===========+===========+===========+===========+===========+===========+
| rt-0.04_mz-241.88396 	| 554   	| 678   	| 674   	| 936   	| 824      	| 940      	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| rt-0.05_mz-143.95911 	| 1364  	| 1340  	| 1692  	| 1948  	| 1928     	| 1956     	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| rt-0.08_mz-553.45905 	| 10972 	| 11190 	| 12172 	| 11820 	| 12026    	| 11604    	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+


write_clean_metabolome_to_csv
*****************************
 
Writes the cleaned metabolome data to the disk as a comma-separated value file.


Usage
-----

write_clean_metabolome_to_csv(self, path_of_cleaned_csv="./data_for_manuals/filtered_metabolome.csv"):


Parameters
----------

path_of_cleaned_csv: str, optional
    The path and filename of the .csv file to save.
    Default to "./data_for_manuals/filtered_metabolome.csv" 


compute_pca_on_metabolites
**************************

Performs a Principal Component Analysis (PCA) on the metabolome data. 

The PCA analysis will return transformed coordinates of the samples in a new space. 
It will also give the percentage of variance explained by each Principal Component. 
Assumes that number of samples < number of features/metabolites
Performs a transpose of the metabolite dataframe if n_samples > n_features (this can be turned off with auto_transpose)


Usage
-----

compute_pca_on_metabolites(
    self, 
    scale=True, 
    n_principal_components=10, 
    auto_transpose=True)


Parameters
----------

scale: `bool`, optional
    Perform scaling (standardize) the metabolite values to zero mean and unit variance. 
    Default is True. 
n_principal_components: int, optional
    number of principal components to keep in the PCA analysis.
    if number of PCs > min(n_samples, n_features) then set to the minimum of (n_samples, n_features)
    Default is to calculate 10 components.
auto_transpose: `bool`, optional. 
    If n_samples > n_features, performs a transpose of the feature matrix.
    Default is True (meaning that transposing will occur if n_samples > n_features).

Returns
-------

self: object
    Object with .exp_variance: dataframe with explained variance per Principal Component
    .metabolome_pca_reduced: dataframe with samples in reduced dimensions
    .pca_performed: `bool`ean set to True


impute_missing_values_with_median
*********************************

Imputes missing values with the median of the column.
This is necessary for PCA to work.


Usage
-----

impute_missing_values_with_median(
    self, 
    missing_value_str='np.nan')


Params
------

missing_value_str: str, optional
    The string that represents missing values in the input dataframe.
    All occurrences of missing_values will be imputed. 
    For pandas’ dataframes with nullable integer dtypes with missing values, missing_values can be set to either np.nan or pd.NA.


Returns
-------

self: object with attribute 'metabolome' updated with imputed values.


create_scree_plot
*****************

Returns a barplot with the explained variance per Principal Component. 
Has to be preceded by perform_pca()


Usage
-----

create_scree_plot(
    self, 
    plot_file_name=None)


Parameters
---------

plot_file_name: string, default='None'
    Path to a file where the plot will be saved.
    For instance 'my_scree_plot.pdf'


Returns
-------

matplotlib Axes
    Returns the Axes object with the scree plot drawn onto it.
    Optionally a saved image of the plot. 


create_sample_score_plot
************************

Returns a sample score plot of the samples on PCx vs PCy. 
Samples are colored based on the grouping variable (e.g. genotype)


Usage
-----

create_sample_score_plot(
    self, 
    pc_x_axis=1, 
    pc_y_axis=2, 
    name_grouping_var='genotype',
    separator_replicates="_",
    show_color_legend=True,
    plot_file_name=None)


Parameters
----------

pc_x_axis: int, optional 
    Principal Component to plot on the x-axis (default is 1 so PC1 will be plotted).
pc_y_axis: int, optional.
    Principal Component to plot on the y-axis (default is 2 so PC2 will be plotted).
name_grouping_var: str, optional
    Name of the variable used to color samples (Default is "genotype"). 
separator_replicates: str, optional.
    String separator that separates grouping factor from biological replicates (default is underscore "_").
show_color_legend: bool, optional.
    Add legend for hue (default is True).
plot_file_name: str, optional 
    A file name and its path to save the sample score plot (default is None).
    For instance "mydir/sample_score_plot.pdf"
    Path is relative to current working directory.


Returns
-------

matplotlib Axes
    Returns the Axes object with the sample score plot drawn onto it.
    Samples are colored by specified grouping variable. 
    Optionally a saved image of the plot. 


compute_metabolome_sparsity
***************************

Determine the sparsity of the metabolome matrix. 
Formula: number of non zero values/number of values * 100
The higher the sparsity, the more zero values 


Usage
-----

compute_metabolome_sparsity(self)


Returns
-------

self: object
    Object with sparsity attribute filled (sparsity is a float).


References
----------

`<https://stackoverflow.com/questions/38708621/how-to-calculate-percentage-of-sparsity-for-a-numpy-array-matrix>`_


plot_features_in_upset_plot
***************************

Visuallises the presence of features per group in an UpSet plot. 
A feature is considered present in a group if the median>0.


Usage
-----

plot_features_in_upset_plot(
    self,
    seperator_replicates="_",
    plot_file_name=None)


Params
------

separator_replicates: string, default="_"
    The separator to split sample names into a grouping variable (e.g. genotype) and the biological replicate number (e.g. 1)
plot_file_name: str, optional 
    A file name and its path to save the sample score plot (default is None).
    For instance "mydir/feature_upset_plot.pdf"
    Path is relative to current working directory.


Returns
-------

Plot:
    UpSet plot with features presence per group.


Examples
-----

Example of an input dataframe

+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| feature_id            | MM_1  	| MM_2  	| MM_3  	| MM_4  	| LA1330_1 	| LA1330_2 	|
+=======================+===========+===========+===========+===========+===========+===========+
| rt-0.04_mz-241.88396 	| 554   	| 678   	| 674   	| 936   	| 824      	| 940      	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| rt-0.05_mz-143.95911 	| 1364  	| 1340  	| 1692  	| 1948  	| 1928     	| 1956     	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+
| rt-0.08_mz-553.45905 	| 10972 	| 11190 	| 12172 	| 11820 	| 12026    	| 11604    	|
+-----------------------+-----------+-----------+-----------+-----------+-----------+-----------+




FeatureSelection
################

A class to perform metabolite feature selection using phenotyping and metabolic data. 

- Perform sanity checks on input dataframes (values above 0, etc.).
- Get a baseline performance of a simple Machine Learning Random Forest ("baseline").
- Perform automated Machine Learning model selection using autosklearn.
    Using metabolite data, train a model to predict phenotypes.
    Yields performance metrics (balanced accuracy, precision, recall) on the selected model.
- Extracts performance metrics from the best ML model. 
- Extracts the best metabolite features based on their feature importance and make plots per sample group. 


Parameters
----------

metabolome_csv: string
    A path to a .csv file with the cleaned up metabolome data (unreliable features filtered out etc.)
    Use the MetabolomeAnalysis class methods. 
    Shape of the dataframe is usually (n_samples, n_features) with n_features >> n_samples
phenotype_csv: string
    A path to a .csv file with the phenotyping data. 
    Should be two columns at least with: 
        - column 1 containing the sample identifiers
        - column 2 containing the phenotypic class e.g. 'resistant' or 'sensitive'
metabolome_feature_id_col: string, default='feature_id'
    The name of the column that contains the feature identifiers.
    Feature identifiers should be unique (=not duplicated).
phenotype_sample_id: string, default='sample_id'
    The name of the column that contains the sample identifiers.
    Sample identifiers should be unique (=not duplicated).


Attributes
----------

metabolome_validated: bool
    Is the metabolome file valid for Machine Learning? (default is False)   

phenotype_validated: bool
    Is the phenotype file valid for Machine Learning? (default is False)

baseline_performance: float 
    The baseline performance computed with get_baseline_performance() i.e. using a simple Random Forest model. 
    Search for the best ML model using search_best_model() should perform better than this baseline performance. 

best_ensemble_models_searched: bool
    Is the search for best ensemble model using auto-sklearn already performed? (default is False)

metabolome: pandas.core.frame.DataFrame
    The validated metabolome dataframe of shape (n_features, n_samples).

phenotype: pandas.core.frame.DataFrame
    A validated phenotype dataframe of shape (n_samples, 1)
    Sample names in the index and one column named 'phenotype' with the sample classes.

baseline_performance: str
    Average balanced accuracy score (-/+ standard deviation) of the basic Random Forest model. 

best_model: sklearn.pipeline.Pipeline
    A scikit-learn pipeline that contains one or more steps.
    It is the best performing pipeline found by TPOT automated ML search.

pc_importances: pandas.core.frame.DataFrame
    A Pandas dataframe that contains Principal Components importances using scikit-learn permutation_importance()
    Mean of PC importance over n_repeats.
    Standard deviation over n_repeats.
    Raw permutation importance scores.

feature_loadings: pandas.core.frame.DataFrame
    A Pandas dataframe that contains feature loadings related to Principal Components


Methods
--------

validate_input_metabolome_df()
    Validates the dataframe read from the 'metabolome_csv' input file.

validate_input_phenotype_df()
    Validates the phenotype dataframe read from the 'phenotype_csv' input file.

get_baseline_performance()
    Fits a basic Random Forest model to get default performance metrics. 

search_best_model_with_tpot_and_get_feature_importances()
    Search for the best ML pipeline using TPOT genetic programming method.
    Computes and output performance metrics from the best pipeline.
    Extracts feature importances using scikit-learn permutation_importance() method. 



Examples
--------

Example of an input metabolome .csv file

+-------------+----------------+----------------+----------------+----------------+
| feature_id  | genotypeA_rep1 | genotypeA_rep2 | genotypeA_rep3 | genotypeA_rep4 |
+=============+================+================+================+================+
| metabolite1 |   1246         | 1245           | 12345          | 12458          |
+-------------+----------------+----------------+----------------+----------------+
| metabolite2 |   0            | 0              | 0              | 0              |
+-------------+----------------+----------------+----------------+----------------+
| metabolite3 |   10           | 0              | 0              | 154            |
+-------------+----------------+----------------+----------------+----------------+


Example of an input phenotype .csv file

+----------------+-----------+
| sample_id      | phenotype | 
+================+===========+
| genotypeA_rep1 | sensitive | 
+----------------+-----------+
| genotypeA_rep2 | sensitive |
+----------------+-----------+   
| genotypeA_rep3 | sensitive |
+----------------+-----------+
| genotypeA_rep4 | sensitive |
+----------------+-----------+ 
| genotypeB_rep1 | resistant |
+----------------+-----------+   
| genotypeB_rep2 | resistant |
+----------------+-----------+


validate_input_metabolome_df
****************************

Validates the dataframe containing the feature identifiers, metabolite values and sample names.
Will place the 'feature_id_col' column as the index of the validated dataframe. 
The validated metabolome dataframe is stored as the 'validated_metabolome' attribute 


Usage
-----

validate_input_metabolome_df(self)


Returns
--------

self: object
    Object with metabolome_validated set to True


Examples
-----
Example of a validated output metabolome dataframe

+-------------+----------------+----------------+----------------+----------------+
| feature_id  | genotypeA_rep1 | genotypeA_rep2 | genotypeA_rep3 | genotypeA_rep4 |
+=============+================+================+================+================+
| metabolite1 |   1246         | 1245           | 12345          | 12458          |
+-------------+----------------+----------------+----------------+----------------+
| metabolite2 |   0            | 0              | 0              | 0              |
+-------------+----------------+----------------+----------------+----------------+
| metabolite3 |   10           | 0              | 0              | 154            |
+-------------+----------------+----------------+----------------+----------------+


validate_input_phenotype_df
***************************

Validates the dataframe containing the phenotype classes and the sample identifiers.


Usage
-----

validate_input_phenotype_df(
    self, 
    phenotype_class_col="phenotype")


Parameters
----------

phenotype_class_col: string, default="phenotype"
    The name of the column to be used 


Returns
-------

self: object
    Object with phenotype_validated set to True


Examples
--------

Example of a validated phenotype dataframe
        
+----------------+-----------+
| sample_id      | phenotype | 
+================+===========+
| genotypeA_rep1 | sensitive | 
+----------------+-----------+
| genotypeA_rep2 | sensitive |
+----------------+-----------+   
| genotypeA_rep3 | sensitive |
+----------------+-----------+
| genotypeA_rep4 | sensitive |
+----------------+-----------+ 
| genotypeB_rep1 | resistant |
+----------------+-----------+   
| genotypeB_rep2 | resistant |
+----------------+-----------+


get_baseline_performance
************************

Takes the phenotype and metabolome dataset and compute a simple Random Forest analysis with default hyperparameters. 
This will give a base performance for a Machine Learning model that has then to be optimised using autosklearn

k-fold cross-validation is performed to mitigate split effects on small datasets. 

get_baseline_performance(
    self, 
    kfold=5, 
    train_size=0.8,
    random_state=123,
    scoring_metric='balanced_accuracy')


Parameters
----------

kfold: int, optional
    Cross-validation strategy. Default is to use a 5-fold cross-validation. 

train_size: float or int, optional
    If float, should be between 0.5 and 1.0 and represent the proportion of the dataset to include in the train split.
    If int, represents the absolute number of train samples. If None, the value is automatically set to the complement of the test size.
    Default is 0.8 (80% of the data used for training).

random_state: int, optional
    Controls both the randomness of the train/test split  samples used when building trees (if bootstrap=True) and the sampling of the features to consider when looking for the best split at each node (if max_features < n_features). See Glossary for details.
    You can change this value several times to see how it affects the best ensemble model performance.
    Default is 123.


scoring_metric: str, optional
    A valid scoring value (default="balanced_accuracy")
    To get a complete list, type:
    >> from sklearn.metrics import SCORERS 
    >> sorted(SCORERS.keys()) 
    balanced accuracy is the average of recall obtained on each class. 


Returns
-------

self: object
    Object with baseline_performance attribute.


search_best_model_with_tpot_and_compute_pc_importances
******************************************************

Search for the best ML model with TPOT genetic programming methodology and extracts best Principal Components.

A characteristic of metabolomic data is to have a high number of features strongly correlated to each other.
This makes it difficult to extract the individual true feature importance. 
Here, this method implements a dimensionality reduction method (PCA) and the importances of each PC is computed. 

A resampling strategy called "cross-validation" will be performed on a subset of the data (training data) to increase 
the model generalisation performance. Finally, the model performance is tested on the unseen test data subset.  

By default, TPOT will make use of a set of preprocessors (e.g. Normalizer, PCA) and algorithms (e.g. RandomForestClassifier)
defined in the default config (classifier.py).
See: `<https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier.py>`_


Usage
-----

search_best_model_with_tpot_and_compute_pc_importances(
    self,
    class_of_interest,
    scoring_metric='balanced_accuracy',
    kfolds=3,
    train_size=0.8,
    max_time_mins=5,
    max_eval_time_mins=1,
    random_state=123,
    n_permutations=10,
    export_best_pipeline=True,
    path_for_saving_pipeline="./best_fitting_pipeline.py")


Parameters
----------

class_of_interest: str
    The name of the class of interest also called "positive class".
    This class will be used to calculate recall_score and precision_score. 
    Recall score = TP / (TP + FN) with TP: true positives and FN: false negatives.
    Precision score = TP / (TP + FP) with TP: true positives and FP: false positives. 

scoring_metric: str, optional
    Function used to evaluate the quality of a given pipeline for the classification problem. 
    Default is 'balanced accuracy'. 
    The following built-in scoring functions can be used:
    'accuracy', 'adjusted_rand_score', 'average_precision', 'balanced_accuracy', 
    'f1', 'f1_macro', 'f1_micro', 'f1_samples', 'f1_weighted', 'neg_log_loss', 
    'precision' etc. (suffixes apply as with ‘f1’), 'recall' etc. (suffixes apply as with ‘f1’), 
    ‘jaccard’ etc. (suffixes apply as with ‘f1’), 'roc_auc', ‘roc_auc_ovr’, ‘roc_auc_ovo’, ‘roc_auc_ovr_weighted’, ‘roc_auc_ovo_weighted’ 

kfolds: int, optional
    Number of folds for the stratified K-Folds cross-validation strategy. Default is 3-fold cross-validation. 
    Has to be comprised between 3 and 10 i.e. 3 <= kfolds =< 10
    See https://scikit-learn.org/stable/modules/cross_validation.html

train_size: float or int, optional
    If float, should be between 0.5 and 1.0 and represent the proportion of the dataset to include in the train split.
    If int, represents the absolute number of train samples. If None, the value is automatically set to the complement of the test size.
    Default is 0.8 (80% of the data used for training).

max_time_mins: int, optional
    How many minutes TPOT has to optimize the pipeline (in total). Default is 5 minutes.
    This setting will allow TPOT to run until max_time_mins minutes elapsed and then stop.
    Try short time intervals (5, 10, 15min) and then see if the model score on the test data improves. 

max_eval_time_mins: float, optional 
    How many minutes TPOT has to evaluate a single pipeline. Default is 1min. 
    This time has to be smaller than the 'max_time_mins' setting.

random_state: int, optional
    Controls both the randomness of the train/test split  samples used when building trees (if bootstrap=True) and the sampling of the features to consider when looking for the best split at each node (if max_features < n_features). See Glossary for details.
    You can change this value several times to see how it affects the best ensemble model performance.
    Default is 123.

n_permutations: int, optional
    Number of permutations used to compute feature importances from the best model using scikit-learn permutation_importance() method.
    Default is 10 permutations.

export_best_pipeline: `bool`, optional
    If True, the best fitting pipeline is exported as .py file. This allows for reuse of the pipeline on new datasets.
    Default is True. 

path_for_saving_pipeline: str, optional
    The path and filename of the best fitting pipeline to save.
    The name must have a '.py' extension. 
    Default to "./best_fitting_pipeline.py"


Returns
------

self: object
    The object with best model searched and feature importances computed. 


.. note:: Principal Component importances are calculated on the training set. Permutation importances can be computed either on the training set or on a held-out testing or validation set. Using a held-out set makes it possible to highlight which features contribute the most to the generalization power of the inspected model. Features that are important on the training set but not on the held-out set might cause the model to overfit. `<https://scikit-learn.org/stable/modules/permutation_importance.html#permutation-importance>`_


get_names_of_top_n_features_from_selected_pc
********************************************

Get the names of features with highest loading scores on selected PC  

Takes the matrix of loading scores of shape (n_samples, n_features) and the metabolome dataframe of shape (n_features, n_samples)
and extract the names of features. 
The loadings matrix is available after running the search_best_model_with_tpot_and_compute_pc_importances() method.


Usage
-----

get_names_of_top_n_features_from_selected_pc(
    self, 
    selected_pc=1, 
    top_n=5)


Parameters
----------

selected_pc: int, optional
    Principal Component to keep. 1-based index (1 selects PC1, 2 selected PC2, etc.)
    Default is 1.
top_n: int, optional
    Number of features to select. 
    The top_n features with the highest absolute loadings will be selected from the selected_pc PC. 
    For instance, the top 5 features from PC1 will be selected with selected_pc=1 and top_n=5.
    Default is 5.


Returns
-------

A list of feature names. 