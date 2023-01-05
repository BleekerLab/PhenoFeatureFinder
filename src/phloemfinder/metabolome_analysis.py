#!/usr/bin/env python3 

import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt

def _extract_samples_to_condition(df, name_grouping_var='genotype', separator_replicates='_'):
    '''
    A utility function to extract the grouping factor (e.g. 'genotype') from sample names. 
    
    Uses melting (wide to long) and split grouping variable from biological replicates using specified separator.
    
    Parameters
    ----------
    df: 'pandas.core.DataFrame'
    name_grouping_var: str, optional 
        Name of the variable used as grouping variable (default is 'genotype').
    separator_replicates: str, optional
        The separator between the grouping variable and the biological replicates ( default is underscore '_')
    
    Returns
    -------
    A dataframe with the correspondence between samples and experimental condition (grouping variable).

    Notes
    -------
    Input dataframe
                        | genotypeA_rep1 | genotypeA_rep2 | genotypeA_rep3 | genotypeA_rep4 |
                        |----------------|----------------|----------------|----------------|
          feature_id
        | metabolite1   |   1246         | 1245           | 12345          | 12458          |
        | metabolite2   |   0            | 0              | 0              | 0              |
        | metabolite3   |   10           | 0              | 0              | 154            |
    
    Output dataframe
        
        | sample_id          | genotype       | replicate      |
        |--------------------|----------------|----------------|
        | genotypeA_rep1     |   genotypeA    | rep1           |
        | genotypeA_rep2     |   genotypeA    | rep2           |
        | genotypeA_rep3     |   genotypeA    | rep3           |
        | genotypeA_rep4     |   genotypeA    | rep4           |
        | etc.
    '''
    melted_df = pd.melt(df.reset_index(), id_vars='feature_id', var_name="sample")
    melted_df[[name_grouping_var, 'rep']] = melted_df["sample"].str.split(pat=separator_replicates, expand=True)
    melted_df_parsed = melted_df.drop(["feature_id", "value"], axis=1)
    melted_df_parsed_dedup = melted_df_parsed.drop_duplicates()
    return melted_df_parsed_dedup

###################
## Class definition 
###################

class MetaboliteAnalysis:
    '''
    A class to streamline the filtering and exploration of a metabolome dataset.   


    Parameters
    ----------
    metabolome_csv: `str`
        A path to a .csv file with the metabolome data (scaled or unscaled).
        Shape of the dataframe is usually (n_samples, n_features) with n_features >> n_samples
        
    metabolome_feature_id_col: `str`, optional
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


    Methods
    -------
    validate_input_metabolome_df
      Check if the provided metabolome file is suitable. 
      Turns attribute metabolome_validated to True. 
    discard_features_detected_in_blanks
      Removes features only detected in blank samples. 
    filter_out_unreliable_features()
      Filter out features not reliably detectable in replicates of the same grouping factor. 
      For instance, if a feature is detected less than 4 times within 4 biological replicates.  
    write_clean_metabolome_to_csv()
       

    Examples
    --------
    >>> from metabolome_analysis import MetaboliteAnalysis
    >>> met = MetaboliteAnalysis(
        metabolome_csv="raw_data_processed_for_machine_learning.csv",
        metabolome_feature_id_col="feature_id")
    >>> met.
    
    Notes
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
    | rt-0.08_mz-343.22434 | 3078    | 2660    | 2786    | 3096    | 3084  | 3190  | 3072  | 2940  | 3620     | 3200     | 3154     | 3146     |
    +----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
    | rt-0.08_mz-287.16171 | 1506    | 1510    | 1546    | 1420    | 1582  | 1786  | 1512  | 1288  | 1606     | 1824     | 1404     | 1356     |
    +----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
    | rt-0.08_mz-629.84892 | 1306    | 2088    | 0       | 1750    | 1514  | 1546  | 0     | 0     | 1292     | 1818     | 1378     | 0        |
    +----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
    | rt-0.08_mz-498.39989 | 1042    | 1274    | 1138    | 1016    | 964   | 1166  | 1076  | 1052  | 1348     | 984      | 1114     | 1042     |
    +----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
    | rt-0.08_mz-497.39637 | 2670    | 2414    | 2646    | 2656    | 2462  | 2836  | 2744  | 2854  | 3188     | 2872     | 2516     | 2946     |
    +----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
    | rt-0.08_mz-318.98083 | 946     | 1068    | 1190    | 1058    | 834   | 858   | 852   | 968   | 784      | 1006     | 1102     | 936      |
    +----------------------+---------+---------+---------+---------+-------+-------+-------+-------+----------+----------+----------+----------+
    
    See also
    --------
    scikit PCA: https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
    
    '''

    ##########################
    # Class constructor method
    ##########################
    def __init__(
        self, 
        metabolome_csv, 
        metabolome_feature_id_col='feature_id'):
        """
        Constructor method. 
        Returns a Python instance of class MetabolomeAnalysis 
        """
        # By default the metabolome and phenotype data imported from .csv files will have to be validated
        # By default all filters have not been executed (blank filtering, etc.)
        self.metabolome_validated = False
        self.blank_features_filtered = False
        self.unreliable_features_filtered = False
        self.pca_performed = False
        
        # Import metabolome dataframe and verify presence of feature id column
        self.metabolome = pd.read_csv(metabolome_csv)
        if metabolome_feature_id_col not in self.metabolome.columns:
            raise ValueError("The specified column with feature identifiers {0} is not present in your '{1}' file.".format(metabolome_feature_id_col,os.path.basename(metabolome_csv)))
        else:
            self.metabolome.set_index(metabolome_feature_id_col, inplace=True)
    
    def validate_input_metabolome_df(self, metabolome_feature_id_col='feature_id'):
        '''
        Validates the dataframe containing the feature identifiers, metabolite values and sample names.
        Will place the 'feature_id_col' column as the index of the validated dataframe. 
        The validated metabolome dataframe is stored as the 'validated_metabolome' attribute. 
        
        Parameters
        ----------
        metabolome_feature_id: str, optional 
            The name of the column that contains the feature identifiers (default is 'feature_id').
            Feature identifiers should be unique (=not duplicated).
            
        Returns
        -------
        self: object
          Object with attribute metabolome_validated set to True if tests are passed. 
        
        Notes
        -----
        Example of a valid input metabolome dataframe


                      | genotypeA_rep1 | genotypeA_rep2 | genotypeA_rep3 | genotypeA_rep4 |
                      |----------------|----------------|----------------|----------------|
          feature_id
        | metabolite1 |   1246         | 1245           | 12345          | 12458          |
        | metabolite2 |   0            | 0              | 0              | 0              |
        | metabolite3 |   10           | 0              | 0              | 154            |
        
        '''

        if np.any(self.metabolome.values < 0):
            raise ValueError("Sorry, metabolite values have to be zero or positive integers (>=0)")
        else:
            self.metabolome_validated = True

    #############################################
    ### Filter features detected in blank samples
    #############################################
    def discard_features_detected_in_blanks(
        self, 
        blank_sample_contains='blank'):
        '''
        Removes 
        Steps:
          1. Sum the abundance of each feature in the blank samples.
          2. Makes a list of features to be discarded (features with a positive summed abundance).
          3. Returns a filtered Pandas dataframe with only features not detected in blank samples

        Parameters
        ----------
        blank_sample_contains: str, optional.
            Column names with this name will be considered blank samples.
            Default is='blank'
        
        Returns
        -------
        metabolome: ndarray
            A filtered Pandas dataframe without features detected in blank samples and with the blank samples removed. 
        '''

        if self.metabolome_validated == True:
            pass
        else:
            self.validate_input_metabolome_df()
        
        blank_cols = [col for col in self.metabolome.columns.values.tolist() if blank_sample_contains in col]

        # If the sum of a feature in blank samples is higher than 0 then 
        # this feature should be removed
        # only keep features that are not detectable in blank samples
        self.metabolome["sum_features"] = self.metabolome[blank_cols].sum(axis=1)
        self.metabolome = self.metabolome.query("sum_features == 0")
        
        # Remove columns with blank samples and the sum column used for filtering
        self.metabolome = self.metabolome.drop(blank_cols, axis=1)
        self.metabolome = self.metabolome.drop("sum_features", axis=1)    

    
    #######################################################################################
    ### Filter features that are not reliable (detected less than nb replicates in a group)
    ######################################################################################
    def filter_out_unreliable_features(
        self,
        name_grouping_var="genotype", 
        nb_times_detected=4,
        separator_replicates='_'):
        '''
        Removes features not reliably detectable in multiple biological replicates from the same grouping factor. 

        Takes a dataframe with feature identifiers in index and samples as columns.
        Step 1: First melt and split the sample names to generate the grouping variable
        Step 2: count number of times a metabolite is detected in the groups. 
        If number of times detected in a group = number of biological replicates then it is considered as reliable
        Each feature receives a tag  'reliable' or 'not_reliable'
        Step 3: discard the 'not_reliable' features and keep the filtered dataframe. 

        Params
        ------
        name_grouping_var: str, optional
            The name used when splitting between replicate and main factor.
            For example "genotype" when splitting MM_rep1 into 'MM' and 'rep1'.
            Default is 
        nb_times_detected: int, optionaldefault=4
            Number of times a metabolite should be detected to be considered 'reliable'. 
            Should be equal to the number of biological replicates for a given group of interest (e.g. genotype)
        separator_replicates: string, default="_"
            The separator to split sample names into a grouping variable (e.g. genotype) and the biological replicate number (e.g. 1)
        

        Returns
        -------
        metabolome: ndarray
            A Pandas dataframe with only features considered as reliable, sample names and their values. 
        
        Notes 
        -----
        Input dataframe

                             	| MM_1  	| MM_2  	| MM_3  	| MM_4  	| LA1330_1 	| LA1330_2 	|
                            	|----------	|----------	|----------	|----------	|----------	|----------	|
          feature_id           	 
        | rt-0.04_mz-241.88396 	| 554   	| 678   	| 674   	| 936   	| 824      	| 940      	|
        | rt-0.05_mz-143.95911 	| 1364  	| 1340  	| 1692  	| 1948  	| 1928     	| 1956     	|
        | rt-0.06_mz-124.96631 	| 0      	| 0     	| 0     	| 888   	| 786      	| 668      	|
        | rt-0.08_mz-553.45905 	| 10972 	| 11190 	| 12172 	| 11820 	| 12026    	| 11604    	|

        Output df (rt-0.06_mz-124.96631 is kicked out because 3x0 and 1x888 in MM groups)
        
                             	| MM_1  	| MM_2  	| MM_3  	| MM_4  	| LA1330_1 	| LA1330_2 	|
                            	|----------	|----------	|----------	|----------	|----------	|----------	|
          feature_id           	 
        | rt-0.04_mz-241.88396 	| 554   	| 678   	| 674   	| 936   	| 824      	| 940      	|
        | rt-0.05_mz-143.95911 	| 1364  	| 1340  	| 1692  	| 1948  	| 1928     	| 1956     	|
        | rt-0.08_mz-553.45905 	| 10972 	| 11190 	| 12172 	| 11820 	| 12026    	| 11604    	|


        '''
        df = self.metabolome

        ### Melt (required to tag reliable/not reliable features)
        melted_df = pd.melt(
            df.reset_index(), 
            id_vars='feature_id', 
            var_name="sample")
        melted_df[[name_grouping_var, 'rep']] = melted_df["sample"].str.split(pat=separator_replicates, expand=True)
        sorted_melted_df = melted_df.sort_values(by=['feature_id', name_grouping_var])
        sorted_melted_df.set_index('feature_id', inplace=True)
        sorted_melted_df_parsed = sorted_melted_df.drop('sample', axis=1)
        
        ### Identify features that are reliable
        # Creates a dictionary to have the feature identifier and a 'reliable'/'not_reliable' tag
        reliability_feature_dict = {}
        features = sorted_melted_df_parsed.index.unique().tolist()

        for feature in features:
          # Dataframe containing only one metabolite
          temp_df = sorted_melted_df_parsed.loc[feature,:]
          temp_df = temp_df.drop('rep', axis=1)

          # number of values above 0 across all groups
          min_number_of_values_above_0_across_all_groups = temp_df[temp_df['value'] > 0].groupby(by=name_grouping_var)[name_grouping_var].count().max()

          # If the feature is detected a minimum of times equal to the number of biological replicates
          # This means the feature is reliably detectable in at least one group (e.g. one genotype)
          if min_number_of_values_above_0_across_all_groups >= nb_times_detected:
            reliability_feature_dict[feature] = "reliable"    
          else:
            reliability_feature_dict[feature] = 'not_reliable'

        # Convert to a Pandas to prepare the filtering of the feature/abundance dataframe
        reliability_df = pd.DataFrame.from_dict(
            reliability_feature_dict, orient="index", columns=["reliability"]).reset_index()
        
        reliability_df = reliability_df.rename({"index":"feature_id"}, axis='columns')
        reliability_df = reliability_df[reliability_df["reliability"] == "reliable"]
        features_to_keep = reliability_df.feature_id.tolist()

        df_reliable_features = df.loc[features_to_keep,:]
                
        self.metabolome = df_reliable_features
        self.unreliable_features_filtered = True

    #################################################
    ### Write filtered metabolomoe data to a csv file
    #################################################

    def write_clean_metabolome_to_csv(self, path_of_cleaned_csv="filtered_metabolome.csv"):
        '''
        A function that verify that the metabolome dataset has been cleaned up. 
        Writes the metabolome data as a comma-separated value file on disk

        Parameters
        ----------
        path_of_cleaned_csv: str, optional
            The path and filename of the .csv file to save.
            Default to "./filtered_metabolome.csv" 
        '''
        try:
            self.blank_features_filtered == True
        except:
            raise ValueError("Features in blank should be removed first using the 'discard_features_detected_in_blanks() method.")
        
        try:
            self.unreliable_features_filtered == True
        except:
            raise ValueError("Features not reliably detected within at least one group should be removed first using the 'filter_out_unreliable_features() method.") 
        
        self.metabolome.to_csv(path_or_buf=path_of_cleaned_csv, sep=',')



    #################################
    ## Principal Component Analysis
    #################################

    def compute_pca_on_metabolites(self, scale=True, n_principal_components=10, auto_transpose=True):
        """
        Performs a Principal Component Analysis (PCA) on the metabolome data. 
        
        The PCA analysis will return transformed coordinates of the samples in a new space. 
        It will also give the percentage of variance explained by each Principal Component. 
        Assumes that number of samples < number of features/metabolites
        Performs a transpose of the metabolite dataframe if n_samples > n_features (this can be turned off with auto_transpose)
        
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
        """
        # Verify that samples are in rows and features in columns
        # Usually n_samples << n_features so we should have n_rows << n_cols
        n_rows = self.metabolome.shape[0]
        n_cols = self.metabolome.shape[1]
        
        metabolite_df = self.metabolome
        if n_rows > n_cols:
            # Likely features are in row so transpose to have samples in rows
            metabolite_df = metabolite_df.transpose()
        else:
            pass
    
        if scale == True:
            scaler = StandardScaler(with_mean=True, with_std=True)
            metabolite_df_scaled = scaler.fit_transform(metabolite_df)
            metabolite_df_scaled = pd.DataFrame(metabolite_df_scaled)
            metabolite_df_scaled.columns = metabolite_df.columns
            metabolite_df_scaled.set_index(metabolite_df.index.values)       
        else:
            pass

        if n_principal_components <= np.minimum(n_rows, n_cols):
            pca = PCA(n_components=n_principal_components)
            metabolite_df_scaled_transformed = pca.fit_transform(metabolite_df_scaled)
            exp_variance = pd.DataFrame(pca.explained_variance_ratio_.round(2)*100, columns=["explained_variance"])
        # If n_principal_components > min(n_samples, n_features)
        # then n_principal_components = min(n_samples, n_features)
        else:
            n_principal_components > np.minimum(n_rows, n_cols)
            n_principal_components = np.minimum(n_rows, n_cols)
            pca = PCA(n_components=n_principal_components)
            metabolite_df_scaled_transformed = pca.fit_transform(metabolite_df_scaled)
            exp_variance = pd.DataFrame(pca.explained_variance_ratio_.round(2)*100, columns=["explained_variance"])         

        # The numbering of the components starts by default at 0. 
        # Setting this to 1 to make it more user friendly
        exp_variance.index = exp_variance.index+1
        
        # Store PCA results 
        self.exp_variance = exp_variance
        self.metabolome_pca_reduced = metabolite_df_scaled_transformed
        self.pca_performed = True

    def create_scree_plot(self, plot_file_name=None):
        '''
        Returns a barplot with the explained variance per Principal Component. 
        Has to be preceded by perform_pca()

        Parameters
        ---------
        plot_file_name: string, default='None'
          Path to a file where the plot will be saved.
          For instance 'my_scree_plot.pdf'
        '''
        try:
            self.pca_performed
        except: 
            raise AttributeError("Please compute the PCA first using the compute_pca_on_metabolites() method.") 
        
        sns.barplot(
            data=self.exp_variance, 
            x=self.exp_variance.index, 
            y="explained_variance")
        plt.xlabel("Principal Component")
        plt.ylabel("Explained variance (%)")

        # Optionally save the plot
        if plot_file_name != None:
            plot_dirname = os.path.dirname(plot_file_name)
            if plot_dirname == '': # means file will be saved in current working directory
                plt.savefig(plot_file_name)
            else:
                os.makedirs(plot_dirname, exist_ok=True)
                plt.savefig(plot_file_name)

        plt.show()


    def create_sample_score_plot(
        self, 
        pc_x_axis=1, 
        pc_y_axis=2, 
        name_grouping_var='genotype',
        separator_replicates="_",
        plot_file_name=None):
        '''
        Returns a sample score plot of the samples on PCx vs PCy. 
        Samples are colored based on the grouping variable (e.g. genotype)

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


        '''
        if self.pca_performed:
            pass
        else:
            raise AttributeError("Please compute the PCA first using the compute_pca_on_metabolites() method.") 

        n_features = self.metabolome.shape[0]
        n_samples = self.metabolome.shape[1]
        min_of_samples_and_features = np.minimum(n_samples, n_features)
        
        samples_to_conditions = _extract_samples_to_condition(df=self.metabolome, name_grouping_var=name_grouping_var, separator_replicates=separator_replicates)

        if pc_x_axis == pc_y_axis:
            raise ValueError("Values for Principal Components on x axis and y axis have to be different.")
        if pc_x_axis > min_of_samples_and_features:
            raise ValueError("Your principal component for x axis should be lower than {0}".format(min_of_samples_and_features))
        if pc_y_axis > np.minimum(n_samples, n_features):
            raise ValueError("Your principal component for y axis should be lower than {0}".format(min_of_samples_and_features))
       
        
        if not name_grouping_var in samples_to_conditions.columns:
            raise IndexError("The grouping variable '{0}' is not present in the samples_to_condition dataframe".format(name_grouping_var))
        else:
            # Build the plot
            plt.figure(figsize=(10,7))
            
            self.scatter_plot = sns.scatterplot(
            x=self.metabolome_pca_reduced[:,pc_x_axis-1],
            y=self.metabolome_pca_reduced[:,pc_y_axis-1],
            hue=samples_to_conditions[name_grouping_var],
            s=200)

            plt.xlabel("PC" + str(pc_x_axis) + ": " + str(self.exp_variance.iloc[pc_x_axis,0].round(2)) + "% variance") 
            plt.ylabel("PC" + str(pc_y_axis) + ": " + str(self.exp_variance.iloc[pc_y_axis,0].round(2)) + "% variance")
            plt.title("PC" + str(pc_x_axis) + " vs PC" + str(pc_y_axis))

            # Optionally save the plot
            if plot_file_name != None:
                plot_dirname = os.path.dirname(plot_file_name)
                if plot_dirname == '': # means file will be saved in current working directory
                    plt.savefig(plot_file_name)
                else:
                    os.makedirs(plot_dirname, exist_ok=True)
                    plt.savefig(plot_file_name)

            plt.show()

