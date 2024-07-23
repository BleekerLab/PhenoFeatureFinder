#!/usr/bin/env python3 

#########################################
# Libraries and TPOT AutoML configuration
#########################################
import os
import sys
from warnings import WarningMessage
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.decomposition import PCA
from sklearn.metrics import balanced_accuracy_score
from sklearn.tree import DecisionTreeClassifier

from tpot import TPOTClassifier
from tpot.export_utils import set_param_recursive

from phenofeaturefinder.utils import compute_metrics_classification 


# TPOT automated ML custom configuration dictionary
# These classifiers and pre-processing steps work best on multidimensional omics datasets (correlated features, N samples << p features)
## See https://epistasislab.github.io/tpot/using/#built-in-tpot-configurations
tpot_custom_config = {
   # Models to test 
  'sklearn.tree.DecisionTreeClassifier': {
        'criterion': ["gini", "entropy"],
        'max_depth': range(1, 11),
        'min_samples_split': range(2, 21),
        'min_samples_leaf':  range(1, 21),
    },
  'sklearn.ensemble.RandomForestClassifier': {
        'n_estimators': [1000],
        'criterion': ["gini", "entropy"],
        'max_features': np.arange(0.05, 1.01, 0.05),
        'min_samples_split': range(2, 21),
        'min_samples_leaf':  range(1, 21),
        'bootstrap': [True, False]
    },
    'sklearn.ensemble.GradientBoostingClassifier': {
        'n_estimators': [1000],
        'learning_rate': [1e-3, 1e-2, 1e-1, 0.5, 1.],
        'max_depth': range(1, 11),
        'min_samples_split': range(2, 21),
        'min_samples_leaf': range(1, 21),
        'subsample': np.arange(0.05, 1.01, 0.05),
        'max_features': np.arange(0.05, 1.01, 0.05)
    }, 
    'xgboost.XGBClassifier': {
        'n_estimators': [1000],
        'max_depth': range(1, 11),
        'learning_rate': [1e-3, 1e-2, 1e-1, 0.5, 1.],
        'subsample': np.arange(0.05, 1.01, 0.05),
        'min_child_weight': range(1, 21),
        'n_jobs': [1],
        'verbosity': [0]
    },

    # Preprocesssors
    'sklearn.preprocessing.Binarizer': {
        'threshold': np.arange(0.0, 1.01, 0.05)
    },

    'sklearn.cluster.FeatureAgglomeration': {
        'linkage': ['ward', 'complete', 'average'],
        'affinity': ['euclidean', 'l1', 'l2', 'manhattan', 'cosine']
    },

    'sklearn.preprocessing.MaxAbsScaler': {
    },

    'sklearn.preprocessing.MinMaxScaler': {
    },

    'sklearn.preprocessing.Normalizer': {
        'norm': ['l1', 'l2', 'max']
    },

    'sklearn.kernel_approximation.Nystroem': {
        'kernel': ['rbf', 'cosine', 'chi2', 'laplacian', 'polynomial', 'poly', 'linear', 'additive_chi2', 'sigmoid'],
        'gamma': np.arange(0.0, 1.01, 0.05),
        'n_components': range(1, 11)
    },

    'sklearn.decomposition.PCA': {
        'svd_solver': ['randomized'],
        'iterated_power': range(1, 11)
    },

    'sklearn.preprocessing.PolynomialFeatures': {
        'degree': [2],
        'include_bias': [False],
        'interaction_only': [False]
    },

    'sklearn.kernel_approximation.RBFSampler': {
        'gamma': np.arange(0.0, 1.01, 0.05)
    },

    'sklearn.preprocessing.RobustScaler': {
    },

    'sklearn.preprocessing.StandardScaler': {
    },

    'tpot.builtins.ZeroCount': {
    },

    'tpot.builtins.OneHotEncoder': {
        'minimum_fraction': [0.05, 0.1, 0.15, 0.2, 0.25],
        'sparse': [False],
        'threshold': [10]
    }
}

####################################
# End of library and config sections
####################################

class FeatureSelection:
    '''
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

    

    Notes
    --------

    Example of an input metabolome .csv file

        | feature_id  | genotypeA_rep1 | genotypeA_rep2 | genotypeA_rep3 | genotypeA_rep4 |
        |-------------|----------------|----------------|----------------|----------------|
        | metabolite1 |   1246         | 1245           | 12345          | 12458          |
        | metabolite2 |   0            | 0              | 0              | 0              |
        | metabolite3 |   10           | 0              | 0              | 154            |

    Example of an input phenotype .csv file

        | sample_id      | phenotype | 
        |----------------|-----------|
        | genotypeA_rep1 | sensitive | 
        | genotypeA_rep2 | sensitive |   
        | genotypeA_rep3 | sensitive |
        | genotypeA_rep4 | sensitive | 
        | genotypeB_rep1 | resistant |   
        | genotypeB_rep2 | resistant |
    
    '''
    # Class attribute shared among all instances of the class
    # By default the metabolome and phenotype data imported from .csv files will have to be validated
    # By default all filters have not been executed (blank filtering, etc.)
    # Baseline performance of a simple ML model as well as search of best model are also None/False by default. 
    metabolome_validated=False
    phenotype_validated=False
    baseline_performance=None
    best_ensemble_models_searched=False

    # Class constructor method
    def __init__(
        self, 
        metabolome_csv, 
        phenotype_csv,
        metabolome_feature_id_col='feature_id', 
        phenotype_sample_id='sample_id'):
        
        # Import metabolome dataframe and verify presence of feature id column
        self.metabolome = pd.read_csv(metabolome_csv)
        if metabolome_feature_id_col not in self.metabolome.columns:
            raise ValueError("The specified column with feature identifiers '{0}' is not present in your '{1}' file.".format(metabolome_feature_id_col,os.path.basename(metabolome_csv)))
        else:
            self.metabolome.set_index(metabolome_feature_id_col, inplace=True)

        # Import phenotype dataframe and verify presence of sample id column
        self.phenotype = pd.read_csv(phenotype_csv)
        if phenotype_sample_id not in self.phenotype.columns:
            raise ValueError("The specified column with sample identifiers '{0}' is not present in your '{1}' file.".format(phenotype_sample_id, os.path.basename(phenotype_csv)))
        else:
            try: 
                self.phenotype.set_index(phenotype_sample_id, inplace=True)
            except:
                raise IndexError("Values for sample identifiers have to be unique. Check your ", phenotype_sample_id, " column.")

    ################
    ## Verify inputs
    ################
    def validate_input_metabolome_df(self):
        '''
        Validates the dataframe containing the feature identifiers, metabolite values and sample names.
        Will place the 'feature_id_col' column as the index of the validated dataframe. 
        The validated metabolome dataframe is stored as the 'validated_metabolome' attribute 
        
        
        Returns
        --------
        self: object
          Object with metabolome_validated set to True

        
        Notes
        -----
        Example of a validated output metabolome dataframe

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
            print("Metabolome data validated.")
    
    def validate_input_phenotype_df(self, phenotype_class_col="phenotype"):
        '''
        Validates the dataframe containing the phenotype classes and the sample identifiers

        Params
        ------
        phenotype_class_col: string, default="phenotype"
            The name of the column to be used 

        Returns
        --------
        self: object
          Object with phenotype_validated set to True

        Notes
        --------
        Example of an input phenotype dataframe
        

        | sample_id      | phenotype | 
        |----------------|-----------|
        | genotypeA_rep1 | sensitive | 
        | genotypeA_rep2 | sensitive |   
        | genotypeA_rep3 | sensitive |
        | genotypeA_rep4 | sensitive | 
        | genotypeB_rep1 | resistant |   
        | genotypeB_rep2 | resistant |

        Example of a validated output phenotype dataframe. 

                         | phenotype | 
                         |-----------|
          sample_id      
        | genotypeA_rep1 | sensitive | 
        | genotypeA_rep2 | sensitive |   
        | genotypeA_rep3 | sensitive |
        | genotypeA_rep4 | sensitive | 
        | genotypeB_rep1 | resistant |   
        | genotypeB_rep2 | resistant |

        Example
        -------
        >> fs = FeatureSelection(
        >>        metabolome_csv="clean_metabolome.csv", 
        >>        phenotype_csv="phenotypes_test_data.csv", 
        >>        phenotype_sample_id='sample')
        >> fs.validate_input_phenotype_df()

        '''
        n_distinct_classes = self.phenotype[phenotype_class_col].nunique()
        try:
            n_distinct_classes == 2
            self.phenotype_validated = True    
            print("Phenotype data validated.")
        except:
            raise ValueError("The number of distinct phenotypic classes in the {0} column should be exactly 2.".format(phenotype_class_col))
    
    #################
    ## Baseline model
    #################
    def get_baseline_performance(
      self, 
      kfold=5, 
      train_size=0.8,
      random_state=123,
      scoring_metric='balanced_accuracy'):
        '''
        Takes the phenotype and metabolome dataset and compute a simple Random Forest analysis with default hyperparameters. 
        This will give a base performance for a Machine Learning model that has then to be optimised using autosklearn

        k-fold cross-validation is performed to mitigate split effects on small datasets. 

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
        
        Example
        -------
        >>> fs = FeatureSelection(
                   metabolome_csv="../tests/clean_metabolome.csv", 
                   phenotype_csv="../tests/phenotypes_test_data.csv", 
                   phenotype_sample_id='sample')
            fs.get_baseline_performance()

        '''
        try:
            self.metabolome_validated == True
        except:
            raise ValueError("Please validate metabolome data first using the validate_input_metabolome_df() method.")

        try:
            self.phenotype_validated == True
        except:
            raise ValueError("Please validate phenotype data first using the validate_input_phenotype_df() method.")

        X = self.metabolome.transpose()
        y = self.phenotype.values.ravel() # ravel makes the array contiguous
        X_train, X_test, y_train, y_test = train_test_split(
          X, y, 
          train_size=train_size, 
          random_state=random_state, 
          stratify=y)

        # Train model and assess performance
        clf = RandomForestClassifier(n_estimators=1000, random_state=random_state)
        scores = cross_val_score(clf, X_train, y_train, scoring=scoring_metric, cv=kfold)
        average_scores = np.average(scores).round(3) * 100
        stddev_scores = np.std(scores).round(3) * 100
        print("====== Training a basic Random Forest model =======")
        baseline_performance = "Average {0} score on training data is: {1:.3f} % -/+ {2:.2f}".format(scoring_metric, average_scores, stddev_scores)
        print(baseline_performance)
        print("\n")
        print("====== Performance on test data of the basic Random Forest model =======")
        clf.fit(X_train, y_train)
        predictions = clf.predict(X_test)
        model_balanced_accuracy_score = balanced_accuracy_score(y_true=y_test, y_pred=predictions).round(3) * 100 
        print("Average {0} score on test data is: {1:.3f} %".format(scoring_metric, model_balanced_accuracy_score))
        self.baseline_performance = baseline_performance
    

    # TODO: make decorator function to check arguments     
    # See -> https://www.pythonforthelab.com/blog/how-to-use-decorators-to-validate-input/
    def search_best_model_with_tpot_and_compute_pc_importances(
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
        path_for_saving_pipeline="./best_fitting_pipeline.py"):
      '''
      Search for the best ML model with TPOT genetic programming methodology and extracts best Principal Components.
 
      A characteristic of metabolomic data is to have a high number of features strongly correlated to each other.
      This makes it difficult to extract the individual true feature importance. 
      Here, this method implements a dimensionality reduction method (PCA) and the importances of each PC is computed. 

      A resampling strategy called "cross-validation" will be performed on a subset of the data (training data) to increase 
      the model generalisation performance. Finally, the model performance is tested on the unseen test data subset.  
      
      By default, TPOT will make use of a set of preprocessors (e.g. Normalizer, PCA) and algorithms (e.g. RandomForestClassifier)
      defined in the default config (classifier.py).
      See: https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier.py

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

      See also
      --------

      Notes
      -----
      Principal Component importances are calculated on the training set
      Permutation importances can be computed either on the training set or on a held-out testing or validation set.
      Using a held-out set makes it possible to highlight which features contribute the most to the generalization power of the inspected model. 
      Features that are important on the training set but not on the held-out set might cause the model to overfit.
      https://scikit-learn.org/stable/modules/permutation_importance.html#permutation-importance 
          
      '''
      X = self.metabolome.transpose().to_numpy(dtype='float64')
      y = self.phenotype.values.ravel()
        
      # Verify input arguments
      if max_time_mins > max_eval_time_mins:
         pass # continue execution
      else:
        print("ValueError: max_time_mins (all pipelines) has to be higher than max_eval_time_mins (single pipeline)")
      
      if 3 <= kfolds <= 10:
        pass
      else:
        print("ValueError: kfolds argument has to be comprised between 3 and 10")
      
      if 0.5 < train_size < 0.9:
        pass
      else: 
        print("ValueError: train_size has to be comprised between 0.5 and 0.9")

      if class_of_interest in set(y):
        pass 
      else:
        print('The class_of_interest value "{0}" has to be in the phenotype labels {1}'.format(class_of_interest, set(y)))
      
      ### Automated search for best model/pipeline
      # First step is a PCA to avoid to work with correlated features 
      # Feature importances become Principal Component importances
      number_of_components = np.min(X.shape) # minimum of (n_samples, n_features)
      pca = PCA(n_components=number_of_components, random_state=random_state)
      X_reduced = pca.fit_transform(X)
      
      X_train, X_test, y_train, y_test = train_test_split(
          X_reduced, y, 
          train_size=train_size, 
          random_state=random_state, 
          stratify=y)    
      tpot = TPOTClassifier(max_time_mins=max_time_mins, max_eval_time_mins=max_eval_time_mins, cv=kfolds, config_dict=tpot_custom_config, random_state=random_state, verbosity=2)
      tpot.fit(X_train, y_train)
      best_pipeline = tpot.fitted_pipeline_
      set_param_recursive(best_pipeline.steps, 'random_state', random_state)
      
      ### Model performance
      predictions = best_pipeline.predict(X_test)
      training_score = best_pipeline.score(X_train, y_train) * 100
      print("============ Performance of ML model on train data =============")
      print("Train {0} score {1:.3f} %".format(scoring_metric, training_score))
      print("\n")
      print("============ Performance of ML model on test data =============")
      print(compute_metrics_classification(y_predictions=predictions, y_trues=y_test, positive_class=class_of_interest))

      ### Compute Principal Components importances
      # Has to be done on the same train/test split. 
      print("\n")
      print("======== Computing Principal Components importances on the training set =======")
      pc_importances_training_set = permutation_importance(
        best_pipeline, 
        X=X_train, 
        y=y_train, 
        scoring=scoring_metric, 
        n_repeats=n_permutations, 
        random_state=random_state)
      mean_imp = pd.DataFrame(pc_importances_training_set.importances_mean, columns=["mean_var_imp"])
      std_imp = pd.DataFrame(pc_importances_training_set.importances_std, columns=["std_var_imp"])
      raw_imp = pd.DataFrame(pc_importances_training_set.importances, columns=["perm" + str(i) for i in range(n_permutations)])
      pc_importances = pd.concat([mean_imp, std_imp, raw_imp], axis=1).sort_values('mean_var_imp', ascending=False)
      pc_importances["pc"] = ["PC" + str(i) for i in pc_importances.index.values]
      pc_importances.set_index("pc", inplace=True)
      # Save results
      self.best_model = best_pipeline
      self.pc_importances = pc_importances
      self.loadings = np.absolute(pca.components_) # required for downstream analyses (extraction of important features based on their loadings)

      if export_best_pipeline is True:
         tpot.export(path_for_saving_pipeline)
      else:
         pass

    
    def get_names_of_top_n_features_from_selected_pc(self, selected_pc=1, top_n=5):
        """
        Get the names of features with highest loading scores on selected PC  

        Takes the matrix of loading scores of shape (n_samples, n_features) and the metabolome dataframe of shape (n_features, n_samples)
        and extract the names of features. 
        The loadings matrix is available after running the search_best_model_with_tpot_and_compute_pc_importances() method.

        Params
        ------
        selected_pc: int, optional
          Principal Component to keep. 1-based index (1 selects PC1, 2 selected PC2, etc.)
          Default is 1.
        top_n: int, optional
          Number of features to select. 
          The top_n features with the highest absolute loadings will be selected from the selected_pc PC. 
          For instance, the top 5 features from PC1 will be selected with selected_pc=1 and top_n=5.
          Default is 5.

        Returns:
          A list of feature names. 
        """
        metabolite_df = self.metabolome.transpose() # to have shape (n_samples, n_features) same as loadings 
        loadings = self.loadings
        if loadings is None:
          print("Please run the search_best_model_with_tpot_and_compute_pc_importances() method first.")

        assert isinstance(selected_pc, int), "Please select an integer higher or equal to 1 for the selected_pc argument"
        assert isinstance(top_n, int), "Please select an integer higher or equal to 1 for the top_n argument"
        assert selected_pc > 0, "Please select a PC equal or higher than 1"
        assert top_n > 0, "Please select a number of top features to return equal or higher than 1"


        zero_off_selected_pc = selected_pc - 1 # avoid 1-off error

        loadings_of_selected_pc = self.loadings[zero_off_selected_pc]
        loadings_indices_top_n_of_selected_pc = np.argsort(loadings_of_selected_pc)[::-1][:top_n] # argsort returns the indices
        loadings_values_top_n_of_selected_pc = loadings_of_selected_pc[loadings_indices_top_n_of_selected_pc]
        top_features_selected_pc = metabolite_df.iloc[:, loadings_indices_top_n_of_selected_pc].columns.tolist()
        names_loadings_top_features = pd.DataFrame({'feature_name': top_features_selected_pc, 
                                                    'loading': list(loadings_values_top_n_of_selected_pc)},
                                                    columns = ['feature_name', 'loading'])
        
        print("Here are the metabolite names with the top {0} absolute loadings on PC{1}".format(top_n, selected_pc))
        return names_loadings_top_features
