#!/usr/bin/env python3 

from enum import auto
import os
from warnings import WarningMessage
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys

from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance

from autosklearn.classification import AutoSklearnClassifier
from autosklearn.metrics import balanced_accuracy, precision, recall, f1

from utils import compute_performance_metrics

class MetabolitePhenotypeFeatureSelection:
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


    Class attributes
    ----------
    metabolome_validated: boolean, default=False
      Is the metabolome file valid for Machine Learning?
    
    phenotype_validated: boolean, default=False
      Is the phenotype file valid for Machine Learning?

    baseline_performance: float, default=None
      The baseline performance computed with get_baseline_performance() i.e. using a simple Random Forest model. 
      Search for the best ML model using search_best_model() should perform better than this baseline performance. 

    best_ensemble_models_searched: boolean, default=False
      Is the search for best ensemble model using auto-sklearn already performed?

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
                raise IndexError("Values for sample identifiers have to be unique. Check your",phenotype_sample_id,"column.")

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

        Examples
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
    def get_baseline_performance(self, kfold=5, scoring_metric='balanced_accuracy'):
        '''
        Takes the phenotype and metabolome dataset and compute a simplistic Random Forest analysis. 
        This will give a base performance for a Machine Learning model that has then to be optimised using autosklearn

        k-fold cross-validation is performed to accomodate for split effects on small datasets. 

        Parameters
        ----------
        fkold: int, default=5
          Cross-validation strategy. 
          Default is to use a 5-fold cross-validation. 
        scoring_metric: string, default="balanced_accuracy"
          A valid scoring value. To get a complete list, type:
          >> from sklearn.metrics import SCORERS 
          >> sorted(SCORERS.keys()) 
          balanced accuracy is the average of recall obtained on each class. 

        Returns
        -------
        self: object
          Object with baseline_performance attribute
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

        clf = RandomForestClassifier(n_estimators=1000, random_state=42)
        scores = cross_val_score(clf, X, y, scoring=scoring_metric, cv=kfold)
        average_scores = np.average(scores).round(3)

        print("Performing a simple Random Forest model training")
        print("N samples: {0}".format(str(X.shape[0])))
        print("N features: {0}".format(str(X.shape[1])))
        print("Average {0} score of the default model is: {1} %".format(scoring_metric, average_scores))

        self.baseline_performance = average_scores


    #FIXME
    # def get_base_confusion_matrix(self):
        '''
        Plot the confusion matrix coming from the simplistic baseline Random Forest model.
        
        A confusion matrix 
        Returns
        -------
        The confusion matrix plot. 
        '''

    #TODO: make decorator function to check arguments
    #See -> https://www.pythonforthelab.com/blog/how-to-use-decorators-to-validate-input/
    def search_best_model_and_calculate_feature_importances(
        self, 
        time_left_for_this_task=300,  
        kfolds=5, 
        class_of_interest='resistant',
        train_size=0.7,
        n_permutations=10, 
        random_state=123):
        '''
        Wrapper function around the AutoSklearnClassifier() autosklearn function. 
        Helpful since the best model will be stored as an class attribute. 
        The best model is actually an ensemble of models 

        Input function arguments have to comply with both names and value type of the 
        original AutoSklearnClassifier() function 
        e.g. 'time_left_for_this_task' has to be named as such and takes integer values for time in seconds 
        https://automl.github.io/auto-sklearn/master/api.html

        The "balanced_accuracy" from autosklearn.metrics will be used as scorer to find the best model.
        https://automl.github.io/auto-sklearn/master/api.html#built-in-metrics

        A resampling strategy called "cross-validation" will be performed to increase 
        the model generalisation performance. 
        Depending on the number of samples, different data subsets will be produced:
          - if n_samples > 30, then a train and a test set are produced. 
            The training set will be subjected to k-fold cv which means that, for each k-fold: 
              part of the training data will be used to train a model. 
              the left out fold will be used as a validation set to compute a model score. 
            Finally, the complete cv model will be evaluated on the never seen test set. 
            This avoids "test data leakage" during the training step of the model.
          - if n_samples =< 30, then only a train set will be produced. 
            Scores will be produced during cv. 
            This is prone to overfitting (overestimated performance on the model) since
            test data "leaks" during model training.  

        Parameters
        ----------
        time_left_for_this_task: int, default=300
          Time (in seconds) allocated to search for the best model. 
        
        metric: string, default='balanced_accuracy'
          Scorer metric. One value such as 'accuracy' or 'balanced_accuracy'
          See https://automl.github.io/auto-sklearn/master/api.html#built-in-metrics
        
        kfolds: integer, default=5
          Number of folds for the stratified K-Folds cross-validation strategy
          Has to be comprised between 3 and 10 i.e. 3 <= kfolds =< 10
          See https://scikit-learn.org/stable/modules/cross_validation.html
        
        class_of_interest: str, default='resistant'
          The name of the class of interest also called "positive class".
          This class will be used to calculate recall_score and precision_score. 
          Recall score = TP / (TP + FN) with TP: true positives and FN: false negatives.
          Precision score = TP / (TP + FP) with TP: true positives and FP: false positives. 
        
        train_size: float, default=0.7
          If float, should be between 0.5 and 1.0 and represent the proportion of the dataset to include in the train split. If int, represents the absolute number of train samples. If None, the value is automatically set to the complement of the test size.
        
        random_state: int, default=123
          Controls both the randomness of the train/test split  samples used when building trees (if bootstrap=True) and the sampling of the features to consider when looking for the best split at each node (if max_features < n_features). See Glossary for details.
          You can change this value several times to see how it affects the best ensemble model performance

        Returns
        ------
        self: object
          The object with a new attribute called 'automated_ml_obj' that contains the fitted AutoSklearnClassifier 

        Example
        -------
        >> fs = MetabolitePhenotypeFeatureSelection(
        >>        metabolome_csv="clean_metabolome.csv", 
        >>        phenotype_csv="phenotypes_test_data.csv", 
        >>        phenotype_sample_id='sample')
        >> fs.validate_input_phenotype_df()
        >> fs.validate_input_metabolome_df()
        >> fs.get_baseline_performance(scoring_metric='accuracy')
        >> fs.search_best_model(time_left_for_this_task=120)
        '''

        X = self.metabolome.transpose().to_numpy(dtype='float64')
        y = self.phenotype.values.ravel()
        
        # Verify input arguments
        try:
          3 <= kfolds <= 10
        except:
          raise ValueError('kfolds argument has to be comprised between 3 and 10')
        
        try:
          0.5 < train_size < 0.9
        except:
          raise ValueError('train_size has to be comprised between 0.5 and 0.9')

        try:
          class_of_interest in set(y)
        except:
          raise ValueError('The class_of_interest value "{0}" has to be in the phenotype labels {1}'.format(class_of_interest, set(y)))
    


        ### CV resampling strategy depends on number of samples ###
        n_samples = int(len(y))
        # Only a train and test set if n_samples < 30 
        # Less than 30 samples is not enough to create distinct train/validation/test sets
        if n_samples > 30:
          X_train, X_test, y_train, y_test = train_test_split(
            X, y, 
            train_size=train_size, 
            random_state=random_state, 
            stratify=y)
          automl = AutoSklearnClassifier(
              time_left_for_this_task=time_left_for_this_task, 
              metric=balanced_accuracy,
              scoring_functions=[precision, recall, f1],
              seed=random_state,
              resampling_strategy='cv', 
              resampling_strategy_arguments={'train_size': train_size, 'folds': kfolds}
              )
          automl.fit(X_train, y_train, dataset_name="metabopheno")
          automl.refit(X_train, y_train) # to fit on the whole train dataset (not only on individual k-folds)
          predictions = automl.predict(X_test)
          compute_performance_metrics(predictions=predictions, y_test=y_test, positive_label=class_of_interest)
        else:
          # All data is used to train the model. 
          # This is prone to overfitting
          X_train = X
          y_train = y 
          automl = AutoSklearnClassifier(
              time_left_for_this_task=time_left_for_this_task, 
              metric=balanced_accuracy,
              seed=random_state,
              scoring_functions=[precision, recall, f1],
              resampling_strategy='cv', 
              resampling_strategy_arguments={'train_size': train_size, 'folds': kfolds}
              )
          automl.fit(X, y, dataset_name="metabopheno")
          predictions = automl.predict(X_test)
          compute_performance_metrics(predictions=predictions, y_test=y_test, positive_label=class_of_interest)
          print("============ Warning =========")
          print("This is prone to overfitting since no external test dataset can be made from less than 30 samples")
        

        ### Get feature importances ###
        # This has to be done from the same test/train split if n_samples > 30
        # See https://scikit-learn.org/stable/modules/permutation_importance.html 
        if n_samples > 30:
            feature_importances = permutation_importance(
            automl, 
            X=X_train, 
            y=y_train,
            scoring=balanced_accuracy,
            n_repeats=n_permutations, 
            random_state=random_state)
        feature_importances_df = pd.DataFrame.from_dict(feature_importances["importances"])
        feature_importances_df.set_index(self.metabolome.index.values, inplace=True)
        feature_importances_df.columns = ["permutation_" + str(i+1) for i in range(n_permutations)]

        self.automated_ml_obj = automl
        self.best_ensemble_models_searched = True
        self.feature_importances = feature_importances_df



    def extract_feature_importances(self, ):
      '''
      Extracts feature importances from a fitted ensemble model obtained using the search_best_model() method

      Parameters
      ---------- 
      '''
      if self.best_ensemble_models_searched == False:
        raise ValueError("Please search best ML model using the search_best")
        

        

