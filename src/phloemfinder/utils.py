#!/usr/bin/env python3 

import os
from warnings import WarningMessage
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import balanced_accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, ConfusionMatrixDisplay

def compute_metrics_classification(y_predictions, y_trues, positive_class):
    '''
    Compute a series of metrics for classification tasks

    Util function designed to work downstream of the search for the best model. 
    Will compute the following metrics:
      - balanced accuracy
      - precision
      - recall
      - f1 score

    Parameters
    ----------
    y_predictions: list
      List of class predictions. 
    y_trues: list
      List of the true values (from the test set)
    positive_class: str
      The name of the positive class for calculation of true positives, true negatives, etc. 

    Returns
    -------
    model_metrics_df: `pandas.core.frame.DataFrame`
      Dataframe with the balanced accuracy, precision, recall and f1 score calculated. 

    See also
    --------
    https://scikit-learn.org/stable/modules/model_evaluation.html
    '''
    
    balanced_accuracy = balanced_accuracy_score(y_pred=y_predictions, y_true=y_trues)
    precision = precision_score(y_pred=y_predictions, y_true=y_trues, pos_label=positive_class)
    recall = recall_score(y_pred=y_predictions, y_true=y_trues, pos_label=positive_class)
    f1 = f1_score(y_true=y_trues, y_pred=y_predictions, pos_label=positive_class)
    
    model_metrics_dict = {"balanced_accuracy": balanced_accuracy, "precision": precision, "recall": recall, "f1 score": f1}
    model_metrics_df = pd.DataFrame.from_dict(model_metrics_dict, orient="index", columns=["value"])
    model_metrics_df_rounded = model_metrics_df.round(3)

    return model_metrics_df_rounded 

def plot_confusion_matrix(y_predictions, y_trues):
    '''
    Plot confusion matrix

    Parameters
    ----------
    y_predictions: list
      List of class predictions. 
    y_trues: list
      List of the true values (from the test set)
    positive_class: str
      The name of the positive class for calculation of true positives, true negatives, etc. 

    Returns
    -------
    model_metrics_df: `pandas.core.frame.DataFrame`
      Dataframe with the balanced accuracy, precision, recall and f1 score calculated. 

    See also
    --------
    https://scikit-learn.org/stable/modules/model_evaluation.html
    '''
    cm = confusion_matrix(y_true=y_trues, y_pred=y_predictions)
    disp = ConfusionMatrixDisplay(
        confusion_matrix=cm)
    disp.plot()
    plt.show()


