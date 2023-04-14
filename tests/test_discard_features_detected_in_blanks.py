from phloemfinder.metabolome_analysis import MetaboliteAnalysis


def test_discard_features_detected_in_blanks():
    """
    Test the discard_features_detected_in_blanks() function
    """
    met = MetaboliteAnalysis(
        metabolome_csv="tests/metabolome_test_datasets/metabolome_filter_blanks_test_data.csv", # contains two metabolites
        metabolome_feature_id_col='feature_id')
    
    # Apply function and get number of metabolites
    met.validate_input_metabolome_df()
    met.discard_features_detected_in_blanks(blank_sample_contains='blank')
    actual_nb_of_metabolites = met.metabolome.shape[1]

    assert actual_nb_of_metabolites == 1, "Number of discarded metabolites in test file is incorrect and should be equal to 1. Please review the discard_features_detected_in_blanks() function."
    