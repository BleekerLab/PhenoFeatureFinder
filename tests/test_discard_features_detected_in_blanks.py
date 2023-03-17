from phloemfinder.metabolome_analysis import MetaboliteAnalysis



def test_discard_features_detected_in_blanks():
    """
    Test the discard_features_detected_in_blanks() function
    """
    met = MetaboliteAnalysis(
        metabolome_csv="metabolome_test_data.csv", 
        metabolome_feature_id_col='feature_id')
    met.validate_input_metabolome_df()
    met.discard_features_detected_in_blanks(blank_sample_contains='blank')
    actual_nb_of_metabolites = met.metabolome.shape[1]

    assert actual_nb_of_metabolites == expected_nb_of_metabolites, "Number of discarded metabolites in blanks is incorrect. Please review the discard_features_detected_in_blanks() function."
    