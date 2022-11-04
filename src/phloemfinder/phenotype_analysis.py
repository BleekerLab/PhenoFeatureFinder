#!/usr/bin/env python3 

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class PhenotypeAnalysis:
    '''
    A class to analyse data from developmental bioassays and group the samples in distict phenotypic classes. 
    - What does this class do? 

    Parameters
    ----------
    The constructor method __init__ takes one arguments:
    bioassay_csv: string
        A path to a .csv file with the bioassay count data.
        Shape of the dataframe is usually ...

    '''
    
    def __init__(
        self, 
        bioassay_csv):

        # Import bioassay dataframe 
        self.bioassay = pd.read_csv(bioassay_csv)



    def reshape_to_wide(
        self,
        sample_id='sample_id',
        grouping_variable='genotype',
        developmental_stages='stage',
        count_values='number',
        time='day'):
        
        '''
        Reshapes the dataframe from a long to a wide format 
        with the counts of each developmental stage in a seperate columns.
        
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

        | sample_id | genotype  | day   | stage         | count |
        |-----------|-----------|-------|---------------|-------|
        | mm_1      |   mm      | 5     | eggs          | 45    |
        | mm_1      |   mm      | 5     | first_instar  | 0     |
        | mm_1      |   mm      | 5     | second_instar | 0     |
        
        
        Example of a reshaped output dataframe

        | sample_id | genotype  | day   | eggs  | first_instar  | second instar | third_instar  |
        |-----------|-----------|-------|-------|---------------|---------------|---------------|
        | mm_1      |   mm      | 5     | 45    | 0             | 0             | 0             |
        | mm_1      |   mm      | 9     | NA    | 10            | 5             | 0             |
        | mm_1      |   mm      | 11    | NA    | 15            | 17            | 4             |
        '''

        
        # check if specified columns exist in dataframe
        if sample_id not in self.bioassay.columns:
            raise ValueError("The specified column with sample identifiers {0} is not present in your '{1}' file.".format(sample_id,os.path.basename(bioassay_csv)))
        else:
            self.sample_id = sample_id

        if grouping_variable not in self.bioassay.columns:
            raise ValueError("The specified column with grouping variable names {0} is not present in your '{1}' file.".format(grouping_variable,os.path.basename(bioassay_csv)))
        else:
            pass

        if developmental_stages not in self.bioassay.columns:
            raise ValueError("The specified column with developmental stages {0} is not present in your '{1}' file.".format(developmental_stages,os.path.basename(bioassay_csv)))
        else:
            pass

        if count_values not in self.bioassay.columns:
            raise ValueError("The specified column with values {0} is not present in your '{1}' file.".format(count_values,os.path.basename(bioassay_csv)))
        else:
            pass

        if time not in self.bioassay.columns:
            raise ValueError("The specified column with time values (e.g. days after infection) {0} is not present in your '{1}' file.".format(time,os.path.basename(bioassay_csv)))
        else:
            pass

        # reshape the dataframe to a wide format with one developmental stage per column
        self.bioassay = self.bioassay.pivot(index=[sample_id, grouping_variable, time], columns=developmental_stages, values=count_values)
        
    

    def total_final_stage(
        self,
        exuviea='exuviea',
        late_last_stage='late_fourth_instar',
        early_last_stage='early_fourth_instar',
        new_last_stage='fourth_instar',
        seperate_exuviea=True,
        late_last_stage_removed=True,
        early_last_stage_kept=True,
        remove_individual_stage_columns=True):
        
        '''
        Calculates the total number of nymphs developed to the final developmental stage per sample on each timepoint.
        This is used when nymphs in the (late) final nymph stage were removed after each counting moment and/or
        when exuviea and last instar stage nymphs were counted seperately.
        
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

        | sample_id | genotype  | day   | eggs  | ... | third_instar  | exuviea   | early_fourth_instar | late_fourth_instar |
        |-----------|-----------|-------|-------|-----|---------------|-----------|---------------------|--------------------|
        | mm_1      |   mm      | 5     | 45    | ... | 0             | 0         | 0                   | 0                  |
        | mm_1      |   mm      | 9     | NA    | ... | 0             | 1         | 5                   | 0                  |
        | mm_1      |   mm      | 11    | NA    | ... | 4             | 0         | 7                   | 4                  |
        
        
        Example of an output dataframe

        | sample_id | genotype  | day   | eggs  | first_instar  | second instar | third_instar  | fourth_instar |
        |-----------|-----------|-------|-------|---------------|---------------|---------------|---------------|
        | mm_1      |   mm      | 5     | 45    | 0             | 0             | 0             | 0             |
        | mm_1      |   mm      | 9     | NA    | 10            | 5             | 0             | 6
        | mm_1      |   mm      | 11    | NA    | 15            | 17            | 4             | 12
        '''

        # if the exuviea and late last stage nymphs were counted individualy but treated simmilarly,
        # the exuviea and late last stage nymphs should be summed before going further
        if seperate_exuviea == True:
            
            # check if specified columns with exuviea and late last stage are present
            if exuviea not in self.bioassay.columns:
                raise ValueError("The specified column with exuviea counts {0} is not present in your file.".format(exuviea))
            else:
                pass

            if late_last_stage not in self.bioassay.columns:
                raise ValueError("The specified column with late last stage counts {0} is not present in your file.".format(late_last_stage))
            else:
                pass

            # if specified columns are present, sum exuviea and late last stage
            self.bioassay['late_exuviea'] = self.bioassay[[exuviea, late_last_stage]].sum(axis=1)
        else:
            if late_last_stage not in self.bioassay.columns:
                raise ValueError("The specified column with late last stage counts {0} is not present in your file.".format(late_last_stage))
            else:
                pass
            self.bioassay['late_exuviea'] = self.bioassay[late_last_stage]

        # if the late last stage nymphs were removed after each count, the cumulative number should be used
        # when analysing the development to this stage
        if late_last_stage_removed == True:
            self.bioassay['late_exuviea'] = self.bioassay.groupby([self.sample_id])['late_exuviea'].cumsum()
        else:
            pass

        # if the early last stage nymphs were kept for further development after counting, the early and late last stage
        # nymphs should be combined for the total number of last stage nymphs
        if early_last_stage_kept == True:

            # check if specified column with early last stage counts is present 
            if early_last_stage not in self.bioassay.columns:
                raise ValueError("The specified column with early last stage counts {0} is not present in your file.".format(early_last_stage))
            else:
                pass

            # check if specified column name for total last instar numbers is not yet present 
            if new_last_stage in self.bioassay.columns:
                raise ValueError("The specified column name for total last instar numbers {0} already exists in your file.".format(new_last_stage))
            else:
                pass

            # if only the early last stage column is present, calculate total last stage nymphs
            self.bioassay[new_last_stage] = self.bioassay[['late_exuviea', early_last_stage]].sum(axis=1)
        else:

            # check if specified column name for total last instar numbers is not yet present 
            if new_last_stage in self.bioassay.columns:
                raise ValueError("The specified column name for total last instar numbers {0} already exists in your file.".format(new_last_stage))
            else:
                pass
            self.bioassay[new_last_stage] = self.bioassay['late_exuviea']

        # cleaning up unwanted columns
        if remove_individual_stage_columns == True:
            self.bioassay = self.bioassay.drop(columns=['late_exuviea'])
            self.bioassay = self.bioassay.drop([exuviea, late_last_stage, early_last_stage], axis=1)
        else:
            pass


    def cumulative_counts(
        self,
        n_developmental_stages=4,
        first_stage='first_instar',
        second_stage='second_instar',
        third_stage='third_instar',
        fourth_stage='fourth_instar',
        fifth_stage='fifth_instar',
        sixth_stage='sixth_instar'):
        
        '''
        Calculates the total number of nymphs developed to or past each stage on each timepoint.
        If nymphs in the (late) final nymph stage were removed after each counting moment and/or
        when exuviea and/or early and late last instar stage nymphs were counted seperately, 
        total_last_stage() should be used first.
        
        Parameters
        ----------
        n_developmental_stages: integer, default=4
            The number of developmental stage swhich were recorded seperately. 
            Can range from 2 till 6.
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
        '''
        
        if n_developmental_stages == 1:
            raise ValueError("n_developmental_stages is set to 1. If only 1 developmental stage was counted and the cumulative number is needed, try using total_final_stage() instead.")
            
        
        elif n_developmental_stages == 2:

            # check if specified columns with counts per stage are present 
            if first_stage not in self.bioassay.columns:
                raise ValueError("The specified column with first stage counts {0} is not present in your file.".format(first_stage))
            else:
                pass
            if second_stage not in self.bioassay.columns:
                raise ValueError("The specified column with second stage counts {0} is not present in your file.".format(second_stage))
            else:
                pass
            

            # if all specified columns are present:
            # calculate the number of nymphs developed to or past first instar stage per timepoint
            self.bioassay[first_stage] = self.bioassay[[first_stage, second_stage]].sum(axis=1)

        
        elif n_developmental_stages == 3:
            
            # check if specified columns with counts per stage are present 
            if first_stage not in self.bioassay.columns:
                raise ValueError("The specified column with first stage counts {0} is not present in your file.".format(first_stage))
            else:
                pass
            if second_stage not in self.bioassay.columns:
                raise ValueError("The specified column with second stage counts {0} is not present in your file.".format(second_stage))
            else:
                pass
            if third_stage not in self.bioassay.columns:
                raise ValueError("The specified column with third stage counts {0} is not present in your file.".format(third_stage))
            else:
                pass

            # if all specified columns are present:
            # calculate the number of nymphs developed to or past first instar stage per timepoint
            self.bioassay[first_stage] = self.bioassay[[first_stage, second_stage, third_stage]].sum(axis=1)

            # calculate the number of nymphs developed to or past second instar stage per timepoint
            self.bioassay[second_stage] = self.bioassay[[second_stage, third_stage]].sum(axis=1)
        
        
        elif n_developmental_stages == 4:

            # check if specified columns with counts per stage are present 
            if first_stage not in self.bioassay.columns:
                raise ValueError("The specified column with first stage counts {0} is not present in your file.".format(first_stage))
            else:
                pass
            if second_stage not in self.bioassay.columns:
                raise ValueError("The specified column with second stage counts {0} is not present in your file.".format(second_stage))
            else:
                pass
            if third_stage not in self.bioassay.columns:
                raise ValueError("The specified column with third stage counts {0} is not present in your file.".format(third_stage))
            else:
                pass
            if fourth_stage not in self.bioassay.columns:
                raise ValueError("The specified column with fourth stage counts {0} is not present in your file.".format(fourth_stage))
            else:
                pass

            # if all specified columns are present:
            # calculate the number of nymphs developed to or past first instar stage per timepoint
            self.bioassay[first_stage] = self.bioassay[[first_stage, second_stage, third_stage, fourth_stage]].sum(axis=1)

            # calculate the number of nymphs developed to or past second instar stage per timepoint
            self.bioassay[second_stage] = self.bioassay[[second_stage, third_stage, fourth_stage]].sum(axis=1)
            
            # calculate the number of nymphs developed to or past third instar stage per timepoint
            self.bioassay[third_stage] = self.bioassay[[third_stage, fourth_stage]].sum(axis=1)


        elif n_developmental_stages == 5:

            # check if specified columns with counts per stage are present 
            if first_stage not in self.bioassay.columns:
                raise ValueError("The specified column with first stage counts {0} is not present in your file.".format(first_stage))
            else:
                pass
            if second_stage not in self.bioassay.columns:
                raise ValueError("The specified column with second stage counts {0} is not present in your file.".format(second_stage))
            else:
                pass
            if third_stage not in self.bioassay.columns:
                raise ValueError("The specified column with third stage counts {0} is not present in your file.".format(third_stage))
            else:
                pass
            if fourth_stage not in self.bioassay.columns:
                raise ValueError("The specified column with fourth stage counts {0} is not present in your file.".format(fourth_stage))
            else:
                pass
            if fifth_stage not in self.bioassay.columns:
                raise ValueError("The specified column with fifth stage counts {0} is not present in your file.".format(fifth_stage))
            else:
                pass

            # if all specified columns are present:
            # calculate the number of nymphs developed to or past first instar stage per timepoint
            self.bioassay[first_stage] = self.bioassay[[first_stage, second_stage, third_stage, fourth_stage, fifth_stage]].sum(axis=1)

            # calculate the number of nymphs developed to or past second instar stage per timepoint
            self.bioassay[second_stage] = self.bioassay[[second_stage, third_stage, fourth_stage, fifth_stage]].sum(axis=1)
            
            # calculate the number of nymphs developed to or past third instar stage per timepoint
            self.bioassay[third_stage] = self.bioassay[[third_stage, fourth_stage, fifth_stage]].sum(axis=1)
            
            # calculate the number of nymphs developed to or past fourth instar stage per timepoint
            self.bioassay[fourth_stage] = self.bioassay[[fourth_stage, fifth_stage]].sum(axis=1)


        elif n_developmental_stages == 6:

            # check if specified columns with counts per stage are present 
            if first_stage not in self.bioassay.columns:
                raise ValueError("The specified column with first stage counts {0} is not present in your file.".format(first_stage))
            else:
                pass
            if second_stage not in self.bioassay.columns:
                raise ValueError("The specified column with second stage counts {0} is not present in your file.".format(second_stage))
            else:
                pass
            if third_stage not in self.bioassay.columns:
                raise ValueError("The specified column with third stage counts {0} is not present in your file.".format(third_stage))
            else:
                pass
            if fourth_stage not in self.bioassay.columns:
                raise ValueError("The specified column with fourth stage counts {0} is not present in your file.".format(fourth_stage))
            else:
                pass
            if fifth_stage not in self.bioassay.columns:
                raise ValueError("The specified column with fifth stage counts {0} is not present in your file.".format(fifth_stage))
            else:
                pass
            if sixth_stage not in self.bioassay.columns:
                raise ValueError("The specified column with sixth stage counts {0} is not present in your file.".format(sixth_stage))
            else:
                pass

            # if all specified columns are present:
            # calculate the number of nymphs developed to or past first instar stage per timepoint
            self.bioassay[first_stage] = self.bioassay[[first_stage, second_stage, third_stage, fourth_stage, fifth_stage, sixth_stage]].sum(axis=1)

            # calculate the number of nymphs developed to or past second instar stage per timepoint
            self.bioassay[second_stage] = self.bioassay[[second_stage, third_stage, fourth_stage, fifth_stage, sixth_stage]].sum(axis=1)
            
            # calculate the number of nymphs developed to or past third instar stage per timepoint
            self.bioassay[third_stage] = self.bioassay[[third_stage, fourth_stage, fifth_stage, sixth_stage]].sum(axis=1)
            
            # calculate the number of nymphs developed to or past fourth instar stage per timepoint
            self.bioassay[fourth_stage] = self.bioassay[[fourth_stage, fifth_stage, sixth_stage]].sum(axis=1)
            
            # calculate the number of nymphs developed to or past fifth instar stage per timepoint
            self.bioassay[fifth_stage] = self.bioassay[[fifth_stage, sixth_stage]].sum(axis=1)