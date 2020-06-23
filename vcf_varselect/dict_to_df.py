#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'Danyang Li'
# __date__ = 2020-04-02

import pandas as pd

def dict_to_df(var_dict):
    """
    Change variant dictionary from all samples to a dataframe
    Argument:
        var_dict: variant dictionary
    Return:
        df with sample ID as row and annotation information as columns

    """
    new_dict = {}
    for key in var_dict:
        for var in var_dict[key]:
            key_var = '_'.join([key, var])
            new_dict[key_var] = {}
            for ID in var_dict[key][var]:
                if ID != 'CSQ':
                    new_dict[key_var][ID] = var_dict[key][var][ID]
                    if type(new_dict[key_var][ID]) == list:
                        new_dict[key_var][ID] = ','.join(new_dict[key_var][ID])

                else:
                    for i in var_dict[key][var][ID]:
                        csq_id = ':'.join([ID, i])
                        new_dict[key_var][csq_id] = var_dict[key][var][ID][i]
                        if type(new_dict[key_var][csq_id]) == list:
                            new_dict[key_var][csq_id] = ','.join(new_dict[key_var][csq_id])

    df = pd.DataFrame.from_dict({i: new_dict[i] for i in new_dict.keys()}, orient='index')
    return df

