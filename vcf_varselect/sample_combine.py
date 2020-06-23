#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'Danyang Li'
# __date__ = 2020-04-02

import os
from vcf_varselect.variant_selection import VariantSelection
from vcf_varselect.match_gene import match_gene
from vcf_varselect.dict_to_df import dict_to_df


def sample_combine(dir=None, innerfreqfile=None, genefile=None, genderfile=None,
                   FILTER=None, DP=None, QD=None, MQ=None,
                   KG=None, EXAC=None, GNOMAD=None, SWEGEN=None,
                   criteria=[]
                   ):
    """
    collect all samples' selected rare damaging variants to a dataframe
    Args:
        dir: directory where all vcf file are
        innerfreqfile: json file of variant inner-freq from all samples
        genefile: disorder related gene list
        genderfile: male samples list

    return:
        df with all sample IDs as rows and selected ndd variants annotation information as columns
    """
    total = {}
    if not dir:
        raise IOError("Please input file directory.")
    os.chdir(dir)
    for filename in os.listdir(dir):
        file_name, file_extension = os.path.splitext(filename)
        if file_extension == 'vcf.gz' or file_extension == '.vcf':
            vcf = VariantSelection(infile=filename)
            damage, lof, mis = vcf.comb_selection(FILTER=FILTER, DP=DP, QD=QD, MQ=MQ,
                                        innerfreqfile=innerfreqfile, KG=KG, EXAC=EXAC, GNOMAD=GNOMAD, SWEGEN=SWEGEN,
                                        criteria=criteria)
            gene_var = match_gene(damage, genefile=genefile, genderfile=genderfile)
            total.update(gene_var)

    df_total = dict_to_df(total)

    return df_total
