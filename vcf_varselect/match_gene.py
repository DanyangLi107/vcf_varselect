#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'Danyang Li'
# __date__ = 2020-04-02

def match_gene(var_dict, genefile=None, genderfile=None):
    """
    Select variants of disorders related genes
    Arguments:
        var_dict: selected variants dictionary;
        genefile: disorder related gene list;
        genderfile: male sample ID
    Return:
        nested dictionary of disorder related variants

    """

    if not genefile:
        raise IOError("Please input disorders gene list.")
    if not genderfile:
        raise IOError("Please input individual gender information")

    # genes in both x and y chromosome
    psudolist = [
        'ENSG00000197976', 'ENSG00000196433', 'ENSG00000169093', 'ENSG00000002586', 'ENSG00000205755',
        'ENSG00000198223', 'ENSG00000169084', 'ENSG00000178605', 'ENSG00000185291', 'ENSG00000182162',
        'ENSG00000182378', 'ENSG00000167393', 'ENSG00000185960', 'ENSG00000169100', 'ENSG00000124343',
        'ENSG00000214717', 'ENSG00000124334', 'ENSG00000168939', 'ENSG00000124333', 'ENSG00000182484'
    ]

    #####generate genefile dictionary######
    gene_file = open(genefile, mode='r')
    gene_dict = {}
    for line in gene_file:
        gene_dict[line.split(',')[1]] = [
            line.strip().split(',')[0], line.strip().split(',')[2]
        ]

    #####generate gender list#####
    malefile = open(genderfile, mode='r')
    male_list = [line.strip() for line in malefile]

    ##divided candidate disorder genes into different groups based on inheritance
    gene_group = {'gene_ar': [], 'gene_xr': [], 'gene_ad': [], 'gene_xd': []}
    for key in gene_dict:
        if 'AR' in gene_dict[key][1] and 'AD' not in gene_dict[key][1]:
            gene_group['gene_ar'].append(key)
        if 'XR' in gene_dict[key][1] and 'XD' not in gene_dict[key][1]:
            gene_group['gene_xr'].append(key)
        if 'AD' in gene_dict[key][1]:
            gene_group['gene_ad'].append(key)
        if 'XD' in gene_dict[key][1]:
            gene_group['gene_xd'].append(key)

            ######select variants#####
    sel_list = []
    genevar_dict = {}
    for sample in var_dict:
        for key in var_dict[sample]:

            #####autosome variant#####
            if key.split(':')[0] != 'X' and key.split(':')[0] != 'Y':
                ##autosome dominant variant
                if var_dict[sample][key]['GT'] == '0/1' or var_dict[sample][key]['GT'] == '1/0':
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if i in gene_group['gene_ad']:
                            sel_list.append(key)

                # autosome recessive variant
                elif var_dict[sample][key]['GT'] == '1/1':
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if i in gene_group['gene_ad'] or i in gene_group['gene_ar']:
                            sel_list.append(key)

            #####x chromosome variant#####
            elif key.split(':')[0] == 'X':

                ##first on female individual
                # dominant
                if sample not in male_list and \
                        (var_dict[sample][key]['GT'] == '0/1' or var_dict[sample][key]['GT'] == '1/0'):
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if i in gene_group['gene_xd']:
                            sel_list.append(key)

                # recessive
                elif sample not in male_list and \
                        var_dict[sample][key]['GT'] == '1/1':
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if i in gene_group['gene_xd'] or i in gene_group['gene_xr']:
                            sel_list.append(key)

                ##next male individual
                # dominant psudogenes variant
                elif sample in male_list and \
                        (var_dict[sample][key]['GT'] == '0/1' or var_dict[sample][key]['GT'] == '1/0'):
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if i in gene_group['gene_xd'] and i in psudolist:
                            sel_list.append(key)

                # recessive psudogenes variant
                elif sample in male_list and \
                        var_dict[sample][key]['GT'] == '1/1':
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if (i in gene_group['gene_xd'] or i in gene_group['gene_xr']) and i in psudolist:
                            sel_list.append(key)

                # non-psudogene variant
                elif sample in male_list and \
                        (var_dict[sample][key]['GT'] == './1' or var_dict[sample][key]['GT'] == '1/.'):
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if (i in gene_group['gene_xd'] or i in gene_group['gene_xr']) and i not in psudolist:
                            sel_list.append(key)

            #####y chromosome#####
            elif key.split(':')[0] == 'Y':
                # dominant psudogene variant
                if sample in male_list and \
                        (var_dict[sample][key]['GT'] == '0/1' or var_dict[sample][key]['GT'] == '1/0'):

                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if i in gene_group['gene_xd'] and i in psudolist:
                            sel_list.append(key)

                # recessive psudogene variant
                elif sample in male_list and \
                        var_dict[sample][key]['GT'] == '1/1':
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if (i in gene_group['gene_xd'] or i in gene_group['gene_xr']) and i in psudolist:
                            sel_list.append(key)

                # non-psudogene variant
                elif sample in male_list and \
                        (var_dict[sample][key]['GT'] == './1' or var_dict[sample][key]['GT'] == '1/.'):
                    for i in var_dict[sample][key]['CSQ']['Gene']:
                        if (i in gene_group['gene_xd'] or i in gene_group['gene_xr']) and i not in psudolist:
                            sel_list.append(key)

        ##generate dictionary
        genevar_dict[sample] = {i: var_dict[sample][i] for i in set(sel_list)}

    return genevar_dict
