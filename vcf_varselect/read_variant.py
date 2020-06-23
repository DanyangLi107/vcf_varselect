#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'Danyang Li'
# __date__ = 2020-04-02

def read_variant(line, parser):
    """
    Yield the variant in the right format.
    Arguments:
        line (str): A string representing a variant line in the vcf
        parser: A MetadataParser object
    Return:
        variant (dict): A dictionary with the variant information.
        dictionary key: variant information, format: chromosome:position:rsID:reference:alternative;
        value: {'QUAL': [], 'FILTER':[], 'GT': [], 'info_ID1':[], 'info_ID2':[],...}
        special parsing of vep annotation: a sub-dictionary with information separated by '|'

    """

    vcf_header = parser.header
    variant = {}

    variant_line = line.rstrip().split('\t')

    if len(vcf_header) != len(variant_line):
        raise SyntaxError("One of the variant lines is malformed: {0}".format(
            line
        ))

    # key information: chromosome, position, rsID, reference allele, alternative allele
    key = ':'.join(
        [variant_line[0], variant_line[1], variant_line[2], variant_line[3], variant_line[4]]
    )
    variant[key] = {}

    variant[key]['QUAL'] = variant_line[5]
    variant[key]['FILTER'] = variant_line[6]
    for i in range(len(variant_line[8].split(':'))):
        if variant_line[8].split(':')[i] == 'GT':
            variant[key]['GT'] = variant_line[9].split(':')[i]

    ##### INFO information #####
    for info in variant_line[7].split(';'):
        info = info.split('=')
        if len(info) > 1:
            if ',' in info[1]:
                variant[key][info[0]] = info[1].split(',')
            else:
                variant[key][info[0]] = [info[1]]
        else:
            variant[key][info[0]] = []

    ##### VEP ANNOTATIONS #####
    if 'CSQ' in variant[key]:
        vep_columns = parser.vep_columns
        vep_total = variant[key]['CSQ']
        variant[key]['CSQ'] = {i: [] for i in vep_columns}
        for vep in vep_total:
            vep_list = [[i] for i in vep.split('|')]
            vep_dict = dict(zip(vep_columns, vep_list))
            variant[key]['CSQ'] = {i: variant[key]['CSQ'][i] + vep_dict[i] for i in vep_dict.keys()}

    return variant