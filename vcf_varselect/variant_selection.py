#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'Danyang Li'
# __date__ = 2020-04-02

import os
import gzip
from codecs import open, getreader
import json

from vcf_varselect.metadata_parser import MetadataParser
from vcf_varselect.read_variant import read_variant


class VariantSelection(object):
    """
    Change vcf file to dictionary, and select variants such as good quality variants, rare variants,
    and damging variants including loss-of-function and missense.
    Argument:
        vcf file
    Return:
        nested variant dictionary {sample:{variant1:{'QUALITY':"", 'FILTER':"", 'GT':"", 'infoID1':[], 'infoID2':[],...}}}

    """

    def __init__(self, infile=None):
        super(VariantSelection, self).__init__()

        if infile == None:
            raise IOError("Please input a file.")
        else:

            file_name, file_extension = os.path.splitext(infile)

            if file_extension == '.gz':
                self.vcf = getreader('utf-8')(gzip.open(infile), errors='replace')
            elif file_extension == '.vcf':
                self.vcf = open(infile, mode='r', encoding='utf-8', errors='replace')
            else:
                raise IOError("File is not in a supported format!\n"
                              "Please use correct ending(.vcf or .vcf.gz)")
            self.next_line = self.vcf.readline().rstrip()
            self.metadata = MetadataParser()

            while self.next_line.startswith('#'):
                if self.next_line.startswith('##fileformat') or self.next_line.startswith('##FILTER') \
                        or self.next_line.startswith('##FORMAT') or self.next_line.startswith('##INFO'):
                    self.metadata.read_metadata(self.next_line)
                elif self.next_line.startswith('#CHROM'):
                    self.metadata.read_header(self.next_line)
                    self.sample = self.next_line.split('\t')[9]
                self.next_line = self.vcf.readline().rstrip()

            self.variant = {}
            self.variant[self.sample] = {}

            while not self.next_line.startswith('#') and len(self.next_line.split('\t')) == 10:
                self.variant[self.sample].update(read_variant(
                    line=self.next_line,
                    parser=self.metadata
                ))
                self.next_line = self.vcf.readline().rstrip()

            self.header = self.metadata.header

            self.id_dict = self.metadata.id_dict
            self.vep_columns = self.metadata.vep_columns

    def __iter__(self):
        return iter(self.__dict__.items())

    def __getitem__(self, item):
        return getattr(self, str(item))

    def quality_selection(self, FILTER=None, DP=None, QD=None, MQ=None):
        """
        Select variants with good quality based on selected criteria (FILTER, DP, QD, MQ)
        Arguments:
            FILTER, DP, QD, MQ: variants quality threshold
        Return:
            nested dictionary of selected variants

        """
        quality = {}
        if FILTER:
            filt = [
                key for key in self.variant[self.sample]
                if 'FILTER' in self.variant[self.sample][key] and self.variant[self.sample][key]['FILTER'] != FILTER
            ]  # PASS
        else:
            filt = []

        if DP:
            dp = [
                key for key in self.variant[self.sample]
                if
                'DP' in self.variant[self.sample][key] and float(''.join(self.variant[self.sample][key]['DP'])) < float(DP)
            ]  # 10.0
        else:
            dp = []

        if QD:
            qd = [
                key for key in self.variant[self.sample]
                if
                'QD' in self.variant[self.sample][key] and float(''.join(self.variant[self.sample][key]['QD'])) < float(QD)
            ]  # 2.0
        else:
            qd = []

        if MQ:
            mq = [
                key for key in self.variant[self.sample]
                if 'MQ' in self.variant[self.sample][key] and float(''.join(self.variant[self.sample][key]['MQ'])) < float(MQ)
            ]  # 40.0
        else:
            mq = []

        quality[self.sample] = {
            i: self.variant[self.sample][i]
            for i in (set(self.variant[self.sample].keys()) - set(filt + dp + qd + mq))
        }

        return quality

    def freq_selection(self, innerfreqfile=None, KG=None, EXAC=None, GNOMAD=None, SWEGEN=None):
        """
        Select rare variants based on selected databases (1000G, EXAC, GNOMAD, SWEGEN)
        Arguments:
            KG, EXAC, GNOMAD, SWEGEN: rare variants frequency threshold;
            innerfreqfile: json file of variant frequency in the samples

        """
        freq = {}
        if KG:
            kg = [
                key for key in self.variant[self.sample]
                if '1000GAF' in self.variant[self.sample][key] and
                   float(''.join(self.variant[self.sample][key]['1000GAF'])) > float(KG)
            ]
        else:
            kg = []

        if EXAC:
            exac = [
                key for key in self.variant[self.sample]
                if 'EXACAF' in self.variant[self.sample][key] and
                   float(''.join(self.variant[self.sample][key]['EXACAF'])) > float(EXAC)
            ]
        else:
            exac = []

        if GNOMAD:
            gnomad = [
                key for key in self.variant[self.sample]
                for i in self.variant[self.sample][key]['CSQ']['gnomAD_AF']
                if i != "" and float(i) > float(GNOMAD)
            ]
        else:
            gnomad = []

        if SWEGEN:
            swegen = [
                key for key in self.variant[self.sample]
                if 'SWEGENAF' in self.variant[self.sample][key] and
                   float(''.join(self.variant[self.sample][key]['SWEGENAF'])) > float(SWEGEN)
            ]
        else:
            swegen = []

        if innerfreqfile:
            innerfreq_file = open(innerfreqfile)
            innerfreq = json.load(innerfreq_file)
            inner = [
                key for key in self.variant[self.sample]
                if key in innerfreq and innerfreq[key] > float(0.01)
            ]
        else:
            inner = []

        freq[self.sample] = {
            i: self.variant[self.sample][i]
            for i in (set(self.variant[self.sample].keys()) - set(kg + exac + gnomad + swegen + inner))
        }

        return freq

    def damaging_selection(self, criteria=[]):
        """
        Select damaging variants including loss-of-function variants (frameshift_variant, stop_gained,
        splice_acceptor_variant, splice_donor_variant, stop_lost, start_lost), and damaging missense
        variants based on vep annotation (missense_variant) and missense algorithms prediction
        (SIFT, POLYPHEN, MPC, CADD, SPIDEX, PHYLOP)
        Argument:
            criteria: selected missense algorithm list

        """
        lof = {}
        mis = {}
        mis_damage = {}
        damaging = {}
        lof[self.sample] = {}
        mis[self.sample] = {}
        mis_damage[self.sample] = {}
        for key in self.variant[self.sample]:
            if self.variant[self.sample][key]['set'] != ['freebayes'] and \
                    self.variant[self.sample][key]['set'] != ['gatk'] and \
                    self.variant[self.sample][key]['set'] != ['samtools']:
                for i in self.variant[self.sample][key]['CSQ']['Consequence']:
                    ###lof####
                    if 'frameshift_variant' in i or \
                            'stop_gained' in i or \
                            'splice_acceptor_variant' in i or \
                            'splice_donor_variant' in i or \
                            'stop_lost' in i or \
                            'start_lost' in i:
                        lof[self.sample][key] = self.variant[self.sample][key]
                    ##missense####
                    if 'missense_variant' in i:
                        mis[self.sample][key] = self.variant[self.sample][key]

        if 'SIFT' in criteria:
            sift = [
                key for key in self.variant[self.sample]
                for i in self.variant[self.sample][key]['CSQ']['SIFT']
                if 'deleterious' in i or 'deleterious_low_confidence' in i
            ]
        else:
            sift = []

        if 'POLYPHEN' in criteria:
            polyphen = [
                key for key in self.variant[self.sample]
                for i in self.variant[self.sample][key]['CSQ']['PolyPhen']
                if 'possibly_damaging' in i or 'probably_damaging' in i
            ]
        else:
            polyphen = []

        if 'MPC' in criteria:
            mpc = []
            for key in self.variant[self.sample]:
                for i in self.variant[self.sample][key]['CSQ']['MPC']:
                    if i != '':
                        if i != 'NA':
                            if float(i) >= 2.0:
                                mpc.append(key)
        else:
            mpc = []

        if 'CADD' in criteria:
            cadd = [
                key for key in self.variant[self.sample]
                if 'CADD' in self.variant[self.sample][key] and
                   float(''.join(self.variant[self.sample][key]['CADD'])) >= 20.0
            ]
        else:
            cadd = []

        if 'SPIDEX' in criteria:
            spidex = [
                key for key in self.variant[self.sample]
                if 'SPIDEX' in self.variant[self.sample][key] and
                   abs(float(''.join(self.variant[self.sample][key]['SPIDEX']))) >= 2.0
            ]
        else:
            spidex = []

        if 'PHYLOP' in criteria:
            phylop = []
            for key in self.variant[self.sample]:
                if 'dbNSFP_phyloP100way_vertebrate' in self.variant[self.sample][key]:
                    for i in self.variant[self.sample][key]['dbNSFP_phyloP100way_vertebrate']:
                        if float(i) >= 2.0:
                            phylop.append(key)
        else:
            phylop = []

        # integrate variants together and sum how many variants are from sift, polyphen, mpc....., respectively
        sum_dict = {}
        cridict = {'sift': sift, 'polyphen': polyphen, 'mpc': mpc, 'cadd': cadd, 'spidex': spidex, 'phylop': phylop}
        for i in cridict:
            for var in cridict[i]:
                if var not in sum_dict:
                    sum_dict[var] = {}
                    sum_dict[var][i] = 1
                else:
                    sum_dict[var][i] = 1

        # select missense variants fulfilling damaging prediction
        for var in sum_dict:
            if sum(sum_dict[var].values()) > 0.5 * (len(criteria)) and var in mis[self.sample]:
                mis_damage[self.sample][var] = mis[self.sample][var]

        damaging[self.sample] = {**lof[self.sample], **mis_damage[self.sample]}

        return damaging, lof, mis_damage

    def comb_selection(self, FILTER=None, DP=None, QD=None, MQ=None,
             innerfreqfile=None, KG=None, EXAC=None, GNOMAD=None, SWEGEN=None,
             criteria=[]):
        """
        Select damaging variants with rare frequency and good quality

        """
        damaging_sel = {}
        lof_sel = {}
        mis_sel = {}
        quality = {}
        freq = {}

        damaging, lof, mis_damage = self.damaging_selection(criteria=criteria)

        if FILTER or DP or QD or MQ:
            quality = self.quality_selection(FILTER=FILTER, DP=DP, QD=QD, MQ=MQ)
        else:
            quality[self.sample] = self.variant[self.sample]

        if innerfreqfile or KG or EXAC or GNOMAD or SWEGEN:
            freq = self.freq_selection(innerfreqfile=innerfreqfile, KG=KG, EXAC=EXAC, GNOMAD=GNOMAD, SWEGEN=SWEGEN)
        else:
            freq[self.sample] = self.variant[self.sample]

        damaging_sel[self.sample] = {
            key: self.variant[self.sample][key]
            for key in self.variant[self.sample]
            if key in quality[self.sample] and key in freq[self.sample] and key in damaging[self.sample]
        }
        lof_sel[self.sample] = {
            key: self.variant[self.sample][key]
            for key in self.variant[self.sample]
            if key in quality[self.sample] and key in freq[self.sample] and key in lof[self.sample]
        }
        mis_sel[self.sample] = {
            key: self.variant[self.sample][key]
            for key in self.variant[self.sample]
            if key in quality[self.sample] and key in freq[self.sample] and key in mis_damage[self.sample]
        }

        return damaging_sel, lof_sel, mis_sel
