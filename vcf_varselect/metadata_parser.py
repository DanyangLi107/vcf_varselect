#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'Danyang Li'
# __date__ = 2020-04-02

import re
class MetadataParser(object):
    """
    Parse vcf metadata

    """

    def __init__(self):
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AGR]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AGR]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.id_dict = {'INFO': [], 'FORMAT': [], 'FILTER': []}
        self.vep_columns = []
        self.header = []
        self.info_keys = ['ID', 'Number', 'Type', 'Description']

    def __iter__(self):
        return iter(self.__dict__.items())

    def read_metadata(self, line):
        """
        Parse metadata to dict

        """
        line = line.rstrip()
        line_info = line[2:].split('=')
        match = False

        if line_info[0] == 'fileformat':
            try:
                self.fileformat = line_info[1]
            except IndexError:
                raise SyntaxError("fileformat must have a value")
        elif line_info[0] == 'INFO':
            match = self.info_pattern.match(line)
            if not match:
                raise SyntaxError("One of the INFO lines is malformed:{0}".format(line))
            matches = [match.group('id'), match.group('number'),
                       match.group('type'), match.group('desc')]

            info_line = dict(list(zip(self.info_keys, matches)))

            if info_line['ID'] == 'CSQ':
                info_line['Format'] = [
                    info.strip() for info in info_line['Description'].split('Format:')
                ][-1]
                self.vep_columns = info_line.get('Format', '').split('|')

            self.id_dict['INFO'].append(match.group('id'))

        elif line_info[0] == 'FILTER':
            match = self.filter_pattern.match(line)
            if not match:
                raise SyntaxError("One of the FILTER lines is malformed: {0}".format(line))
            self.id_dict['FILTER'].append(match.group('id'))

        elif line_info[0] == 'FORMAT':
            match = self.format_pattern.match(line)
            if not match:
                raise SyntaxError("One of the FORMAT lines is malformed: {0}".format(line))
            self.id_dict['FORMAT'].append(match.group('id'))

    def read_header(self, line):
        """
        Parse header

        """
        self.header = line[1:].rstrip().split('\t')
