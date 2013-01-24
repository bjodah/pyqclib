#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os
from itertools import product

import numpy as np

from pyqclib.parsing import QCResult
from pyqclib.defs import UNITS
from pyqclib.tabular_extractor import QCTabularExtractor

bond_lengths = [('C', 'Cl'), ('C', 'Br')]
log_files = ['sub_reactant_a.log', 'reactant.log',
             'ts.log', 'product.log', 'sub_product_a.log']


SRC_DIR = './example_data_01/'

class BondLengthTabularExtractor(QCTabularExtractor):

    prop = 'bond_length'

    def _path_resolver(self, row_var, col_var):
        return os.path.join(SRC_DIR, col_var)

    def mk_prop_getter_args(self, row_var, col_var):
        return [self.qcres(row_var, col_var).get_atom_indices1_of_element(x) for x in row_var]


if __name__ == '__main__':
    blte = BondLengthTabularExtractor(bond_lengths, log_files,
                                      'bond_lengths.dat', UNITS['length']['angstrom'])
    blte.extract()

