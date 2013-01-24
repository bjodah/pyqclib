#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os
from itertools import product

import numpy as np

from pyqclib.parsing import QCResult
from pyqclib.defs import UNITS
from pyqclib.tabular_extractor import QCTabularExtractor

charge_on_atoms = ['C', 'Cl', 'Br']
log_files = ['sub_reactant_a.log', 'sub_reactant_b.log', 'reactant.log',
             'ts.log', 'product.log', 'sub_product_a.log', 'sub_product_b.log']


SRC_DIR = './example_data_01/'

class BondLengthTabularExtractor(QCTabularExtractor):

    prop = 'g09_mulliken_chg'

    def _path_resolver(self, row_var, col_var):
        return os.path.join(SRC_DIR, col_var)

    def mk_prop_getter_args(self, row_var, col_var):
        indices = self.qcres(row_var, col_var).get_atom_indices1_of_element(row_var)
        assert len(indices) <= 1
        return indices


if __name__ == '__main__':
    blte = BondLengthTabularExtractor(charge_on_atoms, log_files,
                                      'charges.dat', UNITS['charge']['elementary_charge'])
    blte.extract()

