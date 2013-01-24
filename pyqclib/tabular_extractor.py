#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os
from itertools import product

import numpy as np

from pyqclib.parsing import QCResult
from pyqclib.defs import UNITLESS

class TabularExtractor(object):
    """
    Interface for creating extractor routines
    where two variables are varied (row_vars and col_vars)
    and the data is conveniently put in 2D array structure
    """

    def __init__(self, row_vars, col_vars, output, unit):
        """

        Arguments:
        - `row_vars`:
        - `col_vars`:
        """
        self._row_vars = row_vars
        self._col_vars = col_vars
        self._output = output
        self._unit = unit
        self._data = np.zeros((len(row_vars), len(col_vars)))


    def _populate(self):
        for (ri, r), (ci, c) in product(enumerate(self._row_vars), enumerate(self._col_vars)):
            prop_getter = self.mk_prop_getter(r, c)
            prop_getter_args = self.mk_prop_getter_args(r, c)
            val = prop_getter(*prop_getter_args)
            self._data[ri, ci] = float(val.rescale(self._unit))


    def extract(self):
        self._populate()
        np.savetxt(self._output, self._data)

    def mk_prop_getter(self, row_var, col_var):
        pass

    def mk_prop_getter_args(self, row_var, col_var):
        pass


class QCTabularExtractor(TabularExtractor):

    _qcres = {}
    prop = 'SET_THIS_IN_SUBCLASS'

    def qcres(self, row_var, col_var):
        """ Caches QCResult instances """
        path = self._path_resolver(row_var, col_var)
        if not path in self._qcres:
            self._qcres[path] = QCResult(path)
        return self._qcres[path]

    def _path_resolver(self, row_var, col_var):
        """ To be subclassed """
        pass

    def mk_prop_getter(self, row_var, col_var):
        return getattr(self.qcres(row_var, col_var), self.prop)
