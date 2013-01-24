#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from pyqclib.parse_g09 import get_mulliken_chg
from pyqclib.defs import ELEMENTARY_CHARGE

def test_get_mulliken_chg():
    path = os.path.join(os.path.dirname(__file__), 'h2o.log')
    chgs = get_mulliken_chg(path)
    assert chgs == [x * ELEMENTARY_CHARGE for x in [-0.733237, 0.366619, 0.366619]]

