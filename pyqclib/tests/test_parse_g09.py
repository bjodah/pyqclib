#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from pyqclib.parse_g09 import get_mulliken_chg, get_mulliken_spin
from pyqclib.defs import ELEMENTARY_CHARGE, UNITLESS

def test_get_mulliken_chg():
    path = os.path.join(os.path.dirname(__file__), 'h2o.log')
    chgs = get_mulliken_chg(path)
    assert chgs == [x * ELEMENTARY_CHARGE for x in [-0.733237, 0.366619, 0.366619]]

def test_get_mulliken_spin():
    # Tests both g09 and g03
    for fname in ['ooh.log', 'ooh_freq_g03.log']:
        path = os.path.join(os.path.dirname(__file__), fname)
        spin = get_mulliken_spin(path)
        assert spin == [x * UNITLESS for x in [0.239162, 0.772545, -0.011707]]



