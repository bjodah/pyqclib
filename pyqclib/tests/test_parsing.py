import os

import numpy as np
from pyqclib.defs import ANGSTROM
from pyqclib.parsing import QCResult


def test_QCResult_get_bond_length():
    qcr = QCResult(os.path.join(os.path.dirname(__file__), 'h2o.log'))
    OH = qcr.bond_length(1, 2)
    HH = qcr.bond_length(2, 3)
    assert float(OH.rescale(ANGSTROM)) - 0.980024 < 1e-6
    assert float(HH.rescale(ANGSTROM)) - 1.530135 < 1e-6

def test_QCResult_get_atom_indices1_of_element():
    qcr = QCResult(os.path.join(os.path.dirname(__file__), 'h2o.log'))
    H_index = qcr.get_atom_indices1_of_element('H')
    O_index = qcr.get_atom_indices1_of_element(8)
    assert np.allclose(H_index, [2, 3])
    assert np.allclose(O_index, [1])
