#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os

import quantities as pq

from pyqclib.defs import UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL

# UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL = (pq.hartree * pq.constants.Avogadro_constant).rescale(KILOJOULE / pq.mol)


def get_gaussian_res(path, token):
    """
    Function to parse gaussian log file
    eg:
    s2 = get_gaussian_res(path = '~/qc/h2o.log', token = '\\S2=')
    """
    start_token = ' 1\\1\\'
    end_token = '\\\\@'
    records=[]
    recording=False
    for line in file(path,'rt'):
        if recording:
            records.append(line) # Skip first whitespace
            if end_token in line:
                recording=False; break
        elif line.startswith(start_token):
            recording=True
            records.append(line)
    if records == []:
        raise IOError("Found no token: %s in %s", token, path)
    record=''.join([x.strip(' ').strip('\n') for x in records])
    if not token in record: return None
    res = record.split(token)[1].split('\\')[0]
    return float(res)


def get_gaussian_ccsd(path):
    return get_gaussian_res(path, '\\CCSD=') * UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL


def get_gaussian_ccsd_t(path):
    return get_gaussian_res(path, '\\CCSD(T)=') * UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL


def get_gaussian_scf(path):
    return get_gaussian_res(path, 'HF=') * UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL


def get_gaussian_after_lineleading_token(path, token, ret_none_when_missing = True):
    for line in file(path,'rt'):
        if line.startswith(token):
            return float(line.split(token)[1])
    if ret_none_when_missing: return None
    raise IOError("Found no free energy in %s", path)


def get_gaussian_after_token(path, token, ret_none_when_missing = True):
    for line in file(path,'rt'):
        if token in line:
            return float(line.split(token)[1].replace('D', 'E'))
    if ret_none_when_missing: return None
    raise IOError("Found no free energy in %s", path)


def get_gaussian_zpe(path, ret_none_when_missing = True):
    raise NotImplementedError
    token = ' '
    res = get_gaussian_after_lineleading_token(path, token, ret_none_when_missing)
    if res != None:
        res *= UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL
    return res


def get_gaussian_thermal(path, ret_none_when_missing = True):
    token = ' Sum of electronic and thermal Free Energies='
    res = get_gaussian_after_lineleading_token(path, token, ret_none_when_missing)
    if res != None:
        res *= UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL
    return res


def get_gaussian_eump2(path, ret_none_when_missing = True):
    token = 'EUMP2 = '
    res = get_gaussian_after_token(path, token, ret_none_when_missing)
    if res != None:
        res *= UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL
    return res


def get_line_in_file_containing(path, token):
    for line in file(path, 'rt'):
        if token in line:
            return line
    raise IOError('Token %s not found in %s' % (token, path))


def get_gaussian_scs_eump2(path, ret_none_when_missing = True,
                           ps = 6.0 / 5.0, pt = 1.0 / 3.0):
    """
    Returns spin component scaled UMP2 energy according to Grimme (2003)
    """
    eump2 = get_gaussian_eump2(path, ret_none_when_missing)

    aa_line   = get_line_in_file_containing(path, 'alpha-alpha')
    aa_E2     = float(aa_line.split('=')[2].replace('D', 'E'))
    ab_line   = get_line_in_file_containing(path, 'alpha-beta')
    ab_E2     = float(ab_line.split('=')[2].replace('D', 'E'))
    bb_line   = get_line_in_file_containing(path, 'beta-beta')
    bb_E2     = float(bb_line.split('=')[2].replace('D', 'E'))
    scs_cor = (ps - 1) * ab_E2 + (pt - 1) * (aa_E2 + bb_E2)
    scs_cor *= UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL

    scs_eump2 = eump2 + scs_cor
    return scs_eump2


