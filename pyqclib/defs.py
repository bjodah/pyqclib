#!/usr/bin/env python
# -*- coding: utf-8 -*-

import quantities as pq

# Alias to be able to easily switch dependency upon quantities if neccessary
UNITLESS = 1.0

#Length
ANGSTROM = pq.angstrom

#Charge
COULOMB = pq.coulomb
ELEMENTARY_CHARGE = pq.elementary_charge

# ANGLE
DEGREE = pq.degree
RADIAN = pq.radian

# ENERGY UNITS
KILOJOULE = pq.UnitQuantity('kilojoule',  pq.joule * 1e3,  symbol = 'kJ')
KILOCALORIE = pq.UnitQuantity('kilocalorie',  pq.calorie * 1e3,  symbol = 'kCal')


UNITLESS_IN_ELECTRONVOLT_TO_TYPED_EV = pq.eV
UNITLESS_IN_HARTREE_TO_TYPED_HARTREE = pq.hartree
UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL = \
    (pq.hartree * pq.constants.Avogadro_constant).rescale(KILOJOULE / pq.mol)
UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL = \
    (pq.eV * pq.constants.Avogadro_constant).rescale(KILOJOULE / pq.mol)
UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL = \
    (pq.hartree * pq.constants.Avogadro_constant).rescale(KILOCALORIE / pq.mol)
UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL = \
    (pq.eV * pq.constants.Avogadro_constant).rescale(KILOCALORIE / pq.mol)

KJ_PER_MOLE   = KILOJOULE / pq.mol
KCAL_PER_MOLE = KILOCALORIE / pq.mol
EV = pq.eV
HARTREE = pq.hartree

UNITS = {'energy': {'kJ_per_mole': KJ_PER_MOLE,
                    'kCal_per_mole': KCAL_PER_MOLE,
                    'eV': EV,
                    'hartree': HARTREE},
         'unitless': {'unitless': UNITLESS},
         'length': {'angstrom': ANGSTROM},
         'charge': {'coulomb': COULOMB,
                    'elementary_charge': ELEMENTARY_CHARGE},
         'angle': {'degree': DEGREE,
                    'radian': RADIAN},
         1.0: {1.0: UNITLESS}
         }

class LevelOfTheory(object):
    """
    Electronic Structure Theory Level of Theory
    """
    def __init__(self, elec_corr_method, basis):
        """

        Arguments:
        - `elec_corr_method`:
        - `basis`:
        """
        self._elec_corr_method = elec_corr_method
        self._basis = basis

class ElecCorrMethod(object):
    # HF / DFT / QMC
    pass

class DFT(ElecCorrMethod):
    # Jacobs ladder
    pass

class Environment(object):
    """
    Electronic Structure Theory Level of Theory
    """
    def __init__(self,
                 vacuum = True,
                 method = None,
                 eps = None,
                 eps_inf = None):
        """

        Arguments:
        - `vacuum`:
        - `method`:
        - `eps`:
        - `eps_inf`:
        """
        self._vacuum = vacuum
        self._method = method
        self._eps = eps
        self._eps_inf = eps_inf


#PROPERTIES = ['scfenergy', 'vibfreq', 'G(298K)']

class Property(object):
    pass


class CalculationType(object):
    """
    Electronic Structure Theory Computation
    """

    def __init__(self, lvl_of_theory, environment, considered_props):
        """

        Arguments:
        - `molecule`:
        - `lvl_of_theory`:
        - `environment`:
        - `considered_props`:
        """
        self._lvl_of_theory = lvl_of_theory
        self._environment = environment
        self._considered_props = considered_props
