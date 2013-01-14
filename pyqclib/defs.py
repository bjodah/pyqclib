#!/usr/bin/env python
# -*- coding: utf-8 -*-

import quantities as pq

KILOJOULE = pq.UnitQuantity('kilojoule',  pq.joule * 1e3,  symbol = 'kJ')
KILOCALORIE = pq.UnitQuantity('kilocalorie',  pq.calorie * 1e3,  symbol = 'kCal')

UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL = \
    (pq.hartree * pq.constants.Avogadro_constant).rescale(KILOJOULE / pq.mol)

UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL = \
    (pq.eV * pq.constants.Avogadro_constant).rescale(KILOJOULE / pq.mol)

UNITLESS_IN_HARTREE_TO_TYPED_KILOJOULE_PER_MOL = \
    (pq.hartree * pq.constants.Avogadro_constant).rescale(KILOCALORIE / pq.mol)

UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL = \
    (pq.eV * pq.constants.Avogadro_constant).rescale(KILOCALORIE / pq.mol)

