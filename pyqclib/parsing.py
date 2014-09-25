import logging, os

import numpy as np
import quantities as pq
import periodictable as pt

import parse_g09
from cclib.parser import ccopen

from pyqclib.defs import UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL, UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOCALORIE_PER_MOL, UNITLESS_IN_ELECTRONVOLT_TO_TYPED_EV, KCAL_PER_MOLE, ANGSTROM, ELEMENTARY_CHARGE, UNITLESS, DEGREE

# Conversion factors
CONVERSION_FROM_UNITLESS_EV_TO = {'kCal_per_mole': UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOCALORIE_PER_MOL,
                                  'kJ_per_mole': UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL,
                                  'eV': UNITLESS_IN_ELECTRONVOLT_TO_TYPED_EV}


def _get_parsed_cclib_data(logfile, logging_level = logging.ERROR):
    ccfile = ccopen(logfile)
    ccfile.logger.setLevel(logging.ERROR)
    cclib_data = ccfile.parse()
    return cclib_data

def get_scfenergy(logfile, energy_unit, index = -1):
    cclib_data = _get_parsed_cclib_data(logfile)
    return cclib_data.scfenergies[-1] * CONVERSION_FROM_UNITLESS_EV_TO[energy_unit]

get_scfenergy.units = CONVERSION_FROM_UNITLESS_EV_TO.keys()

class Property(object):
    """
    Property object to represent type of data that
    can be extracted from a QCResult
    """
    def __init__(self, key, units = None):
        if units == None:
            self._units = [None]
        else:
            self._units = units
        self._key = key

    def __str__(self):
        return self._key


class QCResult(object):
    """
    Class to represent a result from a QC computation
    (one - to one correspondence with logfile)
    """

    default_units = {'energy': KCAL_PER_MOLE,
                     'entropy': KCAL_PER_MOLE / pq.kelvin,
                     'length': ANGSTROM,
                     'charge': ELEMENTARY_CHARGE,
                     'angle': DEGREE,
                     UNITLESS: UNITLESS}

    def __init__(self, logfile, name = None, calculationtype = None,
                 T_ref = 298.15 * pq.Kelvin, linear = False):
        if not os.path.isfile(logfile):
            raise ValueError('Non-existant logfile: {}'.format(logfile))
        self._name = name
        self._logfile = logfile
        self._calculationtype = calculationtype
        self._T_ref = T_ref
        self._cclib_data = None
        self._linear = linear


    def scfenergy(self, index = -1):
        if isinstance(index, basestring): index = int(index) # Neglects string slices
        return (self.cclib_data.scfenergies[index] * pq.eV * \
                pq.constants.Avogadro_constant).rescale(
            self.default_units['energy'])

    scfenergy.returns_with_unit = True
    scfenergy.unit_type = 'energy' # key in defs.UNITS


    def vibtemps(self):
        try:
            rl = self.data.vibfreqs/pq.cm # reciprocal lambda
            rl = rl[:(self.data.natom * 3 - {True: 5, False: 6}[self._linear])]
        except AttributeError:
            return None
        return (pq.speed_of_light*rl*pq.constants.Planck_constant/pq.constants.Boltzmann_constant).rescale(pq.Kelvin)


    def real_vibtemps(self):
        vt = self.vibtemps
        if vt == None: return None
        if vt[0] < 0:
            vt = vt[1:]
        if vt[0] < 0: raise ValueError('Hessian order larger than 1')
        return vt


    def E_vib(self, T):
        vt = self.real_vibtemps()
        if vt == None: return None
        NA = pq.constants.Avogadro_constant
        kB = pq.constants.Boltzmann_constant

        if T < 1e-3 * pq.K:
            Evib = NA * kB * np.sum(vt / 2)
        else:
            Evib = NA * kB * np.sum(vt / 2 + vt * np.exp( - vt / T) /(1 - np.exp( - vt / T)))

        return Evib.rescale(self.default_units['energy'])
    E_vib.returns_with_unit = True
    E_vib.unit_type = 'energy'

    # Convenience function, Zero - point vibrational energy
    def ZPVE(self):
        return self.E_vib(0 * pq.K)
    ZPVE.returns_with_unit = E_vib.returns_with_unit
    ZPVE.unit_type = E_vib.unit_type

    def q_vib(self, T):
        vt = self.real_vibtemps()
        if vt == None: return None
        return np.exp(-vt/(2 * T))/(1-np.exp(-vt/T))

    def S_vib(self, T):
        vt = self.real_vibtemps()
        if vt == None: return None
        NA = pq.constants.Avogadro_constant
        kB = pq.constants.Boltzmann_constant

        if T < 1e-3 * pq.K:
            Svib = NA * kB * 0
        else:
            Svib = NA * kB * np.sum(vt / T * np.exp( - vt / T) / (1 - np.exp( - vt / T)) - \
                                np.log(1 - np.exp( - vt / T)))
        return Svib.rescale(self.default_units['entropy'])
    S_vib.returns_with_unit = True
    S_vib.unit_type = 'entropy'

    def S_elec(self, T):
        pass
        # Selec =
        # return Selec.rescale(default_energy_unit / pq.kelvin)
    S_vib.returns_with_unit = True
    S_vib.unit_type = 'entropy'


    def gibbs_energy_elec(self, T = None):
        pass
        #return

    def bond_length(self, first_index1, second_index1, geom_index0 = -1):
        v = np.linalg.norm(self.cclib_data.atomcoords[geom_index0][first_index1 - 1] - \
                              self.cclib_data.atomcoords[geom_index0][second_index1 - 1])
        if v == 0.0:
            return np.nan * ANGSTROM
        else:
            return v * ANGSTROM
    bond_length.returns_with_unit = True
    bond_length.unit_type = 'length' # key in defs.UNITS

    def bond_angle(self, first_index1, second_index1, third_index1, geom_index0 = -1):
        a, b, c = (self.cclib_data.atomcoords[geom_index0][first_index1 - 1],
                   self.cclib_data.atomcoords[geom_index0][second_index1 - 1],
                   self.cclib_data.atomcoords[geom_index0][third_index1 - 1])
        v1 = a - b
        v2 = c - b

        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)

        angle = np.arccos(np.dot(v1, v2))
        if np.isnan(angle):
            if np.all(v1 == v2):
                angle = 0
            else:
                angle = np.pi
        return angle * 180 / np.pi * DEGREE
    bond_angle.returns_with_unit = True
    bond_angle.unit_type = 'angle' # key in defs.UNITS



    def get_atom_indices1_of_element(self, elem):
        if isinstance(elem, int):
            return np.where(self.cclib_data.atomnos == elem)[0] + 1
        else:
            return np.where(self.cclib_data.atomnos == getattr(pt, elem).number)[0] + 1

    def g09_thermal(self):
        return parse_g09.get_gaussian_thermal(self._logfile).rescale(self.default_units['energy'])
    g09_thermal.returns_with_unit = True
    g09_thermal.unit_type = 'energy' # key in defs.UNITS

    def g09_w_zpve(self):
        return parse_g09.get_g09_w_zpve(self._logfile).rescale(self.default_units['energy'])
    g09_w_zpve.returns_with_unit = True
    g09_w_zpve.unit_type = 'energy' # key in defs.UNITS

    def g09_thermal_enthalpy(self):
        return parse_g09.get_g09_thermal_enthalpy(self._logfile).rescale(self.default_units['energy'])
    g09_thermal_enthalpy.returns_with_unit = True
    g09_thermal_enthalpy.unit_type = 'energy' # key in defs.UNITS


    def g09_mulliken_chg(self, atom_index1 = -1):
        if atom_index1 == -1: return np.nan * self.default_units['charge']
        return parse_g09.get_mulliken_chg(self._logfile)[atom_index1 - 1].rescale(
            self.default_units['charge'])
    g09_mulliken_chg.returns_with_unit = True
    g09_mulliken_chg.unit_type = 'charge' # key in defs.UNITS

    def g09_mulliken_spin(self, atom_index1 = -1):
        if atom_index1 == -1: return np.nan * UNITLESS
        return parse_g09.get_mulliken_spin(self._logfile)[atom_index1 - 1]
    g09_mulliken_spin.returns_with_unit = False
    g09_mulliken_spin.unit_type = UNITLESS # key in defs.UNITS


    qc_properties = [scfenergy, g09_thermal, g09_mulliken_chg, g09_mulliken_spin]

    @property
    def cclib_data(self):
        if self._cclib_data == None:
            self._cclib_data = _get_parsed_cclib_data(self._logfile)
        return self._cclib_data


# class NonExistentQCResult(QCResult):

#     def __init__(self, collect_dir = None):
#         self._collect_dir = collect_dir
#         return None

#     def __getitem__(self, key):
#         return NonExistentStructure(self._collect_dir)

#     def __getattr__(self, attr):
#         return None

#     @property
#     def _logfile(self):
#         return self._collect_dir
