import logging, os

import quantities as pq
import parse_g09
from cclib.parser import ccopen

from pyqclib.defs import UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL, UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL, UNITLESS_IN_ELECTRONVOLT_TO_TYPED_EV, KCAL_PER_MOLE

# Conversion factors
CONVERSION_FROM_UNITLESS_EV_TO = {'kCal_per_mole': UNITLESS_IN_ELECTRONVOLT_TO_TYPED_KILOJOULE_PER_MOL,
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
                     'entropy': KCAL_PER_MOLE / pq.kelvin}

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

    def g09_thermal(self):
        return parse_g09.get_gaussian_thermal(self._logfile).rescale(self.default_units['energy'])
    g09_thermal.returns_with_unit = True
    g09_thermal.unit_type = 'energy' # key in defs.UNITS

    qc_properties = [scfenergy, g09_thermal]

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

