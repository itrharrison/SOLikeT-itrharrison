r"""Class to calculate bias models for haloes and galaxies as cobaya Theory classes

"""

import numpy as np
from typing import Sequence, Union
from cobaya.theory import Theory


class Bias(Theory):

    kmax: float = 10.  # Maximum k (1/Mpc units) for Pk, or zero if not needed
    nonlinear: bool = False  # whether to get non-linear Pk from CAMB/Class
    z: Union[Sequence, np.ndarray] = []  # redshift sampling

    extra_args: dict = {}  # extra (non-parameter) arguments passed to ccl.Cosmology()

    _logz = np.linspace(-3, np.log10(1100), 150)
    _default_z_sampling = 10**_logz
    _default_z_sampling[0] = 0

    def initialize(self):

        self._var_pairs = set()

    def get_requirements(self):
        return {}

    def must_provide(self, **requirements):

        options = requirements.get('linear_bias') or {}

        self.kmax = max(self.kmax, options.get('kmax', self.kmax))
        self.z = np.unique(np.concatenate(
                            (np.atleast_1d(options.get("z", self._default_z_sampling)),
                            np.atleast_1d(self.z))))

        # Dictionary of the things needed from CAMB/CLASS
        needs = {}

        self.nonlinear = self.nonlinear or options.get('nonlinear', False)
        self._var_pairs.update(
            set((x, y) for x, y in
                options.get('vars_pairs', [('delta_tot', 'delta_tot')])))

        needs['Pk_grid'] = {
                'vars_pairs': self._var_pairs or [('delta_tot', 'delta_tot')],
                'nonlinear': (True, False) if self.nonlinear else False,
                'z': self.z,
                'k_max': self.kmax
            }

        return needs

    def _get_Pk_mm(self):

        for pair in self._var_pairs:

            if self.nonlinear:
                k, z, Pk_nonlin = self.provider.get_Pk_grid(var_pair=pair,
                                                            nonlinear=True)
                Pk_mm = np.flip(Pk_nonlin, axis=0)
            else:
                k, z, Pk_lin = self.provider.get_Pk_grid(var_pair=pair,
                                                         nonlinear=False)
                Pk_mm = np.flip(Pk_lin, axis=0)

        return Pk_mm

    # def calculate(self, state, want_derived=True, **params_values_dict):

    def get_Pkgg(self):
        return self._current_state['Pk_gg']

    def get_Pkgm(self):
        return self._current_state['Pk_gm']


class Linear_bias(Bias):

    params = {'b_lin': None}

    def calculate(self, state, want_derived=True, **params_values_dict):

        Pk_mm = self._get_Pk_mm()

        state['Pk_gg'] = params_values_dict['b_lin']**2. * Pk_mm
        state['Pk_gm'] = params_values_dict['b_lin'] * Pk_mm