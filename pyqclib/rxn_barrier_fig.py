from __future__ import division

import matplotlib.pyplot as plt
import matplotlib.patches

import matplotlib.axes
from matplotlib import rcParams

import numpy as np

import scipy.interpolate

import quantities as pq


class BarriersSerie(object):
    """
    Serie to hold data of barrier for plotting
    """

    default_color = 'black'
    default_barrier_lw = 2


    def __init__(self, barriers, zeropoint = 0.0, unit = None, color = None,
                 barrier_lw = None, annotations = None):
        """

        Arguments:
        - `barriers`:
        - `zeropoint`:
        - `unit`:
        - `color`: Default is 'black'
        - `barrier_lw`: Default is 2
        - `annotations`: Dict mapping index to TeX formated string
        """
        self._barriers = barriers
        self._zeropoint = zeropoint
        self._unit = unit
        self._color = color if color != None else self.default_color
        self._barrier_lw = barrier_lw if barrier_lw != None else self.default_barrier_lw
        self._annotations = annotations


    def get_element_args_kwargs(self):
        barrier_line_elems   = self._get_barrier_lines()
        interconn_line_elems = self._get_interconnecting_lines()
        return barrier_line_elems + interconn_line_elems


    @property
    def n(self):
        return len(self._barriers)

    @property
    def x(self):
        return np.linspace(0, 1, 2 * self.n + 2)

    def x_center_barrier_bar(self, idx):
        return (self.x[idx * 2 + 1] + self.x[idx * 2 + 2]) / 2

    def _get_barrier_lines(self):
        element_args_kwargs = []
        for i, yval in enumerate(self._barriers):
            xval = self.x[i * 2 + 1], self.x[i * 2 + 2]
            # Plot barrier bars
            element_args_kwargs.append((
                (xval, [yval] * 2),
                 {'color': self._color, 'linewidth': self._barrier_lw}
                ))
        return element_args_kwargs

    def _get_interconnecting_lines(self):
        element_args_kwargs = []
        for i, yval in enumerate(self._barriers):
            if i < len(self._barriers) - 1:
                interp_xlim = self.x[i * 2 + 2], self.x[i * 2 + 3]
                interp_ylim = self._barriers[i], self._barriers[i + 1]
                interp_deltax = interp_xlim[1] - interp_xlim[0]
                interp_deltay = interp_ylim[1] - interp_ylim[0]
                interp_x = np.linspace(*interp_xlim)
                interp_y = scipy.interpolate.PiecewisePolynomial(
                    [interp_xlim[0], (interp_xlim[0] + interp_xlim[1]) / 2, interp_xlim[1]],
                    [
                        [interp_ylim[0], 0],
                        [(interp_ylim[1] + interp_ylim[0]) / 2, 3 * interp_deltay / interp_deltax],
                        [interp_ylim[1], 0]
                    ],
                    3)(interp_x)
                element_args_kwargs.append((
                    (interp_x, interp_y),
                    {'color': self._color, 'linewidth': self._barrier_lw / 2,
                     'linestyle': 'dashed'}
                    ))
        return element_args_kwargs


class BarrierPlotter(object):
    """
    Convenience class for making plots of reaction
    barriers of results from Quantum Chemical computations
    """

    ax_l, ax_b, ax_w, ax_h = 0.15, 0.15, 0.7, 0.7
    _occ_boxes = []

    def __init__(self, xticks = None, xlabel = None, ylabel = None,
                 outfile = None):
        """

        Arguments:
        - `xticks`:
        - `xlabel`:
        - `ylabel`:
        """
        self._xticks = xticks if xticks != None else []
        self._xlabel = xlabel
        self._ylabel = ylabel
        self._outfile = outfile

    def __call__(self, *args):
        fig, ax = self._pre_plotting()
        for barrier_serie in args:
            for element_args_kwargs in barrier_serie.get_element_args_kwargs():
                elem_args, elem_kwargs = element_args_kwargs
                ax.plot(*elem_args, **elem_kwargs)
        for barrier_serie in args:
            # Now all series are in place, lets annotate
            self._annotate(barrier_serie, ax)
        self._post_plotting()

    def _pre_plotting(self):
        rcParams.update({'text.usetex': True,
                         'savefig.dpi': 600,
                         })
        width_inches = 4.8 # Preferably get fig_with_pt from LaTeX using \showthe\columnwidth
        golden_ratio = 0.6180339
        fig = plt.figure(figsize = [width_inches, width_inches * golden_ratio])
        ax = fig.add_axes([self.ax_l, self.ax_b, self.ax_w, self.ax_h])
        return fig, ax

    def _assign_to_free_space(self, ax, xy, fontsize):
        x, y = xy
        bbox = ax.get_window_extent()
        bounds = bbox.bounds
        yspan = bounds[3] - bounds[1]
        fontsize_in_frac = fontsize / yspan
        return x, y + fontsize_in_frac / 2

    def _annotate(self, barriers_serie, ax):

        for bi, anno_txt in barriers_serie._annotations.iteritems():
            ymin, ymax = ax.get_ylim()
            yspan = ymax - ymin
            norm_y = (barriers_serie._barriers[bi] - ymin) / yspan
            xy = (barriers_serie.x_center_barrier_bar(bi), norm_y)
            xytext = self._assign_to_free_space(ax, xy, rcParams['font.size'])
            ax.annotate(anno_txt, xy, xycoords = 'axes fraction', xytext = xytext,
                         horizontalalignment = 'center', verticalalignment = 'center')


    def _post_plotting(self):
        plt.xticks(self._xticks)
        plt.xlim([0, 1])
        if self._xlabel != None: plt.xlabel(self._xlabel)
        if self._ylabel != None: plt.ylabel(self._ylabel)
        if self._outfile != None:
            plt.savefig(self._outfile)
        else:
            plt.show()

