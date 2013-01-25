from __future__ import division

import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib import rc
import matplotlib.axes
from matplotlib import rcParams

import numpy as np

import scipy.interpolate

import quantities as pq

from dahlgren_pyutils.simple_plot import entex

class BarriersSerie(object):
    """
    Serie to hold data of barrier for plotting
    """

    default_color = 'black'
    default_barrier_lw = 2


    def __init__(self, barriers, zeropoint = 0.0, unit = None, color = None,
                 barrier_lw = None, annotations = None, legend_lbl = None):
        """

        Arguments:
        - `barriers`:
        - `zeropoint`:
        - `unit`:
        - `color`: Default is 'black'
        - `barrier_lw`: Default is 2
        - `annotations`: Dict mapping index to TeX formated string
        - `legend_lbl`: Label in legend
        """
        self._barriers = barriers
        self._zeropoint = zeropoint
        self._unit = unit
        self._color = color if color != None else self.default_color
        self._barrier_lw = barrier_lw if barrier_lw != None else self.default_barrier_lw
        self._annotations = annotations if annotations != None else {}
        self._legend_lbl = legend_lbl


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
            barrier_kwargs = {'color': self._color, 'linewidth': self._barrier_lw}
            if i == 0:
                barrier_kwargs.update({'label': entex(self._legend_lbl)})
            element_args_kwargs.append((
                (xval, [yval] * 2),
                barrier_kwargs
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

    ax_l, ax_b, ax_w, ax_h = 0.15, 0.10, 0.7, 0.7
    _occ_boxes = []

    _fig = None
    _ax = None

    def __init__(self, barrier_series, xticks = None, xlabel = None, ylabel = None,
                 outfile = None, global_annotations = None):
        """

        Arguments:
        - `barrier_series`: list of BarrierSerie instances
        - `xticks`:
        - `xlabel`:
        - `ylabel`:
        - `outfile`:
        - `global_annotations`:
        """
        self._barrier_series = barrier_series
        self._xticks = xticks if xticks != None else []
        self._xlabel = xlabel
        self._ylabel = ylabel
        self._outfile = outfile
        self._global_annotations = global_annotations if global_annotations != None else {}
        self._init_plotting()


    def add_barriers(self, *args):
        self._barrier_series.extend(args)

    def plot(self):
        for barrier_serie in self._barrier_series:
            for element_args_kwargs in barrier_serie.get_element_args_kwargs():
                elem_args, elem_kwargs = element_args_kwargs
                self._ax.plot(*elem_args, **elem_kwargs)
        self._annotate_barrier_series()
        self._global_annotate()
        self._post_plotting()

    def add_legend(self):
        self._ax.legend(loc = 'upper center',  bbox_to_anchor = (0.5, 1.25),
          fancybox = True,  shadow = True,  ncol = len(self._barrier_series))


    def _init_plotting(self):
        rcParams.update({'text.usetex': True,
                         'savefig.dpi': 600,
                         'legend.fontsize': 10,
                         })
        matplotlib.rc('text.latex',  preamble = '\usepackage{color}')
        width_inches = 6 # Preferably get fig_with_pt from LaTeX using \showthe\columnwidth
        golden_ratio = 0.6180339
        self._fig = plt.figure(figsize = [width_inches, width_inches * golden_ratio])
        self._ax = self._fig.add_axes([self.ax_l, self.ax_b, self.ax_w, self.ax_h])


    def _assign_to_free_space(self, xy, fontsize):
        x, y = xy
        bbox = self._ax.get_window_extent()
        bounds = bbox.bounds
        yspan = bounds[3] - bounds[1]
        fontsize_in_frac = fontsize / yspan
        return x, y - fontsize_in_frac / 1.7 # <--- This needs to be IMPROVED!!

    def _get_norm_y_of_barriers_serie(self, barriers_serie, barrier_idx):
        ymin, ymax = self._ax.get_ylim()
        yspan = ymax - ymin
        norm_y = (barriers_serie._barriers[barrier_idx] - ymin) / yspan
        return norm_y

    def _annotate_barrier_series(self):
        for bs in self._barrier_series:
            for bi, anno_txt in bs._annotations.iteritems():
                norm_y = self._get_norm_y_of_barriers_serie(bs, bi)
                xy = (barriers_serie.x_center_barrier_bar(bi), norm_y)
                xytext = self._assign_to_free_space(xy, rcParams['font.size'])
                self._ax.annotate(anno_txt, xy, xycoords = 'axes fraction', xytext = xytext,
                                  horizontalalignment = 'center', verticalalignment = 'center')


    def _global_annotate(self):
        for i, anno_txt in self._global_annotations.iteritems():
            xy = (self._barrier_series[0].x_center_barrier_bar(i),
                  min([self._get_norm_y_of_barriers_serie(barrier_serie, i) for\
                       barrier_serie in self._barrier_series]))
            xytext = self._assign_to_free_space(xy, rcParams['font.size'])
            self._ax.annotate(anno_txt, xy, xycoords = 'axes fraction', xytext = xytext,
                         horizontalalignment = 'center', verticalalignment = 'center')


    def _post_plotting(self):
        plt.xticks(self._xticks)
        plt.xlim([0, 1])
        if self._xlabel != None: plt.xlabel(self._xlabel)
        if self._ylabel != None: plt.ylabel(self._ylabel)

