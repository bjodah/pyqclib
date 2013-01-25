import matplotlib.pyplot as plt
from pyqclib.rxn_barrier_fig import BarriersSerie, BarrierPlotter

def main(barriers, labels = [], outfile = None):
    """
    """
    barrier_serie1 = BarriersSerie(
        barriers, annotations = dict(enumerate(labels)))
    barrier_plotter = BarrierPlotter(ylabel = '$\mathrm{G~/~kcal\cdot mol^{-1}}$',
                                     outfile = outfile)
    barrier_plotter.plot(barrier_serie1)

    #plt.savefig(self._outfile)
    plt.show()


if __name__ == '__main__':
    test_data = [0.0, 4.7, -6.3]
    labels = ['reac', 'TS', 'prod']
    main(barriers = test_data, labels = labels, outfile = 'barriers_example.png')
