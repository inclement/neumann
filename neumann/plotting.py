import numpy as n
import matplotlib.pyplot as plt


def plot_histogram(results):
    energies = sorted(results.keys())
    fig, ax = plt.subplots()

    legend = []

    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
    colours = ['blue', 'green', 'red', 'purple', 'cyan', 'orange']

    i = 0
    for energy in energies:
        ax.hist(results[energy][2],
                linestyle=linestyles[(i % len(linestyles))],
                color = colours[i],
                linewidth=3,
                bins=23,
                normed=True,
                histtype='step',
                )
        legend.append(str(energy))
        i += 1

    l = ax.legend(legend)
    l.set_title('energy')

    ax.set_xlabel('$\\rho$')
    ax.set_ylabel('PDF')

    return fig, ax

def plot_histogram_normalised(results):
    energies = sorted(results.keys())
    fig, ax = plt.subplots()

    legend = []

    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
    colours = ['blue', 'green', 'red', 'purple', 'cyan', 'orange']

    i = 0
    for energy in energies:
        # results[energy][-1] == 2*pi/(c*lambda), so we need to
        # multiply by sqrt(lambda) to complete the normalisation
        ax.hist(results[energy][2] * results[energy][-1] * n.sqrt(energy),
                linestyle=linestyles[(i % len(linestyles))],
                color = colours[i],
                linewidth=3,
                bins=23,
                normed=True,
                histtype='step',
                )
        legend.append(str(energy))
        i += 1

    l = ax.legend(legend)
    l.set_title('energy')

    ax.set_xlabel('$\\rho$')
    ax.set_ylabel('PDF')

    return fig, ax
