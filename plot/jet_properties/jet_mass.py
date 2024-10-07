import numpy as np

from histogram import HistogramData1D
from plotter import property_data_dir, property_figure_dir
from plotter import combine_plotters


# =====================================
# Flags
# =====================================
# Whether to make OD plots
opendata = True

# Which pythia data to plot
pythia = []
# pythia.append('qcd')
pythia.append('w')
pythia.append('top')
npyth  = '10k_200bins'


# =====================================
# Plot parameters
# =====================================
# 1D plots
plot_1d_params = {
    key: {
        'axes.labelsize': 20,
        'ylim': (1e-5, 8e0),
        # 'ylim': (0, 7e-1),
        'xlim': (1e-5, 1e3),
        # 'xlim': (1e-4, 4e0),
        # 'xlim': (0, 3e0),
        'xlabel': r'$R_1$',
        'ylabel':
        r'$\frac{\text{d}\Sigma}{\text{d}\log_{10}m}$',
        # 'y_scale': 'log',
        'y_scale': 'lin',
        'x_scale': 'log',
        # 'x_scale': 'lin',
    }
    for key in ['opendata', 'qcd', 'w', 'top']
}


def stamp_1d(ax, **metadata):
    # Getting jet radius (more complicated than probably needs to be)
    jet_info =  r'$\!|\eta^\text{jet}| <$ 1.9'
    jet_info += r', $p_T^\text{jet} \in\,$[500, 550] GeV'

    # Jet information
    ax.text(0.025, 0.87, jet_info, fontsize=12,
            transform=ax.transAxes)

    # Dataset information
    if metadata['level'] == 'data':
        ax.text(0.025, 0.93,
                r'$\textbf{CMS Open Data}:$ '\
                '2011A Jet Primary Dataset',
                fontsize=14, transform=ax.transAxes)
    else:
        process = metadata['outstate_str']
        if process is None or not process in ['qcd', 'w', 'top']:
            raise ValueError(f"Must specify a pythia process in "
                             "[qcd, w, top] when opendata=False.")
        process_str  = r'$\texttt{Pythia 8.310},\,\,pp\to\,\,$'
        process_str += 'hadrons, ' if process == 'qcd' else \
                       (r'$W^+\,W^-\!,$ ' if process == 'w'
                        else r'$t\overline{t},$ ')
        process_str += r'$\sqrt{s}=\,$14 TeV'
        ax.text(0.025, 0.93, process_str,
                fontsize=14, transform=ax.transAxes)


def plot_jet_mass(file_name=None, hist_data=None, **kwargs):
    # Loading and validating histogram
    hist1d = HistogramData1D(file_name=file_name,
                           hist_data=hist_data)
    print(hist1d)

    # Plotting
    if kwargs:
        # Preparing file name for figure
        save_name = kwargs.pop('save', None)

        # Making plot
        hist1d.make_plot(**kwargs)

        # Saving
        if save_name:
            hist1d.plot.savefig(save_name+'.pdf',
                                property_figure_dir)

    return hist1d



# =====================================
# Main
# =====================================
if __name__ == "__main__":
    # =====================================
    # Opendata Plots
    # =====================================
    if opendata:
        pass

    # =====================================
    # Pythia Plots
    # =====================================
    for process in pythia:
        if process is None:
            continue

        pythia_plotters = []

        pythia_0 = plot_jet_mass(
                file_name=property_data_dir/
                        f'{process}_AKT6_{npyth}_mass.py',
                plot_kwargs=plot_1d_params,
                color='mediumseagreen',
                label=r'AKT6',
                save=None,
                **plot_1d_params[process],
        )
        pythia_plotters.append(pythia_0.plot)

        pythia_0 = plot_jet_mass(
                file_name=property_data_dir/
                        f'{process}_AKT8_{npyth}_mass.py',
                plot_kwargs=plot_1d_params,
                color='cornflowerblue',
                label=r'AKT8',
                save=None,
                **plot_1d_params[process],
        )
        pythia_plotters.append(pythia_0.plot)

        pythia_0 = plot_jet_mass(
                file_name=property_data_dir/
                        f'{process}_AKT10_{npyth}_mass.py',
                plot_kwargs=plot_1d_params,
                color='mediumorchid',
                label=r'AKT10',
                save=None,
                **plot_1d_params[process],
        )
        pythia_plotters.append(pythia_0.plot)

        pythia_0 = plot_jet_mass(
                file_name=property_data_dir/
                        f'{process}_AKT12_{npyth}_mass.py',
                plot_kwargs=plot_1d_params,
                color='lightcoral',
                label=r'AKT12',
                save=None,
                **plot_1d_params[process],
        )
        pythia_plotters.append(pythia_0.plot)

        # Combining plots
        cplotter = combine_plotters(pythia_plotters)
        stamp_1d(cplotter.axes[0], **pythia_0.metadata)
        cplotter.axes[0].legend(loc=(0.03, 0.42))
        cplotter.fig.tight_layout()
        cplotter.savefig(f'{process}_jet_mass.pdf',
                         property_figure_dir)
