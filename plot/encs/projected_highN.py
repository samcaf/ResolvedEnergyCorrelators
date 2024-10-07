import numpy as np

from re import search

from encs.plots import plot_1d_enc

from plotter import enc_data_dir, enc_figure_dir, combine_plotters

# =====================================
# Flags
# =====================================
# Whether to make OD plots
opendata = True

# Which pythia data to plot
pythia = ['qcd']
npyth  = '1M_500bins'


# =====================================
# Plot parameters
# =====================================
# 1D plots
plot_1d_params = {
    key: {
        'axes.labelsize': 20,
        'ylim': (1e-5, 8e0),
        # 'ylim': (0, 7e-1),
        'xlim': (1e-5, 2e0),
        # 'xlim': (1e-4, 4e0),
        # 'xlim': (0, 3e0),
        'xlabel': r'$R_1$',
        'ylabel':
        r'$\frac{\text{d}\Sigma}{\text{d}\log_{10}R_1}$',
        # 'y_scale': 'log',
        'y_scale': 'lin',
        'x_scale': 'log',
        # 'x_scale': 'lin',
    }
    for key in ['opendata', 'qcd', 'w', 'top']
}

plot_1d_params['opendata']['xlim'] = (1e-5, 2e0)
# plot_1d_params['opendata']['xlim'] = (1e-5, 2e0)
# plot_1d_params['opendata']['xlim'] = (0, 1.25)


def stamp_1d(ax, **metadata):
    # Getting jet radius (more complicated than probably needs to be)
    jet_rad = 10*float(metadata['jet_rad'])
    rnd_digit = len(search("\.(0*)", str(jet_rad)).group(1))
    jet_rad = round(jet_rad, rnd_digit)
    jet_info = metadata['jet_alg'].upper()
    if jet_info in ['ANTI-KT', 'ANTIKT', 'AKT']:
        jet_info = 'AK'

    jet_info += str(jet_rad)[:-2] if str(jet_rad)[-2:] == '.0' \
                   else str(jet_rad)
    jet_info += r' Jets, $\!|\eta^\text{jet}| <$ 1.9'
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


# =====================================
# Main
# =====================================
if __name__ == "__main__":
    # =====================================
    # Opendata Plots
    # =====================================
    if opendata:
        od_plotters = []
        od_0 = plot_1d_enc(
                file_name=enc_data_dir/
                        '2particle_od_100k_200bins_nu9-00.py',
                variable_order=['theta1'],
                color='mediumseagreen',
                label=r'$N=10$',
                save=None,
                **plot_1d_params['opendata'],
        )
        od_plotters.append(od_0.plot)
        od_1 = plot_1d_enc(
                file_name=enc_data_dir/
                        '2particle_od_100k_200bins_nu19-00.py',
                plot_kwargs=plot_1d_params,
                variable_order=['theta1'],
                color='cornflowerblue',
                label=r'$N=20$',
                save=None,
                **plot_1d_params['opendata'],
        )
        od_plotters.append(od_1.plot)
        od_2 = plot_1d_enc(
                file_name=enc_data_dir/
                        '2particle_od_100k_200bins_nu49-00.py',
                plot_kwargs=plot_1d_params,
                variable_order=['theta1'],
                color='mediumorchid',
                label=r'$N=50$',
                save=None,
                **plot_1d_params['opendata'],
        )
        od_plotters.append(od_2.plot)
        od_3 = plot_1d_enc(
                file_name=enc_data_dir/
                        '2particle_od_100k_200bins_nu99-00.py',
                plot_kwargs=plot_1d_params,
                variable_order=['theta1'],
                color='lightcoral',
                label=r'$N=100$',
                save=None,
                **plot_1d_params['opendata'],
        )
        od_plotters.append(od_3.plot)
        od_4 = plot_1d_enc(
                file_name=enc_data_dir/
                        '2particle_od_100k_200bins_nu199-00.py',
                plot_kwargs=plot_1d_params,
                variable_order=['theta1'],
                color='sandybrown',
                label=r'$N=200$',
                save=None,
                **plot_1d_params['opendata'],
        )
        od_plotters.append(od_4.plot)

        # Combining plots
        cplotter = combine_plotters(od_plotters)
        stamp_1d(cplotter.axes[0], **od_0.metadata)
        cplotter.axes[0].legend(loc=(0.03, 0.43))
        # cplotter.axes[0].legend(loc='center right')
        cplotter.fig.tight_layout()
        cplotter.savefig(
            f'od_highN_1d.pdf',
            enc_figure_dir/'supplementary/2particle/')


    # =====================================
    # Pythia Plots
    # =====================================
    for process in pythia:
        if process is None:
            continue

        # pythia_plotters = []

        # pythia_0 = plot_1d_enc(
        #         file_name=enc_data_dir/
        #                 f'2particle_{process}_{npyth}_nu9-00.py',
        #         plot_kwargs=plot_1d_params,
        #         variable_order=['theta1'],
        #         color='mediumseagreen',
        #         label=r'$N=10$',
        #         save=None,
        #         **plot_1d_params[process],
        # )
        # pythia_plotters.append(pythia_0.plot)
        # pythia_1 = plot_1d_enc(
        #         file_name=enc_data_dir/
        #                 f'2particle_{process}_{npyth}_nu99-00.py',
        #         plot_kwargs=plot_1d_params,
        #         variable_order=['theta1'],
        #         color='cornflowerblue',
        #         label=r'$N=10^2$',
        #         save=None,
        #         **plot_1d_params[process],
        # )
        # pythia_plotters.append(pythia_1.plot)
        # pythia_2 = plot_1d_enc(
        #         file_name=enc_data_dir/
        #                 f'2particle_{process}_{npyth}_nu999-00.py',
        #         plot_kwargs=plot_1d_params,
        #         variable_order=['theta1'],
        #         color='mediumorchid',
        #         label=r'$N=10^3$',
        #         save=None,
        #         **plot_1d_params[process],
        # )
        # pythia_plotters.append(pythia_2.plot)
        # pythia_3 = plot_1d_enc(
        #         file_name=enc_data_dir/
        #                 f'2particle_{process}_{npyth}_nu9999-00.py',
        #         plot_kwargs=plot_1d_params,
        #         variable_order=['theta1'],
        #         color='lightcoral',
        #         label=r'$N=10^4$',
        #         save=None,
        #         **plot_1d_params[process],
        # )
        # pythia_plotters.append(pythia_3.plot)
        # pythia_4 = plot_1d_enc(
        #         file_name=enc_data_dir/
        #                 f'2particle_{process}_{npyth}_nu99999-00.py',
        #         plot_kwargs=plot_1d_params,
        #         variable_order=['theta1'],
        #         color='sandybrown',
        #         save=None,
        #         label=r'$N=10^5$',
        #         **plot_1d_params[process],
        # )
        # pythia_plotters.append(pythia_4.plot)

        # # Combining plots
        # cplotter = combine_plotters(pythia_plotters)
        # stamp_1d(cplotter.axes[0], **pythia_0.metadata)
        # cplotter.axes[0].legend(loc=(0.03, 0.42))
        # # cplotter.axes[0].legend(loc='center right')
        # cplotter.fig.tight_layout()
        # cplotter.savefig(
        #     f'{process}_highN_1d.pdf',
        #     enc_figure_dir/'supplementary/2particle/')
