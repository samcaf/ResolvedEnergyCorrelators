import numpy as np
from scipy.optimize import curve_fit

from encs.plots import plot_runtime
from plotter import enc_data_dir, enc_figure_dir, combine_plotters

from utils.plot_utils import adjust_lightness

# =====================================
# Flags
# =====================================
# Whether to make OD plots
opendata = True

# Which Pythia plots to make
pythia = ''
# pythia = 'qcd'
# npyth  = '1M'


# =====================================
# Plot parameters
# =====================================
# 1D plots
runtime_plot_params = {
    'axes.labelsize': 20,
    'ylim': (1e-4, 5e3),
    'xlim': (5, 120),
    'y_scale': 'log',
    'x_scale': 'log'
}
fill_between = True

# Linear near 0
# runtime_plot_params = {
#     'axes.labelsize': 20,
#     'ylim': (0, 1e1),
#     'xlim': (0, 20),
#     'y_scale': 'lin',
#     'x_scale': 'lin'
# }

# For fill_between
color_1d = 'cornflowerblue'
t_ind_1d = 5
color_3d = 'firebrick'
t_ind_3d = 35
color_5d = 'mediumorchid'
t_ind_5d = 40
fill_alpha = 0.2

# For estimating scaling
def fit_func(x, a, b, c):
    return a*pow(x, b) + c
# TODO: Include log?

# =====================================
# Main
# =====================================
if __name__ == "__main__":
    # =====================================
    # Opendata Plots
    # =====================================
    if opendata:
        # ---------------------------------
        # 2 particles:
        # ---------------------------------
        hist_1d, od_1d = plot_runtime(
                file_name=enc_data_dir/
                        '2particle_od_100k_200bins_nu1-00.py',
                color=color_1d,
                label='Projected ENC (1 Particle)',
                save=None,
                **runtime_plot_params,
        )
        # Preparing for fill_between
        means_1d = np.array(hist_1d.metadata['runtime_means'])/1000
        stds_1d  = np.array(hist_1d.metadata['runtime_stds'])/1000
        nans     = np.isnan(means_1d) * np.isnan(stds_1d)
        # Removing nans
        nums_1d  = np.arange(0, len(means_1d))[~nans]
        means_1d = means_1d[~nans]
        stds_1d  = stds_1d[~nans]
        # Getting slope
        fit_1d, _ = curve_fit(fit_func, nums_1d, means_1d)

        # ---------------------------------
        # 3 particles:
        # ---------------------------------
        hist_3d, od_3d = plot_runtime(
                file_name=enc_data_dir/
                        '3particle_od_100k_150bins_nus_1-00_1-00.py',
                color=color_3d,
                label='Projected ENC (2 Particles)',
                save=None,
                **runtime_plot_params,
        )
        # Preparing for fill_between
        means_3d = np.array(hist_3d.metadata['runtime_means'])/1000
        stds_3d  = np.array(hist_3d.metadata['runtime_stds'])/1000
        nans     = np.isnan(means_3d) * np.isnan(stds_3d)
        # Removing nans
        nums_3d  =np.arange(0, len(means_3d))[~nans]
        means_3d = means_3d[~nans]
        stds_3d  = stds_3d[~nans]
        # Getting slope
        fit_3d, _ = curve_fit(fit_func, nums_3d, means_3d)

        # ---------------------------------
        # 4 particles:
        # ---------------------------------
        hist_5d, od_5d = plot_runtime(
                file_name=enc_data_dir/
                        '4particle_od_10k_10bins_nus_1-00_1-00_1-00.py',
                color=color_5d,
                alpha=0.2,
                label='Projected ENC (3 Particles)',
                save=None,
                **runtime_plot_params,
        )
        # Preparing for fill_between
        means_5d = np.array(hist_5d.metadata['runtime_means'])/1000
        stds_5d  = np.array(hist_5d.metadata['runtime_stds'])/1000
        nans     = np.isnan(means_5d) * np.isnan(stds_5d)
        # Removing nans
        nums_5d  =np.arange(0, len(means_5d))[~nans]
        means_5d = means_5d[~nans]
        stds_5d  = stds_5d[~nans]
        # Getting slope
        fit_5d, _ = curve_fit(fit_func, nums_5d, means_5d)

        # ---------------------------------
        # Plotting:
        # ---------------------------------
        # Combining plots
        od_plotters = [od_1d, od_3d, od_5d]

        cplotter = combine_plotters(od_plotters)

        # Adding 0,0:
        # nums_1d, means_1d, stds_1d = np.array([0, *nums_1d]), \
        #             np.array([0, *means_1d]), np.array([0, *stds_1d])
        # nums_3d, means_3d, stds_3d = np.array([0, *nums_3d]), \
        #             np.array([0, *means_3d]), np.array([0, *stds_3d])
        # nums_5d, means_5d, stds_5d = np.array([0, *nums_5d]), \
        #             np.array([0, *means_5d]), np.array([0, *stds_5d])
        # Filling standard deviations
        if fill_between:
            cplotter.axes[0].fill_between(nums_1d,
                    means_1d-stds_1d, means_1d+stds_1d,
                    color=color_1d, alpha=fill_alpha)
            cplotter.axes[0].fill_between(nums_3d,
                    means_3d-stds_3d, means_3d+stds_3d,
                    color=color_3d, alpha=fill_alpha)
            cplotter.axes[0].fill_between(nums_5d,
                    means_5d-stds_5d, means_5d+stds_5d,
                    color=color_5d, alpha=fill_alpha)

        # Adding slope information
        poly_1d = fit_1d[1]
        poly_3d = fit_3d[1]
        poly_5d = fit_5d[1]
        cplotter.axes[0].text(nums_1d[t_ind_1d], means_1d[t_ind_1d]/3,
                              rf"$t \sim M^{{{str(poly_1d)[:4]}}}$",
                              color=adjust_lightness(color_1d,0.5),
                              size=18)
        cplotter.axes[0].text(nums_3d[t_ind_3d], means_3d[t_ind_3d]/3,
                              rf"$t \sim M^{{{str(poly_3d)[:4]}}}$",
                              color=adjust_lightness(color_3d,0.5),
                              size=18)
        cplotter.axes[0].text(nums_5d[t_ind_5d], means_5d[t_ind_5d]*8,
                              rf"$t \sim M^{{{str(poly_5d)[:4]}}}$",
                              color=adjust_lightness(color_5d,0.5),
                              size=18)

        # Formatting and saving
        cplotter.axes[0].set_xlabel(r'Number of Particles ($M$)')
        cplotter.axes[0].legend(loc=(0.38, 0.02))
        cplotter.fig.tight_layout()
        cplotter.savefig('opendata_runtimes.pdf',
                         enc_figure_dir/'supplementary/')

    if pythia:
        # =====================================
        # Pythia Plots
        # =====================================
        pythia_1d = plot_1d_enc(
                file_name=enc_data_dir/
                        f'2particle_{pythia}_{npyth}_nu1-00.py',
                color='cornflowerblue',
                alpha=0.2,
                label='$\nu=1$',
                save=None,
                **runtime_plot_params,
        )

        # Combining plots
        # pythia_plotters = [pythia_1d, pythia_2d, pythia_3d]

        # cplotter = combine_plotters(pythia_plotters)
        # cplotter.axes[0].legend(loc=(0.02, 0.3))
        # cplotter.fig.tight_layout()
        # cplotter.savefig(f'pythia_runtimes.pdf',
        #                  enc_figure_dir/'supplementary/')
