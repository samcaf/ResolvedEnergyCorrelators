import numpy as np

from re import search

from histogram import HistogramData, HistogramData2D

from plotter import enc_data_dir, enc_figure_dir
from encs.plots import plot_2d_density
from encs.plots import density_colormap

from qcd.qcd_basics import alpha_s


# =====================================
# Flags
# =====================================
# Whether to plot OD
opendata = True
draw_np = True

# Which pythia data to plot
pythia = ['qcd', 'w', 'top']
pythia = [None]
npyth  = '1M_150bins'


# =====================================
# Plot parameters
# ====================================
# 2D density plots
new_density = {
    'axes.labelsize': 20,
    'ylim': (0, 1.0),
    'xlim': (5e-3, 8e-1),
    'xlabel': r'$R_1$',
    'ylabel': r'$R_2/R_1$',
    'x_scale': 'log',
    'cbar_cmap': density_colormap['opendata']
}
old_density = {
    'axes.labelsize': 20,
    'ylim': (0, 1.0),
    'xlim': (5e-3, 8e-1),
    'xlabel': r'$R_L$',
    'ylabel': r'$R_S/R_L$',
    'x_scale': 'log',
    'cbar_cmap': density_colormap['opendata']
}


def stamp_density(ax, **metadata):
    # Dataset information
    if metadata['level'] == 'data':
        ax.text(-0.00, 1.10,
                r'$\textbf{CMS Open Data}$: '\
                '2011A Jet Primary Dataset',
                fontsize=13, transform=ax.transAxes)
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
        ax.text(0.0, 1.10, process_str,
                fontsize=13, transform=ax.transAxes)

    # Getting jet radius (more complicated than probably needs to be)
    jet_rad = 10*float(metadata['jet_rad'])
    rnd_digit = len(search("\.(0*)", str(jet_rad)).group(1))
    jet_rad = round(jet_rad, rnd_digit)
    jet_rad_str = metadata['jet_alg'].upper()
    if jet_rad_str in ['ANTI-KT', 'ANTIKT', 'AKT']:
        jet_rad_str = 'AK'

    jet_rad_str += str(jet_rad)[:-2] if str(jet_rad)[-2:] == '.0' \
                   else str(jet_rad)
    # Jet information
    ax.text(0.0, 1.04, jet_rad_str+' Jets,  '\
            r'$\!|\eta^\text{jet}| <$ 1.9, '\
            r'$p_T^\text{jet} \in\,$[500, 550] GeV',
            fontsize=11, transform=ax.transAxes)

    # Axes labels
    ax.text(1.02, 1.105, r'$\phi_2$-integrated', fontsize=14,
            transform=ax.transAxes)
    ax.text(1.08, 1.025, r'$\textbf{RE3C}$', fontsize=20,
            transform=ax.transAxes)



if __name__ == "__main__":
    # ==========================================
    # Open Data
    # ==========================================
    if opendata:
        # =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        # Plots with new variables
        # =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        new_hist3d = HistogramData(
            file_name=enc_data_dir/
                  '3particle_od_100k_150bins_nus_1-00_1-00.py'
        )
        new_hist2d = plot_2d_density(
            hist_data=new_hist3d,
            vmax=1, log_colorbar=False,
            save=None, **new_density)
        ax = new_hist2d.density.axes[0]
        # Stamp
        stamp_density(ax, **new_hist2d.metadata)
        # Non-perturbative portion:
        pt = 500
        # Lambda_high = 10/pt
        # Lambda_low  = 5/pt
        # Lambda_mid  = 10**(.7*np.log10(Lambda_high)+
        #                    .3*np.log10(Lambda_low))
        # theta1s = np.logspace(np.log10(Lambda_mid), np.log10(1/2),
        #                       100000)  # many points to shade fully
        Lambda_mid  = 5/pt
        theta1s = np.logspace(-3, np.log10(1),
                              100)

        # Fudged version of naive physics because of some phase
        # space constraints that I coded up in Mathematica after
        # pulling them out of my booty
        # Basically, the place where i2 approaches within LambdaQCD of
        # i1 is not exactly where the NP effects kick in -- you
        # need to take into account when they become comparable
        # to the other physics (when R_2-R_1 is equal to Lambda,
        # the NP physics starts but the associated region of phase
        # space is measure zero. As R_2 increases even further,
        # _then_ the NP effects begin to become more important)
        if draw_np:
            ax.plot(theta1s, Lambda_mid/theta1s,
                    lw=5, ls='dashed',
                    color='darkgoldenrod')
            ax.plot(theta1s, 1-Lambda_mid/2/theta1s,
                    lw=5, ls='dashed',
                    color='mediumseagreen')
            # ax.fill_between(theta1s,
            #             Lambda_low/theta1s,
            #             [Lambda_high/t1 if Lambda_high/t1 < 1-Lambda_high/2/t1
            #              else 1-Lambda_low/2/t1 for t1 in theta1s],
            #             alpha=0.25, color='firebrick',
            #             lw=0)
            # np_high_inds = [i for i, t1 in enumerate(theta1s) if
            #             Lambda_high/t1 < 1-Lambda_high/2/t1]
            # ax.fill_between(theta1s[np_high_inds],
            #             1-Lambda_high/2/theta1s[np_high_inds],
            #             1-Lambda_low/2/theta1s[np_high_inds],
            #             alpha=0.25, color='firebrick',
            #             lw=0)

        new_hist2d.density.fig.tight_layout()
        new_hist2d.density.savefig(
            str(enc_figure_dir/'od_newdef_density.pdf'))

        # Trying something analytic
        """
        new_hist3d = HistogramData(
            file_name=enc_data_dir/
                  '3particle_od_100k_150bins_nus_1-00_1-00.py'
        )
        Rs, Ts, Phis = np.meshgrid(new_hist3d.centers['theta1'],
                            new_hist3d.centers['theta2_over_theta1'],
                            new_hist3d.centers['phi'])

        def law_of_cosines(adj1, adj2, angle):
            c2 = adj1**2. + adj2**2. - 2.*adj1*adj2*np.cos(angle)
            return np.sqrt(c2)

        minR = np.minimum(Rs*Ts, law_of_cosines(Rs, Rs*Ts, Phis))

        pertdensity = 1 \
            * alpha_s(Rs*500, freeze_at=.12) \
            * alpha_s(minR*500, freeze_at=.12) \
            * 1/(minR) \
            * (4/3)**2. / (4*np.pi)**2. / 10000
        # pertdensity = minR

        new_hist3d.hist = pertdensity
        new_hist3d.variable_order = ['theta1', 'theta2_over_theta1',
                                     'phi']
        new_hist2d = HistogramData2D(
                    hist_data=new_hist3d.get_sub_histogram('phi', 0))
        new_hist2d.make_plot('density')
        # new_hist2d = plot_2d_density(
        #     hist_data=new_hist3d,
        #     vmin=0, vmax=1, log_colorbar=False,
        #     save=None, **new_density)
        ax = new_hist2d.density.axes[0]
        # Stamp
        stamp_density(ax, **new_hist2d.metadata)

        new_hist2d.density.fig.tight_layout()
        new_hist2d.density.savefig(
            str(enc_figure_dir/'od_newdef_density.pdf'))
        """


        # =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        # Plots with old variables
        # =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        # Basic
        old_hist3d = HistogramData(
            file_name=enc_data_dir/
                      'old_3particle_od_100k_150bins.py'
        )
        old_hist2d = plot_2d_density(
            hist_data=old_hist3d,
            vmax=2, log_colorbar=False,
            save=None, **old_density)
        ax = old_hist2d.density.axes[0]
        stamp_density(ax, **old_hist2d.metadata)
        old_hist2d.density.fig.tight_layout()
        old_hist2d.density.savefig(
            str(enc_figure_dir/'od_olddef_density.pdf'))


    # ==========================================
    # Pythia
    # ==========================================
    for process in np.atleast_1d(pythia):
        if process is None:
            continue
        pythia_new_density = new_density.copy()
        pythia_new_density['cbar_cmap'] = density_colormap[process]
        pythia_old_density = old_density.copy()
        pythia_old_density['cbar_cmap'] = density_colormap[process]



        # =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        # Plots with new variables
        # =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        new_hist3d = HistogramData(
            file_name=enc_data_dir/
                  f'3particle_{process}_{npyth}_nus_1-00_1-00.py'
        )
        new_hist2d = plot_2d_density(
            hist_data=new_hist3d,
            vmax=1, log_colorbar=False,
            save=None,
            **pythia_new_density)
        ax = new_hist2d.density.axes[0]
        stamp_density(ax, **new_hist2d.metadata)
        new_hist2d.density.fig.tight_layout()
        new_hist2d.density.savefig(
            str(enc_figure_dir/f'supplementary/density/{process}_newdef_density.pdf'))

        # # =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        # # Plots with old variables
        # # =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
        # # Basic
        # old_pythiahist3d = HistogramData(
        #     file_name=enc_data_dir/
        #           f'old_3particle_{process}_{npyth}.py'
        # )

        # old_pythiahist2d = plot_2d_density(
        #     hist_data=old_pythiahist3d,
        #     vmax=2, log_colorbar=False,
        #     save=None,
        #     **pythia_old_density)
        # ax = old_pythiahist2d.density.axes[0]
        # stamp_density(ax, **old_pythiahist2d.metadata)
        # old_pythiahist2d.density.fig.tight_layout()
        # old_pythiahist2d.density.savefig(
        #     str(enc_figure_dir/f'supplementary/{process}_olddef_density.pdf'))

        # Symmetrized
        # symm_old_hist3d = HistogramData(hist_data=old_hist3d)
        # symm_old_hist3d.hist = (old_hist3d.hist[:,:,:] +
        #                         old_hist3d.hist[:,::-1,:])/2

        # symm_old_hist2d = plot_2d_density(
        #     hist_data=symm_old_hist3d,
        #     vmax=1, log_colorbar=False,
        #     save=str(enc_figure_dir/'supplementary/olddef_symm'),
        #     **old_density_color)

        # # Symmetrized - new
        # symm_old_hist3d.hist = new_hist3d.hist - symm_old_hist3d.hist
        # diff_hist2d = plot_2d_density(
        #     hist_data=symm_old_hist3d,
        #     vmin=-0.5, vmax=0.5,
        #     log_colorbar=False,
        #     save=str(enc_figure_dir/'supplementary/new_minus_oldsymm'),
        #     **old_density_color)


# TODO:
# Indicate np region
