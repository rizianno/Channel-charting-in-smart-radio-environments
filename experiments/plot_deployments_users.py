"""Per-user figures of the EMS deployment comparison (results_aoa/serving_deployments/<noise>/per_user.csv).

Separate per-EMS charts; full-NLOS users; training informed (angle prior) or coarse-to-fine ("nested"
in the CSV), K = 1, 3, 6 training acquisitions. Scored users only (anchors excluded). For each training
and surface, in figures/per_user/<training>/ as PDF and PNG:
  map_<surface>_K<K>   : |AoA error| of each user at its position (0-30 deg, darker = worse), EMS1/EMS2
                         rows, EMS trained alone / joint training columns; EMSs and BS marked
  truevsest_<surface>  : estimated vs true azimuth, K = 1, 3, 6 columns, EMS1/EMS2 rows (joint training,
                         trained alone in grey)
  cdf_aoa_<surface>    : CDF of |AoA error| for K = 1, 3, 6 (solid: trained alone, dashed: joint training)
  cdf_snr_<surface>_K<K> : CDF of the SNR gain vs zero phase at K and 2K sweeps: t-SNE (K sweeps) and
                         t-SNE + verify (2K) vs greedy / hierarchical search over the clusters (K and 2K)
Usage: python experiments/plot_deployments_users.py results_aoa/serving_deployments/est_G60dB
"""
import csv, os, sys, time, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

D = sys.argv[1]
INK, MUTED, GRID = '#222220', '#6b6a63', '#e4e3dc'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.8})
BLUES = LinearSegmentedColormap.from_list('blues', ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
GEO = {'EMS1': (-28.85, 4.0), 'EMS2': (-71.9, 4.0), 'BS': (-40.14, 40.5)}   # results_aoa/full_11x11/per_user_errors_geometry.csv
SURF = {'M60_raytraced': 'real 60x60', 'M60': '60x60 (synthetic)', 'M120': '120x120 (synthetic)',
        'M180': '180x180 (synthetic)', 'M240': '240x240 (synthetic)'}
TRN = {'informed': ('1_informed_prior', 'informed training (angle prior)'),
       'nested': ('2_no_prior_coarse_to_fine', 'coarse-to-fine training (no prior)')}
KC = {1: '#eb6834', 3: '#2a78d6', 6: '#1baf7a'}

R = collections.defaultdict(list)                       # (surface, training, deployment, view, K) -> rows
for x in csv.DictReader(open(os.path.join(D, 'per_user.csv'))):
    if x['scored'] in ('1', 'true'):
        R[(x['surface'], x['training'], x['deployment'], x['view'], int(x['K']))].append(x)
col = lambda rows, c: np.array([float(r[c]) for r in rows])


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    for ext in ('png', 'pdf'):
        for attempt in range(5):                       # OneDrive may briefly lock the file
            try:
                fig.savefig(os.path.join(out, f'{name}.{ext}'), dpi=130, bbox_inches='tight')
                break
            except OSError:
                time.sleep(2)
        else:
            raise OSError('could not save ' + name)
    plt.close(fig)


def deps(view):
    return [(f'{view} alone', f'{view} trained alone'), ('joint', 'joint training (both EMSs present)')]


def cdf(a, v, **kw):
    v = np.sort(v[np.isfinite(v)])
    if v.size:
        a.plot(v, np.arange(1, v.size + 1) / v.size, **kw)


for tr, (folder, tname) in TRN.items():
    out = os.path.join(D, 'figures', 'per_user', folder)
    for surf, sname in SURF.items():
        if not any(k[0] == surf and k[1] == tr for k in R):
            continue
        # ---- spatial error maps
        for K in (1, 3, 6):
            fig, ax = plt.subplots(2, 2, figsize=(11, 9), sharex=True, sharey=True)
            for r, view in enumerate(('EMS1', 'EMS2')):
                for c, (dep, dname) in enumerate(deps(view)):
                    a = ax[r, c];  rows = R.get((surf, tr, dep, view, K), [])
                    if rows:
                        e = np.clip(col(rows, 'aoa_err_deg'), 0, 30)
                        sc = a.scatter(col(rows, 'x_m'), col(rows, 'y_m'), c=e, cmap=BLUES, vmin=0, vmax=30, s=9, lw=0)
                        a.set_title(f'{view}, {dname}\nmedian {np.median(col(rows, "aoa_err_deg")):.1f} deg, '
                                    f'RMSE {np.sqrt(np.mean(col(rows, "aoa_err_deg")**2)):.1f} deg', color=INK, fontsize=9)
                    for k, (gx, gy) in GEO.items():
                        a.plot(gx, gy, 's' if k != 'BS' else '^', color=INK, ms=7)
                        a.annotate(k, (gx, gy), xytext=(4, 4), textcoords='offset points', fontsize=8, color=INK)
                    a.set_aspect('equal');  a.grid(False)
                    if r == 1: a.set_xlabel('x [m]')
                    if c == 0: a.set_ylabel('y [m]')
            cb = fig.colorbar(sc, ax=ax, shrink=0.8, label='|AoA error| [deg] (clipped at 30)')
            fig.suptitle(f'Full-NLOS users, surface {sname}, {tname}, K = {K} acquisitions', color=INK, fontsize=12)
            save(fig, out, f'map_{surf}_K{K}')

        # ---- estimated vs true azimuth
        fig, ax = plt.subplots(2, 3, figsize=(13, 8.5))
        for r, view in enumerate(('EMS1', 'EMS2')):
            for c, K in enumerate((1, 3, 6)):
                a = ax[r, c]
                for dep, colr, lab in ((f'{view} alone', '#c3c2b7', 'trained alone'), ('joint', '#2a78d6', 'joint training')):
                    rows = R.get((surf, tr, dep, view, K), [])
                    if rows:
                        a.scatter(col(rows, 'az_true_deg'), col(rows, 'az_est_deg'), s=5, color=colr, lw=0, label=lab, alpha=0.8)
                a.plot([-90, 90], [-90, 90], ':', color=INK, lw=1)
                a.set_xlim(-90, 90);  a.set_ylim(-90, 90);  a.set_aspect('equal')
                a.set_title(f'{view}, K = {K}', color=INK, fontsize=10)
                if r == 1: a.set_xlabel('true azimuth [deg]')
                if c == 0: a.set_ylabel('estimated azimuth [deg]')
        ax[0, 0].legend(frameon=False, fontsize=8, markerscale=2)
        fig.suptitle(f'Full-NLOS users, surface {sname}, {tname}: estimated vs true azimuth', color=INK, fontsize=12)
        save(fig, out, f'truevsest_{surf}')

        # ---- CDF of the AoA error
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
        for c, view in enumerate(('EMS1', 'EMS2')):
            a = ax[c]
            for K in (1, 3, 6):
                for dep, ls in ((f'{view} alone', '-'), ('joint', '--')):
                    rows = R.get((surf, tr, dep, view, K), [])
                    if rows:
                        cdf(a, col(rows, 'aoa_err_deg'), color=KC[K], ls=ls, lw=1.8)
            a.set_xscale('log');  a.set_xlim(0.1, 180);  a.set_title(view, color=INK, fontsize=10)
            a.set_xlabel('|AoA error| [deg] (log scale)')
            if c == 0: a.set_ylabel('CDF')
        h = [Line2D([], [], color=KC[K], lw=1.8, label=f'K = {K}') for K in (1, 3, 6)] + \
            [Line2D([], [], color=MUTED, ls='-', label='trained alone'), Line2D([], [], color=MUTED, ls='--', label='joint training')]
        fig.legend(handles=h, loc='upper center', ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.02))
        fig.suptitle(f'Full-NLOS users, surface {sname}, {tname}: CDF of the AoA error', color=INK, fontsize=12, y=1.08)
        save(fig, out, f'cdf_aoa_{surf}')

        # ---- CDF of the SNR gain at equal sweeps
        for K in (1, 3, 6):
            fig, ax = plt.subplots(2, 2, figsize=(12, 8.5), sharex=True, sharey=True)
            for r, (deplab) in enumerate(('alone', 'joint')):
                for c, view in enumerate(('EMS1', 'EMS2')):
                    a = ax[r, c];  dep = f'{view} alone' if deplab == 'alone' else 'joint'
                    rows = R.get((surf, tr, dep, view, K), [])
                    if not rows:
                        continue
                    for cname, colr, ls, lab in (('dsnr_lut_dB', '#2a78d6', '-', f't-SNE ({K} sweeps)'),
                                                 ('dsnr_greedy_cl_K_dB', '#008300', '-', f'greedy, clusters ({K})'),
                                                 ('dsnr_hier_cl_K_dB', '#6250d6', '-', f'hierarchical, clusters ({K})'),
                                                 ('dsnr_verify_dB', '#2a78d6', '--', f't-SNE + verify ({2*K} sweeps)'),
                                                 ('dsnr_greedy_cl_2K_dB', '#008300', '--', f'greedy, clusters ({2*K})'),
                                                 ('dsnr_hier_cl_2K_dB', '#6250d6', '--', f'hierarchical, clusters ({2*K})')):
                        cdf(a, col(rows, cname), color=colr, ls=ls, lw=1.7, label=lab)
                    a.set_title(f'{view}, ' + ('trained alone' if deplab == 'alone' else 'joint training'), color=INK, fontsize=10)
                    if r == 1: a.set_xlabel('SNR gain vs zero phase [dB]')
                    if c == 0: a.set_ylabel('CDF')
            ax[0, 0].legend(frameon=False, fontsize=8, loc='upper left')
            fig.suptitle(f'Full-NLOS users, surface {sname}, {tname}: SNR gain at {K} and {2*K} sweeps', color=INK, fontsize=12)
            save(fig, out, f'cdf_snr_{surf}_K{K}')
    print(folder)
