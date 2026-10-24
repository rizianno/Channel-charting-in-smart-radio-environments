"""CDF of the t-SNE AoA error per EMS, every K (results_aoa/serving_deployments/<noise>/per_user_aoa.csv).

Each EMS trained alone, full-NLOS users it reaches, scored users only (anchors excluded). One figure
per EMS and surface: uninformed training (no prior: the fixed sweeping order, "nested" in the CSV)
next to informed training (angle prior from the anchors); one curve per K, coloured by K. Dotted:
one beamwidth at broadside. Surfaces are named by their size only (60x60 is the ray-traced one; the
model-based 60x60 is a check and is not drawn).
Output: <results>/figures/aoa_cdf/cdf_aoa_<EMS>_<surface>.png/.pdf
Usage: python experiments/plot_aoa_cdf.py results_aoa/serving_deployments/est_G60dB
"""
import csv, os, sys, time, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

D = sys.argv[1]
OUT = os.path.join(D, 'figures', 'aoa_cdf');  os.makedirs(OUT, exist_ok=True)
INK, MUTED, GRID = '#222220', '#6b6a63', '#e4e3dc'
plt.rcParams.update({'font.size': 10, 'font.family': 'Calibri', 'mathtext.fontset': 'cm', 'axes.edgecolor': MUTED,
                     'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False})
SURF = {'M60_raytraced': (60, '60x60 EMS'), 'M120': (120, '120x120 EMS'), 'M180': (180, '180x180 EMS'), 'M240': (240, '240x240 EMS')}
TRN = [('nested', 'uninformed training'), ('informed', 'informed training (angle prior)')]
KS = list(range(1, 16))
CMAP = ListedColormap(plt.cm.viridis(np.linspace(0.0, 0.88, len(KS))))
NORM = BoundaryNorm(np.arange(0.5, len(KS) + 1), CMAP.N)

E = collections.defaultdict(list)                       # (surface, view, training, K) -> errors of scored users
for r in csv.DictReader(open(os.path.join(D, 'per_user_aoa.csv'))):
    if r['scored'] in ('1', 'true'):
        E[(r['surface'], r['view'], r['training'], int(r['K']))].append(float(r['aoa_err_deg']))

for surf, (M, sname) in SURF.items():
    bw = np.degrees(np.arcsin(0.886 * 2 / M))          # one beamwidth at broadside [deg]
    for view in ('EMS1', 'EMS2'):
        if not any(k[0] == surf and k[1] == view for k in E):
            continue
        fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.9), sharey=True)
        for c, (tr, tname) in enumerate(TRN):
            a = ax[c]
            for i, K in enumerate(KS):
                v = np.sort(np.array(E.get((surf, view, tr, K), [])))
                if v.size:
                    v = np.maximum(v, 1e-2)
                    a.step(v, np.arange(1, v.size + 1) / v.size, where='post', color=CMAP(i), lw=1.6)
            a.axvline(bw, color=INK, ls=':', lw=1)
            a.text(bw * 1.08, 0.04, f'1 beamwidth\n({bw:.2f}°)', fontsize=8.5, color=MUTED, va='bottom')
            a.set_xscale('log');  a.set_xlim(0.05, 180);  a.set_ylim(0, 1.0)
            a.set_xticks([0.1, 1, 10, 100]);  a.set_xticklabels(['0.1', '1', '10', '100'])
            a.grid(True, color=GRID, lw=0.7, which='major')
            a.set_title(tname, fontsize=11, color=INK)
            a.set_xlabel('|AoA error| [deg]', fontsize=10)
            if c == 0:
                a.set_ylabel('CDF over the users', fontsize=10)
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=NORM, cmap=CMAP), ax=ax, ticks=[1, 3, 5, 7, 9, 11, 13, 15],
                          fraction=0.025, pad=0.02)
        cb.set_label('K (training slots)', fontsize=10);  cb.ax.tick_params(labelsize=9)
        n = len(E.get((surf, view, 'informed', 1), []))
        fig.suptitle(f'{view}, {sname}: CDF of the t-SNE AoA error ({n} users)', fontsize=12, color=INK, x=0.45, y=1.02)
        for ext in ('png', 'pdf'):
            for attempt in range(5):                       # OneDrive may briefly lock the file
                try:
                    fig.savefig(os.path.join(OUT, f'cdf_aoa_{view}_{surf}.{ext}'), dpi=200, bbox_inches='tight')
                    break
                except OSError:
                    time.sleep(2)
            else:
                raise OSError(f'could not save cdf_aoa_{view}_{surf}.{ext}')
        plt.close(fig)
        print(view, surf)
