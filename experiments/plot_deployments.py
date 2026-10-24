"""Figures of the EMS deployment comparison (results_aoa/serving_deployments/<noise>/summary.csv).

One figure per user subset (full NLOS, reflection only) and surface. Rows: median and RMSE of the
pointing error of the chosen beam, median SNR gain. Columns: EMS1 alone, EMS2 alone, joint training
seen from EMS1, from EMS2, and both EMSs serving together (angles pooled over both EMSs).
Figures go to one subfolder per training assumption (solid: t-SNE pure and hybrid; dashed: searches):
  0_all_combined            : informed t-SNE + all searches, plus pure random / coarse-to-fine t-SNE
  1_informed_prior          : user angle clusters known (anchors): informed training, hybrids and
                              searches restricted to the clusters
  2_no_prior_coarse_to_fine : no prior: coarse-to-fine training (one loop over N/2, N, N/4, 3N/4, ...;
                              called "nested" in the CSV), hybrids and searches over the full codebook
  3_random_training         : no prior: random training codewords (mean of 10 draws), hybrids and
                              searches over the full codebook
Dotted blue (angle rows): the error of the t-SNE AoA estimate itself.
Usage: python experiments/plot_deployments.py results_aoa/serving_deployments/est_G60dB
"""
import csv, os, sys, time, collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

D = sys.argv[1]
INK, MUTED, GRID = '#222220', '#6b6a63', '#e4e3dc'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.8})

T = collections.defaultdict(dict)
for x in csv.DictReader(open(os.path.join(D, 'summary.csv'))):
    k = (x['population'], x['surface'], x['deployment'], x['view'], x['group'])
    T[k][(x['training'], x['scheme'], int(x['K']), int(x['R']), int(x['sweeps']))] = x

SW = list(range(1, 16))


def tsne_series(tr, sfx, name):
    """t-SNE pure and hybrids for one training; sfx = '[clusters]' (prior known) or '' (full codebook)."""
    return [(f't-SNE (LUT beam)',                 tr, 'tsne',                    None, 0, '#2a78d6', 'o', '-', 2.0),
            (f't-SNE + 2 neighbours',             tr, 'tsne+refine' + sfx,       None, 2, '#eb6834', 's', '-', 1.6),
            (f't-SNE + 4 neighbours',             tr, 'tsne+refine' + sfx,       None, 4, '#1baf7a', '^', '-', 1.6),
            (f't-SNE + 6 neighbours',             tr, 'tsne+refine' + sfx,       None, 6, '#eda100', 'D', '-', 1.6),
            (f't-SNE (K=3) + hierarchical',       tr, 'tsne+hierarchical' + sfx, 3,    0, '#e87ba4', 'v', '-', 1.6)]


SEARCH_CL = [('greedy search, clusters only',       '-', 'greedy[clusters]',       None, 0, '#008300', 'P', '--', 1.6),
             ('hierarchical search, clusters only', '-', 'hierarchical[clusters]', None, 0, '#6250d6', 'X', '--', 1.6)]
SEARCH_FULL = [('greedy search, full codebook',       '-', 'greedy',               None, 0, '#008300', 'P', '--', 1.6),
               ('hierarchical search, full codebook', '-', 'hierarchical',         None, 0, '#e34948', 'x', '--', 1.6)]
# (folder, title, main training, series)
VARIANTS = [
    ('0_all_combined', 'informed t-SNE (solid) vs searches (dashed), grey = pure t-SNE with other trainings', 'informed',
     tsne_series('informed', '[clusters]', 'informed') + SEARCH_CL + [SEARCH_FULL[1]] +
     [('t-SNE, random training',         'random', 'tsne', None, 0, '#9a998f', 's', '-', 1.0),
      ('t-SNE, coarse-to-fine training', 'nested', 'tsne', None, 0, '#c3c2b7', '^', '-', 1.0)]),
    ('1_informed_prior', 'angle prior known: informed t-SNE training, hybrids and searches within the user clusters',
     'informed', tsne_series('informed', '[clusters]', 'informed') + SEARCH_CL),
    ('2_no_prior_coarse_to_fine', 'no prior: coarse-to-fine t-SNE training, hybrids and searches over the full codebook',
     'nested', tsne_series('nested', '', 'coarse-to-fine') + SEARCH_FULL),
    ('3_random_training', 'no prior: random t-SNE training codewords (mean of 10 draws), full-codebook hybrids and searches',
     'random', tsne_series('random', '', 'random') + SEARCH_FULL),
]
COLS = [('EMS1 alone', 'EMS1', 'EMS1 trained alone'), ('EMS2 alone', 'EMS2', 'EMS2 trained alone'),
        ('joint', 'EMS1', 'joint training, EMS1 view'), ('joint', 'EMS2', 'joint training, EMS2 view'),
        ('joint', 'joint', 'joint training, both EMSs serving\n(angles pooled over EMS1 and EMS2)')]
ROWS = [('beam_err_median_deg', 'median angle error [deg]'), ('beam_err_rmse_deg', 'RMSE angle error [deg]'),
        ('dsnr_median_dB', 'median SNR gain [dB]')]
SURF = {'M60_raytraced': 'real 60x60', 'M60': '60x60 (synthetic)', 'M120': '120x120 (synthetic)',
        'M180': '180x180 (synthetic)', 'M240': '240x240 (synthetic)'}


def curve(tab, tr, sc, Kf, R, met):
    pts = {}
    for (t, s, K, r, w), x in tab.items():
        if t == tr and s == sc and r == R and (Kf is None or K == Kf) and 1 <= w <= max(SW):
            pts[w] = float(x[met])
    return sorted(pts), [pts[w] for w in sorted(pts)]


def exhaustive(tab, met):
    v = [float(x[met]) for (t, s, K, r, w), x in tab.items() if s == 'exhaustive']
    return v[0] if v else None


for folder, vtitle, main, series in VARIANTS:
    out = os.path.join(D, 'figures', folder)
    os.makedirs(out, exist_ok=True)
    for pop, grp, gname in [('fullnlos', 'full NLOS', 'full-NLOS users'), ('nolos', 'reflection only', 'reflection-only (orange) users')]:
        for surf, sname in SURF.items():
            if not any(k[0] == pop and k[1] == surf and k[4] == grp for k in T):
                continue
            fig, ax = plt.subplots(len(ROWS), len(COLS), figsize=(19, 10), sharex=True, sharey='row')
            for c, (dep, view, ctitle) in enumerate(COLS):
                tab = T.get((pop, surf, dep, view, grp), {})
                for r, (met, ylab) in enumerate(ROWS):
                    a = ax[r, c]
                    if r == 0:
                        a.set_title(ctitle, color=INK, fontsize=10)
                    for lab, tr, sc, Kf, R, col, mk, ls, lw in series:
                        xs, ys = curve(tab, tr, sc, Kf, R, met)
                        if xs:
                            a.plot(xs, ys, ls, color=col, marker=mk, ms=4, lw=lw)
                    if not curve(tab, main, 'tsne', None, 0, met)[0]:          # e.g. joint training never run
                        a.text(0.5, 0.92, 't-SNE not computed for this training\n(searches only: they need no training)',
                               ha='center', va='top', transform=a.transAxes, color=MUTED, fontsize=8)
                    if met.startswith('beam'):                                 # the t-SNE AoA estimate itself
                        xs, ys = curve(tab, main, 'tsne', None, 0, met.replace('beam', 'tsne'))
                        if xs:
                            a.plot(xs, ys, ':', color='#2a78d6', lw=1.4)
                    else:
                        ex = exhaustive(tab, met)
                        if ex is not None:
                            a.axhline(ex, color=INK, lw=1, ls=':')
                            a.text(15, ex, 'exhaustive', ha='right', va='bottom', color=MUTED, fontsize=8)
                    if c == 0:
                        a.set_ylabel(ylab)
                    if r == len(ROWS) - 1:
                        a.set_xlabel('total sweeps')
                    a.set_xticks(SW[::2])
            h = [Line2D([], [], color=col, marker=mk, ls=ls, lw=lw, label=lab) for lab, tr, sc, Kf, R, col, mk, ls, lw in series]
            h.append(Line2D([], [], color='#2a78d6', ls=':', lw=1.4, label='t-SNE AoA estimate (angle rows)'))
            fig.legend(handles=h, loc='upper center', ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.0), fontsize=9)
            fig.suptitle(f'{gname}, surface {sname}: {vtitle}', y=1.05, color=INK, fontsize=12)
            fig.tight_layout(rect=(0, 0, 1, 0.93))
            fn = f"{'fullNLOS' if pop == 'fullnlos' else 'reflectionOnly'}_{surf}.png"
            for attempt in range(5):                               # OneDrive may briefly lock the file
                try:
                    fig.savefig(os.path.join(out, fn), dpi=100, bbox_inches='tight')
                    break
                except OSError:
                    time.sleep(2)
            else:
                raise OSError('could not save ' + fn)
            plt.close(fig)
        print(folder)
