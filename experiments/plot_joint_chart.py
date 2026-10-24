"""Figures of the joint 2-D chart experiment (results_aoa/serving_joint_chart/<noise>/summary.csv).

One figure per surface and training list. Rows: median and RMSE of the angle error (pointing error of
the chosen beam; dotted: error of the joint t-SNE AoA estimate itself), median SNR gain. Columns:
EMS1 (served by EMS1 alone), EMS2, pair (both EMSs serving; angle error sqrt(e1^2 + e2^2) over the
users reached by both). Solid: joint t-SNE chart with real K x K training (K^2 sweeps) and its
hybrids; dashed: searches. x axis: total sweeps.
  1_informed_prior          : informed K x K lists; hybrids and searches within the user clusters
  2_no_prior_coarse_to_fine : coarse-to-fine K x K lists; hybrids and searches over the full codebook
Usage: python experiments/plot_joint_chart.py results_aoa/serving_joint_chart/est_G60dB
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
    T[(x['surface'], x['view'])][(x['training'], x['scheme'], int(x['K']), int(x['R']), int(x['sweeps']))] = x
XMAX = max(k[4] for t in T.values() for k in t)


def series(tr, sfx, Kh):
    return [('joint t-SNE chart (LUT beam)',    tr, 'tsne',                    None, 0, '#2a78d6', 'o', '-', 2.0),
            ('+ 2 neighbours',                  tr, 'tsne+refine' + sfx,       None, 2, '#eb6834', 's', '-', 1.6),
            ('+ 4 neighbours',                  tr, 'tsne+refine' + sfx,       None, 4, '#1baf7a', '^', '-', 1.6),
            ('+ 6 neighbours',                  tr, 'tsne+refine' + sfx,       None, 6, '#eda100', 'D', '-', 1.6),
            (f'K={Kh} ({Kh*Kh} acq.) + hierarchical', tr, 'tsne+hierarchical' + sfx, Kh, 0, '#e87ba4', 'v', '-', 1.6)]


VARIANTS = [
    ('1_informed_prior', 'angle prior known: informed K x K training, hybrids and searches within the user clusters', 'informed',
     series('informed', '[clusters]', 2) +
     [('greedy search, clusters only',       '-', 'greedy[clusters]',       None, 0, '#008300', 'P', '--', 1.6),
      ('hierarchical search, clusters only', '-', 'hierarchical[clusters]', None, 0, '#6250d6', 'X', '--', 1.6)]),
    ('2_no_prior_coarse_to_fine', 'no prior: coarse-to-fine K x K training, hybrids and searches over the full codebook', 'c2f',
     series('c2f', '', 2) +
     [('greedy search, full codebook',       '-', 'greedy',                 None, 0, '#008300', 'P', '--', 1.6),
      ('hierarchical search, full codebook', '-', 'hierarchical',           None, 0, '#e34948', 'x', '--', 1.6)]),
]
COLS = [('EMS1', 'EMS1 (served by EMS1)'), ('EMS2', 'EMS2 (served by EMS2)'),
        ('pair', 'both EMSs serving\n(angle: sqrt(e1^2 + e2^2), users reached by both)')]
ROWS = [('beam_err_median_deg', 'median angle error [deg]'), ('beam_err_rmse_deg', 'RMSE angle error [deg]'),
        ('dsnr_median_dB', 'median SNR gain [dB]')]
SURF = {'M60_raytraced': 'real 60x60', 'M60': '60x60 (synthetic)', 'M120': '120x120 (synthetic)',
        'M180': '180x180 (synthetic)', 'M240': '240x240 (synthetic)'}


def curve(tab, tr, sc, Kf, R, met):
    pts = {}
    for (t, s, K, r, w), x in tab.items():
        if t == tr and s == sc and r == R and (Kf is None or K == Kf) and w >= 1:
            pts[w] = float(x[met])
    return sorted(pts), [pts[w] for w in sorted(pts)]


for folder, vtitle, main, ser in VARIANTS:
    out = os.path.join(D, 'figures', folder)
    os.makedirs(out, exist_ok=True)
    for surf, sname in SURF.items():
        if not any(k[0] == surf for k in T):
            continue
        fig, ax = plt.subplots(len(ROWS), len(COLS), figsize=(16, 10), sharex=True, sharey='row')
        for c, (view, ctitle) in enumerate(COLS):
            tab = T.get((surf, view), {})
            for r, (met, ylab) in enumerate(ROWS):
                a = ax[r, c]
                if r == 0:
                    a.set_title(ctitle, color=INK, fontsize=10)
                for lab, tr, sc, Kf, R, col, mk, ls, lw in ser:
                    xs, ys = curve(tab, tr, sc, Kf, R, met)
                    if xs:
                        a.plot(xs, ys, ls, color=col, marker=mk, ms=4, lw=lw)
                if met.startswith('beam'):
                    xs, ys = curve(tab, main, 'tsne', None, 0, met.replace('beam', 'tsne'))
                    if xs:
                        a.plot(xs, ys, ':', color='#2a78d6', lw=1.4)
                else:
                    ex = [float(x[met]) for k, x in tab.items() if k[1] == 'exhaustive']
                    if ex:
                        a.axhline(ex[0], color=INK, lw=1, ls=':')
                        a.text(XMAX, ex[0], 'exhaustive', ha='right', va='bottom', color=MUTED, fontsize=8)
                if c == 0:
                    a.set_ylabel(ylab)
                if r == len(ROWS) - 1:
                    a.set_xlabel('total sweeps (pure joint t-SNE at K^2 = 1, 4, 9, ...)')
                a.set_xticks([1, 4, 9, 16, 25, 36, 49, 64, XMAX])
        h = [Line2D([], [], color=col, marker=mk, ls=ls, lw=lw, label=lab) for lab, tr, sc, Kf, R, col, mk, ls, lw in ser]
        h.append(Line2D([], [], color='#2a78d6', ls=':', lw=1.4, label='joint t-SNE AoA estimate (angle rows)'))
        fig.legend(handles=h, loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0), fontsize=9)
        fig.suptitle(f'full-NLOS users, surface {sname}, joint 2-D t-SNE chart: {vtitle}', y=1.05, color=INK, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        fn = f'fullNLOS_{surf}.png'
        for attempt in range(5):                                   # OneDrive may briefly lock the file
            try:
                fig.savefig(os.path.join(out, fn), dpi=100, bbox_inches='tight')
                break
            except OSError:
                time.sleep(2)
        else:
            raise OSError('could not save ' + fn)
        plt.close(fig)
    print(folder)
