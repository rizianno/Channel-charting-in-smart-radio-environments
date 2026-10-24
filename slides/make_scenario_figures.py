"""Scenario figures for the slides, from the dataset users (slides/assets/users_scenario.csv, exported
from data/HEMS_dataset.mat: position, link class, reachability and azimuth at each EMS).

  scenario_link_classes : every user coloured by its link class (LoS / reflection only / full NLOS)
  coverage_EMS1, coverage_EMS2 : full-NLOS users that EMS reaches, coloured by their azimuth at that EMS
Usage: python slides/make_scenario_figures.py
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(HERE, 'assets')
INK, MUTED, GRID = '#222220', '#6b6a63', '#e4e3dc'
plt.rcParams.update({'font.size': 11, 'font.family': 'Calibri', 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False})
U = list(csv.DictReader(open(os.path.join(A, 'users_scenario.csv'))))
N = {r['node']: (float(r['x_m']), float(r['y_m'])) for r in csv.DictReader(open(os.path.join(A, 'nodes_scenario.csv')))}
x = np.array([float(r['x_m']) for r in U]);  y = np.array([float(r['y_m']) for r in U])
cls = np.array([r['link_class'] for r in U])
reach = {e: np.array([r[f'reach_{e}'] == '1' for r in U]) for e in ('EMS1', 'EMS2')}
az = {e: np.array([float(r[f'az_{e}_deg']) for r in U]) for e in ('EMS1', 'EMS2')}
NORMAL = {'EMS1': -1, 'EMS2': +1}                   # EMS1 faces -x, EMS2 faces +x (surface along y)
CLS = [('LoS', '#9fb7d3', 'line of sight'), ('reflection only', '#eda100', 'reflection only (wall path)'),
       ('full NLOS', '#e34948', 'full NLOS (no BS-UE channel)')]


def nodes(a, highlight=None):
    bx, by = N['BS']
    a.plot(bx, by, marker='^', ms=13, color=INK, zorder=5)
    a.annotate('BS', (bx, by), xytext=(8, 4), textcoords='offset points', fontsize=12, weight='bold', color=INK)
    for e in ('EMS1', 'EMS2'):
        ex, ey = N[e];  lw = 6 if (highlight is None or highlight == e) else 3
        col = INK if (highlight is None or highlight == e) else '#b9b7b0'
        a.plot([ex, ex], [ey - 3, ey + 3], color=col, lw=lw, solid_capstyle='butt', zorder=5)
        a.annotate('', xy=(ex + 7 * NORMAL[e], ey), xytext=(ex, ey), arrowprops=dict(arrowstyle='->', color=col, lw=1.6), zorder=5)
        a.annotate(e, (ex, ey), xytext=(-34 if NORMAL[e] < 0 else 10, 12), textcoords='offset points', fontsize=12, weight='bold', color=col)


def finish(a):
    a.set_aspect('equal');  a.set_xlabel('x [m]');  a.set_ylabel('y [m]')
    a.grid(True, color=GRID, lw=0.6);  a.set_axisbelow(True)


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(A, f'{name}.{ext}'), dpi=200, bbox_inches='tight')
    plt.close(fig)


# ---- every user by link class
fig, a = plt.subplots(figsize=(7.2, 6.4))
for c, col, lab in CLS:
    m = cls == c
    a.scatter(x[m], y[m], s=7, color=col, lw=0, label=f'{lab}: {m.sum()} ({100*m.mean():.0f}%)', zorder=2)
nodes(a);  finish(a)
a.legend(loc='upper left', bbox_to_anchor=(0, -0.11), frameon=False, ncol=1, markerscale=2.5, fontsize=11)
save(fig, 'scenario_link_classes')

# ---- full-NLOS users each EMS reaches, by azimuth at that EMS
DIV = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#86b6ef', '#e9e8e3', '#f2a07a', '#b8402a'])
for e in ('EMS1', 'EMS2'):
    fig, a = plt.subplots(figsize=(7.2, 6.4))
    other = cls != 'full NLOS'
    a.scatter(x[other], y[other], s=4, color='#e9e8e3', lw=0, zorder=1)
    nr = (cls == 'full NLOS') & ~reach[e]
    a.scatter(x[nr], y[nr], s=6, color='#c9c7c0', lw=0, zorder=2)
    m = (cls == 'full NLOS') & reach[e]
    sc = a.scatter(x[m], y[m], s=9, c=az[e][m], cmap=DIV, vmin=-90, vmax=90, lw=0, zorder=3)
    nodes(a, highlight=e);  finish(a)
    cb = fig.colorbar(sc, ax=a, shrink=0.75, pad=0.02);  cb.set_label(f'azimuth at {e} [deg]')
    h = [Line2D([], [], ls='', marker='o', ms=6, color='#7f8ea3', label=f'full NLOS, reached by {e}: {m.sum()}'),
         Line2D([], [], ls='', marker='o', ms=6, color='#c9c7c0', label=f'full NLOS, not reached: {nr.sum()}'),
         Line2D([], [], ls='', marker='o', ms=6, color='#e9e8e3', label='LoS or reflection only')]
    a.legend(handles=h, loc='upper left', bbox_to_anchor=(0, -0.11), frameon=False, fontsize=11)
    save(fig, f'coverage_{e}')
print('done')
