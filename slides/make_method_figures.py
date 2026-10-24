"""Figures of the method slides, in slides/assets/method/ (PNG + PDF).

  cc_geometries    : toy example (generated here): user positions -> raw CSI -> covariance features
                     (LED) -> unsupervised t-SNE chart; the four geometries of Studer et al., IEEE
                     Access 2018, Fig. 3
  sre_features     : same toy: three users, the magnitude of their covariances and their LEDs
  aoa_chart_1d     : real data (results_aoa/serving_deployments/est_G60dB/per_user.csv): EMS1 alone,
                     informed training, K = 3, real 60x60; users coloured by their true azimuth, on the
                     map and on the 1-D chart
  codebook_beams   : phase ramp of one codeword on a 60x60 EMS; beam patterns of the 68 codewords
  beam_choice      : probes of each beam-selection scheme for one example user (single-path gain)
  gain_vs_pointing : normalised beam gain vs pointing error in beamwidths (M = 60 and 240)
Toy channel: 30 GHz, BS with an 8-antenna lambda/2 ULA, line of sight plus one reflection on each of
two street walls (image sources),
spherical wavefronts; covariance over 40 snapshots (user moved within 0.3 m) plus the noise floor.
Usage: python slides/make_method_figures.py
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'assets', 'method');  os.makedirs(OUT, exist_ok=True)
INK, MUTED, GRID, PANEL = '#222220', '#6b6a63', '#e4e3dc', '#f7f3ec'
BLUE, ORANGE, GREEN, AMBER, RED, PURPLE, GOLD = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e34948', '#6250d6', '#c9a227'
plt.rcParams.update({'font.size': 10, 'font.family': 'Calibri', 'mathtext.fontset': 'cm', 'axes.edgecolor': MUTED,
                     'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False})
DIV = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#86b6ef', '#e9e8e3', '#f2a07a', '#b8402a'])


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'), dpi=250, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(name)


# ====================================================================== toy channel and t-SNE
LAM = 0.01;  K0 = 2 * np.pi / LAM;  NBS = 8                    # small array: smooth features, readable toy chart
ANT = np.stack([np.zeros(NBS), (np.arange(NBS) - (NBS - 1) / 2) * LAM / 2], 1)     # ULA along y at the origin
WALLS = [(22.0, 0.6), (-22.0, 0.5)]                       # street walls y = const (reflection coefficient)
rng = np.random.default_rng(7)
gx, gy = np.meshgrid(np.linspace(12, 52, 25), np.linspace(-15, 15, 20))
POS = np.stack([gx.ravel(), gy.ravel()], 1) + rng.uniform(-0.6, 0.6, (gx.size, 2))


def channel(p):
    """Uplink channel at the BS antennas, p: (..., 2) user positions -> (..., NBS)."""
    r = np.linalg.norm(p[..., None, :] - ANT, axis=-1)
    h = LAM / (4 * np.pi * r) * np.exp(-1j * K0 * r)
    for yw, g in WALLS:                                    # wall reflection = path from the image source
        im = np.stack([p[..., 0], 2 * yw - p[..., 1]], -1)
        r = np.linalg.norm(im[..., None, :] - ANT, axis=-1)
        h = h + g * LAM / (4 * np.pi * r) * np.exp(-1j * K0 * r)
    return h


T = 40                                                     # snapshots: user moved within 0.3 m
ang = rng.uniform(0, 2 * np.pi, (len(POS), T));  rad = 0.3 * np.sqrt(rng.uniform(0, 1, (len(POS), T)))
H = channel(POS[:, None, :] + np.stack([rad * np.cos(ang), rad * np.sin(ang)], -1))   # (n, T, NBS)
S2 = np.median(np.sum(np.abs(H[:, 0]) ** 2, 1)) / NBS / 10       # noise: 10 dB per antenna (median user)
COV = np.einsum('utn,utm->unm', H, H.conj()) / T + S2 * np.eye(NBS)   # noise floor (exact, no sampling noise)
iu = np.triu_indices(NBS, 1)


def logm_feat(R):
    w, V = np.linalg.eigh(R)
    L = (V * np.log(w)) @ V.conj().T
    return np.concatenate([np.real(np.diag(L)), np.sqrt(2) * np.real(L[iu]), np.sqrt(2) * np.imag(L[iu])])


FEAT = np.array([logm_feat(R) for R in COV])
LED = np.sqrt(np.maximum(np.sum(FEAT ** 2, 1)[:, None] + np.sum(FEAT ** 2, 1)[None] - 2 * FEAT @ FEAT.T, 0))


def d2p(D, perp, tol=1e-5):
    n = len(D);  P = np.zeros((n, n));  D2 = (D / np.median(D)) ** 2;  logU = np.log(perp)
    for i in range(n):
        di = np.delete(D2[i], i);  beta, lo, hi = 1.0, None, None
        for _ in range(200):
            p = np.exp(-(di - di.min()) * beta);  s = p.sum();  Hh = np.log(s) + beta * np.sum((di - di.min()) * p) / s
            if abs(Hh - logU) < tol:
                break
            if Hh > logU:
                lo = beta;  beta = beta * 2 if hi is None else (beta + hi) / 2
            else:
                hi = beta;  beta = beta / 2 if lo is None else (beta + lo) / 2
        P[i, np.arange(n) != i] = p / s
    return P


def tsne(D, dims=2, perp=30, iters=1000, seed=0):
    """Exact t-SNE (van der Maaten's reference scheme: early exaggeration, momentum, gains)."""
    n = len(D);  P = d2p(D, perp);  P = P + P.T;  P = np.maximum(P / P.sum(), 1e-12) * 4
    Y = 1e-4 * np.random.default_rng(seed).standard_normal((n, dims))
    inc = np.zeros_like(Y);  gains = np.ones_like(Y);  mom = 0.5
    for it in range(iters):
        sq = np.sum(Y ** 2, 1);  num = 1 / (1 + sq[:, None] + sq[None] - 2 * Y @ Y.T);  np.fill_diagonal(num, 0)
        Q = np.maximum(num / num.sum(), 1e-12);  L = (P - Q) * num
        G = 4 * ((np.diag(L.sum(0)) - L) @ Y)
        same = np.sign(G) == np.sign(inc)
        gains = np.maximum(np.where(same, gains * 0.8, gains + 0.2), 0.01)
        inc = mom * inc - 200 * gains * G;  Y = Y + inc;  Y -= Y.mean(0)
        if it == 100: P /= 4
        if it == 250: mom = 0.8
    return Y


# bivariate colour of each position (corners: blue, green, amber, magenta)
xn = (POS[:, 0] - POS[:, 0].min()) / np.ptp(POS[:, 0]);  yn = (POS[:, 1] - POS[:, 1].min()) / np.ptp(POS[:, 1])
C00, C10, C01, C11 = (np.array(matplotlib.colors.to_rgb(c)) for c in ('#1c5cab', '#eda100', '#1baf7a', '#c03a8c'))
COL = ((1 - xn) * (1 - yn))[:, None] * C00 + (xn * (1 - yn))[:, None] * C10 + ((1 - xn) * yn)[:, None] * C01 + (xn * yn)[:, None] * C11


def panel(ax, title):
    ax.set_facecolor(PANEL)
    for s in ax.spines.values():
        s.set_visible(True);  s.set_color('#cfc8bb')
    ax.set_xticks([]);  ax.set_yticks([])
    ax.set_title(title, color=INK, fontsize=11.5, pad=5)


def scene(ax, s=11):
    ax.scatter(POS[:, 0], POS[:, 1], c=COL, s=s, lw=0, zorder=3)
    ax.plot(0, 0, marker=(3, 0, -90), ms=12, color=INK, zorder=4);  ax.annotate('BS', (0, 0), xytext=(-4, 9), textcoords='offset points', fontsize=10, weight='bold')
    for yw, _ in WALLS:
        ax.add_patch(Rectangle((-4, yw if yw > 0 else yw - 3), 60, 3, color='#cfc8bb', lw=0, zorder=2))
    ax.text(26, 23.5, 'building wall', ha='center', va='center', fontsize=8.5, color=MUTED, zorder=3)
    ax.set_xlim(-6, 57);  ax.set_ylim(-26, 26);  ax.set_aspect('equal')


def fig_cc_geometries():
    Z = tsne(LED, perp=60)
    h0 = H[:, 0] / np.linalg.norm(H[:, 0], axis=1, keepdims=True)
    X = np.concatenate([np.real(h0), np.imag(h0)], 1);  X = (X - X.mean(0)) / X.std(0)
    Ur = np.linalg.svd(X - X.mean(0), full_matrices=False)[0][:, :2]
    fig = plt.figure(figsize=(7.6, 6.4))
    a1 = fig.add_axes([0.02, 0.56, 0.40, 0.38]);  a2 = fig.add_axes([0.58, 0.56, 0.40, 0.38])
    a3 = fig.add_axes([0.58, 0.04, 0.40, 0.38]);  a4 = fig.add_axes([0.02, 0.04, 0.40, 0.38])
    panel(a1, 'spatial geometry: user positions');  scene(a1)
    panel(a2, 'radio geometry: raw channels (scrambled)')
    a2.scatter(Ur[:, 0], Ur[:, 1], c=COL, s=9, lw=0)
    panel(a3, 'feature geometry: covariances, LED')
    i, j = np.triu_indices(len(POS), 1);  pick = rng.choice(i.size, 9000, replace=False)
    dist = np.linalg.norm(POS[i[pick]] - POS[j[pick]], axis=1);  led = LED[i[pick], j[pick]]
    a3.scatter(dist, led, s=2, color='#9aa7b8', lw=0, alpha=0.5)
    bins = np.linspace(0, 50, 26);  mid = (bins[1:] + bins[:-1]) / 2
    med = [np.median(led[(dist >= b0) & (dist < b1)]) for b0, b1 in zip(bins[:-1], bins[1:])]
    a3.plot(mid, med, color=INK, lw=1.6)
    a3.set_xticks([0, 25, 50]);  a3.set_yticks([]);  a3.tick_params(labelsize=8.5)
    a3.set_xlabel('distance between two users [m]', fontsize=9, labelpad=1);  a3.set_ylabel('LED between them', fontsize=9)
    a3.text(0.97, 0.05, 'close users: small LED', transform=a3.transAxes, ha='right', fontsize=8.5, color=MUTED)
    panel(a4, 'channel chart: t-SNE, no labels')
    a4.scatter(Z[:, 0], Z[:, 1], c=COL, s=9, lw=0)
    a4.margins(0.08)

    def arrow(p0, p1, lab, xy, color=INK, both=False, ha='center'):
        fig.patches.append(FancyArrowPatch(p0, p1, transform=fig.transFigure, arrowstyle='<|-|>' if both else '-|>',
                                           mutation_scale=16, lw=2.2 if both else 1.6, color=color))
        fig.text(*xy, lab, ha=ha, va='center', fontsize=10, color=color, linespacing=1.15)
    arrow((0.44, 0.75), (0.56, 0.75), 'wireless\nchannel', (0.50, 0.82))
    arrow((0.70, 0.54), (0.70, 0.46), 'feature: covariance\ndissimilarity: LED', (0.72, 0.50), ha='left')
    arrow((0.56, 0.23), (0.44, 0.23), 'charting\nfunction\n(t-SNE)', (0.50, 0.32))
    arrow((0.14, 0.54), (0.14, 0.46), 'local geometry\npreserved', (0.16, 0.50), color=RED, both=True, ha='left')
    save(fig, 'cc_geometries')


def fig_sre_features():
    iA = int(np.argmin(np.linalg.norm(POS - [24, -6], axis=1)))
    dAB = np.linalg.norm(POS - POS[iA], axis=1);  dAB[iA] = 99;  iB = int(np.argmin(np.abs(dAB - 1.5)))
    iC = int(np.argmin(np.linalg.norm(POS - [46, 11], axis=1)))
    fig = plt.figure(figsize=(5.4, 5.0))
    a0 = fig.add_axes([0.0, 0.50, 1.0, 0.48]);  scene(a0, s=6)
    a0.set_xticks([]);  a0.set_yticks([]);  a0.set_facecolor(PANEL)
    for s in a0.spines.values(): s.set_visible(True);  s.set_color('#cfc8bb')
    KC = {'A': BLUE, 'B': GREEN, 'C': ORANGE}
    for k, i in zip('ABC', (iA, iB, iC)):
        a0.plot(*POS[i], 'o', ms=10, mfc='none', mec=KC[k], mew=2.2, zorder=5)
        a0.annotate(k, POS[i], xytext=(7, 4) if k != 'B' else (7, -12), textcoords='offset points', fontsize=12, weight='bold', color=KC[k], zorder=5)
    # what the covariance holds: power arriving from each direction, a(psi)^H R a(psi) / N_BS
    psi = np.linspace(-90, 90, 721);  Aps = np.exp(1j * np.pi * np.outer(np.sin(np.radians(psi)), np.arange(NBS)))
    a = fig.add_axes([0.13, 0.08, 0.84, 0.32])
    for k, i in zip('ABC', (iA, iB, iC)):
        Ps = np.real(np.einsum('pn,nm,pm->p', Aps.conj(), COV[i], Aps)) / NBS
        a.plot(psi, 10 * np.log10(Ps / S2), color=KC[k], lw=2.0 if k != 'B' else 1.6, ls='-' if k != 'B' else '--', label=f'user {k}')
    a.set_xlim(-90, 90);  a.set_xticks([-90, -45, 0, 45, 90]);  a.tick_params(labelsize=8.5)
    a.set_xlabel('direction at the BS [deg]', fontsize=9.5, labelpad=1);  a.set_ylabel('power / noise [dB]', fontsize=9.5)
    a.grid(True, color=GRID, lw=0.7);  a.legend(frameon=False, fontsize=9, loc='upper right', ncol=3)
    a.set_title(f'covariance seen as power per direction:  LED $d_{{AB}}$ = {LED[iA, iB]:.1f} (1.5 m apart),  '
                f'$d_{{AC}}$ = {LED[iA, iC]:.1f}', fontsize=10, color=INK, pad=4)
    save(fig, 'sre_features')


# ====================================================================== real 1-D AoA chart
def fig_aoa_chart_1d():
    f = os.path.join(REPO, 'results_aoa', 'serving_deployments', 'est_G60dB', 'per_user.csv')
    rows = [r for r in csv.DictReader(open(f)) if r['surface'] == 'M60_raytraced' and r['deployment'] == 'EMS1 alone'
            and r['view'] == 'EMS1' and r['training'] == 'informed' and r['K'] == '3']
    g = lambda c: np.array([float(r[c]) for r in rows])
    x, y, th, est, err = g('x_m'), g('y_m'), g('az_true_deg'), g('az_est_deg'), g('aoa_err_deg')
    anc = g('anchor') == 1;  sc = g('scored') == 1
    nodes = {r['node']: (float(r['x_m']), float(r['y_m'])) for r in csv.DictReader(open(os.path.join(HERE, 'assets', 'nodes_scenario.csv')))}
    fig = plt.figure(figsize=(5.6, 4.9))
    a0 = fig.add_axes([0.0, 0.42, 0.82, 0.56])
    a0.scatter(x[~anc], y[~anc], c=th[~anc], cmap=DIV, vmin=-90, vmax=90, s=10, lw=0, zorder=2)
    a0.scatter(x[anc], y[anc], c=th[anc], cmap=DIV, vmin=-90, vmax=90, s=22, lw=0.9, edgecolors=INK, zorder=3)
    ex, ey = nodes['EMS1'];  a0.plot([ex, ex], [ey - 3, ey + 3], color=INK, lw=5, solid_capstyle='butt')
    a0.annotate('', xy=(ex - 7, ey), xytext=(ex, ey), arrowprops=dict(arrowstyle='->', color=INK, lw=1.5))
    a0.annotate('EMS1', (ex, ey), xytext=(4, 10), textcoords='offset points', fontsize=11, weight='bold')
    bx, by = nodes['BS'];  a0.plot(bx, by, '^', ms=10, color=INK);  a0.annotate('BS', (bx, by), xytext=(7, 0), textcoords='offset points', fontsize=10, weight='bold')
    a0.set_aspect('equal');  a0.set_xlabel('x [m]', fontsize=9, labelpad=1);  a0.set_ylabel('y [m]', fontsize=9)
    a0.tick_params(labelsize=8.5)
    a0.set_title('users reached by EMS1 (full NLOS), colour = true azimuth $\\theta_u$', fontsize=10, color=INK)
    cax = fig.add_axes([0.85, 0.50, 0.025, 0.40])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(cmap=DIV, norm=matplotlib.colors.Normalize(-90, 90)), cax=cax)
    cb.set_label('azimuth at EMS1 [deg]', fontsize=9);  cb.ax.tick_params(labelsize=8.5)
    a1 = fig.add_axes([0.05, 0.07, 0.77, 0.22])
    jit = np.random.default_rng(3).uniform(-0.35, 0.35, x.size)
    a1.scatter(est[sc], jit[sc], c=th[sc], cmap=DIV, vmin=-90, vmax=90, s=11, lw=0, zorder=3)
    a1.scatter(est[anc], np.ones(anc.sum()) * 0.75, marker='|', s=60, color=INK, lw=1.2, zorder=3)
    xs = np.linspace(-90, 90, 361)
    a1.imshow(xs[None], extent=(-90, 90, -0.75, -0.55), cmap=DIV, vmin=-90, vmax=90, aspect='auto', zorder=2)
    a1.set_xlim(-75, 75);  a1.set_ylim(-0.8, 0.95);  a1.set_yticks([])
    a1.spines['left'].set_visible(False);  a1.tick_params(labelsize=8.5)
    a1.set_xlabel('chart coordinate = estimated azimuth $\\hat\\theta_u$ [deg]', fontsize=9, labelpad=1)
    a1.text(-26, 0.75, '| anchors, clamped at their azimuth', fontsize=8, color=INK, va='center')
    bad = sc & (err > 30)
    a1.annotate(f'{bad.sum()} users of the northern group\nplaced among the southern one\n(large errors: RMSE >> median)',
                (est[bad].mean(), 0.1), xytext=(-12, -0.05), fontsize=8.5, color=INK, va='center',
                arrowprops=dict(arrowstyle='->', color=MUTED, lw=1))
    a1.set_title(f'1-D chart, K = 3 slots: median error {np.median(err[sc]):.1f}°, '
                 f'RMSE {np.sqrt(np.mean(err[sc] ** 2)):.1f}°  ({sc.sum()} scored users)', fontsize=10, color=INK)
    save(fig, 'aoa_chart_1d')


# ====================================================================== codebook and beams
def dirichlet(x, M):
    """Normalised array gain |sum_m exp(j pi m x)|^2 / M^2 of an M-element lambda/2 ramp, x = sin offset."""
    x = np.mod(x + 1, 2) - 1;  num = np.sin(M * np.pi * x / 2);  den = M * np.sin(np.pi * x / 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        g = np.where(np.abs(den) < 1e-12, 1.0, (num / den) ** 2)
    return g


M60 = 60;  BW = 0.886 * 2 / M60;  NC = int(np.ceil(2 / BW));  NC += NC % 2
UN = np.mod((np.arange(1, NC + 1) - NC / 2) * 2 / NC + 1, 2) - 1           # u_n, theta_spec = 0


def fig_codebook_beams():
    fig = plt.figure(figsize=(9.6, 3.6))
    a0 = fig.add_axes([0.0, 0.12, 0.28, 0.76])
    n = 40;  dphi = 2 * np.pi * (n - NC / 2) / NC;  m = np.arange(1, M60 + 1)
    ph = np.mod(-np.outer(np.ones(M60), m) * dphi, 2 * np.pi)
    im = a0.imshow(ph, cmap='twilight', vmin=0, vmax=2 * np.pi, origin='lower', extent=(0.5, 60.5, 0.5, 60.5))
    a0.set_xlabel('element column $m$ (along x)', fontsize=9.5);  a0.set_ylabel('element row', fontsize=9.5)
    a0.set_xticks([1, 20, 40, 60]);  a0.set_yticks([1, 20, 40, 60]);  a0.tick_params(labelsize=8.5)
    a0.set_title(f'phase of codeword $n$ = {n} ($\\Delta\\phi_n$ = {np.degrees(dphi):.0f}° per column)', fontsize=10, color=INK)
    cb = fig.colorbar(im, ax=a0, fraction=0.046, pad=0.03, ticks=[0, np.pi, 2 * np.pi])
    cb.ax.set_yticklabels(['0', '$\\pi$', '$2\\pi$']);  cb.ax.tick_params(labelsize=8.5)
    a1 = fig.add_axes([0.40, 0.14, 0.58, 0.68])
    s = np.linspace(0.06, 0.34, 1400)
    for k in range(NC):
        if 0.0 < UN[k] < 0.40:
            a1.plot(s, 10 * np.log10(np.maximum(dirichlet(s - UN[k], M60), 1e-4)), color='#c3c2b7', lw=0.9)
    for k, c in ((38, BLUE), (39, GREEN), (40, ORANGE)):
        a1.plot(s, 10 * np.log10(np.maximum(dirichlet(s - UN[k], M60), 1e-4)), color=c, lw=1.8)
        a1.text(UN[k], 1.2, f'$n$={k + 1}', ha='center', fontsize=9, color=c)
    a1.axhline(-3, color=MUTED, ls=':', lw=1);  a1.text(0.335, -2.4, '$-3$ dB', ha='right', fontsize=9, color=MUTED)
    a1.annotate('', xy=(UN[38] - BW / 2, -3), xytext=(UN[38] + BW / 2, -3), arrowprops=dict(arrowstyle='<->', color=INK, lw=1))
    a1.text(UN[38], -5.0, '$\\Delta u$', ha='center', fontsize=10)
    est = UN[39] + 0.3 * BW
    a1.axvline(est, color=RED, ls='--', lw=1.3)
    box = dict(boxstyle='round,pad=0.15', fc='white', ec='none')
    a1.text(est + 0.003, -27.5, '$\\sin\\hat\\theta_u$ (estimate)', color=RED, fontsize=10, bbox=box, zorder=6)
    a1.plot(UN[39], 0, '*', ms=14, color=RED, zorder=5)
    a1.text(est + 0.003, -8.5, 'look-up table:\nnearest beam ($\\star$)', color=RED, fontsize=9, va='top', bbox=box, zorder=6)
    a1.set_xlim(0.06, 0.34);  a1.set_ylim(-30, 3.5)
    a1.set_xlabel('$\\sin\\theta$ of the incoming wave (zoom; the 68 beams cover $[-1, 1)$)', fontsize=9.5)
    a1.set_ylabel('beam gain [dB]', fontsize=9.5);  a1.tick_params(labelsize=8.5)
    top = a1.secondary_xaxis('top', functions=(lambda v: np.degrees(np.arcsin(np.clip(v, -1, 1))), lambda d: np.sin(np.radians(d))))
    top.set_xticks([5, 10, 15, 20]);  top.tick_params(labelsize=8.5);  top.set_xlabel('[deg]', fontsize=9, labelpad=2)
    a1.set_title('beams of a 60x60 EMS: 68 codewords, one beamwidth apart', fontsize=10, color=INK, pad=6, loc='left')
    save(fig, 'codebook_beams')


# ====================================================================== beam choice schemes
def nest_order(N):
    f = [1 / 2, 1]
    for lev in range(2, 13):
        f += list(np.arange(1, 2 ** lev, 2) / 2 ** lev)
        if len(set(np.round(np.array(f) * N).astype(int))) >= N:
            break
    o = []
    for v in np.clip(np.round(np.array(f) * N).astype(int), 1, N):
        if v not in o: o.append(int(v))
    return o


def hier(g, N, s):
    """Hierarchical search, as hier1 in Main_AoA_Serving_Deployments.m; returns probes and choice."""
    o = nest_order(N);  pr = o[:min(s, 4)];  b = max(pr, key=lambda c: g[c]);  used = len(pr);  delta = N / 4
    while used < s and delta > 1:
        delta = int(np.ceil(delta / 2))
        for sgn in (-1, 1):
            if used >= s: break
            c = int(np.mod(b - 1 + sgn * delta, N) + 1);  pr.append(c);  used += 1
            if g[c] > g[b]: b = c
    return pr, b


def fig_beam_choice():
    su = 0.43;  g = {n: dirichlet(su - UN[n - 1], M60) for n in range(1, NC + 1)}
    cstar = max(g, key=g.get)
    # illustrative anchor groups (sin theta) for the cluster order and the informed training slots
    anchors = np.concatenate([np.linspace(0.30, 0.46, 30), np.linspace(-0.55, -0.40, 20), np.linspace(0.70, 0.78, 8)])
    near = np.abs(anchors[:, None] - UN[None])
    sup = [n + 1 for n in range(NC) if np.any(near[:, n] <= 2 * BW)]
    order = sorted(sup, key=lambda n: -np.sum(near[:, n - 1] <= BW))
    train = [int(np.argmin(np.abs(UN - c))) + 1 for c in (0.38, -0.475, 0.74)]
    prop = [46, 48, 48]                                    # LUT proposal after slot 1, 2, 3 (illustrative)
    S = 8
    rowsd = []
    rowsd.append(('exhaustive', [(n, 'probe', None) for n in range(1, NC + 1)], cstar, f'{NC} sweeps'))
    gr = order[:S];  rowsd.append(('greedy search\n(informed)', [(n, 'probe', i + 1) for i, n in enumerate(gr)], max(gr, key=g.get), f'{S} sweeps'))
    hp, hb = hier(g, NC, S);  rowsd.append(('hierarchical\nsearch', [(n, 'probe', i + 1) for i, n in enumerate(hp)], hb, f'{S} sweeps'))
    rowsd.append(('t-SNE + look-up\ntable (K = 3)', [(n, 'train', i + 1) for i, n in enumerate(train)] + [(prop[-1], 'lut', None)], prop[-1], '3 sweeps'))
    ver = [(train[0], 'train', 1), (prop[0], 'probe', 2), (train[1], 'train', 3), (prop[1], 'probe', 4), (train[2], 'train', 5), (prop[2], 'probe', 6)]
    bv = max([prop[0], prop[1], prop[2]], key=g.get)
    rowsd.append(('t-SNE + verify', ver, bv, '6 sweeps'))
    nb = [];  b = bv;  off = [1, -1, 2, -2]
    for t, o in enumerate(off):
        c = bv + o;  nb.append((c, 'probe', 7 + t))
        if g[c] > g[b]: b = c
    rowsd.append(('+ neighbour\nsearch', [(bv, 'kept', None)] + nb, b, '10 sweeps'))
    fig, ax = plt.subplots(figsize=(9.4, 4.3))
    ax.axvline(cstar, color=GOLD, lw=6, alpha=0.35, zorder=1)
    ax.text(cstar, len(rowsd) - 0.35, 'best beam $c_u^\\star$', ha='center', fontsize=9.5, color='#8a6d00')
    for r, (lab, probes, chosen, cost) in enumerate(rowsd[::-1]):
        ax.axhline(r, color=GRID, lw=0.8, zorder=0)
        seen = {}
        for n, kind, num in probes:
            y = r + 0.24 * seen.get(n, 0);  seen[n] = seen.get(n, 0) + 1     # a codeword probed again: drawn above
            if kind == 'probe':
                ax.plot(n, y, 'o', ms=10 if num else 4, mfc='white' if num else '#9aa7b8', mec=BLUE if num else '#9aa7b8', mew=1.3, zorder=3)
                if num: ax.text(n, y, str(num), ha='center', va='center', fontsize=7.5, color=BLUE, zorder=4)
            elif kind == 'train':
                ax.plot(n, y, 's', ms=10, mfc='white', mec=GREEN, mew=1.3, zorder=3)
                ax.text(n, y, str(num), ha='center', va='center', fontsize=7.5, color=GREEN, zorder=4)
            elif kind == 'lut':
                ax.plot(n, r + 0.28, 'v', ms=7, color=RED, zorder=4)
            elif kind == 'kept':
                ax.plot(n, r, 'o', ms=10, mfc=BLUE, mec=BLUE, zorder=3)
        ax.plot(chosen, r - 0.3, '*', ms=12, color=RED, zorder=5)
        ax.text(-1, r, lab, ha='right', va='center', fontsize=10, color=INK, linespacing=1.0)
        ax.text(NC + 2, r, cost, ha='left', va='center', fontsize=9.5, color=INK)
    ax.set_xlim(0, NC + 1);  ax.set_ylim(-0.6, len(rowsd) - 0.1);  ax.set_yticks([])
    ax.spines['left'].set_visible(False);  ax.set_xticks([1, 17, 34, 51, 68]);  ax.tick_params(labelsize=9)
    ax.set_xlabel('codeword index $n$  (sin of the steered direction grows with $n$)', fontsize=10)
    h = [Line2D([], [], ls='', marker='s', ms=8, mfc='white', mec=GREEN, label='training slot (shared by all users)'),
         Line2D([], [], ls='', marker='o', ms=8, mfc='white', mec=BLUE, label='measured codeword (number = sweep)'),
         Line2D([], [], ls='', marker='v', ms=7, color=RED, label='look-up proposal, not measured'),
         Line2D([], [], ls='', marker='*', ms=11, color=RED, label='chosen codeword')]
    ax.legend(handles=h, loc='upper center', bbox_to_anchor=(0.45, 1.13), ncol=4, frameon=False, fontsize=9, handletextpad=0.3, columnspacing=1.2)
    save(fig, 'beam_choice')


def fig_gain_vs_pointing():
    fig, ax = plt.subplots(figsize=(5.2, 3.3))
    e = np.linspace(0, 3, 1500)
    for M, c, ls in ((60, BLUE, '-'), (240, ORANGE, '--')):
        bw = 0.886 * 2 / M
        ax.plot(e, 10 * np.log10(np.maximum(dirichlet(e * bw, M), 1e-5)), color=c, ls=ls, lw=1.8, label=f'{M}x{M} EMS')
    ax.axhline(-3, color=MUTED, ls=':', lw=1);  ax.axvline(0.5, color=MUTED, ls=':', lw=1)
    ax.annotate('$-3$ dB at half a beamwidth', (0.5, -3), xytext=(0.75, -7), fontsize=9.5, color=INK,
                arrowprops=dict(arrowstyle='->', color=MUTED, lw=1))
    ax.annotate('first null: 1.13 beamwidths', (1.13, -38), xytext=(1.35, -33), fontsize=9.5, color=INK,
                arrowprops=dict(arrowstyle='->', color=MUTED, lw=1))
    ax.set_xlim(0, 3);  ax.set_ylim(-40, 2)
    ax.set_xlabel('pointing error $\\varepsilon_u$ [beamwidths]', fontsize=10);  ax.set_ylabel('beam gain [dB]', fontsize=10)
    ax.grid(True, color=GRID, lw=0.7);  ax.legend(frameon=False, fontsize=9.5, loc='upper right')
    ax.set_title('same loss for every surface size, in beamwidths', fontsize=10, color=INK)
    save(fig, 'gain_vs_pointing')


def fig_training_slots(view='EMS2', surf='M60_raytraced', K=6):
    """Real data: angle prior of the anchors and the first K training slots of each method."""
    D = os.path.join(REPO, 'results_aoa', 'serving_deployments', 'est_G60dB')
    tl = [r for r in csv.DictReader(open(os.path.join(D, 'training_lists.csv'))) if r['surface'] == surf and r['view'] == view]
    ua = np.sin(np.radians([float(r['az_true_deg']) for r in csv.DictReader(open(os.path.join(D, 'per_user_aoa.csv')))
                            if r['surface'] == surf and r['view'] == view and r['training'] == 'informed'
                            and r['K'] == '1' and r['anchor'] == '1']))
    size = surf.split('_')[0][1:]
    fig = plt.figure(figsize=(6.4, 4.2))
    a0 = fig.add_axes([0.13, 0.42, 0.84, 0.44]);  a1 = fig.add_axes([0.13, 0.12, 0.84, 0.26], sharex=a0)
    a0.hist(ua, bins=np.arange(-1, 1.0001, 0.04), color='#b7cdea', edgecolor='white', lw=0.5)
    g = np.linspace(-1, 1, 801)
    kde = np.exp(-0.5 * ((g[:, None] - ua[None]) / 0.03) ** 2).sum(1) / (ua.size * 0.03 * np.sqrt(2 * np.pi))
    a0.plot(g, kde * ua.size * 0.04, color=BLUE, lw=1.6)
    a0.set_ylabel('anchors', fontsize=9.5);  a0.tick_params(labelsize=8.5, labelbottom=False)
    a0.set_title(f'{view}, {size}x{size} EMS: prior PDF of the user angles (from {ua.size} anchors) '
                 f'and the first {K} slots', fontsize=10, color=INK, pad=24)
    top = a0.secondary_xaxis('top', functions=(lambda v: np.degrees(np.arcsin(np.clip(v, -1, 1))), lambda d: np.sin(np.radians(d))))
    top.set_xticks([-60, -30, 0, 30, 60]);  top.tick_params(labelsize=8.5);  top.set_xlabel('[deg]', fontsize=8.5, labelpad=1)
    for row, (tr, lab, col, mk) in enumerate((('nested', 'uninformed', MUTED, 'o'), ('informed', 'informed', GREEN, 's'))):
        for r in tl:
            if r['training'] == tr and int(r['k']) <= K:
                u = float(r['u'])
                a1.plot(u, row, mk, ms=13, mfc='white', mec=col, mew=1.6, zorder=3)
                a1.text(u, row, r['k'], ha='center', va='center', fontsize=8.5, color=col, zorder=4)
                if tr == 'informed':
                    a0.axvline(u, color=GREEN, lw=0.8, alpha=0.6)
        a1.text(-1.04, row, lab, ha='right', va='center', fontsize=10, color=col)
    a1.set_ylim(-0.6, 1.6);  a1.set_yticks([]);  a1.spines['left'].set_visible(False)
    a1.set_xlim(-1.06, 1.06);  a1.tick_params(labelsize=8.5)
    a1.set_xlabel('$\\sin\\theta$ at the EMS (number = training slot $k$)', fontsize=9.5)
    save(fig, 'training_slots')


if __name__ == '__main__':
    fig_cc_geometries()
    fig_sre_features()
    fig_aoa_chart_1d()
    fig_codebook_beams()
    fig_beam_choice()
    fig_gain_vs_pointing()
    fig_training_slots()
