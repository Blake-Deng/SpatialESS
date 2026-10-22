"""Shared physical-size theme; Matplotlib sizes are points, positions millimetres."""
from pathlib import Path
import json
import os
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LAYOUT = json.loads((ROOT / 'config/layout.json').read_text())
COLORS = {'ESS': '#2166AC', 'reference': '#737B86', 'positive': '#C43B53',
          'negative': '#2166AC', 'nonsignificant': '#BAC3CD', 'ink': '#111111',
          'secondary': '#515C67', 'green': '#16877A', 'purple': '#9C6C9F',
          'pale': '#E8EDF1'}
font_dir = Path(os.environ.get('SPATIALESS_FONT_DIR',
                             str(Path.home() / '.local/share/fonts/arial')))
for path in font_dir.glob('*'):
    if path.suffix.lower() in {'.ttf', '.otf'}:
        font_manager.fontManager.addfont(str(path))
try:
    FONT_PATH = font_manager.findfont('Arial', fallback_to_default=False)
    FONT = 'Arial'
except ValueError:
    FONT = 'Liberation Sans'
    FONT_PATH = font_manager.findfont(FONT)
mpl.rcParams.update({
    'font.family': FONT, 'font.size': 6.3, 'text.color': COLORS['ink'],
    'axes.labelsize': 7, 'axes.titlesize': 7, 'axes.titleweight': 'normal',
    'axes.labelcolor': COLORS['ink'], 'axes.linewidth': .5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.color': COLORS['ink'], 'ytick.color': COLORS['ink'],
    'xtick.major.width': .5, 'ytick.major.width': .5,
    'xtick.major.size': 2.3, 'ytick.major.size': 2.3,
    'xtick.major.pad': 2, 'ytick.major.pad': 2,
    'axes.labelpad': 2.2, 'legend.fontsize': 6, 'legend.title_fontsize': 6,
    'legend.frameon': False, 'legend.handlelength': 1.3, 'legend.borderpad': .2,
    'legend.labelspacing': .25, 'legend.columnspacing': .8,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'svg.hashsalt': 'SpatialESS_Nature_20261009_v2',
    'savefig.facecolor': 'white', 'figure.facecolor': 'white',
    'mathtext.fontset': 'dejavusans', 'lines.linewidth': 1.05,
})
DIVERGING = LinearSegmentedColormap.from_list('SpatialESS_effect',
    [COLORS['negative'], '#FFFFFF', COLORS['positive']])
DENSITY = LinearSegmentedColormap.from_list('SpatialESS_records',
    ['#C4DBEE', COLORS['ESS'], '#123D63'])

def effect_norm(limit):
    return TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)

def figure(name):
    w, h = LAYOUT[name]['size_mm']
    fig = plt.figure(figsize=(w / 25.4, h / 25.4), dpi=160)
    fig._name = name
    return fig

def axes_mm(fig, x, y, w, h, **kwargs):
    fw, fh = fig.get_size_inches() * 25.4
    return fig.add_axes([x/fw, 1-(y+h)/fh, w/fw, h/fh], **kwargs)

def text_mm(fig, x, y, text, **kwargs):
    w, h = fig.get_size_inches() * 25.4
    return fig.text(x/w, 1-y/h, text, va='top', **kwargs)

def heading(fig, letter, title):
    x, y, w, h = LAYOUT[fig._name]['panels'][letter]
    text_mm(fig, x, y, letter, size=8, weight='bold')
    text_mm(fig, x+4.1, y+.15, title, size=7, weight='bold')
    return x, y, w, h

def data_axis(fig, letter, title, left=10, right=2, top=8, bottom=9):
    x, y, w, h = heading(fig, letter, title)
    return axes_mm(fig, x+left, y+top, w-left-right, h-top-bottom)

def cell_label(n):
    return f'{n/1000:g}k' if n >= 1000 else str(int(n))

def log_cells(ax, limits=(200, 550000)):
    ax.set_xscale('log')
    ax.set_xlim(*limits)
    ax.set_xticks([1e3, 1e4, 1e5, 7e5], ['1k', '10k', '100k', '700k'])
    ax.minorticks_off()
    ax.set_xlabel('Cells, N')

def save(fig, name=None):
    name = name or fig._name
    for suffix in ['pdf', 'svg']:
        fig.savefig(ROOT / 'figures' / f'{name}.{suffix}', facecolor='white', dpi=600)
    fig.savefig(ROOT / 'figures' / f'{name}.png', dpi=600, facecolor='white')
    fig.savefig(ROOT / 'preview' / f'{name}.png', dpi=180, facecolor='white')
    # Native panel exports retain text and primitives; the protected A slot is skipped.
    if fig._name in LAYOUT:
        from matplotlib.transforms import Bbox
        height = fig.get_size_inches()[1] * 25.4
        for letter, (x, y, w, h) in LAYOUT[fig._name]['panels'].items():
            if fig._name == 'Figure1' and letter == 'A':
                continue
            box = Bbox.from_extents(x/25.4, (height-y-h)/25.4,
                                   (x+w)/25.4, (height-y)/25.4)
            for suffix in ['pdf', 'svg']:
                fig.savefig(ROOT/'panels'/f'{fig._name}{letter}.{suffix}',
                            bbox_inches=box, pad_inches=0, dpi=600)
    print(f'Exported {name}: {fig.get_size_inches()*25.4} mm', flush=True)

def lr_label(row):
    receptor = str(row['receptor'])
    receptor = '(' + receptor.replace('_', '+') + ')' if '_' in receptor else receptor
    return str(row['ligand']) + '-' + receptor

ABBREVIATIONS = {
    'Proliferating T-cells': 'Prolif. T-cells', 'Proliferating AT2': 'Prolif. AT2',
    'Macrophages - IFN-activated': 'IFN-activated Mph',
    'Alveolar Macrophages': 'Alveolar Mph',
    'Activated Fibrotic FBs': 'Act. fibrotic FBs',
    'Interstitial Macrophages': 'Interst. Mph',
}

def abbrev(value):
    return ABBREVIATIONS.get(value, value)
