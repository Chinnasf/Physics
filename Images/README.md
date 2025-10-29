# My typical settings

```Python

%load_ext autoreload
%autoreload 2

from IPython.core.debugger import Pdb #Pdb().set_trace()
```

## SIMPLE

```Python
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import scienceplots
plt.style.use(['science', 'notebook'])
```

## MORE COMPLETE

```Python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import sys
root_path = Path.cwd().resolve().parent
sys.path.append(str(root_path))

import scienceplots
from matplotlib import animation
from matplotlib.animation import PillowWriter

# Plot settings
plt.style.use(['science', 'notebook'])
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Nimbus Roman', 'Times', 'Times New Roman']
mpl.rcParams['mathtext.fontset'] = 'stix'  # good Times-like math symbols
# Another option that will align all the Greek/math symbols with the same Times-like proportions.
# mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['image.cmap'] = 'jet' 
mpl.rcParams['xtick.direction'] = 'in'  
mpl.rcParams['ytick.direction'] = 'in'   
mpl.rcParams['figure.dpi'] = 150
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['legend.fontsize'] = 12
```

# What journals (e.g., IEEE, Elsevier) prefer for submission

```Python
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---- STYLE (you already have most of this) ----
mpl.rcParams.update({
    # Fonts: open-source Times clone + TrueType embedding
    'font.family': 'serif',
    'font.serif': ['Nimbus Roman', 'Times', 'Times New Roman'],
    'mathtext.fontset': 'stix',
    'pdf.fonttype': 42,   # embed TrueType in PDF
    'ps.fonttype': 42,    # embed TrueType in EPS/PS

    # Legibility
    'axes.titlesize': 10,   # pt
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'lines.linewidth': 1.0,

    # Figure export
    'savefig.dpi': 600,     # for raster exports like TIFF/PNG
})

# ---- SIZE: choose one ----
w_single, h_single = 3.5, 2.2   # IEEE single-column (inches)
w_double, h_double = 7.16, 3.5  # IEEE double-column (inches)

fig, ax = plt.subplots(figsize=(w_single, h_single), layout='constrained')

# ... your plotting code here ...

# ---- EXPORT (vector preferred) ----
fig.savefig("figure1.pdf", bbox_inches="tight", pad_inches=0.01, metadata={
    'Title': 'SXR emissivity and signals',
    'Author': 'Karina Fuentes',
    'Subject': 'Tomography validation figure',
    'Keywords': 'soft X-ray, tomography, emissivity'
})

# If a raster is requested (e.g., Elsevier wants TIFF at 600 dpi):
fig.savefig("figure1.tiff", dpi=600, bbox_inches="tight", pad_inches=0.01)
```
