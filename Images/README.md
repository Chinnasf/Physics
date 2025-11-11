# My typical settings

```Python

%load_ext autoreload
%autoreload 2

from IPython.core.debugger import Pdb #Pdb().set_trace()
```

### Convert from `.ipynb` to `.py`: in bash

```Python
$ jupyter nbconvert --to script your_notebook.ipynb
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
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import scienceplots

import sys
root_path = Path.cwd().resolve().parent
sys.path.append(str(root_path))

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

## Cool Python stuff

This may be needed for the future:

```Python
# check line intersections
def cross2(a, b):
    """
    Expects `a` and `b` to be NumPy arrays (or anything that 
    behaves similar, e.g. JAX arrays, CuPy arrays).

    shape of `a` and `b` = (...,2); where `...` = any dim.
    So valid examples are:

    * Shape (2,) → a single 2D vector
    * Shape (N, 2) → batch of N vectors
    * Shape (M, N, 2) → grid of vectors
    * Shape (M, N, K, 2) → high-dimensional tensor of vectors

    In Python, `...` is    <class 'ellipsis'>.
    see it yourself:       print(type(...))
    This allows you to write extremely general vectorized functions.

    Example: 
    a = np.random.randn(10, 20, 2)  # shape (10, 20, 2)
    a[..., 0].shape  # -> (10, 20)
    a[..., 1].shape  # -> (10, 20)
    """
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]
```
