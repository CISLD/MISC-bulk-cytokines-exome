#!/usr/bin/env python3
"""Figure 6 on one page: panel A (STAT motifs, make_figures_v5.py) above panels B-E (FigR, figr_figure6_v5.R),
both rasterised at 600 dpi.
usage: figure6_composite.py <A.png> <BE.png> <out prefix>
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

a, be, out = sys.argv[1:4]
A, BE = mpimg.imread(a), mpimg.imread(be)
w = 10.0
ha, hb = w * A.shape[0] / A.shape[1], w * BE.shape[0] / BE.shape[1]
fig = plt.figure(figsize=(w, ha + hb))
for img, bottom, h in ((A, hb, ha), (BE, 0, hb)):
    ax = fig.add_axes([0, bottom / (ha + hb), 1, h / (ha + hb)])
    ax.imshow(img, interpolation="lanczos")
    ax.axis("off")
fig.savefig(out + ".pdf", dpi=600)
fig.savefig(out + ".png", dpi=400)
