"""Render the two manuscript figures. Run from any directory."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator
import numpy as np

OUTPUT = Path(__file__).resolve().parent
PDF_META = {"Creator": "Matplotlib", "CreationDate": None, "ModDate": None}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42,
                     "axes.linewidth": .65, "xtick.major.width": .65, "ytick.major.width": .65})
INK, BLUE, RED = "#28343a", "#166579", "#a64023"


def angles(h, t):
    # Parametrize the energy ellipse in the signed Morse coordinates, then invert.
    u = np.sqrt(2*h/np.pi) * np.cos(t)
    v = np.sqrt(2*h) * np.sin(t)

    def inverse(w):
        return np.arctan2(np.sign(w)*np.sqrt(-np.expm1(-w*w)), np.exp(-w*w/2))

    xi, eta = inverse(u), inverse(v)
    error = abs(-np.pi*np.log(np.cos(xi)) - np.log(np.cos(eta)) - h)
    if np.max(error) > 1e-12:
        raise ValueError("energy residual in the sampled oval")
    return xi, eta


def arrow(ax, x, y, index, displacement):
    other = index + displacement
    ax.annotate("", xy=(x[other], y[other]), xytext=(x[index], y[index]),
                arrowprops={"arrowstyle": "-|>", "color": BLUE, "lw": .9,
                            "mutation_scale": 9, "shrinkA": 0, "shrinkB": 0})


def cells():
    fig = plt.figure(figsize=(6.9, 3.05), layout="constrained")
    ax, angle_ax = fig.subplots(1, 2)
    t = np.linspace(0, 2*np.pi, 361)
    roots = np.cos((np.arange(4)+.5)*np.pi/4)
    for root in roots:
        ax.plot([root, root], [roots[-1], roots[0]], color=INK, lw=.8)
        ax.plot([roots[-1], roots[0]], [root, root], color=INK, lw=.8)
    ax.add_patch(Rectangle((roots[2], roots[2]), roots[1]-roots[2], roots[1]-roots[2],
                           facecolor="#ebf4f4", edgecolor="none", zorder=-1))
    for j in range(1, 4):
        for i in range(1, 4):
            ax.plot(np.cos(j*np.pi/4), np.cos(i*np.pi/4), ".", color=INK, ms=3)
            for h in (.35, 1.25):
                xi, eta = angles(h, t)
                x, y = np.cos((j*np.pi+xi)/4), np.cos((i*np.pi+eta)/4)
                ax.plot(x, y, color=BLUE, lw=.85)
                if h == 1.25:
                    arrow(ax, x, y, 52, 12*(-1)**(i+j))
    for x in roots:
        ax.plot([x]*4, roots, "s", ms=2.7, color=INK)
    ax.set(aspect="equal", xlabel="$x$", ylabel="$y$", xlim=(-1, 1), ylim=(-1, 1),
           xticks=[-1, 0, 1], yticks=[-1, 0, 1])
    ax.set_title("Invariant grid", fontsize=10)
    side = np.pi/2
    cmap = LinearSegmentedColormap.from_list("reciprocal", ["#ffffff", "#8ec5bf"])
    norm = Normalize(1, 1.5)
    eta = np.linspace(-side, side, 1025)
    weight = np.broadcast_to((1/(1-np.cos(eta)/3))[:,None], (len(eta), 8))
    angle_ax.imshow(weight, extent=(-side, side, -side, side), origin="lower",
                    cmap=cmap, norm=norm, interpolation="bilinear", zorder=-1)
    angle_ax.plot([0, side], [0, 0], "--", color=INK, lw=.85)
    for h in (.35, 1.25):
        x, y = angles(h, t)
        angle_ax.plot(x, y, color=BLUE, lw=1.15)
        if h == 1.25:
            arrow(angle_ax, x, y, 52, 10)
    angle_ax.plot(0, 0, ".", color=INK)
    angle_ax.set(aspect="equal", xlabel=r"$\xi$", ylabel=r"$\eta$",
                 xlim=(-side, side), ylim=(-side, side),
                 xticks=[-side, 0, side], yticks=[-side, 0, side],
                 xticklabels=[r"$-\pi/2$", "0", r"$\pi/2$"],
                 yticklabels=[r"$-\pi/2$", "0", r"$\pi/2$"])
    angle_ax.set_title("Central cell  $(j,i)=(2,2)$", fontsize=10)
    bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=angle_ax,
                       ticks=[1, 1.25, 1.5], shrink=.70, aspect=18, pad=.08)
    bar.set_label("reciprocal multiplier", fontsize=8)
    fig.savefig(OUTPUT/"cells.pdf", metadata=PDF_META)
    plt.close(fig)


def frequencies():
    fig, (left, right) = plt.subplots(1, 2, figsize=(6.9, 3.25), layout="constrained")
    for ax in (left, right):
        ax.plot([.7, 11.3], [11.3, .7], lw=.75, color="#879398")
        ax.set(xlim=(.2, 12), ylim=(.2, 12), aspect="equal", xlabel="horizontal index $a$",
               ylabel="vertical index $b$")
        ax.xaxis.set_major_locator(FixedLocator(range(1, 12, 2)))
        ax.yaxis.set_major_locator(FixedLocator(range(1, 12, 2)))
        ax.spines[["top", "right"]].set_visible(False)
    for a in range(1, 12, 2):
        for b in range(1, 13-a, 2):
            retained = 2 - int(a > 4) - int(b > 4)
            left.text(a, b, str(retained), color=INK if retained else "#63727a",
                      ha="center", va="center", fontsize=9,
                      bbox=dict(facecolor="white", edgecolor="none", pad=.5))
            if a+b <= 8 and retained:
                right.plot(a, b, "o", mfc=BLUE if retained == 2 else "white",
                           mec=BLUE, ms=6.6, mew=1.1)
    high = [(a+4,b) for a in range(1, 9, 2) for b in (1, 3) if 4 < a+b <= 8]
    right.plot(*zip(*high), ls="", marker="D", color=RED, ms=5.8)
    right.plot([7,7], [1,3], ":", lw=1.4, color=INK, zorder=0)
    right.plot([], [], "o", mfc="white", mec=BLUE, ms=6, label="one border form")
    right.plot([], [], "o", mfc=BLUE, mec=BLUE, ms=6, label="two border forms")
    right.plot([], [], "D", color=RED, ms=5, label="one added form")
    right.legend(loc="upper right", fontsize=7, frameon=False, handlelength=1.1)
    left.set_title("Unweighted periods (rank 22)", fontsize=10)
    right.set_title("Augmented test (14 + 4 forms)", fontsize=10)
    fig.savefig(OUTPUT/"frequencies.pdf", metadata=PDF_META)
    plt.close(fig)


if __name__ == "__main__":
    cells()
    frequencies()
