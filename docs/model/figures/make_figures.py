"""Рисунки статьи docs/model/model.tex из case.json (его пишет app/scripts/model-case.ts).

Запуск:  python3 make_figures.py   (нужны numpy и matplotlib)
"""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "case.json")))

plt.rcParams.update(
    {
        "font.family": "DejaVu Serif",
        "mathtext.fontset": "dejavuserif",
        "font.size": 9,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    }
)
FAST = "#111111"
SLOW = "#0a8f2a"
FAST_C = "#b2182b"


def xt(field, segs, events, name, tmax=None, label_events=True):
    t = np.array(field["t"])
    P = np.array(field["P"])
    xmax = field["xMax"]
    tmax = tmax or t[-1]
    fig, ax = plt.subplots(figsize=(6.3, 3.3))
    lim = np.percentile(np.abs(P), 99.5) or 1e-6
    im = ax.imshow(
        P.T,
        origin="lower",
        aspect="auto",
        extent=[t[0], t[-1], 0, xmax],
        cmap="RdBu_r",
        norm=TwoSlopeNorm(0, -lim, lim),
        interpolation="bilinear",
        rasterized=True,
    )
    for s in segs:
        c = FAST if s["fam"] == "fast" else SLOW
        ls = "-" if s["kind"] == "shock" else (0, (3, 1.5))
        ax.plot([s["t0"], s["t1"]], [s["x0"], s["x1"]], color=c, lw=1.0, ls=ls, solid_capstyle="butt")
    k = 0
    for e in events:
        if e["t"] > tmax:
            continue
        if e["type"] in ("boundary", "collision"):
            k += 1
            ax.plot(e["t"], e["x"], "o", ms=4, mfc="white", mec="black", mew=0.8, zorder=5)
            if label_events:
                ax.annotate(
                    f"$E_{k}$", (e["t"], e["x"]), xytext=(2, 6), textcoords="offset points", fontsize=8, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85),
                )
        elif e["type"] == "untrack":
            ax.plot(e["t"], e["x"], "x", ms=4, color="black", mew=0.8, zorder=5)
    ax.set_xlim(0, tmax)
    ax.set_ylim(0, xmax)
    ax.set_xlabel(r"$t\,c_1/L$")
    ax.set_ylabel(r"$x_1/L$")
    cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.04)
    cb.set_label(r"$p-p_0$")
    cb.formatter.set_powerlimits((-2, 2))
    fig.savefig(os.path.join(HERE, name))
    plt.close(fig)


xt(d["fieldZoom"], d["segmentsZoom"], d["log"], "xt_zoom.pdf", tmax=1.6)
xt(d["fieldFull"], d["segmentsFull"], d["log"], "xt_full.pdf", label_events=False)

# профили
fig, axs = plt.subplots(2, 1, figsize=(6.3, 3.8), sharex=True)
p0, s0 = d["F0"][0], d["F0"][2]
sel = [pr for pr in d["profiles"] if round(pr["t"], 2) in (0.45, 1.0, 1.5)]
for pr, c in zip(sel, ["#2166ac", "#d6604d", "#1b7837"]):
    axs[0].plot(pr["x"], np.array(pr["p"]) * 1e3, color=c, lw=0.9, label=f"$t={pr['t']:.2f}$")
    axs[1].plot(pr["x"], np.array(pr["s"]) * 1e3, color=c, lw=0.9)
for ax in axs:
    ax.axhline(0, color="0.6", lw=0.5)
axs[0].set_ylabel(r"$p\cdot10^{3}$")
axs[1].set_ylabel(r"$s\cdot10^{3}$")
axs[1].set_xlabel(r"$x_1/L$")
axs[0].legend(ncol=3, fontsize=7, frameon=False, loc="upper right")
axs[1].set_xlim(0, 1.8)
fig.savefig(os.path.join(HERE, "profiles.pdf"))
plt.close(fig)

# скорости характеристик как функции p (разномодульность)
sc = d["speedCurve"]
p = np.array([q["p"] for q in sc]) * 1e3
fig, ax = plt.subplots(figsize=(4.2, 2.4))
ax.plot(p, [q["cf"] for q in sc], color=FAST_C, lw=1.1, label=r"$c_f$")
ax.plot(p, [q["cs"] for q in sc], color=SLOW, lw=1.1, label=r"$c_s$")
for x in (p0 * 1e3,):
    ax.axvline(x, color="0.5", lw=0.5, ls=":")
ax.set_xlabel(r"$p\cdot10^{3}$  ($s=s_0$)")
ax.set_ylabel(r"скорость$/c_1$")
ax.legend(frameon=False, fontsize=8)
fig.savefig(os.path.join(HERE, "speeds.pdf"))
plt.close(fig)
print("ok")
