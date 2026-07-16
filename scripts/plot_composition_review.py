#!/usr/bin/env python3
"""Publication-style figure for the biological-accuracy composition review.

Panel A: observed phage fraction (of classified viruses) vs the literature
         expected band, per collection, coloured by verdict.
Panel B: observed unknown-family fraction vs its data-quality ceiling.

Reads validation/composition_scorecard.json. Okabe-Ito palette, DejaVu Sans.
Outputs PNG (150 dpi review) + PDF (vector) to docs/figures/.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.labelsize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# Okabe-Ito
OK = {"ok": "#009E73", "minor": "#56B4E9", "moderate": "#E69F00", "major": "#D55E00"}


def main() -> None:
    sc = json.loads(Path("validation/composition_scorecard.json").read_text())["sites"]
    order = {"major": 0, "moderate": 1, "minor": 2, "ok": 3}
    items = sorted(sc.items(), key=lambda kv: (order.get(kv[1]["site_verdict"], 9), int(kv[0])))

    labels, verdicts = [], []
    ph_obs, ph_lo, ph_hi = [], [], []
    fam_unk, unclassified = [], []
    for cid, s in items:
        labels.append(f"{cid}. {s['name'].split(' - ')[0].split(' (')[0]}")
        verdicts.append(s["site_verdict"])
        pf = next((f for f in s["findings"] if f["metric"] == "phage_fraction_of_classified"), None)
        if pf:
            ph_obs.append(pf["observed"]); ph_lo.append(pf["effective_band"][0]); ph_hi.append(pf["effective_band"][1])
        else:
            ph_obs.append(float("nan")); ph_lo.append(float("nan")); ph_hi.append(float("nan"))
        # Panel B: family-unassigned (ICTV omits family for many phages) vs the
        # genuinely unclassified fraction (no class even in NCBI).
        fam_unk.append(s["observed"]["unknown_family"]["abundance"])
        unclassified.append(s["observed"].get("unclassified", {}).get("abundance", 0.0))

    n = len(labels)
    y = list(range(n))[::-1]  # top = worst verdict

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 6.2), sharey=True)

    # Panel A: phage fraction
    for yi, lo, hi, obs, v in zip(y, ph_lo, ph_hi, ph_obs, verdicts):
        if lo == lo:  # not nan
            axA.plot([lo, hi], [yi, yi], color="#BBBBBB", lw=4, solid_capstyle="round", zorder=1)
        axA.scatter([obs], [yi], color=OK[v], s=34, zorder=3, edgecolor="white", linewidth=0.6)
    axA.set_xlim(-0.02, 1.02)
    axA.set_yticks(y)
    axA.set_yticklabels(labels)
    axA.set_xlabel("Phage fraction (classified, abundance)")
    axA.set_title("A", loc="left", fontsize=12, fontweight="bold")
    axA.text(0.5, n + 0.2, "grey = expected band, dot = observed", ha="center", fontsize=7, color="#555555")

    # Panel B: family-unassigned (light, ICTV omits family for phages) with the
    # genuinely-unclassified fraction (dark, no class in NCBI) overlaid.
    for yi, fu, uc in zip(y, fam_unk, unclassified):
        axB.barh(yi, fu, color="#CCCCCC", height=0.6, zorder=2)
        if uc > 0:
            axB.barh(yi, uc, color="#D55E00", height=0.6, zorder=3)
    axB.set_xlim(0, max(fam_unk) * 1.1 + 0.01)
    axB.set_xlabel("Fraction of community")
    axB.set_title("B", loc="left", fontsize=12, fontweight="bold")
    axB.text(axB.get_xlim()[1] * 0.5, n + 0.2, "grey = no ICTV family (expected for phages)",
             ha="center", fontsize=7, color="#555555")

    handles = [Line2D([0], [0], marker="o", color="white", markerfacecolor=OK[k], markersize=7,
                      label=k.capitalize()) for k in ["ok", "minor", "moderate", "major"]]
    handles.append(Line2D([0], [0], marker="s", color="white", markerfacecolor="#D55E00",
                          markersize=7, label="Unclassified (no class)"))
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, -0.02), title="Panel A: site verdict   |   Panel B bar overlay")

    fig.suptitle("ViroForge collections: observed viral composition vs literature expectation",
                 fontsize=10, y=0.99)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])

    outdir = Path("docs/figures")
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "composition_review.png", dpi=150, bbox_inches="tight")
    fig.savefig(outdir / "composition_review.pdf", bbox_inches="tight")
    print(f"wrote {outdir}/composition_review.png and .pdf")


if __name__ == "__main__":
    main()
