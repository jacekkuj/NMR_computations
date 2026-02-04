from __future__ import annotations

from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

from nmr_app.stats import StatsResult, stats_label


def fig_compare_raw_vs_calibrated(
    exp: Sequence[float],
    calc_raw: Sequence[float],
    calc_cal: Optional[Sequence[float]],
    stats_raw: StatsResult,
    stats_cal: Optional[StatsResult],
):
    exp = np.asarray(exp, dtype=float)
    raw = np.asarray(calc_raw, dtype=float)
    cal = np.asarray(calc_cal, dtype=float) if calc_cal is not None else None

    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)

    # punkty
    ax.scatter(raw, exp, label=stats_label(stats_raw, "Raw"), s=60)
    if cal is not None and stats_cal is not None:
        ax.scatter(cal, exp, label=stats_label(stats_cal, "Calibrated"), s=60)

    # linia idealna y=x (to jest ta “niebieska/pomarańczowa” w Twoich wcześniejszych rysunkach)
    mn = float(min(exp.min(), raw.min(), (cal.min() if cal is not None else raw.min())))
    mx = float(max(exp.max(), raw.max(), (cal.max() if cal is not None else raw.max())))
    ax.plot([mn, mx], [mn, mx], linestyle="--", alpha=0.6, label="Ideal (y = x)")

    ax.set_xlabel("Calculated δ (ppm)")
    ax.set_ylabel("Experimental δ (ppm)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    return fig
