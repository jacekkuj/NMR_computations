from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, List

import numpy as np
from scipy import stats
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


@dataclass(frozen=True)
class StatsResult:
    n: int
    mae: float
    rmse: float
    r2: float
    pearson_r: float
    pearson_p: float
    spearman_r: float
    spearman_p: float
    max_abs_error: float
    mean_error: float

    def to_dict(self) -> Dict[str, float]:
        return asdict(self)


def compute_stats(y_true: List[float], y_pred: List[float]) -> StatsResult:
    y_true_arr = np.asarray(y_true, dtype=float)
    y_pred_arr = np.asarray(y_pred, dtype=float)

    if y_true_arr.shape != y_pred_arr.shape:
        raise ValueError("Experimental and computed arrays must have the same length.")
    if y_true_arr.size < 2:
        raise ValueError("Need at least 2 points to compute correlation / R² reliably.")

    err = y_true_arr - y_pred_arr

    mae = float(mean_absolute_error(y_true_arr, y_pred_arr))
    mse = float(mean_squared_error(y_true_arr, y_pred_arr))
    rmse = float(np.sqrt(mse))
    r2 = float(r2_score(y_true_arr, y_pred_arr))

    pearson_r, pearson_p = stats.pearsonr(y_true_arr, y_pred_arr)
    spearman_r, spearman_p = stats.spearmanr(y_true_arr, y_pred_arr)

    return StatsResult(
        n=int(y_true_arr.size),
        mae=mae,
        rmse=rmse,
        r2=r2,
        pearson_r=float(pearson_r),
        pearson_p=float(pearson_p),
        spearman_r=float(spearman_r),
        spearman_p=float(spearman_p),
        max_abs_error=float(np.max(np.abs(err))),
        mean_error=float(np.mean(err)),
    )


def stats_label(s: StatsResult, name: str) -> str:
    # krótka etykieta do legendy / adnotacji
    return f"{name}: R²={s.r2:.3f}, RMSE={s.rmse:.3f} ppm"
