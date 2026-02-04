from __future__ import annotations

import io
from typing import List, Tuple

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from nmr_app.parser import parse_text_blob
from nmr_app.stats import compute_stats
from nmr_app.plots import fig_compare_raw_vs_calibrated


def _try_read_log(uploaded_file) -> str:
    data = uploaded_file.getvalue()
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("latin-1")


def linear_calibration(calc: List[float], exp: List[float]) -> Tuple[float, float, List[float]]:
    x = np.asarray(calc).reshape(-1, 1)
    y = np.asarray(exp)
    model = LinearRegression()
    model.fit(x, y)
    a = float(model.coef_[0])
    b = float(model.intercept_)
    return a, b, model.predict(x).tolist()


def fig_scatter(exp, calc, a=None, b=None):
    x = np.asarray(calc)
    y = np.asarray(exp)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    ax.scatter(x, y)

    mn = float(min(x.min(), y.min()))
    mx = float(max(x.max(), y.max()))
    ax.plot([mn, mx], [mn, mx])

    if a is not None:
        xx = np.linspace(mn, mx, 200)
        ax.plot(xx, a * xx + b)

    ax.set_xlabel("Calculated δ (ppm)")
    ax.set_ylabel("Experimental δ (ppm)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def fig_residuals(exp, pred):
    y = np.asarray(exp)
    res = y - np.asarray(pred)

    fig, ax = plt.subplots(figsize=(7, 4), dpi=200)
    ax.scatter(y, res)
    ax.axhline(0)
    ax.set_xlabel("Experimental δ (ppm)")
    ax.set_ylabel("Residual (Exp − Pred)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


st.set_page_config("1H NMR assignment & validation", layout="wide")
st.title("¹H NMR – manual assignment, grouping & calibration")

with st.sidebar:
    uploaded = st.file_uploader("Gaussian NMR log (*.log)", type=["log", "txt"])
    h_ref = st.number_input("H_ref (TMS shielding)", value=31.6674, format="%.4f")
    do_cal = st.checkbox("Apply linear calibration", value=True)

if not uploaded:
    st.info("Upload Gaussian NMR log to start")
    st.stop()

# -------- Parse Gaussian log --------
try:
    parsed = parse_text_blob(_try_read_log(uploaded))
except Exception as e:
    st.error(f"Parsing error: {e}")
    st.stop()

computed = parsed.computed_shifts(h_ref)

df_atoms = pd.DataFrame(
    {
        "Atom": parsed.atom_numbers,
        "ComputedShift_ppm": computed,
        "ExperimentalShift_ppm": "",
    }
)

st.subheader("Manual experimental assignment (per proton)")
df_atoms = st.data_editor(
    df_atoms,
    use_container_width=True,
    hide_index=True,
    column_config={
        "ExperimentalShift_ppm": st.column_config.TextColumn(
            "Enter experimental δ (ppm)"
        )
    },
)

# -------- Validation --------
if df_atoms["ExperimentalShift_ppm"].astype(str).str.strip().eq("").any():
    st.info("Fill experimental chemical shift for ALL protons.")
    st.stop()

try:
    df_atoms["ExperimentalShift_ppm"] = (
        df_atoms["ExperimentalShift_ppm"]
        .astype(str)
        .str.replace(",", ".")
        .astype(float)
    )
except Exception:
    st.error("All experimental shifts must be numeric (e.g. 4.23)")
    st.stop()

# -------- Automatic grouping: peaks --------
df_peaks = (
    df_atoms
    .groupby("ExperimentalShift_ppm", as_index=False)
    .agg(
        CalcMean=("ComputedShift_ppm", "mean"),
        N_Assigned=("ComputedShift_ppm", "count"),
    )
    .sort_values("ExperimentalShift_ppm")
    .reset_index(drop=True)
)

exp_vals = df_peaks["ExperimentalShift_ppm"].tolist()
calc_vals = df_peaks["CalcMean"].tolist()

st.subheader("Peak-level grouped data")
st.dataframe(df_peaks, use_container_width=True, hide_index=True)

# -------- Statistics --------
stats_raw = compute_stats(exp_vals, calc_vals)

a = b = None
pred_cal = None
stats_cal = None

if do_cal:
    a, b, pred_cal = linear_calibration(calc_vals, exp_vals)
    stats_cal = compute_stats(exp_vals, pred_cal)

st.subheader("Statistics")

c1, c2 = st.columns(2)
with c1:
    st.markdown("### Raw")
    st.metric("N peaks", stats_raw.n)
    st.metric("MAE", f"{stats_raw.mae:.3f}")
    st.metric("RMSE", f"{stats_raw.rmse:.3f}")
    st.metric("R²", f"{stats_raw.r2:.3f}")

with c2:
    st.markdown("### Calibrated")
    if do_cal:
        st.write(f"δ_exp = {a:.5f} · δ_calc + {b:.5f}")
        st.metric("MAE", f"{stats_cal.mae:.3f}")
        st.metric("RMSE", f"{stats_cal.rmse:.3f}")
        st.metric("R²", f"{stats_cal.r2:.3f}")

# -------- Plots --------
st.subheader("Publication plots")

fig_cmp = fig_compare_raw_vs_calibrated(
    exp=exp_vals,
    calc_raw=calc_vals,
    calc_cal=pred_cal if do_cal else None,
    stats_raw=stats_raw,
    stats_cal=stats_cal if do_cal else None,
)
st.pyplot(fig_cmp)

pred_use = pred_cal if do_cal else calc_vals
fig2 = fig_residuals(exp_vals, pred_use)
st.pyplot(fig2)

# -------- Downloads --------
st.subheader("Downloads")

csv_atoms = df_atoms.to_csv(index=False).encode()
csv_peaks = df_peaks.to_csv(index=False).encode()

st.download_button("Download per-atom assignment CSV", csv_atoms, "assignment_atoms.csv")
st.download_button("Download peak-level results CSV", csv_peaks, "nmr_peaks_results.csv")
