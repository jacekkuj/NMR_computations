# NMR Computations (Gaussian ¹H NMR → assignment, grouping & calibration)

[![CI](https://github.com/jacekkuj/NMR_computations/actions/workflows/ci.yml/badge.svg)](https://github.com/jacekkuj/NMR_computations/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue.svg)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)
[![Streamlit](https://img.shields.io/badge/Streamlit-app-ff4b4b.svg)](#running-the-app)

A lightweight Streamlit app that:
- parses **Gaussian** logs in the section `SCF GIAO Magnetic shielding tensor (ppm):`
- extracts **¹H isotropic shieldings**
- converts them to chemical shifts using `δ(¹H) = H_ref − σ_iso`
- supports manual experimental assignment per proton
- performs peak-level grouping and optional **linear calibration** (Exp vs Calc)
- reports statistics (MAE, RMSE, R², Pearson/Spearman) and publication-ready plots

## Demo (screenshots)

> Add your screenshots here (e.g. `docs/img/ui.png`, `docs/img/plots.png`)

## Installation

### Option A: pip (recommended)
```bash
python -m venv .venv
# Windows PowerShell:
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
