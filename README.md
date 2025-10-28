# Langmuir Probe Data Analysis Functions

This repository provides a set of Python functions and classes for **analyzing Langmuir probe measurements** in low-temperature plasmas. The tools allow you to process I-V characteristics, estimate plasma parameters, calculate the electron energy distribution function (EEDF), and iteratively refine electron density and temperature.

---

## Features

- Read Langmuir probe CSV data.

- Fit I-V curves using polynomials or hyperbolic tangent functions.
- Determine **floating potential** and **plasma potential**.
- Estimate **electron temperature (T_e)** and **electron density (n_e)**.
- Iterative refinement using **Laframboise ion collection model**.
- Compute **Debye length** and probe-radius-to-Debye-length ratio.
- Calculate **electron energy distribution function (EEDF)** from probe data.
- Plot I-V characteristics and fitted curves.