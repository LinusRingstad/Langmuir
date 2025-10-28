# langmuir-LR: Basic Langmuir Probe Data Analysis in Python

Langmuir-LR is a Python package for **analyzing Langmuir probe data**. It provides tools for reading I-V curves, estimating plasma parameters (electron temperature, density, plasma potential, floating potential), and computing the **electron energy distribution function (EEDF)**.

---

## Installation

Install the package locally (after building):

```bash
pip install .
```

or from **TestPyPI**:

```bash
pip install -i https://test.pypi.org/simple/ langmuir-LR
```

---

## Basic Usage

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from langmuir_LR.basic import lp_Settings, read_file, plasma_potential, floating_potential, guess_Te, guess_ne, interp_laframboise
from langmuir_LR.eedf import compute_eedf, EEDFBundle, plot_eedf_vs_maxwellian
```

---

## 1. Probe Settings

`lp_Settings` holds constants and probe parameters:

```python
settings = lp_Settings()
settings.set_probe_r(0.001)     # probe radius in meters
settings.set_probe_l(0.004)     # probe length in meters
settings.set_params(pressure=75.1, power=1000, seqno=1, length=2)
settings.set_skipheader(61) #skips desired num of rows in data file
settings.set_skipfooter(5) 
```

Attributes include physical constants (e, k (J), k (eV), pi), probe geometry, and default file naming (this will be generalized later but for now needs to be changed.

---

## 2. Reading I-V Data

```python
V, I, coeffs = read_file(settings, degree=9)
```

* `V`: Voltage array [V]
* `I`: Current array [A]
* `coeffs`: Polynomial fit coefficients

---

## 3. Estimating Plasma Parameters

```python
Vp, Ip = plasma_potential(V, coeffs)  # Plasma potential & saturation current
Vf = floating_potential(V, coeffs)    # Floating potential
Te_guess = guess_Te(Vp, Vf, settings) # Electron temperature in eV
ne_guess = guess_ne(Ip, Te_guess, settings)  # Electron density in m^-3
```

* `plasma_potential`: Finds `Vp` and ion saturation current `Ip` using polynomial curvature
* `floating_potential`: Finds `Vf` where `I ≈ 0`
* `guess_Te`: Estimates electron temperature
* `guess_ne`: Estimates electron density from ion current

---

## 4. Iterative Electron Density & Temperature

```python
from langmuir_LR.compute import find_ne_te_iterative

Te, ne, Ie = find_ne_te_iterative(settings, max_iterations=1000, tolerance=1e-4)
```

* Iteratively refines `Te` and `ne` using the **Laframboise correction**
* Returns electron temperature (eV), density (m^-3), and electron current `Ie`

---

## 5. Electron Energy Distribution Function (EEDF)

Compute per-particle EEDF:

```python
E_eV, p_eV, f_eV, d2I = compute_eedf(settings, Ie, ne, smooth_span_V=5.0, polyorder=3)
```

* `E_eV`: Energy axis in eV
* `p_eV`: Per-particle EEDF (dimensionless, integrates to 1)
* `f_eV`: Absolute EEDF (m^-3 / eV)
* `d2I`: Second derivative of I-V curve

Normalize EEDF:

```python
from langmuir_LR.eedf import normalize_eedf

p_norm, area = normalize_eedf(E_eV, f_eV, ne)
```

* Ensures integral of EEDF over energy = 1

---

## 6. Bundling EEDF Data

```python
from langmuir_LR.eedf import EEDFBundle

bundle = EEDFBundle(E_eV, p_eV, f_eV, ne, Te*11604)  # Te in K
```

* `EEDFBundle` stores all EEDF data and plasma parameters
* `bundle.maxwellian` gives the Maxwellian PDF for comparison

---

## 7. Plotting EEDF vs Maxwellian

```python
fig, ax = plot_eedf_vs_maxwellian(bundle, calibrate=True)
plt.show()
```

* Shows per-particle EEDF and corresponding Maxwellian distribution
* Automatically normalizes EEDF if `calibrate=True`

---

## Example Workflow

```python
settings = lp_Settings()
settings.set_params(pressure=75.1, power=1000, seqno=1, length=2)

# Read I-V curve
V, I, coeffs = read_file(settings)

# Plasma parameters
Vp, Ip = plasma_potential(V, coeffs)
Vf = floating_potential(V, coeffs)
Te_init = guess_Te(Vp, Vf, settings)
ne_init = guess_ne(Ip, Te_init, settings)

# Refined Te and ne
Te, ne, Ie = find_ne_te_iterative(settings)

# Compute EEDF
E_eV, p_eV, f_eV, d2I = compute_eedf(settings, Ie, ne)

# Normalize and bundle
from langmuir_LR.eedf import normalize_eedf, EEDFBundle
p_norm, _ = normalize_eedf(E_eV, f_eV, ne)
bundle = EEDFBundle(E_eV, p_norm, f_eV, ne, Te*11604)

# Plot
fig, ax = plot_eedf_vs_maxwellian(bundle)
plt.show()
```

---

## References

* Laframboise, J. G., *Theory of Spherical and Cylindrical Langmuir Probes in a Collisionless, Maxwellian Plasma*, University of Toronto, 1966.

---

## License

MIT License
