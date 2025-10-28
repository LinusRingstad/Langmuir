#### ---------- EEDF calculation ----------

import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from langmuir_LR.basic import read_file, plasma_potential
from dataclasses import dataclass
from functools import cached_property
from numba import njit
import matplotlib.pyplot as plt



k_eV = 8.617333262145e-5  # eV/K


@njit(cache=True, fastmath=True)
def maxwellian_pdf(E_eV: np.ndarray, Te_eV: float) -> np.ndarray:
    """Per-particle Maxwellian energy PDF"""
    E = np.maximum(E_eV, 0.0)
    Te = max(Te_eV, 1e-9)
    coeff = 2.0 / np.sqrt(np.pi) / (Te ** 1.5)
    return coeff * np.sqrt(E) * np.exp(-E / Te)



def compute_eedf(
    settings: object,
    I_e: np.ndarray,
    ne_m3: float,
    smooth_span_V: float = 5.0,
    polyorder: int = 3
):
    
    e = settings.e  # C
    m_e = settings.m_e  # kg
    
    A_probe_m2 = settings.Aprob


    V_0,I_0,coeffs = read_file(settings)


    Vp, Ip = plasma_potential(V_0, coeffs)

    VminusVp = V_0 - Vp



    """Compute per-particle EEDF from probe IV curve."""
    VminusVp, I_e = map(np.asarray, (VminusVp, I_e))

    # Retarding region (V <= Vp
    m = VminusVp <= 0
    x, I = VminusVp[m], I_e[m]


    dV = np.median(np.diff(x))
    win = max(7, int(round(smooth_span_V / max(abs(dV), 1e-12))) | 1)
    win = min(win, len(I) - (len(I) + 1) % 2)  # ensure odd and within array

    d2I = savgol_filter(I, window_length=win, polyorder=polyorder,
                        deriv=2, delta=dV, mode='interp')

    # Energy in eV and J
    E_eV = -x
    E_J = e * np.clip(E_eV, 0, None)

    # Druyvesteyn relation (SI)
    pref = 4.0 / (A_probe_m2 * e**2) * np.sqrt(m_e / (2.0 * e))
    f_J = pref * np.sqrt(E_J) * d2I

    # Convert to m^-3/eV and per-particle
    f_eV = f_J / e
    p_eV = f_eV / ne_m3

    order = np.argsort(E_eV)
    return tuple(arr[order] for arr in (E_eV, p_eV, f_eV, d2I))


# --- EEDF normalization helper ---
def normalize_eedf(E_eV: np.ndarray, f_eV: np.ndarray, ne_m3: float):
    """Normalize EEDF"""
    f_eV_pos = np.clip(f_eV, 0, None)
    # print(f_eV_pos)
        
    # integrate
    area = np.trapezoid(f_eV_pos, E_eV)
    
    # normalize to integral = 1
    f_norm = f_eV_pos / area
    
    print('Normalized EEDF integral:', area)
    return f_norm, area


@dataclass(slots=True)
class EEDFBundle:
    E_eV: np.ndarray
    p_eV: np.ndarray
    f_eV: np.ndarray
    ne: float
    Te_K: float

    @property
    def Te_eV(self) -> float:
        return self.Te_K * k_eV

    @property
    def maxwellian(self) -> np.ndarray:
        return maxwellian_pdf(self.E_eV, self.Te_eV)


def plot_eedf_vs_maxwellian(bundle: EEDFBundle, calibrate: bool = True, figsize=(7, 4.5)):
    """Plot EEDF (per-particle) and Maxwellian."""
    E, p, f = bundle.E_eV, bundle.p_eV, bundle.f_eV
    if calibrate:
        f_cal, _ = normalize_eedf(E, f, bundle.ne)
        p = f_cal

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(E, np.clip(p, 0, None), label="EEDF per particle", lw=2)
    ax.plot(E, bundle.maxwellian, "--", label=f"Maxwellian (Te={bundle.Te_eV:.2f} eV)", lw=2)
    ax.set(xlabel="Energy  E = Vp - V  (eV)", ylabel="Relative Probability Density (a.u.)",)
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    return fig, ax
