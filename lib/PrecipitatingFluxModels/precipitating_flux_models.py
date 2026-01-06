"""
Python wrapper for PrecipitatingFluxModels.jl

Provides access to empirical precipitating electron flux models derived from
DMSP-ELFIN conjunction observations.

Setup
-----
1. Install dependencies: `pip install juliacall`
2. The Julia package will be automatically installed on first import via juliapkg.json

Example
-------
>>> from precipitating_flux_models import load_model, get_model, flux, n_flux, e_flux
>>> model = load_model()
>>> subset = get_model(model, mlat=70, ae=200)
>>> flux(subset, 50.0)  # Mean flux at 50 keV
>>> n_flux(subset, 30, 1000)  # Number flux 30-1000 keV
"""

import os as _os
from pathlib import Path as _Path

# Set JULIA_PROJECT to the package directory for juliapkg to find juliapkg.json
_pkg_dir = _Path(__file__).parent
_os.environ.setdefault("JULIA_PROJECT", str(_pkg_dir))

from juliacall import Main as jl

# Initialize Julia and load the package
jl.seval("using PrecipitatingFluxModels")

_PFM = jl.PrecipitatingFluxModels


def load_model(path=None):
    """
    Load an EmpiricalFluxModel from disk.

    Parameters
    ----------
    path : str, optional
        Path to the model file. If None, loads the default model.

    Returns
    -------
    model
        Julia EmpiricalFluxModel object.
    """
    if path is None:
        return _PFM.load_model()
    return _PFM.load_model(path)


def get_model(model, *, mlat=None, mlt=None, ae=None):
    """
    Query subset of models matching spatial bin criteria.

    Parameters
    ----------
    model : EmpiricalFluxModel
        The model to query.
    mlat : float, optional
        Filter by magnetic latitude bin (uses absolute value).
    mlt : float, optional
        Filter by magnetic local time bin.
    ae : float, optional
        Filter by AE index bin.

    Returns
    -------
    EmpiricalFluxModel
        Subset model containing only matching rows.
    """
    return _PFM.get_model(model, mlat=mlat, mlt=mlt, ae=ae)


def flux(model, energy):
    """
    Compute mean differential flux at given energy.

    Parameters
    ----------
    model : EmpiricalFluxModel
        The model (or subset) to evaluate.
    energy : float
        Energy in keV.

    Returns
    -------
    float
        Mean differential flux [#/cm²/s/sr/keV].
    """
    return float(model(energy))


def n_flux(model, Emin, Emax):
    """
    Compute mean number flux integrated over energy range.

    Parameters
    ----------
    model : EmpiricalFluxModel
        The model (or subset) to evaluate.
    Emin : float
        Minimum energy in keV.
    Emax : float
        Maximum energy in keV.

    Returns
    -------
    float
        Mean integrated number flux [#/cm²/s/sr].
    """
    return float(_PFM.n_flux(model, Emin, Emax))


def e_flux(model, Emin, Emax):
    """
    Compute mean energy flux integrated over energy range.

    Parameters
    ----------
    model : EmpiricalFluxModel
        The model (or subset) to evaluate.
    Emin : float
        Minimum energy in keV.
    Emax : float
        Maximum energy in keV.

    Returns
    -------
    float
        Mean integrated energy flux [keV/cm²/s/sr].
    """
    return float(_PFM.e_flux(model, Emin, Emax))


def nrow(model):
    """Return the number of parameter rows in the model."""
    return int(jl.seval("nrow")(model.parameters))


def spatial_bins(model):
    """
    Get the spatial bin edges.

    Returns
    -------
    dict
        Dictionary with 'mlat', 'mlt', 'ae' keys containing bin edge arrays.
    """
    bins = model.spatial_bins
    return {
        "mlat": list(bins.mlat),
        "mlt": list(bins.mlt),
        "ae": list(bins.ae),
    }


def energy_range(model):
    """Get the valid energy range (Emin, Emax) in keV."""
    return tuple(model.energy_range)
