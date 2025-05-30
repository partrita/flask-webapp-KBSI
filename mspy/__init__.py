#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MSPy Package
============

A package for mass spectrometry data analysis.
"""

# Import functions from py_calculations to make them available
# at the mspy package level, e.g., mspy.signal_area()
from .py_calculations import *

# Optionally, define __all__ to specify what `from mspy import *` imports
# This can be useful to avoid polluting the namespace if py_calculations
# has many internal/helper functions not meant for direct user access.
# For now, let's import all, assuming most are intended for use.

__all__ = [
    # List all function and class names from py_calculations intended for public API
    # Example:
    "MBox", "MNoise", "array_print",
    "signal_median", "signal_interpolate_y", "signal_interpolate_x",
    "signal_locate_x", "signal_locate_max_y", "signal_box",
    "signal_intensity", "signal_centroid", "signal_width", "signal_area",
    "signal_noise", "signal_local_maxima",
    "signal_crop", "signal_offset", "signal_multiply", "signal_normalize",
    "signal_smooth_ma",
    "formula_composition", # Though it's a placeholder
    "signal_gaussian", "signal_lorentzian", "signal_gausslorentzian",
    "signal_profile", "signal_profile_to_raster"
]

# Import functions from other modules if they are part of the public API
# For example, if mod_basics, mod_formulator etc. have public functions:
from .mod_basics import delta, mz, md, nominalmass, rdbe, frules
from .mod_formulator import formulator
from .mod_pattern import pattern #, gaussian, lorentzian, gausslorentzian, profile, matchpattern
# Note: gaussian, lorentzian, etc. from mod_pattern are now in py_calculations
# and imported above. If mod_pattern's versions are different and preferred,
# they might need aliasing or selective import. For now, assuming py_calculations'
# versions are the replacements.

from .obj_compound import compound
# from .obj_peak import peak # Assuming this would be defined if it existed
# from .obj_peaklist import peaklist # Assuming this would be defined

# Filter out the functions that are also defined in py_calculations from mod_basics and mod_pattern
# to avoid confusion, if we want py_calculations to be the primary source.
# For now, the last import takes precedence if names clash.

# Add version
__version__ = "0.1.0" # Example version

# It's also good practice to update the __init__.py in the `mspy` directory
# to make it clear that `calculations.c` is no longer part of the package.
# This has been done by deleting the C file and not including it in setup.py for compilation.
# The setup.py would now build py_calculations if it were a C extension, but it's pure Python.
# So, no setup.py changes needed for this specific Python module, but if there was a build process,
# it would need to be updated.

# The formulator function was originally in mspy.py, which was merged into calculator.py
# We need to ensure it's properly exposed if it's part of the mspy package API.
# The `formulator` is already imported from `mod_formulator`.

# The original isocalc.py contained functions like iso_table, draw_plot, molmass, molcharge.
# These seem to be used by calculator.py via `import isocalc_`.
# If these are intended to be part of the `mspy` package's direct API (e.g. `mspy.molmass`),
# they should be imported here. However, the current structure suggests they are used
# via the `isocalc_` alias, which might be an internal convention or a renamed import.
# For now, I will not add them to `mspy/__init__.py` unless explicitly stated they should be part of the top-level API.

# The old `blocks.py` defined element, monomer, enzyme, fragment, modification classes and default dicts.
# And also load/save functions for these.
# These are likely part of the mspy API.
from .blocks import (
    element, monomer, enzyme, fragment, modification,
    elements, monomers, enzymes, fragments, modifications, # The data dicts
    loadMonomers, loadEnzymes, loadModifications,
    saveMonomers, saveEnzymes, saveModifications
)

# Add these to __all__ as well
__all__.extend([
    "delta", "mz", "md", "nominalmass", "rdbe", "frules",
    "formulator",
    "pattern",
    "compound",
    "element", "monomer", "enzyme", "fragment", "modification",
    "elements", "monomers", "enzymes", "fragments", "modifications",
    "loadMonomers", "loadEnzymes", "loadModifications",
    "saveMonomers", "saveEnzymes", "saveModifications"
])

# Remove duplicates from __all__ just in case
__all__ = sorted(list(set(__all__)))
