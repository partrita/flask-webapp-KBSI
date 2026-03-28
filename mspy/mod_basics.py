#!/usr/bin/env python3
# -------------------------------------------------------------------------
#     Copyright (C) 2005-2013 Martin Strohalm <www.mmass.org>
# -------------------------------------------------------------------------

import re

# regex for chemical formula
# Standard: ElementSymbol{Isotope}Count, e.g., C{13}2H6
# Supports negative counts for losses: e.g., O-1C-1H-3
# Also supports empty formulas
FORMULA_PATTERN = re.compile(r"^([A-Z][a-z]?(\{([0-9]+)\})?(-?\d+)?|\(|\)|[0-9])*$")
ELEMENT_PATTERN = re.compile(r"([A-Z][a-z]?)(\{([0-9]+)\})?(-?\d+)?")
GROUP_PATTERN = re.compile(r"\(([^()]+)\)(\d+)?")


# BASIC CALCULATIONS
# ------------------


def delta(exactMass, mz):
    """Calculate error (ppm) between mz and exact mass."""

    if exactMass == 0 or mz == 0:
        return 0.0

    return (mz - exactMass) / exactMass * 10**6


# ----


def mz(exactMass, charge, agentMass=1.00727646677, agentCharge=1):
    """Calculate m/z for given exact mass and charge."""

    if charge == 0:
        return exactMass

    m = exactMass + abs(charge / agentCharge) * agentMass
    return m / abs(charge / agentCharge)


# ----


def md(mz):
    """Calculate mass defect (m - dm) for given mz."""

    return mz - round(mz)


# ----


def nominalmass(mz):
    """Calculate nominal mass for given mz."""

    return round(mz)


# ----


def rdbe(formula):
    """Calculate RDBE value for given formula."""

    # Late import to avoid circular dependency
    from . import blocks, obj_compound

    # get composition
    cmpd = obj_compound.compound(formula)
    composition = cmpd.composition()

    # get valence sum
    rdbe_val = 0.0
    for el, count in composition.items():
        if el not in blocks.elements:
            continue
        rdbe_val += count * (blocks.elements[el].valence - 2)

    return 0.5 * rdbe_val + 1.0


# ----


def frules(formula):
    """
    Check basic formula rules:
    - RDBE > -0.5
    - H/C ratio between 0.2 and 3.0
    - NOPS/C ratio between 0.0 and 2.0
    """

    # Late import
    from . import obj_compound

    # get composition
    cmpd = obj_compound.compound(formula)
    composition = cmpd.composition()

    # check RDBE
    if rdbe(formula) < -0.5:
        return False

    # check ratios
    cCount = float(composition.get("C", 0))
    if cCount > 0:
        hCount = float(composition.get("H", 0))
        nCount = float(composition.get("N", 0))
        oCount = float(composition.get("O", 0))
        pCount = float(composition.get("P", 0))
        sCount = float(composition.get("S", 0))

        if hCount / cCount < 0.2 or hCount / cCount > 3.0:
            return False
        if (nCount + oCount + pCount + sCount) / cCount > 2.0:
            return False

    return True
