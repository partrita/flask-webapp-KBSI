#!/usr/bin/env python3
# -------------------------------------------------------------------------
#     Copyright (C) 2005-2013 Martin Strohalm <www.mmass.org>
# -------------------------------------------------------------------------

# load libs
import numpy

# load stopper
from .mod_stopper import CHECK_FORCE_QUIT

# load blocks
from . import blocks

# load objects
from . import obj_compound

# load modules
from . import py_calculations as calculations


# ISOTOPIC PATTERN FUNCTIONS
# --------------------------


def pattern(
    compound,
    fwhm=0.1,
    threshold=0.01,
    charge=0,
    agentFormula="H",
    agentCharge=1,
    agentMass=1.00727646677,
    real=False,
    model="gaussian",
):
    """Calculate isotopic pattern for given compound.
    compound (mspy.compound, str) - compound or formula to use
    fwhm (float) - peak width
    threshold (float) - signals threshold
    charge (int) - charge of the compound
    agentFormula (str) - ionization agent formula
    agentCharge (int) - charge of the ionization agent
    agentMass (float) - exact mass of the ionization agent
    real (bool) - use centroiding on profile to get real peaks
    model (gaussian, lorentzian, gausslorentzian) - peak shape function
    """

    # check compound type
    if not isinstance(compound, obj_compound.compound):
        compound = obj_compound.compound(compound)

    # get composition
    composition = compound.composition()

    # get isotopic patterns for each element
    patterns = []
    for element, count in composition.items():
        if element not in blocks.elements:
            continue
        p = blocks.elements[element].pattern(count, threshold)
        patterns.append(p)

    # get ionization agent m/z shift
    mzShift = 0.0
    if charge != 0:
        agent = obj_compound.compound(agentFormula)
        mzShift = (
            abs(charge / agentCharge) * agentMass
            + (agent.mass() - agentMass) * abs(charge / agentCharge)
        ) / abs(charge / agentCharge)

    # combine isotopic patterns
    finalPattern = []
    if patterns:
        finalPattern = patterns[0]
        for i in range(1, len(patterns)):
            finalPattern = _combine(finalPattern, patterns[i], threshold)
            CHECK_FORCE_QUIT()

    # apply charge and shift
    if charge != 0:
        for i in range(len(finalPattern)):
            finalPattern[i][0] = (finalPattern[i][0] / abs(charge)) + mzShift

    # get real peaks from profile
    if real:
        prof = profile(finalPattern, fwhm=fwhm, points=100, model=model)
        finalPattern = []
        for isotope in calculations.signal_local_maxima(prof):
            finalPattern.append(isotope)
            centroid = calculations.signal_centroid(prof, isotope[0], isotope[1] * 0.99)
            if abs(centroid - isotope[0]) < fwhm / 100.0:
                finalPattern[-1][0] = centroid

    # normalize pattern
    finalPattern = _normalize(finalPattern)

    return finalPattern


# ----


def profile(
    pattern,
    fwhm=0.1,
    threshold=0.01,
    points=500,
    noise=0.0,
    forceFwhm=False,
    model="gaussian",
    raster=None,
):
    """Make profile spectrum for given isotopic pattern.
    pattern (list of [mz,intens]) - isotopic pattern
    fwhm (float) - peak width
    threshold (float) - signals threshold
    points (int) - number of data points
    noise (float) - noise level
    forceFwhm (bool) - use default fwhm for all peaks
    model (gaussian, lorentzian, gausslorentzian) - peak shape function
    """

    # check raster type
    if raster is not None and not isinstance(raster, numpy.ndarray):
        raster = numpy.array(raster)

    # get peaks
    peaks = []
    for isotope in pattern:
        if isotope[1] > threshold:
            peaks.append([isotope[0], isotope[1], fwhm])

    # get model shape
    shape = 0
    if model == "gaussian":
        shape = 0
    elif model == "lorentzian":
        shape = 1
    elif model == "gausslorentzian":
        shape = 2

    # make profile
    if raster is not None:
        data = calculations.signal_profile_to_raster(
            numpy.array(peaks), raster, float(noise), shape
        )
    else:
        data = calculations.signal_profile(
            numpy.array(peaks), int(points), float(noise), shape
        )

    return data


# ----


def matchpattern(signal, pattern, pickingHeight=0.75, baseline=None):
    """Compare signal with given isotopic pattern.
    signal (numpy array) - signal data points
    pattern (list of [mz,intens]) - theoretical pattern to compare
    pickingHeight (float) - centroiding height
    baseline (numpy array) - signal baseline
    """

    # check signal type
    if not isinstance(signal, numpy.ndarray):
        signal = numpy.array(signal)

    # check pattern type
    if not isinstance(pattern, list):
        pattern = list(pattern)

    # normalize pattern
    pattern = _normalize(pattern)

    # label peaks in signal for each pattern isotope
    peaklist = []
    for isotope in pattern:
        # Simplified: find closest index in signal and get intensity
        idx = (numpy.abs(signal[:, 0] - isotope[0])).argmin()
        peak_mz = signal[idx, 0]
        peak_intensity = signal[idx, 1]
        peaklist.append([peak_mz, peak_intensity])

    # normalize signal peaks
    maxIntensity = 0.0
    for peak in peaklist:
        if peak[1] > maxIntensity:
            maxIntensity = peak[1]

    if maxIntensity > 0:
        for i in range(len(peaklist)):
            peaklist[i][1] = peaklist[i][1] / maxIntensity * 100.0

    # calculate match
    match = 0.0
    for i in range(len(pattern)):
        match += abs(pattern[i][1] - peaklist[i][1])

    return match / len(pattern)


# ----


def _combine(pattern1, pattern2, threshold=0.01):
    """Combine two isotopic patterns."""

    # combine patterns
    combined = {}
    for mz1, intens1 in pattern1:
        for mz2, intens2 in pattern2:
            mz = mz1 + mz2
            intens = intens1 * intens2
            if intens > threshold / 1000.0:
                mz_key = round(mz, 6)
                combined[mz_key] = combined.get(mz_key, 0.0) + intens

    # get final pattern
    final = []
    for mz in sorted(combined.keys()):
        final.append([mz, combined[mz]])

    return _normalize(final)


# ----


def _normalize(pattern):
    """Normalize isotopic pattern to 100%."""

    # get max intensity
    maxIntensity = 0.0
    for mz, intens in pattern:
        if intens > maxIntensity:
            maxIntensity = intens

    # normalize
    if maxIntensity > 0:
        for i in range(len(pattern)):
            pattern[i][1] = pattern[i][1] / maxIntensity * 100.0

    return pattern
