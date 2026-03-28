#!/usr/bin/env python3
import re
import sys
import operator
from itertools import groupby
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")
import mpld3
from mspy.periodic_table import PERIODIC_TABLE

# Configuration constants
DEFAULT_SIGMA = 0.15
DEFAULT_PTS = 500
CUTOFF = 0.0001


def get_mass(segment):
    """Calculates the atomic mass for a formula segment (e.g., 'H2')."""
    atom = re.findall(r"[A-Z][a-z]*", segment)[0]
    number = re.findall(r"[0-9]+", segment)
    multiplier = float(number[0]) if number else 1.0

    masses = np.array(PERIODIC_TABLE[atom][2])
    abundances = np.array(PERIODIC_TABLE[atom][3])
    atomic_mass = float(np.dot(masses, abundances))

    return atomic_mass * multiplier


def get_charge(segment):
    """Calculates the default charge for a formula segment."""
    atom = re.findall(r"[A-Z][a-z]*", segment)[0]
    number = re.findall(r"[0-9]+", segment)
    multiplier = float(number[0]) if number else 1.0
    atomic_charge = float(PERIODIC_TABLE[atom][1])
    return atomic_charge * multiplier


def get_isotope(segment, data_index):
    """Retrieves isotope data (mass or ratio) for a formula segment."""
    atom = re.findall(r"[A-Z][a-z]*", segment)[0]
    number = re.findall(r"[0-9]+", segment)
    multiplier = int(number[0]) if number else 1
    isotope_data = PERIODIC_TABLE[atom][data_index]
    return [isotope_data] * multiplier


def molmass(formula):
    """Calculates total molecular mass."""
    segments = re.findall("[A-Z][a-z]*[0-9]*", formula)
    return sum(get_mass(s) for s in segments)


def molcharge(formula):
    """Calculates total molecular charge based on default oxidation states."""
    segments = re.findall("[A-Z][a-z]*[0-9]*", formula)
    return sum(get_charge(s) for s in segments)


def isotoperatios(formula):
    """Returns a list of isotope abundance ratios for each atom in the formula."""
    segments = re.findall("[A-Z][a-z]*[0-9]*", formula)
    ratios = []
    for s in segments:
        ratios += get_isotope(s, 3)
    return ratios


def isotopemasses(formula):
    """Returns a list of isotope masses for each atom in the formula."""
    segments = re.findall("[A-Z][a-z]*[0-9]*", formula)
    masses = []
    for s in segments:
        masses += get_isotope(s, 2)
    return masses


def formulaExpander(formula):
    """Expands formulas with parentheses, e.g., 'Al2(NO3)4'."""
    while "(" in formula:
        matches = re.finditer(r"\(([^()]*)\)([0-9]+)", formula)
        for match in matches:
            full_match = match.group(0)
            inner_content = match.group(1)
            multiplier = int(match.group(2))

            inner_segments = re.findall("[A-Z][a-z]*[0-9]*", inner_content)
            expanded = ""
            for s in inner_segments:
                atom = re.findall(r"[A-Z][a-z]*", s)[0]
                num_match = re.findall(r"[0-9]+", s)
                count = int(num_match[0]) if num_match else 1
                expanded += f"{atom}{count * multiplier}"

            formula = formula.replace(full_match, expanded)

        # Clean up empty parentheses
        formula = formula.replace("()", "")
        if "(" in formula and not re.search(r"\(([^()]*)\)([0-9]+)", formula):
            formula = formula.replace("(", "").replace(")", "")

    return formula


def trim(ratios, masses):
    """Combines duplicate masses and trims low abundance signals."""
    if len(ratios) > 10:
        pairs = sorted(zip(masses, ratios))
        new_masses = []
        new_ratios = []
        for key, group in groupby(pairs, operator.itemgetter(0)):
            new_masses.append(key)
            new_ratios.append(sum(item[1] for item in group))
        return new_ratios, new_masses
    return ratios, masses


def slowcartesian(rx, mx, ry, my, idx):
    """Recursively calculates the cartesian product of isotope distributions."""
    xp = []
    xs = []
    max_ry = max(ry) if ry else 1.0
    ry_norm = [n / max_ry for n in ry]

    for k in range(len(rx[idx])):
        rk = rx[idx][k]
        mk = mx[idx][k]
        for j in range(len(ry_norm)):
            rj = ry_norm[j]
            mj = my[j]
            comp = rj * rk
            if comp > CUTOFF:
                xp.append(comp)
                xs.append(mj + mk)

    xp, xs = trim(xp, xs)
    idx += 1
    if idx < len(rx) and len(xp) < 1000000:
        return slowcartesian(rx, mx, xp, xs, idx)
    return xp, xs


def isotopes(ratios, masses):
    """Computes the full isotopic distribution."""
    xs, xp = slowcartesian(ratios, masses, ratios[0], masses[0], 1)
    return [round(n, 8) for n in xs], [round(n, 8) for n in xp]


def genDict(masses, ratios, charges):
    """Generates a dictionary of m/z vs relative intensity."""
    if not masses:
        return {}

    pairs = sorted(zip(masses, ratios))
    combined = {}
    for m, r in pairs:
        mz = m / abs(charges) if charges != 0 else m
        combined[mz] = combined.get(mz, 0) + r

    max_val = max(combined.values()) if combined else 1.0
    final = {
        k: (v * 100 / max_val)
        for k, v in combined.items()
        if v * 100 / max_val > CUTOFF
    }
    return final


def gaussian(x, c, s):
    """Gaussian function."""
    return (1.0 / (s * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - c) / s) ** 2)


def genGaussian(final, sigma, pts):
    """Generates gaussian-broadened spectrum data."""
    if not final:
        return [], []

    x = sorted(final.keys())
    y = [final[k] for k in x]

    x_vec = np.linspace(min(x) - 1, max(x) + 1, pts)
    y_vec = np.zeros(pts)

    for i, xv in enumerate(x_vec):
        val = sum(y[j] * gaussian(x[j], xv, sigma) for j in range(len(x)))
        y_vec[i] = val

    max_y = np.max(y_vec) if len(y_vec) > 0 else 1.0
    y_vec = (y_vec / max_y) * 100 if max_y > 0 else y_vec

    return x_vec, y_vec


def iso_table(formula):
    """Main entry point for calculating isotope table."""
    formula = formulaExpander(formula)
    masses = isotopemasses(formula)
    ratios = isotoperatios(formula)
    charge = 1

    if len(masses) > 1:
        xs, xp = isotopes(ratios, masses)
        final = genDict(xs, xp, charge)
    elif masses:
        final = genDict(masses[0], ratios[0], charge)
    else:
        final = {}

    return {k: v for k, v in sorted(final.items()) if v > CUTOFF}


def draw_plot(formula, sigma=DEFAULT_SIGMA, pts=DEFAULT_PTS):
    """Main entry point for generating isotope plot HTML."""
    table = iso_table(formula)
    x_vec, y_vec = genGaussian(table, sigma, pts)

    fig, ax = plt.subplots()
    ax.plot(x_vec, y_vec)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative Intensity (%)")
    ax.set_title(f"Isotopic Distribution for {formula}")

    html = mpld3.fig_to_html(fig)
    plt.close(fig)
    return html


if __name__ == "__main__":
    if len(sys.argv) > 1:
        formula = sys.argv[1]
        print(f"Calculating for: {formula}")
        print(f"Molecular Mass: {molmass(formula)}")
        print(f"Isotope Table: {iso_table(formula)}")
    else:
        print("Usage: ./isocalc.py <formula>")
