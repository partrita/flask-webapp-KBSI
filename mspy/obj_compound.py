#!/usr/bin/env python3
# -------------------------------------------------------------------------
#     Copyright (C) 2005-2013 Martin Strohalm <www.mmass.org>
# -------------------------------------------------------------------------

# load modules
from . import mod_basics


# COMPOUND OBJECT DEFINITION
# --------------------------


class compound:
    """Compound object definition."""

    def __init__(self, expression, **attr):
        # check formula
        self._checkFormula(expression)
        self.expression = expression

        # buffers
        self._composition = None
        self._formula = None
        self._mass = None
        self._nominalmass = None

        # get additional attributes
        self.attributes = {}
        for name, value in attr.items():
            self.attributes[name] = value

    def count(self, item, groupIsotopes=True):
        """Count given atom or group in the compound."""

        # Late import to avoid circular dependency
        from . import blocks

        # get composition
        comp = self.composition()

        # get atoms to count
        atoms = [item]
        if groupIsotopes and item in blocks.elements:
            for massNo in blocks.elements[item].isotopes:
                atom = "%s{%d}" % (item, massNo)
                atoms.append(atom)

        # count atom
        atomsCount = 0
        for atom in atoms:
            if atom in comp:
                atomsCount += comp[atom]

        return atomsCount

    def composition(self, formula=None):
        """Return the molecular composition of the compound."""

        if formula is None:
            if self._composition is not None:
                return self._composition
            formula = self.expression

        # calculate composition
        composition = {}
        while "(" in formula:
            match = mod_basics.GROUP_PATTERN.search(formula)
            if not match:
                break
            group, count = match.groups()
            count = int(count) if count else 1

            groupComposition = self.composition(group)
            for atom, atomCount in groupComposition.items():
                composition[atom] = composition.get(atom, 0) + atomCount * count

            formula = formula.replace(match.group(0), "", 1)

        for symbol, iso_full, iso_val, count in mod_basics.ELEMENT_PATTERN.findall(
            formula
        ):
            atom = symbol
            if iso_val:
                atom = "%s{%s}" % (symbol, iso_val)
            count = int(count) if count else 1
            composition[atom] = composition.get(atom, 0) + count

        if formula == self.expression:
            self._composition = composition

        return composition

    def formula(self):
        """Return the formatted chemical formula."""

        if self._formula is not None:
            return self._formula

        # get composition
        composition = self.composition()

        # format formula
        formula = ""
        elements = sorted(composition.keys())
        if "C" in elements:
            elements.remove("C")
            elements.insert(0, "C")
        if "C" in elements and "H" in elements:
            elements.remove("H")
            elements.insert(1, "H")

        for el in elements:
            count = composition[el]
            if count == 1:
                formula += el
            elif count != 0:
                formula += "%s%d" % (el, count)

        self._formula = formula
        return formula

    def mass(self, type=0):
        """Calculate the exact mass of the compound."""

        # Late import
        from . import blocks

        if self._mass is not None:
            return self._mass[type]

        # calculate mass
        massMo = 0.0
        massAv = 0.0
        composition = self.composition()

        for atom, count in composition.items():
            # check specified isotope and mass
            match = mod_basics.ELEMENT_PATTERN.match(atom)
            if not match:
                continue
            symbol, iso_full, iso_val, count_str = match.groups()
            if iso_val:
                isotope = blocks.elements[symbol].isotopes[int(iso_val)]
                atomMass = (isotope[0], isotope[0])
            else:
                atomMass = blocks.elements[symbol].mass

            # multiply atom mass
            massMo += atomMass[0] * count
            massAv += atomMass[1] * count

        # store mass in buffer
        self._mass = (massMo, massAv)

        return self._mass[type]

    def nominalmass(self):
        """Calculate the nominal mass of the compound."""

        # Late import
        from . import blocks

        if self._nominalmass is not None:
            return self._nominalmass

        # calculate nominal mass
        nominalmass = 0
        composition = self.composition()

        for atom, count in composition.items():
            # check specified isotope and mass
            match = mod_basics.ELEMENT_PATTERN.match(atom)
            if not match:
                continue
            symbol, iso_full, iso_val, count_str = match.groups()
            if iso_val:
                isotope = blocks.elements[symbol].isotopes[int(iso_val)]
                atomMass = isotope[0]
            else:
                atomMass = blocks.elements[symbol].mass[0]

            # multiply atom mass
            nominalmass += round(atomMass) * count

        # store mass in buffer
        self._nominalmass = nominalmass

        return nominalmass

    def _checkFormula(self, formula):
        """Check if given formula is valid."""

        # Late import
        from . import blocks

        if not mod_basics.FORMULA_PATTERN.match(formula):
            raise ValueError("Wrong formula! --> " + formula)

        # check elements and isotopes
        for symbol, iso_full, iso_val, count in mod_basics.ELEMENT_PATTERN.findall(
            formula
        ):
            if symbol not in blocks.elements:
                raise ValueError(
                    "Unknown element in formula! --> " + symbol + " in " + formula
                )
            elif iso_val and int(iso_val) not in blocks.elements[symbol].isotopes:
                raise ValueError(
                    "Unknown isotope in formula! --> "
                    + symbol
                    + "{"
                    + iso_val
                    + "}"
                    + " in "
                    + formula
                )

    def __str__(self):
        return self.formula()

    def __repr__(self):
        return "<mspy.compound: %s>" % self.formula()
