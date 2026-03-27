import os

files = [
    "mspy/mod_basics.py", "mspy/mod_formulator.py", "mspy/mod_pattern.py",
    "mspy/mod_stopper.py", "mspy/obj_compound.py", "mspy/blocks.py",
    "mspy/__init__.py", "mspy/py_calculations.py"
]

for file in files:
    with open(file, "r") as f:
        content = f.read()

    # Apply fixes
    content = content.replace("import blocks", "from . import blocks")
    content = content.replace("import obj_compound", "from . import obj_compound")
    content = content.replace("import mod_basics", "from . import mod_basics")
    content = content.replace("import mod_pattern", "from . import mod_pattern")
    content = content.replace("import elements", "from . import elements")
    content = content.replace("from mod_stopper import CHECK_FORCE_QUIT", "from .mod_stopper import CHECK_FORCE_QUIT")
    content = content.replace("import obj_peaklist", "#import obj_peaklist")
    content = content.replace("import calculations", "from . import py_calculations as calculations")
    content = content.replace("import mod_signal", "#import mod_signal")
    content = content.replace("import mod_peakpicking", "#import mod_peakpicking")

    if file == "mspy/py_calculations.py":
        content = content.replace("num_elements = len(masses)", "num_elements = len(element_masses)")
        content = content.replace("masses[i]", "element_masses[i]")
        content = content.replace("masses[element_idx]", "element_masses[element_idx]")
        content = content.replace("lo_mass_target", "low_mass_target")

    with open(file, "w") as f:
        f.write(content)
