import re

with open("mspy/py_calculations.py", "r") as f:
    content = f.read()

content = content.replace("num_elements = len(element_masses)", "num_elements = len(element_masses)")
content = content.replace("masses", "element_masses")
content = content.replace("lo_mass_target", "low_mass_target")

with open("mspy/py_calculations.py", "w") as f:
    f.write(content)
