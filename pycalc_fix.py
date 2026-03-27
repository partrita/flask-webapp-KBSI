with open("mspy/py_calculations.py", "r") as f:
    content = f.read()

content = content.replace("element_element_masses", "masses")
content = content.replace("element_masses", "masses")
content = content.replace("lo_mass_target", "low_mass_target")

with open("mspy/py_calculations.py", "w") as f:
    f.write(content)
