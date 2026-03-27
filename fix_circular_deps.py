with open("mspy/mod_basics.py", "r") as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    if "from . import blocks" in line or "from . import obj_compound" in line:
        continue
    new_lines.append(line)
    if "ELECTRON_MASS = 0.00054857990924" in line:
        new_lines.append("\nfrom . import blocks\nfrom . import obj_compound\n")

with open("mspy/mod_basics.py", "w") as f:
    f.writelines(new_lines)
