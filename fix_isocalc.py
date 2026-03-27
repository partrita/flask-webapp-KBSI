with open("isocalc.py", "r") as f:
    content = f.read()

orig_getMass = """def getMass(x):
    atom = re.findall(r"[A-Z][a-z]*", x)
    number = re.findall(r"[0-9]+", x)
    multiplier = float(number[0]) if number else 1
    atomic_mass = float(
        matrix(PeriodicTable[atom[0]][2]) * transpose(matrix(PeriodicTable[atom[0]][3]))
    )
    return atomic_mass * multiplier"""

new_getMass = """def getMass(x):
    atom = re.findall(r"[A-Z][a-z]*", x)
    number = re.findall(r"[0-9]+", x)
    multiplier = float(number[0]) if number else 1
    atomic_mass = float(
        (matrix(PeriodicTable[atom[0]][2]) * transpose(matrix(PeriodicTable[atom[0]][3]))).item()
    )
    return atomic_mass * multiplier"""

content = content.replace(orig_getMass, new_getMass)

with open("isocalc.py", "w") as f:
    f.write(content)
