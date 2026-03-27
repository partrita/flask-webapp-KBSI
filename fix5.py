import sys

with open("calculator.py", "r") as f:
    content = f.read()

m2f_orig = """        expression = request.form.get("expression")
        tol = int(request.form.get("tol"))
        formula = float(expression)
        table = mspy.formulator(mz=formula, tolerance=tol)
        return render_template("m2f-table.html", result=table, name=formula)"""

m2f_new = """        expression = request.form.get("expression")
        try:
            tol = int(request.form.get("tol"))
            formula = float(expression)
            table = mspy.formulator(mz=formula, tolerance=tol)
            return render_template("m2f-table.html", result=table, name=formula)
        except Exception:
            return "An error occurred while processing your request.", 400"""

content = content.replace(m2f_orig, m2f_new)

with open("calculator.py", "w") as f:
    f.write(content)
