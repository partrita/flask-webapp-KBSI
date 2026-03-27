import sys
import re

with open("calculator.py", "r") as f:
    content = f.read()

isoc_orig = """        expression = request.form.get("expression")
        pts = int(request.form.get("pts"))
        sigma = float(request.form.get("sigma"))

        # result = str(expression)
        formula = str(expression)
        charge = isocalc_.molcharge(formula)
        mass = isocalc_.molmass(formula)
        table = isocalc_.iso_table(formula)
        plot = isocalc_.draw_plot(formula, sigma, pts)
        return render_template(
            "isoc-table.html",
            result=table,
            mass=mass,
            name=formula,
            plot=plot,
            charge=charge,
        )"""

isoc_new = r"""        expression = request.form.get("expression")
        try:
            if not expression or not re.match(r"^[A-Za-z0-9\(\)]+$", expression):
                return "Invalid expression provided.", 400
            pts = int(request.form.get("pts"))
            sigma = float(request.form.get("sigma"))

            formula = str(expression)
            charge = isocalc.molcharge(formula)
            mass = isocalc.molmass(formula)
            table = isocalc.iso_table(formula)
            plot = isocalc.draw_plot(formula, sigma, pts)
            return render_template(
                "isoc-table.html",
                result=table,
                mass=mass,
                name=formula,
                plot=plot,
                charge=charge,
            )
        except Exception:
            return "An error occurred while processing your request.", 400"""

content = content.replace(isoc_orig, isoc_new)
content = content.replace("from flask import Flask, request, render_template\n", "from flask import Flask, request, render_template\nimport re\n")

with open("calculator.py", "w") as f:
    f.write(content)
