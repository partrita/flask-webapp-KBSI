#!/usr/bin/env python3
from flask import Flask, request, render_template

import isocalc
import mspy

# create app
app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/about")
def about():
    return render_template("about.html")


@app.route("/tools")
def tools():
    return render_template("tools.html")


@app.route("/isoc", methods=["GET", "POST"])
def isoc():
    if request.method == "GET":
        # show html form
        return render_template("isoc-submit.html")
    elif request.method == "POST":
        # calculate result
        expression = request.form.get("expression")
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
        )


@app.route("/m2f", methods=["GET", "POST"])
def m2f():
    if request.method == "GET":
        # show html form
        return render_template("m2f-submit.html")
    elif request.method == "POST":
        # calculate result
        expression = request.form.get("expression")
        tol = int(request.form.get("tol"))
        formula = float(expression)
        table = mspy.formulator(mz=formula, tolerance=tol)
        return render_template("m2f-table.html", result=table, name=formula)


# run app
if __name__ == "__main__":
    # Debug mode should be controlled by FLASK_DEBUG environment variable
    # app.run(debug=True) # Old way
    # For modern Flask, 'flask run' command handles this.
    # The app.run() call is often not needed for development if using 'flask run'.
    # However, to keep it runnable directly with 'python calculator.py',
    # we can check an environment variable.
    import os
    debug_mode = os.environ.get('FLASK_DEBUG', 'false').lower() == 'true'
    app.run(debug=debug_mode, host='0.0.0.0')
