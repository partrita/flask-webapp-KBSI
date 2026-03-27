#!/usr/bin/env python3
from flask import Flask, request, render_template
import re

import isocalc
import mspy

# create app
app = Flask(__name__)


@app.after_request
def add_security_headers(response):
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['X-Frame-Options'] = 'SAMEORIGIN'
    response.headers['Strict-Transport-Security'] = 'max-age=31536000; includeSubDomains'
    return response

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
            return "An error occurred while processing your request.", 400


@app.route("/m2f", methods=["GET", "POST"])
def m2f():
    if request.method == "GET":
        # show html form
        return render_template("m2f-submit.html")
    elif request.method == "POST":
        # calculate result
        expression = request.form.get("expression")
        try:
            tol = int(request.form.get("tol"))
            formula = float(expression)
            table = mspy.formulator(mz=formula, tolerance=tol)
            return render_template("m2f-table.html", result=table, name=formula)
        except Exception:
            return "An error occurred while processing your request.", 400


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
