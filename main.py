#!/usr/bin/env python3
from flask import Flask, request, render_template
import re

import isocalc
import mspy

# create app
app = Flask(__name__)


@app.after_request
def add_security_headers(response):
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    response.headers["Strict-Transport-Security"] = (
        "max-age=31536000; includeSubDomains"
    )
    return response


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/about")
def about():
    return render_template("about.html")


@app.route("/isoc", methods=["GET", "POST"])
def isoc():
    if request.method == "GET":
        # show html form
        return render_template("isoc-submit.html")
    elif request.method == "POST":
        # calculate result
        expression = request.form.get("expression")
        try:
            if not expression or not re.match(
                r"^[A-Za-z0-9\(\)\{\}\[\]\-\+]+$", expression
            ):
                return f"Invalid expression provided: {expression}", 400

            pts_val = request.form.get("pts")
            pts = int(pts_val) if pts_val else 500

            sigma_val = request.form.get("sigma")
            sigma = float(sigma_val) if sigma_val else 0.15

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
    import os

    debug_mode = os.environ.get("FLASK_DEBUG", "false").lower() == "true"
    port = int(os.environ.get("PORT", 8080))
    app.run(debug=debug_mode, host="0.0.0.0", port=port)
