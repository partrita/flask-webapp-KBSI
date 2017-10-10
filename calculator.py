#!/usr/bin/env python
from flask import Flask, request, render_template
#from flask.ext.script import Manager

import isocalc
import mspy

# create app
app = Flask(__name__)

#manager = Manager(app)
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/tools')
def tools():
    return render_template('tools.html')

@app.route('/isoc', methods=['GET', 'POST'])
def isoc():
    if request.method == 'GET':
        # show html form
        return render_template('isoc-submit.html')
    elif request.method == 'POST':
        #calculate result
        expression = request.form.get('expression')
        pts = int(request.form.get('pts'))
        sigma = float(request.form.get('sigma'))

        #result = str(expression)
        formula = str(expression)
        charge = isocalc.molcharge(formula)
        mass = isocalc.molmass(formula)
        table = isocalc.iso_table(formula)
        plot = isocalc.draw_plot(formula,sigma,pts)
        return render_template('isoc-table.html',result=table,mass=mass,name=formula, plot = plot, charge   = charge)

@app.route('/m2f', methods=['GET', 'POST'])
def m2f():
    if request.method == 'GET':
        # show html form
        return render_template('m2f-submit.html')
    elif request.method == 'POST':
        #calculate result
        expression = request.form.get('expression')
        tol = int(request.form.get('tol'))
        formula = float(expression)
        table = mspy.formulator(formula)
        return render_template('m2f-table.html',result=table, name=formula)



# run app
if __name__ == '__main__':
    app.run(debug=True)
#debug option must be off when deployed
