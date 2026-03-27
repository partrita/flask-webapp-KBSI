from flask import Flask, request
import re
import isocalc

expression = "H2O"
pts = 500
sigma = 0.15

formula = str(expression)
charge = isocalc.molcharge(formula)
mass = isocalc.molmass(formula)
table = isocalc.iso_table(formula)
plot = isocalc.draw_plot(formula, sigma, pts)
print("SUCCESS!")
