import mspy

expression = "18.01"
tol = 5

formula = float(expression)
table = mspy.formulator(mz=formula, tolerance=tol)
print("SUCCESS!", len(table))
