with open("calculator.py", "r") as f:
    content = f.read()

headers_code = """
@app.after_request
def add_security_headers(response):
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['X-Frame-Options'] = 'SAMEORIGIN'
    response.headers['Strict-Transport-Security'] = 'max-age=31536000; includeSubDomains'
    return response

@app.route("/")
"""

content = content.replace('\n@app.route("/")\n', headers_code)

with open("calculator.py", "w") as f:
    f.write(content)
