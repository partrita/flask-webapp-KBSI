from setuptools import setup, Extension
import numpy
import sys # Import sys to get Python include path

# Get Python include path dynamically
python_include_path = sys.prefix + "/include/python" + str(sys.version_info.major) + "." + str(sys.version_info.minor)

module = Extension(
    'mspy.calculations',  # Module name as it will be imported
    sources=['mspy/calculations.c'],
    include_dirs=[
        numpy.get_include(),
        python_include_path # Add Python's include directory
    ]
)

setup(
    name='MSPyCalculations',
    version='1.0',
    description='C extension for mspy package',
    ext_modules=[module]
)
