from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'CD45 isoform detection'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="CD45er", 
        version=VERSION,
        author="Nick Haradhvala",
        description=DESCRIPTION,
        packages=find_packages()
)
        
