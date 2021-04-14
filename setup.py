import setuptools

with open("README.md", "r") as readme:
    project_description = readme.read()

setup(
   name='protein_compiler',
   version='1.0',
   url = "https://github.com/rcarcamo-dot/PYTHON-SBI"
   description= project_description,
   author='Patrick James Gohl, Oumout Egkemen Moustafa, Roberto Cárcamo Calvo',
   author_email='patrick.gohl01@estudiant.upf.edu, oumoutegkemen.moustafa01@estudiant.upf.edu, roberto.carcamo01@estudiant.upf',
   packages= setuptools.find_packages()
   
)
