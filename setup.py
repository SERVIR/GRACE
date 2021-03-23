import os
import sys
from setuptools import setup, find_namespace_packages
from tethys_apps.app_installation import find_resource_files

### Apps Definition ###
app_package = 'grace'
release_package = 'tethysapp-' + app_package

# -- Get Resource File -- #
resource_files = find_resource_files('tethysapp/' + app_package + '/templates','tethysapp/' + app_package )
resource_files += find_resource_files('tethysapp/' + app_package + '/public','tethysapp/' + app_package )

### Python Dependencies ###
dependencies = []

setup(
    name=release_package,
    version='0.0.1',
    tags='Hydrology',
    description='View GRACE mission data',
    long_description='',
    keywords='',
    author='Sarva Pulla',
    author_email='spulla@usra.edu',
    url='',
    license='',
    packages=find_namespace_packages(),
    include_package_data=True,
    package_data={'': resource_files},
    zip_safe=False,
    install_requires=dependencies,
)
