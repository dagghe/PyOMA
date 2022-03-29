import datetime

from setuptools import setup, find_packages
from pip._internal.req import parse_requirements


# Function to load requirements from file
def load_requirements(fname):
    return [
        req.requirement
        for req
        in parse_requirements(fname, session=False)
    ]


# Load requirements from .txt
requirements = load_requirements("requirements.txt")

# Run setuptools setup
setup(
    name='PyOMAQT',
    version='1.0.0',
    description='The QT sample gui',
    author='stef',
    author_email='PyOMA.com',
    url='PyOMA',
    license=f'Copyright {datetime.date.today().year} PyOMA, All rights reserved.',
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'PyOMAQT=py_oma.py_oma_gui:run', # Function to start the app
        ],
    },
    packages=find_packages(exclude=['tests']),
    package_data={
        # 'py_oma.sources': ['*.*'], # Declare folder
        'py_oma': ['*.ui', '*/*'] # Declare extra type of files, everything of py_oma
    }
)
