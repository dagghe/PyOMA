from distutils.core import setup

setup(
    name = 'Py-OMA',
    packages = ['PyOMA'],
    install_requires=['numpy','scipy','pandas','matplotlib','seaborn','mplcursors'],
    version = 'v1.1',  # Ideally should be same as your GitHub release tag varsion
    description = 'PyOMA allows the experimental estimation of the modal parameters (natural frequencies, mode shapes, damping ratios) of a structure from measurements of the vibration response in operational condition.',
    author = 'Dag Pasquale Pasca, Angelo Aloisio, Lorenzo De Lauretis',
    author_email = '',
    url = 'https://github.com/dagghe/PyOMA/tree/master/source',
    download_url = 'https://github.com/dagghe/PyOMA/archive/refs/tags/v1.1.tar.gz',
    keywords = ['OMA', 'structure'],
    classifiers = [],
)