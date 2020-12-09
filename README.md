
# PyOMA
This software was created for ................. and various images





## What is PyOMA?
PyOMA is a powerful software that does ....... bla bla

PyOMA includes complex algorithms that allows to calculate the following: 

 1. FDD - Frequency Domain Decomposition
 2. eFDD - Enhanced Frequency Domain Decomposition
 3. SSI-cov - Covariance Driven Stochastic Sub-space Identification





## Installing PyOMA
As a prerequisite to install PyOMA, you need to install [Anaconda](https://docs.anaconda.com/anaconda/install/) first.
You should install a Python version greather equal 3.5 or the software may run in troubles.

To install PyOMA, you need to run the following commands:

 - pip install pyOMA
 - conda install -c conda-forge pyOMA
 
 ### Dependencies
 - numpy (https://numpy.org/)
 - pandas (https://pandas.pydata.org/)
 - scipy -> signal (https://www.scipy.org/)
 - scipy.optimize -> curve_fit (https://www.scipy.org/)
 - matplotlib.pyplot (https://matplotlib.org/)
 - matplotlib.ticker -> [MultipleLocator,FormatStrFormatter] (https://matplotlib.org/)
 - seaborn (https://seaborn.pydata.org/)
 - mplcursors (https://mplcursors.readthedocs.io/en/stable/)

 

# How to use
We parametrized each one of the three main functions that we offer in our project.
## FDD
Our FDD functions has a number of parameters and two stages: in the first stage, it computes the , in the second stage, …....

Parameters for FDD in stage 1:
 - stage: is the stage (int)
 - Fs: is the sampling frequency (int)
 - q: is the decimation factor (to avoid decimation, simply put it to 1) (int)
 - nSec: is the number of segments (int)

How to launch:
python3 PyOMA.py fileName fdd stage Fs q nSec

Example:
python3 PyOMA.py testData.txt fdd 1 1200 2 30

Parameters for FDD in stage 2:
 - stage: is the stage (int)
 - Fs: is the sampling frequency (int)
 - q: is the decimation factor (to avoid decimation, simply put it to 1) (int)
 - nSec: is the number of segments (int)
 - **x**: is the spectral line (float)
 
How to launch:
python3 PyOMA.py fileName fdd stage Fs q nSec x

Example:
python3 PyOMA.py testData.txt fdd 2 1200 2 30 62.2
