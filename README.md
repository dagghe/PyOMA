# Welcome to PyOMA
This project is still in beta phase.
It uses python > 3.5.

In our project we have three main functions:

 1. FDD - Frequency Domain Decomposition
 2. eFDD - Enhanced Frequency Domain Decomposition
 3. SSI-cov - Covariance Driven Stochastic Sub-space Identification
 

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

 



 


