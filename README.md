
# PyOMA
This software was created to perform output-only modal identification (Operational Modal Analysis, OMA).

OMA allows the experimental estimation of the modal parameters (natural frequencies, mode shapes, damping ratios) of a structure from measurements of the vibration response in operational condition.





## What is PyOMA?
PyOMA is a python module that allows to perform OMA on ambient vibration measurments datasets.

PyOMA include the following algorithms:

1. Frequency Domain Decomposition (FDD)

	1a. Original Frequency Domain Decomposition (FDD)
	
	2a. Enhanced Frequency Domain Decomposition (EFDD)
	
	3a. Frequency Spatial Domain Decomposition (FSDD)
	
2. Stochastic Subspace Identification (SSI)

	2a. Covariance-driven Stochastic Subspace Identification (cov-SSI)
	
	2b. Data-driven Stochastic Subspace Identification (dat-SSI)	
	




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


# Workflow

![title](Images/Flowchart.png)

FDD:

1 run FDDsvp

	2a run FDDmodEX to run original FDD
		
	and/or
		
	2b run EFDDmodEX(method='EFDD') to run EFDD
		
	and/or
		
	2c run EFDDmodEX(method='FSDD') to run FSDD

SSI

1.1 run SSIcovStaDiag 
	
	2 run SSImodEX to run cov-SSI
		
1.2 run SSIdatStaDiag 
	
	2 run SSImodEX to run dat-SSI 

# Function Description

# pyOMA.SSIcovStaDiag

**SSIcovStaDiag**(*data, fs, br, [ordmax=[None], lim=[(0.01,0.05,0.02,0.1)],
method=['1']]*)

This function returns a plot of the stabilization diagram, and the results
organised in a dictionary, for the given data calculated according to the
Covariance driven Stochastic Sub-space Identification algorithm.

**Parameters**

-   **data: array of float**

Bidimensional array Nxl where N=number of sampled data, and l=number of output
channels. Each channel is a time series of measurement values.

-   **fs: float**

Sampling frequency of the time series.

-   **br: int**

Number of block rows.

-   **ordmax: int, optional**

Specifies the maximum order of the model for the construction of the
stabilisation diagram. If *None* the maximum allowable model order is used
(ordmax = br\*l, br=number of block rows, l=number of channels). If specified,
ordmax, must be lower to the maximum allowable model order.

-   **lim: tuple of float, optional**

Limits used in the construction of the stabilisation diagrams. The poles
associated to a given model order are compared with those obtained from a
one-order lower model. Only the poles that fulfil the (first 3) stability
requirements: \|f(n)-f(n+1)\|/f(n) \< 0.01, \|x(n)- x (n+1)\|/ x (n) \< 0.05,
\|f(n)- f (n+1)\|/ f (n) \< 0.02 are labelled as stable in the plot. The last
limit (0.1 by Default) is used to remove all the poles that have a damping above
the specified limit.

-   **method: {‘1’, ’2’}, optional**

Method to estimate the *state matrix* [A]. Method ‘1’ is based on the
decomposition property of the one-lag shifted Toeplitz matrix, and is basically
equivalent to the “NExT-ERA” method. In Method ‘2’ the so-called Balanced
Realization (BR) Cov-SSI procedure is adopted: the state matrix [A] is estimated
by exploiting the shift structure of the observability matrix, and using
identity matrices as weighting matrices.

**Returns**

-   **fig1: Figure**

A matplotlib Figure with a plot of the stabilisation diagram.

-   **Results: dictionary**

Dictionary that collects the results of the identification algorithm (poles and
mode shapes) and that is passed further to the SSIModEx function.

The dictionary is organised as follows:

-   Results['Data'] = {'Data': data} : raw data

-   Results['Data']['Samp. Freq.'] = fs: sampling frequency

-   Results['All Poles'] = \_df1: Pandas DataFrame that contains the results
    (poles).

-   Results['Reduced Poles'] = df2: Pandas DataFrame that contains the results
    (poles).

-   Results['Modes'] = Ms: matrix that contain the mode shapes associated to the
    poles.

References

Rainieri, Carlo, and Giovanni Fabbrocino. "Operational modal analysis of civil
engineering structures." *Springer, New York* 142 (2014): 143.

Brincker, Rune, and Carlos Ventura. *Introduction to operational modal
analysis*. John Wiley & Sons, 2015.

Peeters, Bart. "System identification and damage detection in civil
engeneering." (2000).

Van Overschee, Peter, and B. L. De Moor. *Subspace identification for linear
systems: Theory—Implementation—Applications*. Springer Science & Business Media,
2012.

# pyOMA.SSIdatStaDiag

**SSIdatStaDiag**(*data, fs, br, [ordmax=[None], lim=[(0.01,0.05,0.02,0.1)],
method=['1']]*)

This function returns a plot of the stabilization diagram, and the results
organised in a dictionary, for the given data calculated according to the Data
driven Stochastic Sub-space Identification algorithm.

**Parameters**

-   **data: array of float**

Bidimensional array Nxl where N=number of sampled data, and l=number of output
channels. Each channel is a time series of measurement values.

-   **fs: float**

Sampling frequency of the time series.

-   **br: int**

Number of block rows.

-   **ordmax: int, optional**

Specifies the maximum order of the model for the construction of the
stabilisation diagram. If *None* the maximum allowable model order is used
(ordmax = br\*l, br=number of block rows, l=number of channels). If specified,
ordmax, must be lower to the maximum allowable model order.

-   **lim: tuple of float, optional**

Limits used in the construction of the stabilisation diagrams. The poles
associated to a given model order are compared with those obtained from a
one-order lower model. Only the poles that fulfil the (first 3) stability
requirements: \|f(n)-f(n+1)\|/f(n) \< 0.01, \|x(n)- x (n+1)\|/ x (n) \< 0.05,
\|f(n)- f (n+1)\|/ f (n) \< 0.02 are labelled as stable in the plot. The last
limit (0.1 by Default) is used to remove all the poles that have a damping above
the specified limit.

-   **method: {‘1’, ’2’}, optional**

Method to estimate the *state matrix* [A]. Method ‘1’ uses the Kalman filter
state sequences [Si] and [Si+1]. On the analogy with Cov-SSI, Method ‘2’ takes
advantage of the shift structure of the observability matrix, and uses identity
matrices as weighting (Unweighted Principal Component, UPC).

**Returns**

-   **fig1: Figure**

A matplotlib Figure with a plot of the stabilisation diagram.

-   **Results: dictionary**

Dictionary that collects the results of the identification algorithm (poles and
mode shapes) and that is needed passed further to the SSIModEx function.

The dictionary is organised as follows:

-   Results['Data'] = {'Data': data}: raw data

-   Results['Data']['Samp. Freq.'] = fs: sampling frequency

-   Results['All Poles'] = \_df1: Pandas DataFrame that contains the results
    (poles).

-   Results['Reduced Poles'] = df2: Pandas DataFrame that contains the results
    (poles).

-   Results['Modes'] = Ms: matrix that contain the mode shapes associated to the
    poles.

References

Rainieri, Carlo, and Giovanni Fabbrocino. "Operational modal analysis of civil
engineering structures." *Springer, New York* 142 (2014): 143.

Brincker, Rune, and Carlos Ventura. *Introduction to operational modal
analysis*. John Wiley & Sons, 2015.

Peeters, Bart. "System identification and damage detection in civil
engeneering." (2000).

Van Overschee, Peter, and B. L. De Moor. *Subspace identification for linear
systems: Theory—Implementation—Applications*. Springer Science & Business Media,
2012.

# pyOMA.SSIModEX

**SSIModEX***(FreQ, Results, [deltaf=[0.05], aMaClim=[0.9]]):*

This function extracts the modal properties (natural frequency, damping, mode
shape) at the frequency lines listed in FreQ.

**Parameters**

-   **FreQ: list/array of float**

Mono dimensional array that contains the frequencies that will be extracted.

-   **Results: dict**

Dictionary of results (returned by either SSIcovStaDiag or SSIdatStaDiag).

-   **deltaf: float, optional**

For each frequency line specified in FreQ the quantities “FreQ[k] + deltaf” and
“FreQ[k] – deltaf” define the limits where the routine search for the reference
modal shape.

-   **aMaClim: float, optional**

Modal Assurance Criterion limit value. The poles which have a MAC value less
than aMaClim are excluded from the calculation of the statistics of the modal
properties.

**Returns**

-   **Results: dictionary**

Dictionary that contains the results in terms of modal properties.

The dictionary is organised as follows:

-   Results['Frequencies'] = Freq: frequencies (undamped).

-   Results['Damping'] = Damp: damping.

-   Results['Mode Shapes'] = Fi.T: mode shapes.

References

Rainieri, Carlo, and Giovanni Fabbrocino. "Operational modal analysis of civil
engineering structures." *Springer, New York* 142 (2014): 143.

Brincker, Rune, and Carlos Ventura. *Introduction to operational modal
analysis*. John Wiley & Sons, 2015.

Peeters, Bart. "System identification and damage detection in civil
engeneering." (2000).

Van Overschee, Peter, and B. L. De Moor. *Subspace identification for linear
systems: Theory—Implementation—Applications*. Springer Science & Business Media,
2012.

# pyOMA.FDDsvp

**FDDsvp***(data, fs, [df=[0.01], pov=[0.5], window=['hann']]):*

This function returns the plot of the singular values of the PSD matrix, and the
results organised in a dictionary, for the given data, calculated according to
Frequency Domain Decomposition algorithm.

**Parameters**

-   **data: array of float**

Bidimensional array Nxl where N=number of sampled data, and l=number of output
channels. Each channel is a time series of measurement values.

-   **fs: float**

Sampling frequency of the time series.

-   **df: float, optional**

Desired frequency resolution of the PSD plot. Default to 0.01[Hz].

-   **pov: float, optional**

It is the (percentage) overlapping between segments. Default 0.5.

-   **window: str or tuple or array like, optional**

Desired window to use. If window is a string or tuple, it is passed to
get_window to generate the window values, which are DFT-even by default. See
get_window for a list of windows and required parameters. If window is
array_like it will be used directly as the window and its length must be
nperseg. Defaults to a hann window (descrizione copiata da scipy).

**Returns**

-   **fig1: Figure**

A matplotlib Figure with a plot of the Singular Value of the PSD matrix.

-   **Results: dictionary**

Dictionary that collects the results of the identification algorithm, and that
is passed further to the FDDModEx function and/or EFDDModEx function.

The dictionary is organised as follows:

-   Results['Data'] = {'Data': data}: raw data

-   Results['Data']['Samp. Freq.'] = fs: sampling frequency

-   Results['Data']['Freq. Resol.'] = df: frequency resolution

-   Results['Singular Values'] = S_val: matrix of the singular values of the PSD
    matrix

-   Results['Singular Vectors'] = S_vecl: matrix of the singular vectors of the
    PSD matrix

-   Results['PSD Matrix'] = PSD_matr: PSD matrix

Notes

The PSD matrix is estimated using “Welch’s method” (through scipy.signal.csd
function). The length of each segment (parameter “nperseg” of the function
scipy.signal.csd) is defined by the frequency resolution parameter “df” through
the relation “nperseg=fs/df”. The window used is passed to scipy.signal.csd
through the window parameter (default “hann”). The overlap between segments can
be passed to scipy.signal.csd through the “pov” parameter (default 0.5), which
is related to the parameter “noverlap” of the scipy.signal.csd function through
the expression “noverlap = nperseg // (1/pov)”.

References

Rainieri, Carlo, and Giovanni Fabbrocino. "Operational modal analysis of civil
engineering structures." *Springer, New York* 142 (2014): 143.

Brincker, Rune, and Carlos Ventura. *Introduction to operational modal
analysis*. John Wiley & Sons, 2015.

Peeters, Bart. "System identification and damage detection in civil
engeneering." (2000).

Van Overschee, Peter, and B. L. De Moor. *Subspace identification for linear
systems: Theory—Implementation—Applications*. Springer Science & Business Media,
2012.

# pyOMA.FDDmodEX

**FDDmodEX***(FreQ, Results, [ndf=[2]]):*

This function extracts the modal properties (natural frequency and mode shape)
at the frequency lines listed in FreQ according to the Frequency Domain
Decomposition algorithm.

**Parameters**

-   **FreQ: list/array of float**

Mono dimensional array that contains the frequencies that will be extracted.

-   **Results: dict**

Dictionary of results (returned by FDDsvp).

-   **ndf: int, optional**

Number of spectral lines near the k-th element of FreQ where the algorithm
search for the maximum difference between the first and the second singular
value.

**Returns**

-   **Results: dictionary**

Dictionary that contains the results in terms of modal properties.

The dictionary is organised as follows:

-   Results['Frequencies'] = Freq: frequencies (undamped).

-   Results['Mode Shapes'] = Fi.T: mode shapes.

References

Rainieri, Carlo, and Giovanni Fabbrocino. "Operational modal analysis of civil
engineering structures." *Springer, New York* 142 (2014): 143.

Brincker, Rune, and Carlos Ventura. *Introduction to operational modal
analysis*. John Wiley & Sons, 2015.

Peeters, Bart. "System identification and damage detection in civil
engeneering." (2000).

Van Overschee, Peter, and B. L. De Moor. *Subspace identification for linear
systems: Theory—Implementation—Applications*. Springer Science & Business Media,
2012.

# pyOMA.EFDDmodEX

**EFDDmodEX***(FreQ, Results, [ndf=[2], MAClim=[0.8], sppk=[1], npmax=[15],*
*method='FSDD']):*

This function extracts the modal properties (natural frequency, damping, mode
shape) at the frequency lines listed in FreQ according to the Enhanced Frequency
Domain Decomposition algorithm.

**Parameters**

-   **FreQ: list/array**

Mono dimensional array that contains the frequencies that will be extracted.

-   **Results: dict**

Dictionary of results (returned by FDDsvp).

-   **ndf: int, optional**

Number of spectral lines near the k-th element of FreQ where the algorithm
search for the maximum difference between the first and the second singular
value.

-   **MaClim: float, optional**

Modal Assurance Criterion limit value for the extraction of the SDOF bell
function.

-   **sppk: int, optional**

Number of peaks to skip in the SDOF auto-correlation function for the
calculations.

-   **npmax: int, optional**

Number of maximum points in the SDOF auto-correlation function to be used in the
calculations.

-   **method: {‘FSDD’, ’EFDD’}, optional**

Method used to estimate the SDOF bell function. Default to ‘FSDD’ uses the
Frequency Spatial Domain Decomposition algorithm. Method ‘EFDD’ uses the
classical Enhanced Frequency Domain Decomposition algorithm.

**Returns**

-   **fig: Figure**

len(FreQ) matplotlib figures that contain four sub-plots: the extracted SDOF
bell function, the autocorrelation function associated to the SDOF bell, the
portion of the autocorrelation function used to calculate the results, and the
results of the regression.

-   **Results: dictionary**

Dictionary that contains the results in terms of modal properties.

The dictionary is organised as follows:

-   Results['Frequencies'] = Freq: frequencies (undamped).

-   Results['Damping'] = Damp: damping.

-   Results['Mode Shapes'] = Fi.T: mode shapes.

References

Rainieri, Carlo, and Giovanni Fabbrocino. "Operational modal analysis of civil
engineering structures." *Springer, New York* 142 (2014): 143.

Brincker, Rune, and Carlos Ventura. *Introduction to operational modal
analysis*. John Wiley & Sons, 2015.

Peeters, Bart. "System identification and damage detection in civil
engeneering." (2000).

Van Overschee, Peter, and B. L. De Moor. *Subspace identification for linear
systems: Theory—Implementation—Applications*. Springer Science & Business Media,
2012.


