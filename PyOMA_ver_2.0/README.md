
# PyOMA python module

## ➤ Update version 2.0, what's new?
* Work in progress...

---

## ➤ Update version 1.5, what's new?
* `PSD_welch` function added
* `PSD_welch1` function added
* `FDDsvp` modified

	### PSD_welch brief description
	```
		This function calculate the Power Spectral Density (PSD) matrix of the 
		signals according to the Periodogram approach (Welch estimator). 
		(N.B. This function uses SciPy's "scipy.signal.csd" function. The function
		 has sometimes crashed when using big datasets)
	```

	### PSD_welch brief description
	```
		This function calculate the Power Spectral Density (PSD) matrix of the 
		signals according to the Periodogram approach (Welch estimator). 
		This function is less efficient than the previous since the calculations
		are performed "manually" and not using the more efficient "scipy.signal.csd" 
		function.
		This function was introduced since, the other version sometime crashed
		when using big dataset. 
	```

	### FDDsvp updates
	This function now perform the Frequency Domain Decomposition algorithm on the Power Spectral Density (PSD) passed as input. The function return the plot of the singular values of the Power Spectral Density (PSD) to perform the peak-peaking procedure. 

---

# How to contact us
If you have any issue, please feel free to contact us at our official e-mail address:

> [supportPyOMA@polito.it](mailto:supportPyOMA@polito.it)

# How to cite
If you use this code, please don't forget to cite this work:

> Dag Pasquale Pasca, Angelo Aloisio, Marco Martino Rosso et al., PyOMA and PyOMA_GUI: A Python module and software for Operational Modal Analysis. SoftwareX (2022) 101216, [https://doi.org/10.1016/j.softx.2022.101216](https://doi.org/10.1016/j.softx.2022.101216).