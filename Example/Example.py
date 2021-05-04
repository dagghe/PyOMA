# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 11:05:29 2021

@author: dagpa
"""
# =============================================================================
# EXAMPLE pyOMA
# =============================================================================
# Import modules
import numpy as np
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt
import pyOMA as oma


# ======== PRE-PROCESSING =====================================================
# To open a .txt file create a variable containing the path to the file
# Path to the txt file
_file = r"C:<Path to the txt file>\Ex_file.txt"

# open the file with pandas and create a dataframe 
data = pd.read_csv(_file, header=0, sep="\t", index_col=False) 
data = np.array(data)

# retrieve the example data from a function 
data, (fex, FI_ex, xi_ex) = oma.Exdata()

# Sampling frequency
fs = 100 # [Hz] Sampling Frequency
q = 5 # Decimation factor

# Using SciPy's signal module we can pre-process our data e.g. performing
# decimation, trend removal and filtering. 
# Detrend and decimate
data = signal.detrend(data, axis=0) # Trend rmoval
data = signal.decimate(data,  q, ftype='fir', axis=0) # Decimation
fs = fs/q # [Hz] Decimated sampling frequency

# Filter
_b, _a = signal.butter(12, (0.3,6.5), fs=fs, btype='bandpass')
filtdata = signal.filtfilt(_b, _a, data,axis=0)


# ======== ANALYSIS ===========================================================
# Run FDD
FDD = oma.FDDsvp(data,  fs)
# FDD = oma.FDDsvp(filtdata,  fs)

# Define list/array with the peaks identified from the plot
FreQ = [0.89, 2.6, 4.1, 5.27, 6] # identified peaks


# Extract the modal properties 
Res_FDD = oma.FDDmodEX(FreQ, FDD[1])
Res_EFDD = oma.EFDDmodEX(FreQ, FDD[1], method='EFDD', npmax = 25, MAClim=0.95, plot=False)
Res_FSDD = oma.EFDDmodEX(FreQ, FDD[1], method='FSDD', npmax = 35, MAClim=0.95, plot=True)

# Run SSI
br = 15
SSIcov= oma.SSIcovStaDiag(data, fs, br, ordmax=60)
SSIdat = oma.SSIdatStaDiag(data, fs, br, ordmax=60) 



# Extract the modal properties
Res_SSIcov = oma.SSIModEX(FreQ, SSIcov[1])
Res_SSIdat= oma.SSIModEX(FreQ, SSIdat[1])

R_SSI_auto = oma.SSIAutoModEX(SSIdat)

# =============================================================================
# 
# =============================================================================

_MS_EFDD = Res_EFDD['Mode Shapes']
_MS_FSDD = Res_FSDD['Mode Shapes']
_MS_SSIcov = Res_SSIcov['Mode Shapes']
_MS_SSIdat = Res_SSIdat['Mode Shapes']
_nch = data.shape[1]

MACmatr = np.reshape(
        [MaC(FI_ex[:,_l],_MS_FSDD[:,_k].real) for _k in range(_nch) for _l in range(_nch)], # (_nch*_nch) list of MAC values 
        (_nch,_nch)) # new (real) shape (_nch x _nch) of the MAC matrix

autoMAC = np.reshape(
        [MaC(_MS_SSIcov[:,_l].real,_MS_SSIdat[:,_k].real) for _k in range(_nch) for _l in range(_nch)], # (_nch*_nch) list of MAC values 
        (_nch,_nch)) # new (real) shape (_nch x _nch) of the MAC matrix

# PLOTTO MATRICE MAC 
meas = ["mode I", "mode II", "mode III", "mode IV", "mode V"]
num = ["mode I", "mode II", "mode III", "mode IV", "mode V"]

fig, ax = plt.subplots()
im, cbar = heatmap(MACmatr*100, num, meas, ax=ax,
                   cmap="jet", cbarlabel="MAC [%]")

texts = annotate_heatmap(im, valfmt="{x:.2f}")

fig.tight_layout()
plt.show()
