---
title: 'PyOma: A Python module for Operational Modal Analysis'
tags:
  - Python
  - operational modal analysis
  - ambient vibration test
  - structural dynamics
  - modal testing
authors:
  - name: Dag Pasquale Pasca
    orcid: 0000-0002-3830-2835
    affiliation: 1
  - name: Angelo Aloisio
    orcid: 0000-0002-6190-0139
    affiliation: 2
  - name: Lorenzo De Lauretis
    orcid: 0000-0002-8796-3672
    affiliation: 3
affiliations:
 - name: Faculty of Science and Technology, Norwegian University of Life Sciences, \AA s, Norway
   index: 1
 - name: Department of Civil, Construction-Architectural and Environmental Engineering, Università degli Studi dell'Aquila, L'Aquila, Italy
   index: 1
 - name: Department of Information Engineering, Computer Science and Mathematics, Università degli Studi dell'Aquila, L'Aquila, Italy
   index: 1
date: 24 March 2021
bibliography: paper.bib
---

# Summary

Operational modal analysis (OMA) is an exciting and relatively new branch of 
modal testing. It allows to estimate the modal properties (natural frequencies, 
mode shapes and damping ratios) of a structure from ambient vibration tests. 
In these tests the structure is under its operating conditions (hence the name), 
the input excitation is neither known or controlled, but assumed to be 
'broadband random'. Ambient vibration tests have both practical and economical 
advantages over free vibration or forced vibration tests (with known input), 
especially for civil engineering structures (e.g., buildings, bridges) 
[@brincker2015introduction], [@rainieri2014operational].

# Statement of Need

By measuring and understanding the modal parameters, engineers can design 
structures, machines and devices that perform better, last longer, and are 
more comfortable for their users or occupants. There are currently some 
commercial software packages, and a Matlab toolbox, that offer OMA procedures. 
However, being proprietary software, these are closed environments that are sold 
at high prices. 

PyOMA is an open source Python module that allows to perform OMA procedures on 
ambient vibration datasets. The module offers the following identification 
algorithms:
- **Frequency Domain Decomposition** [@brincker2001modal];
- **Enhanced Frequency Domain Decomposition** [brincker2001damping];
- **Frequency Spatial Domain Decomposition** [@zhang2010frequency];
- **Covariance driven Stochastic Subspace Identification** [@peeters1999reference];
- **Data driven Stochastic Subspace Identification** [@van2012subspace]


# Acknowledgements


# References