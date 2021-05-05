---
title: 'PyOMA: A Python module for Operational Modal Analysis'
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
The practice of Operational Modal Analysis (OMA) of civil structures and
infrastructures has been growing significantly in the last decades 
[@rainieri2014operational]. OMA allows to estimate the modal properties 
(natural frequencies, mode shapes and damping ratios) of a structure from ambient 
vibration tests, that is while the structure is in its operating conditions. 
It is nearly impossible to measure the input ambient excitations for a civil
engineering structure. Neither it is possible to adequately excite them using 
conventional devices, such as impact hammer and shakers. On the other hand,
output vibration data acquired through accelerometers are relatively easy to get.
In OMA, the deterministic knowledge of the input excitation is replaced by the 
assumption that the input is a realization of astochastic process. Accordingly, 
Operational Modal Analysis and Output-Only dynamic identification are considered synonyms. 
There is a large variety of algorithms for OMA. Among them, the Stochastic Subspace 
Identification (SSI) and the Frequency Domain Decomposition (FDD) have proved to be 
effective and reliable algorithms for Output-Only dynamic identification 
[@rainieri2014operational]. The SSI is a time-domain method that estimates the 
stochastic state-space model from Output-Only data [@peeters1999reference]. The FDD 
is a frequency domain method that estimates the modal parameters using the input/output 
relationship of an n-DOF system stochastic process [@brincker2001modal]. There are
commercial software that offer the two algorithms, the most known probably are
ARTeMIS, by Strctural Vibration Solutions, and MACEC, a Matlab toolbox for modal 
testing and OMA [@reynders2014macec]. However, there is no Python module for OMA.


# Statement of Need
PyOMA is an open-source Python module for OMA. It can estimate the modal parameters
of a structure using five acknowledged techniques derived from SSI and FDD:
- **Frequency Domain Decomposition** [@brincker2001modal];
- **Enhanced Frequency Domain Decomposition** [brincker2001damping];
- **Frequency Spatial Domain Decomposition** [@zhang2010frequency];
- **Covariance driven Stochastic Subspace Identification** [@peeters1999reference];
- **Data driven Stochastic Subspace Identification** [@van2012subspace]

The API for PyOMA provides a set of functions for a quick and simple estimation of the 
natural frequencies, mode shapes and damping using the experimental data recorded 
by the user. Specifically, the user, in addition to the measurments data, need to
specify only a minimal amount of input parameters. For a full description of the 
functionalities please refer to the [documentation page](https://github.com/dagghe/PyOMA/wiki/Function-Description). 
The flowchart in Fig.1 shows the general architecture of PyOMA. PyOMA is designed to be 
used by both researchers and practitioners. It has already been used in several applications,
as proved by several scientific publications: [@alaggio2021two], [@aloisio2020assessment], 
[@aloisio2021assessment], [@aloisio2020bayesian], [@aloisio2020identification], 
[@aloisio2020recorded], [@aloisio2020dynamic],[@aloisio2020time], [@capanna2021operational]

![`PyOMA` flowchart.](Fig1.png)

# Acknowledgements
The authors acknowledge the meaningful contribution of Professor Rocco Alaggio, who
encouraged the authors to study and develop these topics.

# References
