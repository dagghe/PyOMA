---
title: 'PyOMA and PyOMA_GUI: A Python module and software for Operational Modal Analysis'
tags:
  - Python
  - operational modal analysis
  - ambient vibration test
  - structural dynamics
  - modal testing
  - graphical user interface
authors:
  - name: Dag Pasquale Pasca
    orcid: 0000-0002-3830-2835
    affiliation: 1
  - name: Angelo Aloisio
    orcid: 0000-0002-6190-0139
    affiliation: 2
  - name: Marco Martino Rosso
    orcid: 0000-0002-9098-4132
    affiliation: 3
  - name: Stefanos Sotiropoulos
    orcid: 0000-0001-5669-1320
    affiliation: 3
affiliations:
 - name: Norsk Treteknisk Institutt, Børrestuveien 3, 0373, Oslo, Norway
   index: 1
 - name: Department of Civil, Construction-Architectural and Environmental Engineering, Università degli Studi dell'Aquila, L'Aquila, Italy
   index: 2
 - name: Department of Structural, Geotechnical and Building Engineering Politecnico di Torino, Corso Duca degli Abruzzi, 24, 10129, Torino, Italy
   index: 3
date: 18 May 2022
bibliography: paper.bib
---

# Summary
The practice of Operational Modal Analysis (OMA) of civil structures and
infrastructures has been growing significantly in the last decades 
[@rainieri2014operational]. OMA allows one to estimate the modal properties 
(natural frequencies, mode shapes and damping ratios) of a structure from ambient 
vibration tests, while the structure is in its operating conditions. 
It is nearly impossible to measure the input ambient excitations for a civil
engineering structure. Neither is it possible to adequately excite them using 
conventional devices, such as impact hammers and shakers. On the other hand,
output vibration data acquired through accelerometers are relatively easy to get.
In OMA, the deterministic knowledge of the input excitation is replaced by the 
assumption that the input is a realization of a stochastic process. Accordingly, 
OMA and Output-Only dynamic identification are considered synonyms. 
There is a large variety of algorithms for OMA. Among them, the Stochastic Subspace 
Identification (SSI) and the Frequency Domain Decomposition (FDD) have proved to be 
effective and reliable algorithms for Output-Only dynamic identification 
[@rainieri2014operational]. The SSI is a time-domain method that estimates the 
stochastic state-space model from Output-Only data [@peeters1999reference]. The FDD 
is a frequency domain method that estimates the modal parameters using the input/output 
relationship of an n-DOF system stochastic process [@brincker2001modal]. 
So in conclusion, pyOMA is an open source Python module that allows the estimation 
of the modal parameters of a structure using five acknowledged techniques derived from 
SSI and FDD:

- **Frequency Domain Decomposition** [@brincker2001modal];

- **Enhanced Frequency Domain Decomposition** [@brincker2001damping];

- **Frequency Spatial Domain Decomposition** [@zhang2010frequency];

- **Covariance driven Stochastic Subspace Identification** [@peeters1999reference];

- **Data driven Stochastic Subspace Identification** [@van2012subspace]


# Statement of Need
To this date, there are a few commercial software that implement the two aforementioned 
algorithms, the most known probably are ARTeMIS, by Structural Vibration Solutions, and 
MACEC, a Matlab toolbox for modal testing and OMA [@reynders2014macec]. However, to the 
authors best knowledge, there is no Python module, nor any other open source toolbox to 
perform OMA, for this reason we have developed pyOMA.
The API for PyOMA provides a set of functions for a quick and simple estimation of the 
natural frequencies, mode shapes and damping using the experimental data recorded 
by the user. Specifically, the user, in addition to the measurements data, needs to
specify only a minimal amount of input parameters. For a full description of the 
functionalities please refer to the [documentation page](https://github.com/dagghe/PyOMA/wiki/Function-Description). 
The flowchart in Fig.1 shows the general architecture of PyOMA. PyOMA is designed to be 
used by both researchers and practitioners. It has already been used in several applications,
as proved by several scientific publications: @alaggio2021two, @aloisio2020assessment, 
@aloisio2021assessment, @aloisio2020bayesian, @aloisio2020identification, 
@aloisio2020recorded, @aloisio2020dynamic,@aloisio2020time, @capanna2021operational. 

![`PyOMA` flowchart.](Fig1.png)

# Software overview
PyOMA_GUI is a graphical user interface software developed in [PyQt5](https://pypi.org/project/PyQt5/), which implements in a single integrated tool the operational modal analysis of civil structures with output-only measurement data. This software utilises the aforementioned functionalities offered by the [PyOMA](https://github.com/dagghe/PyOMA) python module. Therefore, PyOMA_GUI provides a remarkably user-friendly interface to improve the accessibility of the PyOMA module, ensuring widespread usage both for scientists, researchers, and even for applied civil and structural engineers. The main features PyOMA_GUI provides are listed below:
- Importing data tab;
- Definition of the geometry of the structure and the monitoring system (channels and degrees of freedom, DOFs);
- Preprocessing of signals tool with detrending and decimation options;
- Dynamic identification algorithms with visualization of the results (graphs, modal shapes);
- Post-processing tabs and output export functionalities;

![`PyOMA_GUI` general overview.](Fig2.png)

# Acknowledgements
The authors acknowledge the meaningful contribution of Professor Rocco Alaggio from Università degli Studi dell'Aquila, who encouraged the authors to study and develop these topics. The authors acknowledge the meaningful contribution of Professor Giuseppe Carlo Marano from Politecnico di Torino for promoting the Graphical User Interface programming and coordinating the team activities.

# References
