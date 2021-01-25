![badge-OS](https://img.shields.io/badge/OS-tested%20under%20Windows%2010-brightgreen)

Support this project and keep always updated about recent software releases, bug fixes and major improvements by [following on github](https://github.com/dpscience?tab=followers).

![badge-followers](https://img.shields.io/github/followers/dpscience?style=social)
![badge-stars](https://img.shields.io/github/stars/dpscience/DQuickLTFit?style=social)
![badge-forks](https://img.shields.io/github/forks/dpscience/DQuickLTFit?style=social)

# DQuickLTFit <img src="https://github.com/dpscience/DQuickLTFit/blob/master/Images/IconPNGRounded.png" width="25" height="25">   

![badge-license](https://img.shields.io/badge/OS-Windows-blue)
![badge-language](https://img.shields.io/badge/language-C++-blue)
![badge-license](https://img.shields.io/badge/license-GPL-blue)

Copyright (c) 2016-2021 Danny Petschke (danny.petschke@uni-wuerzburg.de) All rights reserved.<br><br>
<b>DQuickLTFit</b> - A least-square fitting tool for the analysis of positron lifetime spectra consisting of discrete characteristic lifetimes using the Levenberg-Marquardt algorithm.<br>

<br>![DQuickLTFit](/TestData/Software.png)

# Quickstart Guide on Windows OS <img src="https://github.com/dpscience/DQuickLTFit/blob/master/Images/IconPNGRounded.png" width="25" height="25"> 

## ``Option 1: via installer``
* Download the latest installer (<b>installer_DDRS4PALS-v1.14.exe</b>): https://github.com/dpscience/DDRS4PALS/releases
* Run the installer.
* Run the <b>DQuickLTFit</b> executable.

<b>The repository of the DQuickLTFit-installer can be found</b> [here.](https://github.com/dpscience/DQuickLTFit-installer)

## ``Option 2: manual installation``
* Download the latest software release (<b>DQuickLTFit_v4_2.rar</b>): https://github.com/dpscience/DQuickLTFit/releases
* Unzip <b>DQuickLTFit_v4_2.rar</b>.
* Download and install the Visual C++ Redistributable Package (x64) if required: https://www.microsoft.com/de-de/download/details.aspx?id=48145
* Run the <b>DQuickLTFit</b> executable.

# How to cite this Software? <img src="https://github.com/dpscience/DQuickLTFit/blob/master/Images/IconPNGRounded.png" width="25" height="25">

* <b>You must cite the applied version of this software in your study.</b><br>

You can cite all versions by using the <b>DOI 10.5281/zenodo.1168285</b>. This DOI represents all versions, and will always resolve to the latest one.<br>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1168285.svg)](https://doi.org/10.5281/zenodo.1168285)

## ``v4.x``
DQuickLTFit v4.2:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3356830.svg)](https://doi.org/10.5281/zenodo.3356830)<br>

DQuickLTFit v4.1:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3356830.svg)](https://doi.org/10.5281/zenodo.3356830)<br>

DQuickLTFit v4.0:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1414142.svg)](https://doi.org/10.5281/zenodo.1414142)<br>

## ``v3.x``
DQuickLTFit v3.02:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1219482.svg)](https://doi.org/10.5281/zenodo.1219482)<br>

DQuickLTFit v3.01:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1208613.svg)](https://doi.org/10.5281/zenodo.1208613)<br>

DQuickLTFit v3.0:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1202345.svg)](https://doi.org/10.5281/zenodo.1202345)<br>

## ``v2.x``
DQuickLTFit v2.06:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1168286.svg)](https://doi.org/10.5281/zenodo.1168286)<br>

# License (GNU General Public License) <img src="https://github.com/dpscience/DQuickLTFit/blob/master/Images/IconPNGRounded.png" width="25" height="25">

Copyright (c) 2016-2021 Danny Petschke (danny.petschke@uni-wuerzburg.de) All rights reserved.<br>

<p align="justify">This program is free software: you can redistribute it and/or modify<br>
it under the terms of the GNU General Public License as published by<br>
the Free Software Foundation, either version 3 of the License, or<br>
(at your option) any later version.<br><br>

This program is distributed in the hope that it will be useful,<br>
but WITHOUT ANY WARRANTY; without even the implied warranty of<br>
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.<br><br></p>

For more details see [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0)

# Used Third Party Libraries and Licenses <img src="https://github.com/dpscience/DQuickLTFit/blob/master/Images/IconPNGRounded.png" width="25" height="25">

<b>DQuickLTFit</b> is written in C++ using the [Qt-Framework](https://www.qt.io/) licensed under the [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0)

The following 3rd party library is used by <b>DQuickLTFit</b> software.<br>

### ``MPFIT (MINPACK-1)``
* [MPFIT: A MINPACK-1 Least Squares Fitting Library in C](https://www.physics.wisc.edu/~craigm/idl/cmpfit.html)<br>

The disclaimer of MPFIT library can be found [here.](/Fit/mpfit_DISCLAIMER).

# Deploy DQuickLTFit from Sources using QtCreator<img src="https://github.com/dpscience/DQuickLTFit/blob/master/Images/IconPNGRounded.png" width="25" height="25">

* Download the QtCreator and the [Qt-framework](https://www.qt.io/download) (at least v5.x).
* Download and Setup the MS Visual Studio compiler (at least version 2013). It should also work with any other compiler e.g. MinGW but I recommend using VS compiler.
* Open the .pro file in QtCreator. 
* Deploy DQuickLTFit. It should finish without any errors.
* Finished.
