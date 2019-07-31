# DQuickLTFit  <img src="https://github.com/dpscience/DQuickLTFit/blob/master/Images/IconPNGRounded.png" width="25" height="25">   

Copyright (c) 2016-2019 Danny Petschke (danny.petschke@uni-wuerzburg.de) All rights reserved.<br><br>
<b>DQuickLTFit</b> - A Least-Square Fitting Tool for the Analysis of Positron-Lifetime Spectra consisting of discrete specific Lifetimes using the Levenberg-Marquardt Algorithm<br>

#### We recommend the use of DQuickLTFit v4.0 or higher as this is the first stable version.

<br>![DQuickLTFit](/TestData/Software.png)

## How to cite this Software?

You can cite all versions by using the <b>DOI 10.5281/zenodo.1168285</b>. This DOI represents all versions, and will always resolve to the latest one.<br>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1168285.svg)](https://doi.org/10.5281/zenodo.1168285)

## v4.x
DQuickLTFit v4.1:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3356830.svg)](https://doi.org/10.5281/zenodo.3356830)<br>
DQuickLTFit v4.0:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1414142.svg)](https://doi.org/10.5281/zenodo.1414142)<br>

## v3.x
DQuickLTFit v3.02:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1219482.svg)](https://doi.org/10.5281/zenodo.1219482)<br>
DQuickLTFit v3.01:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1208613.svg)](https://doi.org/10.5281/zenodo.1208613)<br>
DQuickLTFit v3.0:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1202345.svg)](https://doi.org/10.5281/zenodo.1202345)<br>

## v2.x
DQuickLTFit v2.06:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1168286.svg)](https://doi.org/10.5281/zenodo.1168286)<br>

# About
<b>DQuickLTFit</b> software is written in C++ ([Qt-framework](https://www.qt.io/)) and has implemented the [MPFIT](https://www.physics.wisc.edu/~craigm/idl/cmpfit.html) C library [1] for solving the non-linear least-square problem using the Levenberg-Marquardt algorithm. 
MPFIT was ported from [MINPACK-1](http://www.netlib.org/minpack/) library [2,3].<br><br>
The disclaimer of MPFIT library can be found [here](/Fit/mpfit_DISCLAIMER).<br>

# License (GNU General Public License)
Copyright (c) 2016-2019 Danny Petschke (danny.petschke@uni-wuerzburg.de) All rights reserved.<br><br>

<p align="justify">This program is free software: you can redistribute it and/or modify<br>
it under the terms of the GNU General Public License as published by<br>
the Free Software Foundation, either version 3 of the License, or<br>
(at your option) any later version.<br><br>

This program is distributed in the hope that it will be useful,<br>
but WITHOUT ANY WARRANTY; without even the implied warranty of<br>
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.<br><br></p>

For more details see [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0)

# Run executable (*.exe)
The executable (and binary files) can be either downloaded from <b>zenodo platform</b> by following the related DOI (see above) or from the corresponding release.

# Setup from source files
1. Download [Qt-framework](https://www.qt.io/download) (at least v5.x).
2. Download and Setup MS Visual Studio Compiler 2013 (or any newer version) .
3. Compile (x86/x64) and Run the code on your Machine.
4. Open the example project ([TestData.dquicklt](/TestData/TestData.dquicklt)) in the folder [TestData](/TestData/).
5. Finished. It should show up the [data](/TestData/TestData_5psPerChannel.dat) and fit(-results) as displayed above.

# References
[1] [C.B. Markwardt. Astron. Data Anal. Softw. Syst. XVIII ASP Conf. Ser. 2009; 411: 251.](https://arxiv.org/abs/0902.2850)<br>
[2] [J.J. Moré. in: Springer, Berlin, Heidelberg; 1978: pp. 105–116. doi:10.1007/BFb0067700](https://link.springer.com/chapter/10.1007/BFb0067700).<br>
[3] [MINPACK-1, Jorge More'](http://www.netlib.org/minpack/)




