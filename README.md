# DQuickLTFit
Copyright (c) 2016-2018 Danny Petschke (danny.petschke@uni-wuerzburg.de)<br><br>
<b>DQuickLTFit</b> - A Least-Square Fitting Tool for the Analysis of Positron-Lifetime Spectra consisting of discrete specific Lifetimes using the Levenberg-Marquardt Algorithm<br>

# How to cite this Software?

[![DOI](https://zenodo.org/badge/120471586.svg)](https://zenodo.org/badge/latestdoi/120471586)<br>

![DQuickLTFit](/TestData/Software.png)

# About
<b>DQuickLTFit</b> software is written in C++ ([Qt-framework](https://www.qt.io/)) and implements the [MPFIT](https://www.physics.wisc.edu/~craigm/idl/cmpfit.html) C library [1] for solving the non-linear least-square problem using the Levenberg-Marquardt algorithm. 
MPFIT was ported from [MINPACK-1](http://www.netlib.org/minpack/) library [2,3].<br><br>

# License (GNU General Public License)
This program is free software: you can redistribute it and/or modify<br>
it under the terms of the GNU General Public License as published by<br>
the Free Software Foundation, either version 3 of the License, or<br>
(at your option) any later version.<br><br>

This program is distributed in the hope that it will be useful,<br>
but WITHOUT ANY WARRANTY; without even the implied warranty of<br>
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.<br><br>

For more details see [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0)

# Setup
1. Download [Qt-framework](https://www.qt.io/download) (at least v5.x).
2. Download and Setup MS Visual Studio Compiler 2013 (or any newer version) .
3. Compile (x86/x64) and Run the code on your Machine.
4. Open the example project ([TestData.dquicklt](/TestData/TestData.dquicklt)) in the folder [TestData](/TestData/).
5. Finished. It should show up the [data](/TestData/TestData_5psPerChannel.dat) and fit(-results) as displayed above.

# References
[1] [C.B. Markwardt. Astron. Data Anal. Softw. Syst. XVIII ASP Conf. Ser. 2009; 411: 251.](https://arxiv.org/abs/0902.2850)<br>
[2] [J.J. Moré. in: Springer, Berlin, Heidelberg; 1978: pp. 105–116. doi:10.1007/BFb0067700](https://link.springer.com/chapter/10.1007/BFb0067700).
[3] [MINPACK-1, Jorge More'](http://www.netlib.org/minpack/)





