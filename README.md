# DQuickLTFit
Copyright (c) 2016-2018 Danny Petschke (danny.petschke@uni-wuerzburg.de)<br><br>
DQuickLTFit - A Least-Square Fitting Tool for the Analysis of Positron-Lifetime Spectra consisting of discrete specific Lifetimes using the Levenberg-Marquardt Algorithm

![DQuickLTFit](/TestData/Software.png)

# About
DQuickLTFit software is written in C++ (Qt-framework) and implements the MPFIT C library [1] for solving the non-linear least-square problem using the Levenberg-Marquardt algorithm. 
MPFIT was ported from MINPACK-1 library [2].

# License (GNU General Public License)
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License<br> 
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.<br><br>

The complete form of this license is available at <ref>http://www.gnu.org/licenses/</ref>

# Setup
1. Download Qt-framework (at least v5.x) at <ref>https://www.qt.io/download</ref><br>
2. Download and Setup MS Visual Studio Compiler 2013 (or any newer version) 
3. Compile (x86/x64) and Run the code
4. Open the example project (TestData.dquicklt) in the folder /TestData/

# References
[1] C.B. Markwardt. Astron. Data Anal. Softw. Syst. XVIII ASP Conf. Ser. 2009; 411: 251.<br>
[2] J.J. Moré. in: Springer, Berlin, Heidelberg; 1978: pp. 105–116. doi:10.1007/BFb0067700.

