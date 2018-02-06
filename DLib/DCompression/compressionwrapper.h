/****************************************************************************
**
**  DQuickLTFit, a software for the analysis of Positron-Lifetime Spectra
**  based on the Least-Square Optimization using the Levenberg-Marquardt
**  Algorithm.
**
**  Copyright (C) 2016-2018 Danny Petschke
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see http://www.gnu.org/licenses/.
**
*****************************************************************************
**
**  @author: Danny Petschke
**  @contact: danny.petschke@uni-wuerzburg.de
**
*****************************************************************************/

#ifndef COMPRESSIONWRAPPER_H
#define COMPRESSIONWRAPPER_H

#include "../DTypes/types.h"

class DCompressor
{
public:
    typedef enum{
        NO_COMPRESSION = 0,
        BEST_SPEED = 1,
        BEST_COMPRESSION = 9,
        UBER_COMPRESSION = 10,
        DEFAULT_LEVEL = 6,
        DEFAULT_COMPRESSION = -1
    }COMPRESSION_LEVEL;

    static bool compressIt(QByteArray *pDest, const QByteArray& pSource, COMPRESSION_LEVEL level = DEFAULT_LEVEL);
    static bool uncompressIt(QByteArray *pDest, const QByteArray& pSource);

    static QByteArray zip(const QByteArray& pSource, COMPRESSION_LEVEL level = DEFAULT_LEVEL);
    static QByteArray unzip(const QByteArray& pSource);

private:
    DCompressor();
    virtual ~DCompressor();
};

#endif // COMPRESSIONWRAPPER_H
