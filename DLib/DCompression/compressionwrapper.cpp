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

#include "compressionwrapper.h"

#include "miniz.c"

bool DCompressor::compressIt(QByteArray *pDest, const QByteArray &pSource, DCompressor::COMPRESSION_LEVEL level)
{
    if ( !pDest )
        return false;


    unsigned char* toCompress = (unsigned char*)pSource.data();
    mz_ulong len = (mz_ulong)pSource.size();

    mz_ulong outputLen = len;
    unsigned char* compressedData = (unsigned char*)malloc((size_t)len);

    const int returnVal = mz_compress2(compressedData, &outputLen, toCompress, len, level);

    pDest->resize(outputLen);
    for ( int i = 0 ; i < outputLen ; ++ i )
        pDest->append(compressedData[i]);

    if ( returnVal == MZ_OK )
        return true;
    else
        return false;
}

bool DCompressor::uncompressIt(QByteArray *pDest, const QByteArray &pSource)
{
    if ( !pDest )
        return false;


    const unsigned char* toUncompress = (const unsigned char*)pSource.data();
    const mz_ulong len = (mz_ulong)pSource.size();

    mz_ulong outputLen = len;
    unsigned char* uncompressedData;

    const int returnVal = mz_uncompress(uncompressedData, &outputLen, toUncompress, len);

    pDest->resize(outputLen);
    for ( int i = 0 ; i < outputLen ; ++ i )
        pDest->append(uncompressedData[i]);

    if ( returnVal == MZ_OK )
        return true;
    else
        return false;
}

QByteArray DCompressor::unzip(const QByteArray &pSource)
{
    QByteArray pDest(0);

    const bool returnVal = uncompressIt(&pDest, pSource);
    DUNUSED_PARAM(returnVal);

    return pDest;
}

QByteArray DCompressor::zip(const QByteArray &pSource, DCompressor::COMPRESSION_LEVEL level)
{
    QByteArray pDest(0);

    const bool returnVal = compressIt(&pDest, pSource, level);
    DUNUSED_PARAM(returnVal);

    return pDest;
}

DCompressor::DCompressor() {}
DCompressor::~DCompressor() {}
