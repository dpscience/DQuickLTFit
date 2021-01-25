/****************************************************************************
**
**  DQuickLTFit, a software for the analysis of Positron-Lifetime Spectra
**  based on the Least-Square Optimization using the Levenberg-Marquardt
**  Algorithm.
**
**  Copyright (C) 2016-2021 Danny Petschke
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

#ifndef LIFETIMEDECAYFIT_H
#define LIFETIMEDECAYFIT_H

#include <cmath>

#include "../Settings/projectmanager.h"
#include "../Settings/settings.h"

#include "mpfit.h"

//#define __FITPARAM_DEBUG
//#define __FITQUEUE_DEBUG

#define __MAX_NUMBER_OF_FIT_RUNS 20

//additional error-enums for the mpfit.h
#define MP_ERR_NULLPTR_DATASTRUCTURE (-60) /*PALSDataStructure = nullptr;*/
#define MP_ERR_NULLPTR_FITSET_DATASET (-61) /*PALSFitSet || PALSDataSet = nullptr;*/
#define MP_ERR_NO_DATA (-62) /*no data to fit*/

typedef enum : int {
    yerror_Weighting = 1 /* assumption: Poisson noise */
} residualWeighting;

typedef struct {
  double *x;
  double *y;
  double *yInitial;
  double *ey;

  int dataCnt;

  double peakValue;
  double startChannel;
  double stopChannel;
  int startChannelIndex;
  int stopChannelIndex;
  int peakChannelIndex;

  int integralCountsInROI;
  double peakToBackgroundRatio;

  int countOfDeviceResolutionParams;

  double chiSquareOrig;

  int weighting;

  int mpfitRuns;

  int niter[__MAX_NUMBER_OF_FIT_RUNS];
  double chiSquareStart[__MAX_NUMBER_OF_FIT_RUNS];
  double chiSquareFinal[__MAX_NUMBER_OF_FIT_RUNS];

} values;

int multiExpDecay(int dataCnt, int ltParam, double *ltFitParamArray, double *dy, double **dvec, void *vars);

class LifeTimeDecayFitEngine : public QObject
{
    Q_OBJECT
public:
    LifeTimeDecayFitEngine();

public slots:
    void init(PALSDataStructure *dataStructure);
    void fit();

public:
    QList<QPointF> getFitPlotPoints() const;

private:
    void updateDataStructureFromResult(PALSDataStructure *dataStructure, mp_result *result, values *v, double *params);
    void createResultString(PALSDataStructure *dataStructure, values *v);

signals:
    void finished();

private:
    QList<QPointF> m_fitPlotSet;
    PALSDataStructure *m_dataStructure;
};

class PALSFitErrorCodeStringBuilder
{
public:
    inline static QString errorString(int errorCode) {
        switch ( errorCode ) {
        case 0:
            return QString("General Input Parameter Error.");
            break;

        case 1:
            return QString("OK. Convergence in &#967;<sup>2</sup>.");
            break;

        case 2:
            return QString("OK. Convergence in Parameter Value.");
            break;

        case 3:
            return QString("OK. Convergence in &#967;<sup>2</sup> & Parameter Value.");
            break;

        case 4:
            return QString("OK. Convergence in Orthogonality.");
            break;

        case 5:
            return QString("OK. Maximum Number of Iterations reached.");
            break;

        case 6:
            return QString("OK. No further Imrovements: Relative &#967;<sup>2</sup>-Convergence Criterium.");
            break;

        case 7:
            return QString("OK. No further Imrovements: Relative Parameter-Convergence Criterium.");
            break;

        case 8:
            return QString("OK. No further Imrovements: Orthogonality-Convergence Criterium.");
            break;

        case -16:
            return QString("Error. User-Function produced non-finite Values.");
            break;

        case -17:
            return QString("Error. No User Function was supplied.");
            break;

        case -18:
            return QString("Error. No User Data-Points were supplied.");
            break;

        case -19:
            return QString("Error. No free Parameters.");
            break;

        case -20:
            return QString("Error. Memory Allocation Error.");
            break;

        case -21:
            return QString("Error. Initial Values inconsistent with Constraints.");
            break;

        case -22:
            return QString("Error. Initial Constraints inconsistent.");
            break;

        case -23:
            return QString("Error. General Input Parameter Error.");
            break;

        case -24:
            return QString("Error. Not enough degrees of freedom.");
            break;

        case -60:
            return QString("Error: Internal Nullptr.");
            break;

        case -61:
            return QString("Error: Internal Nullptr.");
            break;

        case -62:
            return QString("Error: Internal Nullptr.");
            break;

        default:
            return QString("");
            break;
        }

        Q_UNREACHABLE();
    }
};

#endif // LIFETIMEDECAYFIT_H
