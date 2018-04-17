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

#include "lifetimedecayfit.h"

double weightedResidual(double y, double yi, double erryi, int weightingCode)
{
    DUNUSED_PARAM(erryi);

    double residual = 0.0f;
    double weight = 0.0f;

    switch ( weightingCode )
    {
    case 0: //no weighting:
        weight = 1;

        residual = weight*(yi - y);
        break;

    case 1: //y-variance weighting:
        if (qFuzzyCompare(yi, 0.0) || yi < 0 || (int)yi <= 0)
            weight = 1;
        else
            weight = 1.0/yi;

        residual = weight*(yi - y);
        break;

    case 2: //y-error weighting:
        if (qFuzzyCompare(yi, 0.0) || yi < 0 || (int)yi <= 0)
            weight = 1;
        else
            weight = 1.0/sqrt(yi);

        residual = weight*(yi - y);
        break;

    default:
        residual = 0;
        break;
    }

    return residual;
}

/*
 * fit function declarations:
 *---------------------------
 *
 * dataCnt - number of data points
 * ltParam - number of lifetime-components -> means <tau> and <I>
 * ltFitParamArray - array of lifetime-fit parameters (tau and I)
 * gaussianSigma - sigma-value of gaussian distribution
 * gaussianMu - t0 of gaussian distribution
 * bgrd - background
 *
 * dy - array of residuals to be returned
 * vars - private data (struct values *)
 *        which contains the x/y-values and the y-uncertainties
 *
 * returns: 0 for success and 1 for failed!
 */

int lifeTimeDecaySum(int dataCnt, int paramCnt, double *fitParamArray, double *dy, double **dvec, void *vars)
{
        DUNUSED_PARAM(dvec);

        values *v = (values*) vars;

        double *x = v->x;
        double *y = v->y;
        double *ey = v->ey;

        const int cntGaussian = v->countOfDeviceResolutionParams;
        const double bkgrd = fitParamArray[paramCnt-1];


        for ( int i = 0 ; i < dataCnt-1 ; ++ i )
        {
            double f = 0;

            x[i] -= v->startChannel;
            x[i+1] -= v->startChannel;

            for ( int device = paramCnt - cntGaussian - 1 ; device < paramCnt - 1 ; device += 3 )
            {
                const double gaussianSigma = fitParamArray[device]/(2*sqrt(log(2)));
                const double gaussianMu = fitParamArray[device+1];
                const double gaussianIntensity = fitParamArray[device+2];

                for ( int param = 0 ; param <  paramCnt - cntGaussian - 1/*3*/ ; param += 2 ) // first param[0] = tau; second param[1] = Intensity
                {
                    //Eldrup-Formula:
                    const double yji = exp(-(x[i]-gaussianMu-(gaussianSigma*gaussianSigma)/(4*fitParamArray[param]))/fitParamArray[param])*(1-erf((0.5*gaussianSigma/fitParamArray[param])-(x[i]-gaussianMu)/gaussianSigma));
                    const double yji_plus_1 = exp(-(x[i+1]-gaussianMu-(gaussianSigma*gaussianSigma)/(4*fitParamArray[param]))/fitParamArray[param])*(1-erf((0.5*gaussianSigma/fitParamArray[param])-(x[i+1]-gaussianMu)/gaussianSigma));

                    f += 0.5*fitParamArray[param+1]/**fitParamArray[param]*/*(yji-yji_plus_1-erf((x[i]-gaussianMu)/gaussianSigma)+erf((x[i+1]-gaussianMu)/gaussianSigma));
                }

                f *= gaussianIntensity;
            }

            x[i] += v->startChannel;
            x[i+1] += v->startChannel;

            f *= ((double)v->integralCountsInROI)-(v->stopChannel-v->startChannel)*bkgrd;
            f += bkgrd;

            dy[i] = weightedResidual(f, y[i], ey[i], v->weighting);
        }

        return 1;
}


LifeTimeDecayFitEngine::LifeTimeDecayFitEngine() :
    m_dataStructure(nullptr) {}

void LifeTimeDecayFitEngine::init(PALSDataStructure *dataStructure)
{
    m_dataStructure = dataStructure;
}

void LifeTimeDecayFitEngine::fit()
{
    PALSDataStructure *dataStructure = m_dataStructure;

    if ( !dataStructure )
        return;

    if ( !dataStructure->getDataSetPtr() || !dataStructure->getFitSetPtr() )
        return;

    if ( dataStructure->getDataSetPtr()->getLifeTimeData().isEmpty() )
        return;

    //initialize the data-set:
    const int paramCnt = dataStructure->getFitSetPtr()->getComponentsCount() + dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize() + 1;

    double channelResolution = dataStructure->getFitSetPtr()->getChannelResolution(); //[ps/chn]

    double startChannel = dataStructure->getFitSetPtr()->getStartChannel();
    double stopChannel = dataStructure->getFitSetPtr()->getStopChannel();
    double peakChannel = 0;
    double countsInPeak = -LONG_MAX;

     int startChannelIndex = 0;
     int stopChannelIndex = 0;
     int peakChannelIndex = 0;

     int channelCnt = 0;
     int dataCntInRange = stopChannel-startChannel+1;

     double *x = new double[dataCntInRange];
     double *y = new double[dataCntInRange];
     double *ey = new double[dataCntInRange];

     int inRangeCnt = 0;
     int integralCountROI = 0;

     for ( QPointF p : dataStructure->getDataSetPtr()->getLifeTimeData() )
     {
         if ( p.x() >= startChannel && p.x() <= stopChannel )
         {
             x[inRangeCnt] = p.x();
             y[inRangeCnt] = p.y();
             ey[inRangeCnt] = 0;

             integralCountROI += (int)p.y();

             if ( p.x() == startChannel )
                 startChannelIndex = channelCnt;

             if ( p.x() == stopChannel )
                 stopChannelIndex = channelCnt;

             if ( p.y() > countsInPeak )
             {
                 countsInPeak = p.y();
                 peakChannel = p.x();
                 peakChannelIndex = channelCnt;
             }

             inRangeCnt ++;
         }

         channelCnt ++;
     }

     values v;

     v.dataCnt = dataCntInRange;

     v.x = x;
     v.y = y;
     v.yInitial = y;
     v.ey = ey;
     v.peakValue = countsInPeak;
     v.startChannelIndex = startChannelIndex;
     v.stopChannelIndex = stopChannelIndex;
     v.peakChannelIndex = peakChannelIndex;
     v.startChannel = startChannel;
     v.stopChannel = stopChannel;
     v.integralCountsInROI = integralCountROI;

     v.countOfDeviceResolutionParams = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize();

     v.weighting = (dataStructure->getFitSetPtr()->usingYVariance()?residualWeighting::yvariance_Weighting:residualWeighting::no_Weighting);

    //Parameter constraints:
    mp_par *paramContraints = new mp_par[paramCnt];

    for ( int t = 0 ; t < paramCnt ; ++ t ) {
        paramContraints[t] = {0};
    }

    //Parameter initial conditions:
    double *params = new double[paramCnt]; //source->sample->gaussian->bkgrd

    int i = 0;
    for ( i = 0 ; i < dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i+=2 )
    {
        const PALSFitParameter *fitParam = dataStructure->getFitSetPtr()->getSourceParamPtr()->getParameterAt(i);
        params[i] = fitParam->getStartValue()/channelResolution; //tau

#ifdef __FITPARAM_DEBUG
        qDebug() << "";
        qDebug() << "param: " % QString(fitParam->getName()) % " (" % QString(fitParam->getAlias()) % ")";
        qDebug() << "value:" % QVariant(fitParam->getStartValue()).toString();
#endif

        paramContraints[i].deriv_debug = 0;
        paramContraints[i+1].deriv_debug = 0;

        if ( fitParam->isFixed() )
            paramContraints[i].fixed = 1;
        else
            paramContraints[i].fixed = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "fixed?: " % QVariant(fitParam->isFixed()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding?: " % QVariant(fitParam->isLowerBoundingEnabled()).toString();
#endif

        if ( fitParam->isLowerBoundingEnabled() )
        {
            paramContraints[i].limited[0] = 1;
            paramContraints[i].limits[0] = fitParam->getLowerBoundingValue()/channelResolution;

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding-value: " % QVariant(fitParam->getLowerBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i].limited[0] = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding?: " % QVariant(fitParam->isUpperBoundingEnabled()).toString();
#endif

        if ( fitParam->isUpperBoundingEnabled() )
        {
            paramContraints[i].limited[1] = 1;
            paramContraints[i].limits[1] = fitParam->getUpperBoundingValue()/channelResolution;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding-value: " % QVariant(fitParam->getUpperBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i].limited[1] = 0;


        const PALSFitParameter *fitParam2 = dataStructure->getFitSetPtr()->getSourceParamPtr()->getParameterAt(i+1);
        params[i+1] = fitParam2->getStartValue(); //I

#ifdef __FITPARAM_DEBUG
        qDebug() << "";
        qDebug() << "param: " % QString(fitParam2->getName()) % " (" % QString(fitParam2->getAlias()) % ")";
        qDebug() << "value:" % QVariant(fitParam2->getStartValue()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "fixed?: " % QVariant(fitParam2->isFixed()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding?: " % QVariant(fitParam2->isLowerBoundingEnabled()).toString();
#endif

        if ( fitParam2->isFixed() )
            paramContraints[i+1].fixed = 1;
        else
            paramContraints[i+1].fixed = 0;

        if ( fitParam2->isLowerBoundingEnabled() )
        {
            paramContraints[i+1].limited[0] = 1;
            paramContraints[i+1].limits[0] = fitParam2->getLowerBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding-value: " % QVariant(fitParam2->getLowerBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i+1].limited[0] = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding?: " % QVariant(fitParam2->isUpperBoundingEnabled()).toString();
#endif

        if ( fitParam2->isUpperBoundingEnabled() )
        {
            paramContraints[i+1].limited[1] = 1;
            paramContraints[i+1].limits[1] = fitParam2->getUpperBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding-value: " % QVariant(fitParam2->getUpperBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i+1].limited[1] = 0;
    }

    int cnt = 0;
    for ( i = dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i < dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i+=2 )
    {
        const PALSFitParameter *fitParam = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getParameterAt(cnt);
        params[i] = fitParam->getStartValue()/channelResolution; //tau

#ifdef __FITPARAM_DEBUG
        qDebug() << "";
        qDebug() << "param: " % QString(fitParam->getName()) % " (" % QString(fitParam->getAlias()) % ")";
        qDebug() << "value:" % QVariant(fitParam->getStartValue()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "fixed?: " % QVariant(fitParam->isFixed()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding?: " % QVariant(fitParam->isLowerBoundingEnabled()).toString();
#endif

        paramContraints[i].deriv_debug = 0;
        paramContraints[i+1].deriv_debug = 0;

        if ( fitParam->isFixed() )
            paramContraints[i].fixed = 1;
        else
            paramContraints[i].fixed = 0;

        if ( fitParam->isLowerBoundingEnabled() )
        {
            paramContraints[i].limited[0] = 1;
            paramContraints[i].limits[0] = fitParam->getLowerBoundingValue()/channelResolution;

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding-value: " % QVariant(fitParam->getLowerBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i].limited[0] = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding?: " % QVariant(fitParam->isUpperBoundingEnabled()).toString();
#endif

        if ( fitParam->isUpperBoundingEnabled() )
        {
            paramContraints[i].limited[1] = 1;
            paramContraints[i].limits[1] = fitParam->getUpperBoundingValue()/channelResolution;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding-value: " % QVariant(fitParam->getUpperBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i].limited[1] = 0;


        cnt ++;

        const PALSFitParameter *fitParam2 = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getParameterAt(cnt);
        params[i+1] = fitParam2->getStartValue(); //I

#ifdef __FITPARAM_DEBUG
        qDebug() << "";
        qDebug() << "param: " % QString(fitParam2->getName()) % " (" % QString(fitParam2->getAlias()) % ")";
        qDebug() << "value:" % QVariant(fitParam2->getStartValue()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "fixed?: " % QVariant(fitParam2->isFixed()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding?: " % QVariant(fitParam2->isLowerBoundingEnabled()).toString();
#endif

        if ( fitParam2->isFixed() )
            paramContraints[i+1].fixed = 1;
        else
            paramContraints[i+1].fixed = 0;

        if ( fitParam2->isLowerBoundingEnabled() )
        {
            paramContraints[i+1].limited[0] = 1;
            paramContraints[i+1].limits[0] = fitParam2->getLowerBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding-value: " % QVariant(fitParam2->getLowerBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i+1].limited[0] = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding?: " % QVariant(fitParam2->isUpperBoundingEnabled()).toString();
#endif

        if ( fitParam2->isUpperBoundingEnabled() )
        {
            paramContraints[i+1].limited[1] = 1;
            paramContraints[i+1].limits[1] = fitParam2->getUpperBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding-value: " % QVariant(fitParam2->getUpperBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i+1].limited[1] = 0;

        cnt ++;
    }

    /***********************************************************************/

    int cntGaussian= 0;
    for ( i = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i < dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() + dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize() ; i+=3 )
    {
        const PALSFitParameter *fitParam = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        params[i] = fitParam->getStartValue()/channelResolution; //FWHM


#ifdef __FITPARAM_DEBUG
        qDebug() << "";
        qDebug() << "param: " % QString(fitParam->getName()) % " (" % QString(fitParam->getAlias()) % ")";
        qDebug() << "value:" % QVariant(fitParam->getStartValue()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "fixed?: " % QVariant(fitParam->isFixed()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding?: " % QVariant(fitParam->isLowerBoundingEnabled()).toString();
#endif

        paramContraints[i].deriv_debug = 0;
        paramContraints[i+1].deriv_debug = 0;
        paramContraints[i+2].deriv_debug = 0;

        if ( fitParam->isFixed() )
            paramContraints[i].fixed = 1;
        else
            paramContraints[i].fixed = 0;

        if ( fitParam->isLowerBoundingEnabled() )
        {
            paramContraints[i].limited[0] = 1;
            paramContraints[i].limits[0] = fitParam->getLowerBoundingValue()/channelResolution; //FWHM
            paramContraints[i].limits[0] /= paramContraints[i].limits[0]/(2*sqrt(log(2))); //sigma

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding-value: " % QVariant(fitParam->getLowerBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i].limited[0] = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding?: " % QVariant(fitParam->isUpperBoundingEnabled()).toString();
#endif

        if ( fitParam->isUpperBoundingEnabled() )
        {
            paramContraints[i].limited[1] = 1;
            paramContraints[i].limits[1] = fitParam->getUpperBoundingValue()/channelResolution; //FWHM
            paramContraints[i].limits[1] /= paramContraints[i].limits[1]/(2*sqrt(log(2))); //sigma

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding-value: " % QVariant(fitParam->getUpperBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i].limited[1] = 0;


        cntGaussian ++;

        const PALSFitParameter *fitParam2 = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        params[i+1] = fitParam2->getStartValue()/channelResolution; //mu

#ifdef __FITPARAM_DEBUG
        qDebug() << "";
        qDebug() << "param: " % QString(fitParam2->getName()) % " (" % QString(fitParam2->getAlias()) % ")";
        qDebug() << "value:" % QVariant(fitParam2->getStartValue()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "fixed?: " % QVariant(fitParam2->isFixed()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding?: " % QVariant(fitParam2->isLowerBoundingEnabled()).toString();
#endif

        if ( fitParam2->isFixed() )
            paramContraints[i+1].fixed = 1;
        else
            paramContraints[i+1].fixed = 0;

        if ( fitParam2->isLowerBoundingEnabled() )
        {
            paramContraints[i+1].limited[0] = 1;
            paramContraints[i+1].limits[0] = fitParam2->getLowerBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding-value: " % QVariant(fitParam2->getLowerBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i+1].limited[0] = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding?: " % QVariant(fitParam2->isUpperBoundingEnabled()).toString();
#endif

        if ( fitParam2->isUpperBoundingEnabled() )
        {
            paramContraints[i+1].limited[1] = 1;
            paramContraints[i+1].limits[1] = fitParam2->getUpperBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding-value: " % QVariant(fitParam2->getUpperBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i+1].limited[1] = 0;

        cntGaussian ++;


        const PALSFitParameter *fitParam3 = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        params[i+2] = fitParam3->getStartValue(); //I

#ifdef __FITPARAM_DEBUG
        qDebug() << "";
        qDebug() << "param: " % QString(fitParam3->getName()) % " (" % QString(fitParam3->getAlias()) % ")";
        qDebug() << "value:" % QVariant(fitParam3->getStartValue()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "fixed?: " % QVariant(fitParam3->isFixed()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding?: " % QVariant(fitParam3->isLowerBoundingEnabled()).toString();
#endif

        if ( fitParam3->isFixed() )
            paramContraints[i+2].fixed = 1;
        else
            paramContraints[i+2].fixed = 0;

        if ( fitParam3->isLowerBoundingEnabled() )
        {
            paramContraints[i+2].limited[0] = 1;
            paramContraints[i+2].limits[0] = fitParam3->getLowerBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding-value: " % QVariant(fitParam3->getLowerBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i+2].limited[0] = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding?: " % QVariant(fitParam3->isUpperBoundingEnabled()).toString();
#endif

        if ( fitParam3->isUpperBoundingEnabled() )
        {
            paramContraints[i+2].limited[1] = 1;
            paramContraints[i+2].limits[1] = fitParam3->getUpperBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding-value: " % QVariant(fitParam3->getUpperBoundingValue()).toString();
#endif
        }
        else
            paramContraints[i+2].limited[1] = 0;

        cntGaussian ++;
    }

    /***********************************************************************/

    //Gaussian-sigma, Gaussian-mu, Background:
    PALSFitParameter *bkgrd = dataStructure->getFitSetPtr()->getBackgroundParamPtr()->getParameter();

    //at the moment there is no variation of the background = const.!
    bkgrd->setAsFixed(true);
    bkgrd->setLowerBoundingEnabled(false);
    bkgrd->setUpperBoundingEnabled(false);

    const int bkgrdIndex = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() + dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize();//+2;

    params[bkgrdIndex] = bkgrd->getStartValue();

    v.peakToBackgroundRatio = (double)(countsInPeak-bkgrd->getStartValue())/bkgrd->getStartValue();

#ifdef __FITPARAM_DEBUG
        qDebug() << "";
        qDebug() << "param: " % QString(bkgrd->getName()) % " (" % QString(bkgrd->getAlias()) % ")";
        qDebug() << "value:" % QVariant(bkgrd->getStartValue()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "fixed?: " % QVariant(bkgrd->isFixed()).toString();
#endif

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding?: " % QVariant(bkgrd->isLowerBoundingEnabled()).toString();
#endif

    paramContraints[bkgrdIndex].deriv_debug = 0;

    if ( bkgrd->isFixed() )
        paramContraints[bkgrdIndex].fixed = 1;
    else
        paramContraints[bkgrdIndex].fixed = 0;

    if ( bkgrd->isLowerBoundingEnabled() )
    {
        paramContraints[bkgrdIndex].limited[0] = 1;
        paramContraints[bkgrdIndex].limits[0] = bkgrd->getLowerBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "lower-bounding-value: " % QVariant(bkgrd->getLowerBoundingValue()).toString();
#endif
    }
    else
        paramContraints[bkgrdIndex].limited[0] = 0;

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding?: " % QVariant(bkgrd->isUpperBoundingEnabled()).toString();
#endif

    if ( bkgrd->isUpperBoundingEnabled() )
    {
        paramContraints[bkgrdIndex].limited[1] = 1;
        paramContraints[bkgrdIndex].limits[1] = bkgrd->getUpperBoundingValue();

#ifdef __FITPARAM_DEBUG
            qDebug() << "upper-bounding-value: " % QVariant(bkgrd->getUpperBoundingValue()).toString();
#endif
    }
    else
        paramContraints[bkgrdIndex].limited[1] = 0;

    ///Returned parameter errors:
    double *paramErrors = new double[paramCnt];
    double *finalResiduals = new double[dataCntInRange];

    mp_result result;
    memset(&result,0,sizeof(result));

    result.xerror = paramErrors;
    result.resid = finalResiduals;

    mp_config config;
    memset(&config,0,sizeof(config));

    config.maxiter = dataStructure->getFitSetPtr()->getMaximumIterations();

    //run the fit:
    int stat = mpfit(lifeTimeDecaySum,
                              dataCntInRange,
                              paramCnt,
                              params,
                              paramContraints,
                              &config,
                              (void*) &v,
                              &result);

    DUNUSED_PARAM(stat);

    updateDataStructureFromResult(dataStructure, &result, &v, params);

    delete [] x;
    delete [] y;
    delete [] ey;

    delete [] params;
    delete [] paramContraints;

    delete [] paramErrors;
    delete [] finalResiduals;

    emit finished();
}

QList<QPointF> LifeTimeDecayFitEngine::getFitPlotPoints() const
{
    return m_fitPlotSet;
}

void LifeTimeDecayFitEngine::updateDataStructureFromResult(PALSDataStructure *dataStructure, mp_result *result, values *v, double *params)
{
    m_fitPlotSet.clear();

    dataStructure->getFitSetPtr()->setNeededIterations((unsigned int)result->niter);
    dataStructure->getFitSetPtr()->setCountsInRange(v->integralCountsInROI);
    dataStructure->getFitSetPtr()->setPeakToBackgroundRatio(v->peakToBackgroundRatio);
    dataStructure->getFitSetPtr()->setTimeStampOfLastFitResult(QDateTime::currentDateTime().toString());
    dataStructure->getFitSetPtr()->setFitFinishCodeValue(result->status);
    dataStructure->getFitSetPtr()->setFitFinishCode(PALSFitErrorCodeStringBuilder::errorString(result->status));

    const double channelResolution = dataStructure->getFitSetPtr()->getChannelResolution();
    const int paramCnt = dataStructure->getFitSetPtr()->getComponentsCount() + dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize() + 1;

    double sumOfIntensities = 0.0f;

    int i = 0;
    for ( i = 0 ; i < dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i+=2 )
    {
        PALSFitParameter *param_tau = dataStructure->getFitSetPtr()->getSourceParamPtr()->getParameterAt(i);
        param_tau->setFitValue(params[i]*channelResolution);
        param_tau->setFitValueError(result->xerror[i]*channelResolution);

        PALSFitParameter *param_I = dataStructure->getFitSetPtr()->getSourceParamPtr()->getParameterAt(i+1);
        param_I->setFitValue(params[i+1]);
        param_I->setFitValueError(result->xerror[i+1]);

        sumOfIntensities += param_I->getFitValue();
    }

    double tauAverage = 0.0f;
    double tauAverageError = 0.0f;

    int cnt = 0;
    for ( i = dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i < dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize()+dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize() ; i+=2 )
    {
        PALSFitParameter *param_tau = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getParameterAt(cnt);
        param_tau->setFitValue(params[i]*channelResolution);
        param_tau->setFitValueError(result->xerror[i]*channelResolution);

        tauAverage += param_tau->getFitValue();
        tauAverageError += param_tau->getFitValueError()*param_tau->getFitValueError();

        cnt ++;

        PALSFitParameter *param_I = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getParameterAt(cnt);
        param_I->setFitValue(params[i+1]);
        param_I->setFitValueError(result->xerror[i+1]);

        sumOfIntensities += param_I->getFitValue();

        cnt ++;
    }

    tauAverage /= (double)(cnt/2);
    tauAverageError = sqrtf(tauAverageError);

    dataStructure->getFitSetPtr()->setAverageLifeTime(tauAverage);
    dataStructure->getFitSetPtr()->setAverageLifeTimeError(tauAverageError);
    dataStructure->getFitSetPtr()->setSumOfIntensities(sumOfIntensities);


    /********************************************************************************************/

    int cntGaussian= 0;
    for ( i = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i < dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() + dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize() ; i+=3 )
    {
        PALSFitParameter *param_sigma = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        param_sigma->setFitValue(params[i]*channelResolution);
        param_sigma->setFitValueError(result->xerror[i]*channelResolution);

        cntGaussian ++;

        PALSFitParameter *param_mu = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        param_mu->setFitValue(params[i+1]*channelResolution);
        param_mu->setFitValueError(result->xerror[i+1]*channelResolution);

        cntGaussian ++;

        PALSFitParameter *param_I = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        param_I->setFitValue(params[i+2]);
        param_I->setFitValueError(result->xerror[i+2]);

        cntGaussian ++;
    }

    /********************************************************************************************/

    const int bkgrdIndex = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() + dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize();//+2;

    PALSFitParameter *bkgrd = dataStructure->getFitSetPtr()->getBackgroundParamPtr()->getParameter();
    bkgrd->setFitValue(params[bkgrdIndex]);
    bkgrd->setFitValueError(result->xerror[bkgrdIndex]);


    double chiSquare = 0.0f;

    QList<QPointF> residuals;

    const double bkgrdVal = params[bkgrdIndex];
    double intergralCountsWithoutBkgrd = ((double)v->integralCountsInROI)-(v->stopChannel-v->startChannel)*bkgrdVal;
    double tZeroChannel = 0;
    int tZeroIndex = 0;
    double maxf = -1;

    for ( int i = 0 ; i < v->dataCnt-1 ; ++ i )
    {
        double f = 0;
        double x = v->x[i]-v->startChannel;
        double x_plus_1 = v->x[i+1]-v->startChannel;

        for ( int device = paramCnt - cntGaussian - 1 ; device < paramCnt - 1 ; device += 3 )
        {
            const double gaussianSigmaVal = params[device]/(2*sqrt(log(2)));
            const double gaussianMuVal = params[device+1];
            const double gaussianIntensity = params[device+2];

            for ( int param = 0 ; param <  paramCnt - cntGaussian - 1/*3*/ ; param += 2 ) // first param[0] = tau; second param[1] = Intensity
            {
                //Eldrup-Formula:
                const double yji = exp(-(x-gaussianMuVal-(gaussianSigmaVal*gaussianSigmaVal)/(4*params[param]))/params[param])*(1-erf((0.5*gaussianSigmaVal/params[param])-(x-gaussianMuVal)/gaussianSigmaVal));
                const double yji_plus_1 = exp(-(x_plus_1-gaussianMuVal-(gaussianSigmaVal*gaussianSigmaVal)/(4*params[param]))/params[param])*(1-erf((0.5*gaussianSigmaVal/params[param])-(x_plus_1-gaussianMuVal)/gaussianSigmaVal));

                f += 0.5*params[param+1]/**params[param]*/*(yji-yji_plus_1-erf((x-gaussianMuVal)/gaussianSigmaVal)+erf((x_plus_1-gaussianMuVal)/gaussianSigmaVal));
            }

            f *= gaussianIntensity;
        }

        x += v->startChannel;

        f *=  intergralCountsWithoutBkgrd;
        f += bkgrdVal;

        if (f > maxf) {
            maxf = f;
            tZeroChannel = x;
            tZeroIndex = i;
        }

        if ( !qFuzzyCompare(v->y[i], 0.0) )
            chiSquare += (f-v->y[i])*(f-v->y[i])/v->y[i];

        const double res = (f-v->y[i])/sqrt(v->y[i]);

        m_fitPlotSet.append(QPointF(x, f));
        residuals.append(QPointF(x, res));
    }

    //center of mass (spectral centroid):
    double tCenter = 0;
    double sumOfCounts = 0;
    for ( int i = tZeroIndex ; i < m_fitPlotSet.size()-1 ; ++ i ) {
        const double time = ((m_fitPlotSet.at(i).x()-tZeroChannel) + 0.5)*channelResolution;
        const double counts = 0.5*(m_fitPlotSet.at(i).y()+m_fitPlotSet.at(i+1).y());

        tCenter += time*counts;
        sumOfCounts += counts;
    }

    tCenter /= sumOfCounts;

    chiSquare /= (double)v->dataCnt;

    dataStructure->getFitSetPtr()->setChiSquareAfterFit(chiSquare);

    dataStructure->getFitSetPtr()->setTZeroSpectralCentroid((tZeroChannel-v->startChannel)*channelResolution);
    dataStructure->getFitSetPtr()->setSpectralCentroid(tCenter);

    dataStructure->getDataSetPtr()->setResiduals(residuals);
    dataStructure->getDataSetPtr()->setFitData(m_fitPlotSet);

    createResultString(dataStructure);
}

void LifeTimeDecayFitEngine::createResultString(PALSDataStructure *dataStructure)
{
    if ( !dataStructure )
        return;


    const PALSFitSet *fitSet = dataStructure->getFitSetPtr();

    const QString lineBreak("<br>");

    const QString tableStart("<table>");
    const QString tableEnd("</table>");

    const QString tableBorderStart("<table border=\"1\" style=\"width:100%\">");
    const QString tableBorderEnd("</table>");

    const QString startRow("<tr>");
    const QString finishRow("</tr>");

    const QString startHeader("<th>");
    const QString endHeader("</th>");

    const QString startContent("<td>");
    const QString finishContent("</td>");

    const QString alignCenterStart("<div align=\"center\">");
    const QString alignCenterEnd("</div>");

    const QString spacer("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;");

    const QString alertHtml = "<font color=\"DeepPink\">";
    const QString notifyHtml = "<font color=\"Lime\">";
    const QString infoHtml = "<font color=\"Aqua\">";
    const QString endHtml = "</font>";

    const QString projectName("<nobr><b>Project:</b></nobr>");
    const QString asciiFileName("<nobr><b>Raw-Data:</b></nobr>");

    const QString fitFinishCode("<nobr><b>Finish-Code:</b></nobr>");

    QString fitFinishCodeVal = "";
    if ( fitSet->getFitFinishCodeValue() == 1 )
        fitFinishCodeVal = QString("<nobr><b>" % notifyHtml % fitSet->getFitFinishCode() % " </b>" % endHtml % "[" % fitSet->getTimeStampOfLastFitResult() % "]</nobr>");
    else if ( fitSet->getFitFinishCodeValue() > 1 )
        fitFinishCodeVal = QString("<nobr><b>" % infoHtml % fitSet->getFitFinishCode() % " </b>" % endHtml % "[" % fitSet->getTimeStampOfLastFitResult() % "]</nobr>");
    if ( fitSet->getFitFinishCodeValue() <= 0 )
        fitFinishCodeVal = QString("<nobr><b>" % alertHtml % fitSet->getFitFinishCode() % " </b>" % endHtml % "[" % fitSet->getTimeStampOfLastFitResult() % "]</nobr>");

    const QString chiSquare("<nobr><b>&#935;<sup>2</sup>:</b></nobr>");
    const QString chiSquareVal("<nobr><b>" % QString::number(fitSet->getChiSquareAfterFit(), 'g', 4) % "</b> (" % QString::number(fitSet->getChiSquareOnStart(), 'g', 4) % " @ start)" % "</nobr>");

    const QString fitWeighting("<nobr><b>Fit-Weighting:</b></nobr>");
    const QString fitWeightingVal("<nobr><b>" % QString((!fitSet->getUsingYVariance()?QString("no weighting"):QString("variance"))) % "</b></nobr>");


    const QString iterations("<nobr><b>Iterations:</b></nobr>");
    QString iterationsVal = "";

    if ( fitSet->getNeededIterations() == fitSet->getMaximumIterations() || fitSet->getNeededIterations() == 0 )
        iterationsVal = QString("<nobr><b>" % alertHtml % QVariant((int)fitSet->getNeededIterations()).toString() % "/" % QVariant((int)fitSet->getMaximumIterations()).toString() % endHtml % "</b></nobr>");
    else
        iterationsVal = QString("<nobr><b>" % QVariant((int)fitSet->getNeededIterations()).toString() % "/" % QVariant((int)fitSet->getMaximumIterations()).toString() % "</b></nobr>");

    const QString binFac("<nobr>Bin-Factor:</nobr>");
    const QString binFacVal("<nobr><b>" % QVariant(dataStructure->getDataSetPtr()->getBinFactor()).toString() % " </b></nobr>");

    const QString channelRange("<nobr>ROI:</nobr>");
    const QString channelRangeVal("<nobr><b>" % QVariant(fitSet->getStartChannel()).toString() % ":" % QVariant(fitSet->getStopChannel()).toString() % "</b> [" % QVariant(PALSProjectManager::sharedInstance()->getMinChannel()).toString()  % ":" % QVariant(PALSProjectManager::sharedInstance()->getMaxChannel()).toString() % "]</nobr>");

    const QString channelResolution("<nobr>Channel-Resolution:</nobr>");
    const QString channelResolutionVal("<nobr><b>" % QVariant(fitSet->getChannelResolution()).toString() % " </b>ps</nobr>");

    const QString backgroundCounts("<nobr>Background:</nobr>");
    const QString backgroundCountsVal("<nobr><b>" % QString::number(fitSet->getBackgroundParamPtr()->getParameter()->getStartValue(), 'f', 4) % "</b></nobr>");

    const QString countsInRange("<nobr>Integral Counts in ROI:</nobr>");
    const QString countsInRangeVal("<nobr><b>" % QVariant(fitSet->getCountsInRange()).toString() % "</b></nobr>");

    const QString peakToBackgroundRatio("<nobr>Peak-to-Background Ratio:</nobr>");
    const QString peakToBackgroundRatioVal("<nobr><b>" % QString::number(fitSet->getPeakToBackgroundRation(), 'f', 3) % "</b></nobr>");

    const QString centerOfMass("<nobr>Center of Mass:</nobr>");
    const QString centerOfMassVal("<nobr><b>" %  QString::number(fitSet->getSpectralCentroid(), 'f', 4) % " </b>ps (estimated t<sub>0</sub>: <b>" % QString::number(fitSet->getT0SpectralCentroid(), 'f', 4) % "</b> ps) - ROI: [" % QVariant(fitSet->getStartChannel()).toString() % ":" % QVariant(fitSet->getStopChannel()).toString() % "]</nobr>");


    const QString fitParamCount("<nobr>Fit-Parameter Count:</nobr>");
    const QString fitParamCountVal("<nobr><b>" % QVariant(fitSet->getComponentsCount()+fitSet->getDeviceResolutionParamPtr()->getSize()).toString() % "</b></nobr>");


    const QString sumOfIntensities("<nobr>Sum of Component's Intensities:       </nobr>");
    QString sumOfIntensitiesVal;
    if ( fitSet->getSumOfIntensities() > 1.0f || fitSet->getSumOfIntensities() < 0.0f )
        sumOfIntensitiesVal = QString("<nobr><b>" % alertHtml % QString::number(fitSet->getSumOfIntensities(), 'f', 3) % endHtml % "</b></nobr>");
    else
        sumOfIntensitiesVal = QString("<nobr><b>" % notifyHtml % QString::number(fitSet->getSumOfIntensities(), 'f', 3) % endHtml % "</b></nobr>");


    const QString tauAverage("<nobr><b>&#964;<sub>average</sub>:</b></nobr>");
    const QString tauAverageVal("<nobr><b>( " % QString::number(fitSet->getAverageLifeTime(), 'f', 4) % " &plusmn; " % QString::number(fitSet->getAverageLifeTimeError(), 'f', 4) % " ) </b>ps</nobr>");


    QString resultString = "";

    resultString = resultString % tableStart;

    /*project-name:*/   resultString = resultString % startRow % startContent % projectName % finishContent % startContent % PALSProjectManager::sharedInstance()->getFileName() % finishContent % finishRow;
    /*ascii-file-name:*/   resultString = resultString % startRow % startContent % asciiFileName % finishContent % startContent % ((PALSProjectManager::sharedInstance()->getASCIIDataName()==QString("unknown"))?QString("unknown source"):PALSProjectManager::sharedInstance()->getASCIIDataName()) % finishContent % finishRow % lineBreak;

    /*finish code and time/date:*/resultString = resultString % startRow % startContent % fitFinishCode % finishContent % startContent % fitFinishCodeVal % finishContent % finishRow % lineBreak;

    /*chi-square:*/resultString = resultString % startRow % startContent % chiSquare % finishContent % startContent % chiSquareVal % finishContent % finishRow % lineBreak;
    /*fit-weighting:*/resultString = resultString % startRow % startContent % fitWeighting % finishContent % startContent % fitWeightingVal % finishContent % finishRow % lineBreak;

    /*iterations:*/resultString = resultString % startRow % startContent % iterations % finishContent % startContent % iterationsVal % finishContent % finishRow % lineBreak;


    /*Channel-Range:*/resultString = resultString % startRow % startContent % channelRange % finishContent % startContent % channelRangeVal % finishContent % finishRow;
    /*Channel-Resolution:*/resultString = resultString % startRow % startContent % channelResolution % finishContent % startContent % channelResolutionVal % finishContent % finishRow;
    /*Channel-Resolution:*/resultString = resultString % startRow % startContent % binFac % finishContent % startContent % binFacVal % finishContent % finishRow % lineBreak;

    /*Background-Counts:*/resultString = resultString % startRow % startContent % backgroundCounts % finishContent % startContent % backgroundCountsVal % finishContent % finishRow;
    /*Counts in Range:*/resultString = resultString % startRow % startContent % countsInRange % finishContent % startContent % countsInRangeVal % finishContent % finishRow;
    /*Peak-to-Background Ratio:*/resultString = resultString % startRow % startContent % peakToBackgroundRatio % finishContent % startContent % peakToBackgroundRatioVal % finishContent % finishRow % lineBreak;

    /*Center-of-Mass:*/resultString = resultString % startRow % startContent % centerOfMass % finishContent % startContent % centerOfMassVal % finishContent % finishRow % lineBreak;

    /*Fit-Parameter-Count:*/resultString = resultString % startRow % startContent % fitParamCount % finishContent % startContent % fitParamCountVal % finishContent % finishRow % lineBreak;


    /*Sum of Intensities:*/resultString = resultString % startRow % startContent % sumOfIntensities % finishContent % startContent % sumOfIntensitiesVal % finishContent % finishRow % lineBreak;

    /*tau average:*/resultString = resultString % startRow % startContent % tauAverage % finishContent % startContent % tauAverageVal % finishContent % finishRow % lineBreak;

    resultString = resultString % tableEnd % lineBreak;


    /*Sample-Components:*/
    resultString = resultString % "<nobr><b><big>" % "Sample-Components [" % QVariant(fitSet->getLifeTimeParamPtr()->getSize()).toString() % "/" % QVariant(fitSet->getComponentsCount()+fitSet->getDeviceResolutionParamPtr()->getSize()).toString() % "]</b></big></nobr>";

    resultString = resultString % tableBorderStart;

    //header:
    resultString = resultString % startRow % startHeader % spacer % "   name   " % spacer % endHeader % startHeader % spacer % "   fit-value   " % spacer % endHeader % startHeader % spacer % "   start-value   " % spacer % endHeader % startHeader % spacer %  "   lower-limit  \n reached? " % spacer % endHeader % startHeader % spacer %  "   upper-limit  \n reached? " % spacer % endHeader % startHeader % spacer %  "   fixed?   " % spacer % endHeader % finishRow;

    for ( int i = 0 ; i < fitSet->getLifeTimeParamPtr()->getSize() ; i += 2 )
    {
        const QString nameAndAliasTau("<nobr><b>" % spacer % QString(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getAlias()) % "</b> (" % QString(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getName()) % ")" % spacer %"</nobr>");
        const QString nameAndAliasIntensity("<nobr><b>"  % spacer % QString(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getAlias()) % "</b> (" % QString(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getName()) % ")" % spacer % "</nobr>");

        const QString tau("<nobr><b>" % spacer % "( " % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getFitValue(), 'f', 4) % " &plusmn; " % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getFitValueError(), 'f', 4) % " )</b> ps" % spacer % "</nobr>");
        const QString intensity("<nobr><b>" % spacer % "( " % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 4) % " &plusmn; " % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getFitValueError(), 'f', 4) % " )" % spacer % "</b></nobr>");

        const QString tauStart("<nobr>" % spacer % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getStartValue(), 'f', 4) % " ps" % spacer % "</nobr>");
        const QString IntensityStart("<nobr>" % spacer % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getStartValue(), 'f', 4) % spacer % "</nobr>");

        QString lowerLimitTau = "";
        if ( fitSet->getLifeTimeParamPtr()->getParameterAt(i)->isLowerBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getFitValue(), fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getLowerBoundingValue()) )
                lowerLimitTau = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getFitValue(), 'f', 4) % " ps)" % spacer % "</nobr>");
        }

        QString lowerLimitIntensity = "";
        if ( fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->isLowerBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getFitValue(), fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getLowerBoundingValue()) )
                lowerLimitIntensity = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 4) % ")" % spacer % "</nobr>");
        }

        QString upperLimitTau = "";
        if ( fitSet->getLifeTimeParamPtr()->getParameterAt(i)->isUpperBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getFitValue(), fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getUpperBoundingValue()) )
                upperLimitTau = QString("<nobr><b>" % alertHtml % spacer %  "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i)->getFitValue(), 'f', 4) % " ps)" % spacer % "</nobr>");
        }

        QString upperLimitIntensity = "";
        if ( fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->isUpperBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getFitValue(), fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getUpperBoundingValue()) )
                upperLimitIntensity = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 3) % ")" % spacer % "</nobr>");
        }

        QString fixedTau = "";
        if ( fitSet->getLifeTimeParamPtr()->getParameterAt(i)->isFixed() )
            fixedTau = QString("<nobr><b>" % infoHtml % spacer % "&#8226;" % spacer % endHtml % "</b></nobr>");

        QString fixedIntensity = "";
        if ( fitSet->getLifeTimeParamPtr()->getParameterAt(i+1)->isFixed() )
            fixedIntensity = QString("<nobr><b>" % infoHtml % spacer % "&#8226;" % spacer % endHtml % "</b></nobr>");


        const QString startContentAligned = startContent % alignCenterStart;
        const QString finishContentAligned = alignCenterEnd % finishContent;

        resultString = resultString % startRow % startContentAligned % nameAndAliasTau % finishContentAligned % startContentAligned % tau % finishContentAligned % startContentAligned % tauStart % finishContentAligned % startContentAligned % lowerLimitTau % finishContentAligned % startContentAligned % upperLimitTau % finishContentAligned % startContentAligned % fixedTau % finishContentAligned % finishRow;
        resultString = resultString % startRow % startContentAligned % nameAndAliasIntensity % finishContentAligned % startContentAligned % intensity % finishContentAligned % startContentAligned % IntensityStart % finishContentAligned % startContentAligned % lowerLimitIntensity % finishContentAligned % startContentAligned % upperLimitIntensity % finishContentAligned % startContentAligned % fixedIntensity % finishContentAligned % finishRow;
    }

    resultString = resultString % tableBorderEnd;
    resultString = resultString % lineBreak % lineBreak;



    /*Source-Components:*/
    resultString = resultString % "<nobr><b><big>" % "Source-Components [" % QVariant(fitSet->getSourceParamPtr()->getSize()).toString() % "/" % QVariant(fitSet->getComponentsCount()+fitSet->getDeviceResolutionParamPtr()->getSize()).toString() % "]</b></big></nobr>";

    resultString = resultString % tableBorderStart;

    //header:
    resultString = resultString % startRow % startHeader % spacer % "   name   " % spacer % endHeader % startHeader % spacer % "   fit-value   " % spacer % endHeader % startHeader % spacer % "   start-value   " % spacer % endHeader % startHeader % spacer %  "   lower-limit  \n reached? " % spacer % endHeader % startHeader % spacer %  "   upper-limit  \n reached? " % spacer % endHeader % startHeader % spacer %  "   fixed?   " % spacer % endHeader % finishRow;

    for ( int i = 0 ; i < fitSet->getSourceParamPtr()->getSize() ; i += 2 )
    {
        const QString nameAndAliasTau("<nobr><b>" % spacer % QString(fitSet->getSourceParamPtr()->getParameterAt(i)->getAlias()) % "</b> (" % QString(fitSet->getSourceParamPtr()->getParameterAt(i)->getName()) % ")" % spacer %"</nobr>");
        const QString nameAndAliasIntensity("<nobr><b>"  % spacer % QString(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getAlias()) % "</b> (" % QString(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getName()) % ")" % spacer % "</nobr>");

        const QString tau("<nobr><b>" % spacer % "( " % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i)->getFitValue(), 'f', 4) % " &plusmn; " % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i)->getFitValueError(), 'f', 4) % " )</b> ps" % spacer % "</nobr>");
        const QString intensity("<nobr><b>" % spacer % "( " % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 4) % " &plusmn; " % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getFitValueError(), 'f', 4) % " )" % spacer % "</b></nobr>");

        const QString tauStart("<nobr>" % spacer % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i)->getStartValue(), 'f', 4) % " ps" % spacer % "</nobr>");
        const QString IntensityStart("<nobr>" % spacer % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getStartValue(), 'f', 4) % spacer % "</nobr>");

        QString lowerLimitTau = "";
        if ( fitSet->getSourceParamPtr()->getParameterAt(i)->isLowerBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getSourceParamPtr()->getParameterAt(i)->getFitValue(), fitSet->getSourceParamPtr()->getParameterAt(i)->getLowerBoundingValue()) )
                lowerLimitTau = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i)->getFitValue(), 'f', 3) % " ps)" % spacer % "</nobr>");
        }

        QString lowerLimitIntensity = "";
        if ( fitSet->getSourceParamPtr()->getParameterAt(i+1)->isLowerBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getFitValue(), fitSet->getSourceParamPtr()->getParameterAt(i+1)->getLowerBoundingValue()) )
                lowerLimitIntensity = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 4) % ")" % spacer % "</nobr>");
        }

        QString upperLimitTau = "";
        if ( fitSet->getSourceParamPtr()->getParameterAt(i)->isUpperBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getSourceParamPtr()->getParameterAt(i)->getFitValue(), fitSet->getSourceParamPtr()->getParameterAt(i)->getUpperBoundingValue()) )
                upperLimitTau = QString("<nobr><b>" % alertHtml % spacer %  "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i)->getFitValue(), 'f', 4) % " ps)" % spacer % "</nobr>");
        }

        QString upperLimitIntensity = "";
        if ( fitSet->getSourceParamPtr()->getParameterAt(i+1)->isUpperBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getFitValue(), fitSet->getSourceParamPtr()->getParameterAt(i+1)->getUpperBoundingValue()) )
                upperLimitIntensity = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getSourceParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 4) % ")" % spacer % "</nobr>");
        }

        QString fixedTau = "";
        if ( fitSet->getSourceParamPtr()->getParameterAt(i)->isFixed() )
            fixedTau = QString("<nobr><b>" % infoHtml % spacer % "&#8226;" % spacer % endHtml % "</b></nobr>");

        QString fixedIntensity = "";
        if ( fitSet->getSourceParamPtr()->getParameterAt(i+1)->isFixed() )
            fixedIntensity = QString("<nobr><b>" % infoHtml % spacer % "&#8226;" % spacer % endHtml % "</b></nobr>");


        const QString startContentAligned = startContent % alignCenterStart;
        const QString finishContentAligned = alignCenterEnd % finishContent;

        resultString = resultString % startRow % startContentAligned % nameAndAliasTau % finishContentAligned % startContentAligned % tau % finishContentAligned % startContentAligned % tauStart % finishContentAligned % startContentAligned % lowerLimitTau % finishContentAligned % startContentAligned % upperLimitTau % finishContentAligned % startContentAligned % fixedTau % finishContentAligned % finishRow;
        resultString = resultString % startRow % startContentAligned % nameAndAliasIntensity % finishContentAligned % startContentAligned % intensity % finishContentAligned % startContentAligned % IntensityStart % finishContentAligned % startContentAligned % lowerLimitIntensity % finishContentAligned % startContentAligned % upperLimitIntensity % finishContentAligned % startContentAligned % fixedIntensity % finishContentAligned % finishRow;
    }

    resultString = resultString % tableBorderEnd;
    resultString = resultString % lineBreak % lineBreak;


    /*IRF-Components:*/
    resultString = resultString % "<nobr><b><big>" % "IRF-Components [" % "2/" % QVariant(fitSet->getComponentsCount()+fitSet->getDeviceResolutionParamPtr()->getSize()).toString() % "]</b></big></nobr>";

    resultString = resultString % tableBorderStart;

    //header:
    resultString = resultString % startRow % startHeader % spacer % "   name   " % spacer % endHeader % startHeader % spacer % "   fit-value   " % spacer % endHeader % startHeader % spacer % "   start-value   " % spacer % endHeader % startHeader % spacer %  "   lower-limit  \n reached? " % spacer % endHeader % startHeader % spacer %  "   upper-limit  \n reached? " % spacer % endHeader % startHeader % spacer %  "   fixed?   " % spacer % endHeader % finishRow;

    for ( int i = 0 ; i < fitSet->getDeviceResolutionParamPtr()->getSize() ; i += 3 )
    {
        const QString nameAndAliasSigma("<nobr><b>" % spacer % QString(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getAlias()) % "</b> (" % QString(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getName()) % ")" % spacer %"</nobr>");
        const QString nameAndAliasMu("<nobr><b>"  % spacer % QString(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getAlias()) % "</b> (" % QString(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getName()) % ")" % spacer % "</nobr>");
        const QString nameAndAliasIGauss("<nobr><b>"  % spacer % QString(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getAlias()) % "</b> (" % QString(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getName()) % ")" % spacer % "</nobr>");

        const QString sigma("<nobr><b>" % spacer % "( " % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getFitValue(), 'f', 4) % " &plusmn; " % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getFitValueError(), 'f', 4) % " )</b> ps" % spacer % "</nobr>");
        const QString mu("<nobr><b>" % spacer % "( " % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 4) % " &plusmn; " % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getFitValueError(), 'f', 4) % " )</b> ps" % spacer % "</nobr>");
        const QString IGauss("<nobr><b>" % spacer % "( " % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getFitValue(), 'f', 4) % " &plusmn; " % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getFitValueError(), 'f', 4) % " )</b>" % spacer % "</nobr>");

        const QString sigmaStart("<nobr>" % spacer % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getStartValue(), 'f', 4) % " ps" % spacer % "</nobr>");
        const QString muStart("<nobr>" % spacer % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getStartValue(), 'f', 4) % " ps" % spacer % "</nobr>");
        const QString IGaussStart("<nobr>" % spacer % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getStartValue(), 'f', 4) % "" % spacer % "</nobr>");

        QString lowerLimitSigma = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->isLowerBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getFitValue(), fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getLowerBoundingValue()) )
                lowerLimitSigma = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getFitValue(), 'f', 4) % " ps)" % spacer % "</nobr>");
        }

        QString lowerLimitMu = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->isLowerBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getFitValue(), fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getLowerBoundingValue()) )
                lowerLimitMu = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 4) % " ps)" % spacer % "</nobr>");
        }

        QString lowerLimitIGauss = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->isLowerBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getFitValue(), fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getLowerBoundingValue()) )
                lowerLimitIGauss = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getFitValue(), 'f', 4) % ")" % spacer % "</nobr>");
        }

        QString upperLimitSigma = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->isUpperBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getFitValue(), fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getUpperBoundingValue()) )
                upperLimitSigma = QString("<nobr><b>" % alertHtml % spacer %  "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->getFitValue(), 'f', 4) % " ps)" % spacer % "</nobr>");
        }

        QString upperLimitMu = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->isUpperBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getFitValue(), fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getUpperBoundingValue()) )
                upperLimitMu = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->getFitValue(), 'f', 4) % " ps)" % spacer % "</nobr>");
        }

        QString upperLimitIGauss = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->isUpperBoundingEnabled() )
        {
            if ( qFuzzyCompare(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getFitValue(), fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getUpperBoundingValue()) )
                upperLimitIGauss = QString("<nobr><b>" % alertHtml % spacer % "&#8226;" % endHtml % "</b>" % "   (" % QString::number(fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->getFitValue(), 'f', 4) % ")" % spacer % "</nobr>");
        }

        QString fixedSigma = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->isFixed() )
            fixedSigma = QString("<nobr><b>" % infoHtml % spacer % "&#8226;" % spacer % endHtml % "</b></nobr>");

        QString fixedMu = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->isFixed() )
            fixedMu = QString("<nobr><b>" % infoHtml % spacer % "&#8226;" % spacer % endHtml % "</b></nobr>");

        QString fixedIGauss = "";
        if ( fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->isFixed() )
            fixedIGauss = QString("<nobr><b>" % infoHtml % spacer % "&#8226;" % spacer % endHtml % "</b></nobr>");


        const QString startContentAligned = startContent % alignCenterStart;
        const QString finishContentAligned = alignCenterEnd % finishContent;

        resultString = resultString % startRow % startContentAligned % nameAndAliasSigma % finishContentAligned % startContentAligned % sigma % finishContentAligned % startContentAligned % sigmaStart % finishContentAligned % startContentAligned % lowerLimitSigma % finishContentAligned % startContentAligned % upperLimitSigma % finishContentAligned % startContentAligned % fixedSigma % finishContentAligned % finishRow;
        resultString = resultString % startRow % startContentAligned % nameAndAliasMu % finishContentAligned % startContentAligned % mu % finishContentAligned % startContentAligned % muStart % finishContentAligned % startContentAligned % lowerLimitMu % finishContentAligned % startContentAligned % upperLimitMu % finishContentAligned % startContentAligned % fixedMu % finishContentAligned % finishRow;
        resultString = resultString % startRow % startContentAligned % nameAndAliasIGauss % finishContentAligned % startContentAligned % IGauss % finishContentAligned % startContentAligned % IGaussStart % finishContentAligned % startContentAligned % lowerLimitIGauss % finishContentAligned % startContentAligned % upperLimitIGauss % finishContentAligned % startContentAligned % fixedIGauss % finishContentAligned % finishRow;
    }

    resultString = resultString % tableBorderEnd;

    PALSResult *result = new PALSResult(fitSet->getResultHistoriePtr());
    DUNUSED_PARAM(result);

    result->setResultText(resultString);
}
