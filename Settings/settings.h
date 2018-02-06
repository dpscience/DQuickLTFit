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

#ifndef PALSSETTINGS_H
#define PALSSETTINGS_H

#include <QDateTime>
#include <QBuffer>

#include "../DLib/DLib.h"

#define SETTINGS_READ              /*Reading and Loading*/
#define SETTINGS_WRITE             /*Writing and Saving */
#define LOAD_CONSTRUCTOR      /*Load Constructor*/
#define SAVE_CONSTRUCTOR       /*Save Constructor*/
#define CREATE_FROM_OUTSIDE  /*default constructor*/

///File-Hierarchy:
class PALSProject;
    class PALSDataStructure;
        class PALSDataSet;
        class PALSFitSet;
            class PALSSourceParameter;
                class PALSFitParameter;
            class PALSDeviceResolutionParameter;
                class PALSFitParameter;
            class PALSBackgroundParameter;
                class PALSFitParameter;
            class PALSLifeTimeParameter;
                class PALSFitParameter;
            class PALSResultHistorie;
                class PALSResult;

class PALSProject
{
    DSimpleXMLNode *m_rootNode;
    DSimpleXMLNode *m_parentNode;
    DSimpleXMLNode *m_lastSaveTimeNode;
    DSimpleXMLNode *m_projectNameNode;

    QList<PALSDataStructure*> m_dataStructureList;

public:
    PALSProject();

    virtual ~PALSProject();

    bool save(const DString& projectPath);
    bool load(const DString& projectPath);

    unsigned int getSize() const;
    DSimpleXMLNode *getParent() const;

SETTINGS_WRITE
    void setLastSaveTime(const QDateTime& dateTime);
    void setName(const DString& name);

    void addDataStructure(PALSDataStructure* dataStructure);
    void removeDataStructure(unsigned int index);
    void removeDataStructure(PALSDataStructure* dataStructure);

SETTINGS_READ
    PALSDataStructure* operator[](unsigned int index) const;
    PALSDataStructure* getDataStructureAt(unsigned int index) const;
    QDateTime getLastSaveTime() const;
    DString getName() const;

private:
    void clear();
};


class PALSDataStructure
{
    DSimpleXMLNode *m_parentNode;
    DSimpleXMLNode *m_nameNode;

    PALSDataSet *m_dataSet;
    PALSFitSet *m_fitSet;

public:
    CREATE_FROM_OUTSIDE SAVE_CONSTRUCTOR PALSDataStructure(PALSProject *parent);
    LOAD_CONSTRUCTOR PALSDataStructure(PALSProject* parent, const DSimpleXMLTag& tag, const DString& name);

    virtual ~PALSDataStructure();

    DSimpleXMLNode *getParent() const;

    PALSDataSet* getDataSetPtr() const;
    PALSFitSet* getFitSetPtr() const;

SETTINGS_WRITE
    void setName(const DString& name);

SETTINGS_READ
    DString getName() const;
};


class PALSDataSet
{
    DSimpleXMLNode *m_parentNode;
    DSimpleXMLNode *m_xyDataNode;
    DSimpleXMLNode *m_fitDataNode;
    DSimpleXMLNode *m_residualNode;
    DSimpleXMLNode *m_colorResidualsNode;
    DSimpleXMLNode *m_colorDataNode;

    QList<QPointF> m_xyData;
    QList<QPointF> m_fitData;
    QList<QPointF> m_residualData;

public:
    SAVE_CONSTRUCTOR PALSDataSet(PALSDataStructure *parent);
    LOAD_CONSTRUCTOR PALSDataSet(PALSDataStructure *parent, const DSimpleXMLTag& tag);

    virtual ~PALSDataSet();

    DSimpleXMLNode* getParent() const;

    void clearFitData();
    void clearResidualData();

SETTINGS_WRITE
    void setLifeTimeData(const QList<QPointF>& dataSet);
    void setFitData(const QList<QPointF>& dataSet);
    void setResiduals(const QList<QPointF>& residuals);
    void setLifeTimeDataColor(const DColor& color);
    void setResidualsColor(const DColor& color);

SETTINGS_READ
    QList<QPointF> getLifeTimeData() const;
    QList<QPointF> getFitData() const;
    QList<QPointF> getResiduals() const;

    DColor getLifeTimeDataColor() const;
    DColor getResidualsColor() const;
};


class PALSFitSet
{
    DSimpleXMLNode *m_parentNode;

    DSimpleXMLNode *m_maxIterationsNode;
    DSimpleXMLNode *m_neededIterationsNode;

    DSimpleXMLNode *m_usingYVarianceNode;

    DSimpleXMLNode *m_chiSquareOnStart;
    DSimpleXMLNode *m_chiSquareAfterFit;

    DSimpleXMLNode *m_channelResolutionNode;

    DSimpleXMLNode *m_startChannelNode;
    DSimpleXMLNode *m_stopChannelNode;

    DSimpleXMLNode *m_dateTimeOfLastFitResultsNode;
    DSimpleXMLNode *m_fitFinishCodeNode;
    DSimpleXMLNode *m_fitFinishCodeValueNode;

    DSimpleXMLNode *m_countsInRangeNode;

    DSimpleXMLNode *m_averageLifeTimeNode;
    DSimpleXMLNode *m_averageLifeTimeErrorNode;

    DSimpleXMLNode *m_peakToBackgroundRatioNode;
    DSimpleXMLNode *m_sumOfIntensitiesNode;

    DSimpleXMLNode *m_dataPlotImageNode;
    DSimpleXMLNode *m_residualPlotImageNode;


    PALSSourceParameter *m_sourceParams;
    PALSDeviceResolutionParameter *m_deviceResolutionParams;
    PALSBackgroundParameter *m_bgParam;
    PALSLifeTimeParameter *m_lifeTimeParams;

    PALSResultHistorie *m_resultHistorie;

public:
    SAVE_CONSTRUCTOR PALSFitSet(PALSDataStructure *parent);
    LOAD_CONSTRUCTOR PALSFitSet(PALSDataStructure *parent, const DSimpleXMLTag& tag);

    virtual ~PALSFitSet();

    DSimpleXMLNode *getParent() const;

    PALSSourceParameter *getSourceParamPtr() const;
    PALSDeviceResolutionParameter *getDeviceResolutionParamPtr() const;
    PALSLifeTimeParameter *getLifeTimeParamPtr() const;
    PALSBackgroundParameter *getBackgroundParamPtr() const;

    PALSResultHistorie *getResultHistoriePtr() const;

    int getComponentsCount() const; //without background and device

SETTINGS_WRITE
    void setMaximumIterations(unsigned int iterations);
    void setNeededIterations(unsigned int iterations);
    void setUsingYVariance(bool usingYVar);
    void setChiSquareOnStart(double value);
    void setChiSquareAfterFit(double value);
    void setChannelResolution(double resolution);
    void setStartChannel(int startChannel);
    void setStopChannel(int stopChannel);
    void setAverageLifeTime(double avgLifeTime);
    void setAverageLifeTimeError(double avgLtError);
    void setCountsInRange(int countsInRange);
    void setFitFinishCode(const QString& finishCode);
    void setFitFinishCodeValue(int value);
    void setTimeStampOfLastFitResult(const QString& timeStamp);
    void setPeakToBackgroundRatio(double ratio);
    void setSumOfIntensities(double sum);
    void setDataPlotImage(const QImage& image);
    void setResidualPlotImage(const QImage& image);

SETTINGS_READ
    unsigned int getMaximumIterations() const;
    unsigned int getNeededIterations() const;
    bool usingYVariance() const;
    double getChiSquareOnStart() const;
    double getChiSquareAfterFit() const;
    bool getUsingYVariance() const;
    double getChannelResolution() const;
    int getStartChannel() const;
    int getStopChannel() const;
    double getAverageLifeTime() const;
    double getAverageLifeTimeError() const;
    int getCountsInRange() const;
    QString getFitFinishCode() const;
    int getFitFinishCodeValue() const;
    QString getTimeStampOfLastFitResult() const;
    double getPeakToBackgroundRation() const;
    double getSumOfIntensities() const;
    QImage getDataPlotImage() const;
    QImage getResidualPlotImage() const;
};

class PALSResultHistorie
{
    DSimpleXMLNode *m_parentNode;

    QList<PALSResult*> m_resultList;

public:
    SAVE_CONSTRUCTOR PALSResultHistorie(PALSFitSet *parent);
    LOAD_CONSTRUCTOR PALSResultHistorie(PALSFitSet *parent, const DSimpleXMLTag& tag);

    virtual ~PALSResultHistorie();

    DSimpleXMLNode *getParent() const;
    unsigned int getSize() const;

SETTINGS_WRITE
    void addResult(PALSResult *result);
    void removeResult(unsigned int index);

SETTINGS_READ
    PALSResult* operator[](unsigned int index) const;
    PALSResult* getResultAt(unsigned int index) const;
};


class PALSSourceParameter
{
    DSimpleXMLNode *m_parentNode;

    QList<PALSFitParameter* > m_fitParamSettings;

public:
    SAVE_CONSTRUCTOR PALSSourceParameter(PALSFitSet *parent);
    LOAD_CONSTRUCTOR PALSSourceParameter(PALSFitSet *parent, const DSimpleXMLTag& tag);

    virtual ~PALSSourceParameter();

    DSimpleXMLNode* getParent() const;
    unsigned int getSize() const;

SETTINGS_WRITE
    void addParameter(PALSFitParameter *param);
    void removeParameter(unsigned int index);
    void removeParameter(PALSFitParameter *param);

SETTINGS_READ
    PALSFitParameter* operator[](unsigned int index) const;
    PALSFitParameter* getParameterAt(unsigned int index) const;
};


class PALSDeviceResolutionParameter
{
    DSimpleXMLNode *m_parentNode;

    QList<PALSFitParameter* > m_fitParamSettings;

public:
    SAVE_CONSTRUCTOR PALSDeviceResolutionParameter(PALSFitSet *parent);
    LOAD_CONSTRUCTOR PALSDeviceResolutionParameter(PALSFitSet *parent, const DSimpleXMLTag& tag);

    virtual ~PALSDeviceResolutionParameter();

    DSimpleXMLNode* getParent() const;
    unsigned int getSize() const;

SETTINGS_WRITE
    void addParameter(PALSFitParameter *param);
    void removeParameter(unsigned int index);
    void removeParameter(PALSFitParameter *param);

SETTINGS_READ
    PALSFitParameter* operator[](unsigned int index) const;
    PALSFitParameter* getParameterAt(unsigned int index) const;
};


class PALSBackgroundParameter
{
    DSimpleXMLNode *m_parentNode;

    PALSFitParameter *m_fitParamSettings;

public:
    SAVE_CONSTRUCTOR PALSBackgroundParameter(PALSFitSet *parent);
    LOAD_CONSTRUCTOR PALSBackgroundParameter(PALSFitSet *parent, const DSimpleXMLTag& tag);

    virtual ~PALSBackgroundParameter();

    DSimpleXMLNode* getParent() const;

SETTINGS_WRITE
    void addParameter(PALSFitParameter *param);

SETTINGS_READ
    PALSFitParameter* getParameter() const;
};


class PALSLifeTimeParameter
{
    DSimpleXMLNode *m_parentNode;

    QList<PALSFitParameter* > m_fitParamSettings;

public:
    SAVE_CONSTRUCTOR PALSLifeTimeParameter(PALSFitSet *parent);
    LOAD_CONSTRUCTOR PALSLifeTimeParameter(PALSFitSet *parent, const DSimpleXMLTag& tag);

    virtual ~PALSLifeTimeParameter();

    DSimpleXMLNode* getParent() const;
    unsigned int  getSize() const;

SETTINGS_WRITE
    void addParameter(PALSFitParameter *param);
    void removeParameter(unsigned int index);
    void removeParameter(PALSFitParameter *param);

SETTINGS_READ
    PALSFitParameter* operator[](unsigned int index) const;
    PALSFitParameter* getParameterAt(unsigned int index) const;
};


class PALSResult
{
    DSimpleXMLNode *m_parentNode;

    DSimpleXMLNode *m_htmlTextNode;

public:
    SAVE_CONSTRUCTOR PALSResult(PALSResultHistorie *parent);
    LOAD_CONSTRUCTOR PALSResult(PALSResultHistorie *parent, const DSimpleXMLTag &tag, const DString &name);

    virtual ~PALSResult();

    DSimpleXMLNode* getParent() const;

SETTINGS_WRITE
    void setResultText(const QString& resultText);

SETTINGS_READ
    QString getResultText() const;
};

class PALSFitParameter
{
    DSimpleXMLNode *m_parentNode;

    DSimpleXMLNode *m_active;
    DSimpleXMLNode *m_name;
    DSimpleXMLNode *m_alias;
    DSimpleXMLNode *m_startValue;
    DSimpleXMLNode *m_upperBounding;
    DSimpleXMLNode *m_bupperBounding;
    DSimpleXMLNode *m_lowerBounding;
    DSimpleXMLNode *m_blowerBounding;

    DSimpleXMLNode *m_fixed;

    DSimpleXMLNode *m_fitValue;
    DSimpleXMLNode *m_fitValueError;

public:
    SAVE_CONSTRUCTOR PALSFitParameter(PALSSourceParameter *parent);
    SAVE_CONSTRUCTOR PALSFitParameter(PALSDeviceResolutionParameter *parent);
    SAVE_CONSTRUCTOR PALSFitParameter(PALSBackgroundParameter *parent);
    SAVE_CONSTRUCTOR PALSFitParameter(PALSLifeTimeParameter *parent);
    LOAD_CONSTRUCTOR PALSFitParameter(PALSSourceParameter *parent, const DSimpleXMLTag &tag, const DString &name);
    LOAD_CONSTRUCTOR PALSFitParameter(PALSDeviceResolutionParameter* parent, const DSimpleXMLTag &tag, const DString &name);
    LOAD_CONSTRUCTOR PALSFitParameter(PALSBackgroundParameter *parent, const DSimpleXMLTag &tag, const DString &name);
    LOAD_CONSTRUCTOR PALSFitParameter(PALSLifeTimeParameter *parent, const DSimpleXMLTag &tag, const DString &name);

    virtual ~PALSFitParameter();

    DSimpleXMLNode* getParent() const;

SETTINGS_WRITE
    void setActive(bool active);
    void setName(const DString& name);
    void setAlias(const DString& name);
    void setStartValue(double value);
    void setUpperBoundingValue(double value);
    void setUpperBoundingEnabled(bool enabled);
    void setLowerBoundingValue(double value);
    void setLowerBoundingEnabled(bool enabled);
    void setAsFixed(bool fixed);
    void setFitValue(double value);
    void setFitValueError(double error);

SETTINGS_READ
    bool isActive() const;
    DString getName() const;
    DString getAlias() const;
    double  getStartValue() const;
    bool isUpperBoundingEnabled() const;
    double getUpperBoundingValue() const;
    double getLowerBoundingValue() const;
    bool isLowerBoundingEnabled() const;
    bool isFixed() const;
    double getFitValue() const;
    double getFitValueError() const;
};

#endif // PALSSETTINGS_H
