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

#include "ltparameterlistview.h"
#include "ui_ltparameterlistview.h"

#define MIN_COMPONENTS      2 //tau and I is separate
#define MAX_COMPONENTS   24 //--

#define WINDOWS_FONT(__pointSize__)  QFont("Arial", __pointSize__)

ParameterListView::ParameterListView(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ParameterListView),
    m_fitSet(nullptr)
{
    ui->setupUi(this);

    ui->tableWidget_Source->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidget_Sample->setSelectionBehavior(QAbstractItemView::SelectRows);

    QStringList headerList_1;
    headerList_1.append(""); //alias-placeholder
    headerList_1.append("Description/Alias");
    headerList_1.append("Start-Value");
    headerList_1.append("Lower-Limit");
    headerList_1.append("Upper-Limit");
    headerList_1.append("Fixed?");

    QStringList headerList_2;
    headerList_2.append("Start-Value");
    headerList_2.append("Lower-Limit");
    headerList_2.append("Upper-Limit");
    headerList_2.append("Fixed?");

    QStringList headerList_3;
    headerList_3.append(""); //alias-placeholder
    headerList_3.append("Description/Alias");
    headerList_3.append("Start-Value");
    headerList_3.append("Lower-Limit");
    headerList_3.append("Upper-Limit");
    headerList_3.append("Fixed?");

    ui->tableWidget_Source->setColumnCount(6);
    ui->tableWidget_Sample->setColumnCount(6);
    ui->tableWidget_Device->setColumnCount(6);

    ui->tableWidget_Device->verticalHeader()->setVisible(false);

    ui->tableWidget_Source->setHorizontalHeaderLabels(headerList_1);
    ui->tableWidget_Sample->setHorizontalHeaderLabels(headerList_1);
    ui->tableWidget_Device->setHorizontalHeaderLabels(headerList_3);

    ui->tableWidget_Source->setFrameStyle(QFrame::NoFrame);
    ui->tableWidget_Sample->setFrameStyle(QFrame::NoFrame);
    ui->tableWidget_Device->setFrameStyle(QFrame::NoFrame);

    ui->tableWidget_Source->setCornerButtonEnabled(false);
    ui->tableWidget_Sample->setCornerButtonEnabled(false);
    ui->tableWidget_Device->setCornerButtonEnabled(false);

    ui->tableWidget_Sample->verticalHeader()->setVisible(false);
    ui->tableWidget_Source->verticalHeader()->setVisible(false);
    ui->tableWidget_Device->verticalHeader()->setVisible(false);

    QFont widgetFont("Helvetica", 12);
    widgetFont.setBold(true);

    ui->tableWidget_Source->setFont(widgetFont);
    ui->tableWidget_Sample->setFont(widgetFont);
    ui->tableWidget_Device->setFont(widgetFont);

    ui->spinBox_backgroundChannel->setRange(2, 20000);

    ui->checkBox_FirstChannelBkgrd->setChecked(false);

    ui->pushButtonAdd_Sample->setLiteralSVG(":/localImages/Images/add");
    ui->pushButtonRemove_Sample->setLiteralSVG(":/localImages/Images/remove");

    ui->pushButtonAdd_Source->setLiteralSVG(":/localImages/Images/add");
    ui->pushButtonRemove_Source->setLiteralSVG(":/localImages/Images/remove");

    ui->pushButtonAdd_Device->setLiteralSVG(":/localImages/Images/add");
    ui->pushButtonRemove_Device->setLiteralSVG(":/localImages/Images/remove");

    ui->pushButton_Background->setLiteralSVG(":/localImages/Images/arrowRight");

    connect(ui->pushButtonAdd_Source, SIGNAL(clicked()), SLOT(addSourceComponent()));
    connect(ui->pushButtonRemove_Source, SIGNAL(clicked()), SLOT(removeSourceComponent()));
    connect(ui->pushButtonAdd_Sample, SIGNAL(clicked()), SLOT(addSampleComponent()));
    connect(ui->pushButtonRemove_Sample, SIGNAL(clicked()), SLOT(removeSampleComponent()));
    connect(ui->pushButtonAdd_Device, SIGNAL(clicked()), SLOT(addDeviceResolutionComponent()));
    connect(ui->pushButtonRemove_Device, SIGNAL(clicked()), SLOT(removeDeviceResolutionComponent()));

    connect(ui->pushButton_Background, SIGNAL(clicked()), this, SLOT(updateBackgroundValue()));

    connect(ui->tableWidget_Source, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));
    connect(ui->tableWidget_Sample, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));
    connect(ui->tableWidget_Device, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    connect(ui->spinBox_backgroundChannel, SIGNAL(valueChanged(int)), this, SLOT(saveBackgroundChannelRanges(int)));
    connect(ui->checkBox_FirstChannelBkgrd, SIGNAL(clicked()), this, SLOT(setUsingFirstChannelsForBkgrdCalc()));

    ui->pushButtonAdd_Sample->setToolTip("Add Component...");
    ui->pushButtonAdd_Source->setToolTip("Add Component...");
    ui->pushButtonAdd_Device->setToolTip("Add Component...");
    ui->pushButtonRemove_Sample->setToolTip("Remove Selected Component...");
    ui->pushButtonRemove_Source->setToolTip("Remove Selected Component...");
    ui->pushButtonRemove_Device->setToolTip("Remove Selected Component...");

    ui->spinBox_backgroundChannel->setToolTip("<nobr>Select the Channel Count of the upper ROI to calculate the Background</nobr>");
    ui->pushButton_Background->setToolTip("<nobr>Calculate the Background from the selected Channel Count of ROI</nobr>");
    ui->checkBox_FirstChannelBkgrd->setToolTip("<nobr>Using first Channels of ROI for Background Calculation?</nobr>");

    ui->spinBox_iterations->setToolTip("<nobr>Maximum Count of Iterations used to converge in &#967;<sup>2</sup></nobr>");
    ui->doubleSpinBox_channelResolution->setToolTip("<nobr>Type here the Channel Resolution [ps]</nobr>");

    ui->doubleSpinBox_background->setToolTip("<nobr>Type here the Background Counts or calculate it.</nobr>");

    ui->widget->setToolTip("<nobr>Select the Region of Interest (ROI).<br>Data outside ROI will be ignored by the Fit.</nobr>");

#if defined(Q_OS_WIN)
    ui->labelBGCounts->setFont(WINDOWS_FONT(9));
    ui->labelChnResolution->setFont(WINDOWS_FONT(9));
    ui->labelLstChnOfROI->setFont(WINDOWS_FONT(9));
    ui->labelMaxIterations->setFont(WINDOWS_FONT(9));
    ui->checkBox_FirstChannelBkgrd->setFont(WINDOWS_FONT(9));

    ui->groupBox->setFont(WINDOWS_FONT(11));
    ui->groupBox_2->setFont(WINDOWS_FONT(11));
    ui->groupBox_3->setFont(WINDOWS_FONT(11));

    ui->doubleSpinBox_background->setFont(WINDOWS_FONT(9));
    ui->doubleSpinBox_channelResolution->setFont(WINDOWS_FONT(9));
    ui->spinBox_backgroundChannel->setFont(WINDOWS_FONT(9));
    ui->spinBox_iterations->setFont(WINDOWS_FONT(9));

    ui->tabWidget->setFont(WINDOWS_FONT(9));
#endif
}

ParameterListView::~ParameterListView()
{
    while ( m_sourceWidgetCollection.size() > 0  )
    {
        PALSSourceTableWidgetItemCollector *item = m_sourceWidgetCollection.takeFirst();
        DDELETE_SAFETY(item);
    }

    while ( m_sampleWidgetCollection.size() > 0 )
    {
        PALSSampleTableWidgetItemCollector *item = m_sampleWidgetCollection.takeFirst();
        DDELETE_SAFETY(item);
    }

    while ( m_deviceWidgetCollection.size() > 0 )
    {
        PALSDeviceTableWidgetItemCollector *item = m_deviceWidgetCollection.takeFirst();
        DDELETE_SAFETY(item);
    }

    DDELETE_SAFETY(m_fitSet);
    DDELETE_SAFETY(ui);
}

void ParameterListView::updateParamterList()
{
    if ( !PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr() )
        return;


    disconnect(ui->widget, SIGNAL(rangeChanged(double, double)), this, SLOT(updateChannelRange(double, double)));

    m_fitSet = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr();

    disconnect(ui->tableWidget_Source, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));
    disconnect(ui->tableWidget_Sample, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));
    disconnect(ui->tableWidget_Device, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    while ( ui->tableWidget_Source->rowCount() > 0 )
        ui->tableWidget_Source->removeRow(0);

    while ( ui->tableWidget_Sample->rowCount() > 0 )
        ui->tableWidget_Sample->removeRow(0);

    while ( ui->tableWidget_Device->rowCount() > 0 )
        ui->tableWidget_Device->removeRow(0);

    initializeSourceTableWidget();
    initializeSampleTableWidget();
    initializeDeviceTableWidget();

    ui->doubleSpinBox_channelResolution->setRange(0, 2000);
    ui->doubleSpinBox_channelResolution->setDecimals(3);
    ui->doubleSpinBox_channelResolution->setSingleStep(0.05);
    ui->doubleSpinBox_channelResolution->setValue(m_fitSet->getChannelResolution());

    ui->spinBox_iterations->setRange(0, 100000);
    ui->spinBox_iterations->setSingleStep(1);
    ui->spinBox_iterations->setValue(m_fitSet->getMaximumIterations());

    ui->doubleSpinBox_background->setRange(0, 1000000000);
    ui->doubleSpinBox_background->setDecimals(3);
    ui->doubleSpinBox_background->setSingleStep(0.001);
    ui->doubleSpinBox_background->setValue(m_fitSet->getBackgroundParamPtr()->getParameter()->getStartValue());

    ui->widget->setLimits((double)PALSProjectManager::sharedInstance()->getMinChannel(), (double)PALSProjectManager::sharedInstance()->getMaxChannel());
    ui->widget->setLowerLevel((double)m_fitSet->getStartChannel());
    ui->widget->setUpperLevel((double)m_fitSet->getStopChannel());

    emit fitRangeChanged(m_fitSet->getStartChannel(), m_fitSet->getStopChannel());

    connect(ui->tableWidget_Source, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));
    connect(ui->tableWidget_Sample, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));
    connect(ui->tableWidget_Device, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    connect(ui->doubleSpinBox_channelResolution, SIGNAL(valueChanged(double)), this, SLOT(updateChannelResolution(double)));
    connect(ui->spinBox_iterations, SIGNAL(valueChanged(int)), this, SLOT(updateIterations(int)));
    connect(ui->doubleSpinBox_background, SIGNAL(valueChanged(double)), this, SLOT(updateBackground(double)));
    connect(ui->widget, SIGNAL(rangeChanged(double, double)), this, SLOT(updateChannelRange(double, double)));

    emit dataChanged();
}

void ParameterListView::setEnabled(bool enable)
{
    ui->pushButtonAdd_Source->enableWidget(enable);
    ui->pushButtonAdd_Sample->enableWidget(enable);
    ui->pushButtonAdd_Device->enableWidget(enable);

    ui->pushButtonRemove_Source->enableWidget(enable);
    ui->pushButtonRemove_Sample->enableWidget(enable);
    ui->pushButtonRemove_Device->enableWidget(enable);

    ui->pushButton_Background->enableWidget(enable);

    QWidget::setEnabled(enable);
}

void ParameterListView::initializeSourceTableWidget()
{
    if ( !m_fitSet )
        return;

    if ( !m_fitSet->getSourceParamPtr() )
        return;

    while ( m_sourceWidgetCollection.size() > 0 )
    {
        PALSSourceTableWidgetItemCollector *item = m_sourceWidgetCollection.takeFirst();
        DDELETE_SAFETY(item);
    }

    int cnt = 0;
    for ( int i = 0 ; i < m_fitSet->getSourceParamPtr()->getSize() ; i += 2 )
    {
        PALSSourceTableWidgetItemCollector *sourceTableRowWidget_tau = new PALSSourceTableWidgetItemCollector(m_fitSet->getSourceParamPtr()->getParameterAt(i), ui->tableWidget_Source, i);
        m_sourceWidgetCollection.append(sourceTableRowWidget_tau);

        PALSSourceTableWidgetItemCollector *sourceTableRowWidget_I = new PALSSourceTableWidgetItemCollector(m_fitSet->getSourceParamPtr()->getParameterAt(i+1), ui->tableWidget_Source, i+1);
        m_sourceWidgetCollection.append(sourceTableRowWidget_I);

        cnt ++;

        if ( cnt % 2 )
        {
            sourceTableRowWidget_tau->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
            sourceTableRowWidget_I->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
        }
        else
        {
            sourceTableRowWidget_tau->setBackgroundColor(/*Qt::lightGray*/QColor(102, 134, 145));
            sourceTableRowWidget_I->setBackgroundColor(/*Qt::lightGray*/QColor(102, 134, 145));
        }
    }
}

void ParameterListView::initializeSampleTableWidget()
{
    if ( !m_fitSet )
        return;

    if ( !m_fitSet->getLifeTimeParamPtr() )
        return;

    while ( m_sampleWidgetCollection.size() > 0 )
    {
        PALSSampleTableWidgetItemCollector *item = m_sampleWidgetCollection.takeFirst();
        DDELETE_SAFETY(item);
    }

    int cnt = 0;
    for ( int i = 0 ; i < m_fitSet->getLifeTimeParamPtr()->getSize() ; i += 2 )
    {
        PALSSampleTableWidgetItemCollector *sampleTableRowWidget_tau = new PALSSampleTableWidgetItemCollector(m_fitSet->getLifeTimeParamPtr()->getParameterAt(i), ui->tableWidget_Sample, i);
        m_sampleWidgetCollection.append(sampleTableRowWidget_tau);

        PALSSampleTableWidgetItemCollector *sampleTableRowWidget_I = new PALSSampleTableWidgetItemCollector(m_fitSet->getLifeTimeParamPtr()->getParameterAt(i+1), ui->tableWidget_Sample, i+1);
        m_sampleWidgetCollection.append(sampleTableRowWidget_I);

        cnt ++;

        if ( cnt % 2 )
        {
            sampleTableRowWidget_tau->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
            sampleTableRowWidget_I->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
        }
        else
        {
            sampleTableRowWidget_tau->setBackgroundColor(/*Qt::lightGray*/QColor(102, 134, 145));
            sampleTableRowWidget_I->setBackgroundColor(/*Qt::lightGray*/QColor(102, 134, 145));
        }
    }
}

void ParameterListView::initializeDeviceTableWidget()
{
    if ( !m_fitSet )
        return;

    if ( !m_fitSet->getDeviceResolutionParamPtr() )
        return;

    if ( m_fitSet->getDeviceResolutionParamPtr()->getSize() % 3 )
        return;

    while ( m_deviceWidgetCollection.size() > 0 )
    {
        PALSDeviceTableWidgetItemCollector *item = m_deviceWidgetCollection.takeFirst();
        DDELETE_SAFETY(item);
    }

    int cnt = 0;
    for ( int i = 0 ; i < m_fitSet->getDeviceResolutionParamPtr()->getSize() ; i += 3 )
    {
        const DString sigmaText(QString("<b>FWHM<sub>") % QVariant(cnt+1).toString() % QString("</sub></b> [ps]"));
        const DString muText(QString("<b>t0<sub>") % QVariant(cnt+1).toString() % QString("</sub></b> [ps]"));
        const DString IText(QString("<b>I<sub>") % QVariant(cnt+1).toString() % QString("</sub></b> [ps]"));

        m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(i)->setAlias(sigmaText);

        PALSDeviceTableWidgetItemCollector *item_sigma = new PALSDeviceTableWidgetItemCollector(m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(i), ui->tableWidget_Device, i);
        m_deviceWidgetCollection.append(item_sigma);

        m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1)->setAlias(muText);

        PALSDeviceTableWidgetItemCollector *item_mu = new PALSDeviceTableWidgetItemCollector(m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+1), ui->tableWidget_Device, i+1);
        m_deviceWidgetCollection.append(item_mu);

        m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2)->setAlias(IText);

        PALSDeviceTableWidgetItemCollector *item_I = new PALSDeviceTableWidgetItemCollector(m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(i+2), ui->tableWidget_Device, i+2);
        m_deviceWidgetCollection.append(item_I);

        cnt ++;

        if ( cnt % 3 )
        {
            item_sigma->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
            item_mu->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
            item_I->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
        }
        else
        {
            item_sigma->setBackgroundColor(/*Qt::lightGray*/QColor(102, 134, 145));
            item_mu->setBackgroundColor(/*Qt::lightGray*/QColor(102, 134, 145));
            item_I->setBackgroundColor(/*Qt::lightGray*/QColor(102, 134, 145));
        }
    }
}

QCheckBox *ParameterListView::fixedBackgroundCheckBox() const
{
    return ui->checkBox_backgroundFixed;
}

void ParameterListView::updateChannelResolution(double value)
{
    m_fitSet->setChannelResolution(value);
    emit dataChanged();
}

void ParameterListView::updateIterations(int iter)
{
    m_fitSet->setMaximumIterations(iter);
}

void ParameterListView::updateBackground(double value)
{
    if ( !m_fitSet )
        return;

    if ( !m_fitSet->getBackgroundParamPtr() )
        return;

    if ( !m_fitSet->getBackgroundParamPtr()->getParameter() )
        return;

    m_fitSet->getBackgroundParamPtr()->getParameter()->setStartValue(value);
    emit dataChanged();
}

void ParameterListView::updateChannelRange(double lower, double upper)
{
    m_fitSet->setStartChannel((int)lower);
    m_fitSet->setStopChannel((int)upper);

    emit fitRangeChanged((int)lower, (int)upper);
}

void ParameterListView::sendToInstantPreview(int row, int col)
{
    DUNUSED_PARAM(row);
    DUNUSED_PARAM(col);

    emit dataChanged();
}

void ParameterListView::updateBackgroundValue()
{
    const int channels = ui->spinBox_backgroundChannel->value();

    PALSDataStructure *dataStructure = PALSProjectManager::sharedInstance()->getDataStructure();

    if ( !dataStructure )
        return;

    if ( dataStructure->getDataSetPtr()->getLifeTimeData().isEmpty() )
    {
        DMSGBOX("<nobr>No data available. Please import any data before.</nobr>")
        return;
    }

    const int channelMin = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel();
    const int channelMax = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel();
    const int channelDifferenz = channelMax-channelMin;

    if ( channels > channelDifferenz )
    {
        DMSGBOX("<nobr>The range of channels is larger than the ROI.</nobr>");
        return;
    }

    int startChannelIndex = -1;
    int stopChannelIndex = -1;

    int cnt = 0;
    for ( QPointF p : dataStructure->getDataSetPtr()->getLifeTimeData() )
    {
        if ( qFuzzyCompare(p.x(), channelMin) )
            startChannelIndex = cnt;

        if ( qFuzzyCompare(p.x(), channelMax) )
            stopChannelIndex = cnt;

         cnt ++;

         if ( startChannelIndex != -1 && stopChannelIndex != -1 )
             break;
    }

    if ( startChannelIndex == -1 || stopChannelIndex == -1 )
    {
        DMSGBOX("<nobr>Sorry, an unknown error occurred while calculating the background.</nobr>");
    }

    double average = 0.0f;
    int channelDiff = dataStructure->getDataSetPtr()->getLifeTimeData().last().x() -dataStructure->getDataSetPtr()->getLifeTimeData().at(stopChannelIndex).x();
    int indexOffset = dataStructure->getDataSetPtr()->getLifeTimeData().size()-channels-channelDiff-1;

    if ( PALSProjectSettingsManager::sharedInstance()->getBackgroundCalculationFromFirstChannels() )
        indexOffset = 0;

    for ( int i = 0 ; i < channels ; ++ i )
    {
        const int index = indexOffset+i;

        average += dataStructure->getDataSetPtr()->getLifeTimeData().at(index).y();
    }

    average /= channels;

    ui->doubleSpinBox_background->setValue(average);
}

void ParameterListView::saveBackgroundChannelRanges(int channels)
{
    PALSProjectSettingsManager::sharedInstance()->setLastBackgroundChannelRange(channels);
}

void ParameterListView::setUsingFirstChannelsForBkgrdCalc()
{
    PALSProjectSettingsManager::sharedInstance()->setBackgroundCalculationFromFirstChannels(ui->checkBox_FirstChannelBkgrd->isChecked());
}

void ParameterListView::addSourceComponent()
{
    if ( m_fitSet->getComponentsCount() == MAX_COMPONENTS )
    {
        DMSGBOX("<nobr>Sorry, the maximum count of lifetime components is reached.</nobr>");
        return;
    }

    disconnect(ui->tableWidget_Source, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    const int paramCnt = m_fitSet->getSourceParamPtr()->getSize()/2;

    const DString tauText(QString("<b>&#964;<sub>") % QVariant(paramCnt+1).toString() % QString("</sub></b> [ps]"));
    const DString iText(QString("<b>I<sub>") % QVariant(paramCnt+1).toString() % QString("</sub></b>"));

    PALSFitParameter *tau = new PALSFitParameter(m_fitSet->getSourceParamPtr());
    tau->setAlias(tauText);

    tau->setStartValue(120.0f);

    PALSFitParameter *intensity = new PALSFitParameter(m_fitSet->getSourceParamPtr());
    intensity->setAlias(iText);

    intensity->setStartValue(0.1f);

    //add new Param to list:
    PALSSourceTableWidgetItemCollector *item_tau = new PALSSourceTableWidgetItemCollector(tau, ui->tableWidget_Source, ui->tableWidget_Source->rowCount());
    m_sourceWidgetCollection.append(item_tau);

    PALSSourceTableWidgetItemCollector *item_I = new PALSSourceTableWidgetItemCollector(intensity, ui->tableWidget_Source, ui->tableWidget_Source->rowCount());
    m_sourceWidgetCollection.append(item_I);

    updateSourceComponentNames();

    connect(ui->tableWidget_Source, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    emit dataChanged();
}

void ParameterListView::addSampleComponent()
{
    if ( m_fitSet->getComponentsCount() == MAX_COMPONENTS )
    {
        DMSGBOX("<nobr>Sorry, the maximum count of lifetime components is reached.</nobr>");
        return;
    }

    disconnect(ui->tableWidget_Sample, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    const int paramCnt = m_fitSet->getLifeTimeParamPtr()->getSize()/2;

    const DString tauText(QString("<b>&#964;<sub>") % QVariant(paramCnt+1).toString() % QString("</sub></b> [ps]"));
    const DString iText(QString("<b>I<sub>") % QVariant(paramCnt+1).toString() % QString("</sub></b>"));

    PALSFitParameter *tau = new PALSFitParameter(m_fitSet->getLifeTimeParamPtr());
    tau->setAlias(tauText);

    tau->setStartValue(120.0f);

    PALSFitParameter *intensity = new PALSFitParameter(m_fitSet->getLifeTimeParamPtr());
    intensity->setAlias(iText);

    intensity->setStartValue(0.1f);

    //add new Param to list:
    PALSSampleTableWidgetItemCollector *item_tau = new PALSSampleTableWidgetItemCollector(tau, ui->tableWidget_Sample, ui->tableWidget_Sample->rowCount());
    m_sampleWidgetCollection.append(item_tau);

    PALSSampleTableWidgetItemCollector *item_I = new PALSSampleTableWidgetItemCollector(intensity, ui->tableWidget_Sample, ui->tableWidget_Sample->rowCount());
    m_sampleWidgetCollection.append(item_I);

    updateSampleComponentNames();

    connect(ui->tableWidget_Sample, SIGNAL(cellChanged(int, int)), this, SLOT(sendToInstantPreview(int, int)));

    emit dataChanged();
}

void ParameterListView::addDeviceResolutionComponent()
{
    disconnect(ui->tableWidget_Device, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    const int paramCnt = m_fitSet->getDeviceResolutionParamPtr()->getSize()/3;

    const DString sigmaText(QString("<b>FWHM<sub>") % QVariant(paramCnt+1).toString() % QString("</sub></b> [ps]"));
    const DString muText(QString("<b>t0<sub>") % QVariant(paramCnt+1).toString() % QString("</sub></b> [ps]"));
    const DString IText(QString("<b>I<sub>") % QVariant(paramCnt+1).toString() % QString("</sub></b> [ps]"));

    PALSFitParameter *sigma = new PALSFitParameter(m_fitSet->getDeviceResolutionParamPtr());
    sigma->setAlias(sigmaText);

    sigma->setStartValue(220.0f);

    PALSFitParameter *mu = new PALSFitParameter(m_fitSet->getDeviceResolutionParamPtr());
    mu->setAlias(muText);

    mu->setStartValue(1.0f);

    PALSFitParameter *I = new PALSFitParameter(m_fitSet->getDeviceResolutionParamPtr());
    I->setAlias(IText);

    I->setStartValue(0.0f);

    //add new Param to list:
    PALSDeviceTableWidgetItemCollector *item_sigma = new PALSDeviceTableWidgetItemCollector(sigma, ui->tableWidget_Device, ui->tableWidget_Device->rowCount());
    m_deviceWidgetCollection.append(item_sigma);

    PALSDeviceTableWidgetItemCollector *item_mu = new PALSDeviceTableWidgetItemCollector(mu, ui->tableWidget_Device, ui->tableWidget_Device->rowCount());
    m_deviceWidgetCollection.append(item_mu);

    PALSDeviceTableWidgetItemCollector *item_I = new PALSDeviceTableWidgetItemCollector(I, ui->tableWidget_Device, ui->tableWidget_Device->rowCount());
    m_deviceWidgetCollection.append(item_I);

    updateDeviceComponentNames();

    connect(ui->tableWidget_Device, SIGNAL(cellChanged(int, int)), this, SLOT(sendToInstantPreview(int, int)));

    emit dataChanged();
}

void ParameterListView::removeSourceComponent()
{
    if ( m_fitSet->getComponentsCount() == MIN_COMPONENTS )
    {
        DMSGBOX("<nobr>Sorry, no components available.</nobr>");
        return;
    }

    const QModelIndexList indexList = ui->tableWidget_Source->selectionModel()->selectedIndexes();

    int row = -1;

    foreach ( QModelIndex index, indexList )
        row = index.row();

    if ( row == -1 )
    {
        DMSGBOX("Please select the component you want to delete!");
        return;
    }

    disconnect(ui->tableWidget_Source, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    //index == even: tau-value; odd: intensity-value
    if ( row % 2 ) //intensity-row?
    {
        PALSFitParameter *tauParam = m_fitSet->getSourceParamPtr()->getParameterAt(row-1);
        PALSFitParameter *intensityParam = m_fitSet->getSourceParamPtr()->getParameterAt(row);

        m_fitSet->getSourceParamPtr()->removeParameter(tauParam);
        m_fitSet->getSourceParamPtr()->removeParameter(intensityParam);

        ui->tableWidget_Source->removeRow(row-1);
        ui->tableWidget_Source->removeRow(row-1);

        PALSSourceTableWidgetItemCollector* collector_tau = m_sourceWidgetCollection.takeAt(row-1);
        DDELETE_SAFETY(collector_tau);

        PALSSourceTableWidgetItemCollector* collector_I = m_sourceWidgetCollection.takeAt(row-1);
        DDELETE_SAFETY(collector_I);
    }
    else //tau-row?
    {
        PALSFitParameter *tauParam = m_fitSet->getSourceParamPtr()->getParameterAt(row);
        PALSFitParameter *intensityParam = m_fitSet->getSourceParamPtr()->getParameterAt(row+1);

        m_fitSet->getSourceParamPtr()->removeParameter(tauParam);
        m_fitSet->getSourceParamPtr()->removeParameter(intensityParam);

        ui->tableWidget_Source->removeRow(row+1);
        ui->tableWidget_Source->removeRow(row);

        PALSSourceTableWidgetItemCollector* collector_tau = m_sourceWidgetCollection.takeAt(row+1);
        DDELETE_SAFETY(collector_tau);

        PALSSourceTableWidgetItemCollector* collector_I = m_sourceWidgetCollection.takeAt(row);
        DDELETE_SAFETY(collector_I);
    }

    updateSourceComponentNames();

    ui->tableWidget_Source->selectRow(ui->tableWidget_Source->rowCount()-1);

    connect(ui->tableWidget_Source, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    emit dataChanged();
}

void ParameterListView::removeSampleComponent()
{
    if ( m_fitSet->getComponentsCount() == MIN_COMPONENTS )
    {
        DMSGBOX("<nobr>Sorry, no components available.</nobr>");
        return;
    }

    const QModelIndexList indexList = ui->tableWidget_Sample->selectionModel()->selectedIndexes();

    int row = -1;

    foreach ( QModelIndex index, indexList )
        row = index.row();

    if ( row == -1 )
    {
        DMSGBOX("Please select the component you want to delete!");
        return;
    }

    disconnect(ui->tableWidget_Sample, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    //index == even: tau-value; odd: intensity-value
    if ( row % 2 ) //intensity-row?
    {
        PALSFitParameter *tauParam = m_fitSet->getLifeTimeParamPtr()->getParameterAt(row-1);
        PALSFitParameter *intensityParam = m_fitSet->getLifeTimeParamPtr()->getParameterAt(row);

        m_fitSet->getLifeTimeParamPtr()->removeParameter(tauParam);
        m_fitSet->getLifeTimeParamPtr()->removeParameter(intensityParam);

        ui->tableWidget_Sample->removeRow(row-1);
        ui->tableWidget_Sample->removeRow(row-1);

        PALSSampleTableWidgetItemCollector* collector_tau = m_sampleWidgetCollection.takeAt(row-1);
        DDELETE_SAFETY(collector_tau);

        PALSSampleTableWidgetItemCollector* collector_I = m_sampleWidgetCollection.takeAt(row-1);
        DDELETE_SAFETY(collector_I);
    }
    else //tau-row?
    {
        PALSFitParameter *tauParam = m_fitSet->getLifeTimeParamPtr()->getParameterAt(row);
        PALSFitParameter *intensityParam = m_fitSet->getLifeTimeParamPtr()->getParameterAt(row+1);

        m_fitSet->getLifeTimeParamPtr()->removeParameter(tauParam);
        m_fitSet->getLifeTimeParamPtr()->removeParameter(intensityParam);

        ui->tableWidget_Sample->removeRow(row+1);
        ui->tableWidget_Sample->removeRow(row);

        PALSSampleTableWidgetItemCollector* collector_tau = m_sampleWidgetCollection.takeAt(row+1);
        DDELETE_SAFETY(collector_tau);

        PALSSampleTableWidgetItemCollector* collector_I = m_sampleWidgetCollection.takeAt(row);
        DDELETE_SAFETY(collector_I);
    }

    updateSampleComponentNames();

    ui->tableWidget_Sample->selectRow(ui->tableWidget_Sample->rowCount()-1);

    connect(ui->tableWidget_Sample, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    emit dataChanged();
}

void ParameterListView::removeDeviceResolutionComponent()
{
    if ( m_fitSet->getDeviceResolutionParamPtr()->getSize() == MIN_COMPONENTS + 1 )
    {
        DMSGBOX("<nobr>Sorry, at least 1 component is required.</nobr>");
        return;
    }

    const QModelIndexList indexList = ui->tableWidget_Device->selectionModel()->selectedIndexes();

    int row = -1;

    foreach ( QModelIndex index, indexList )
        row = index.row();

    if ( row == -1 )
    {
        DMSGBOX("Please select the component you want to delete!");
        return;
    }

    disconnect(ui->tableWidget_Device, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));


    if ( !((row+1) % 3) ) //I-row?
    {
        PALSFitParameter *sigmaParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row-2);
        PALSFitParameter *muParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row-1);
        PALSFitParameter *IParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row);

        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(sigmaParam);
        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(muParam);
        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(IParam);

        ui->tableWidget_Device->removeRow(row-2);
        ui->tableWidget_Device->removeRow(row-2);
        ui->tableWidget_Device->removeRow(row-2);

        PALSDeviceTableWidgetItemCollector* collector_sigma = m_deviceWidgetCollection.takeAt(row-2);
        DDELETE_SAFETY(collector_sigma);

        PALSDeviceTableWidgetItemCollector* collector_mu = m_deviceWidgetCollection.takeAt(row-2);
        DDELETE_SAFETY(collector_mu);

        PALSDeviceTableWidgetItemCollector* collector_I = m_deviceWidgetCollection.takeAt(row-2);
        DDELETE_SAFETY(collector_I);
    }
    else if ( !((row+1) % 2) ) //Mu-row?
    {
        PALSFitParameter *sigmaParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row-1);
        PALSFitParameter *muParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row);
        PALSFitParameter *IParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row+1);

        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(sigmaParam);
        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(muParam);
        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(IParam);

        ui->tableWidget_Device->removeRow(row-1);
        ui->tableWidget_Device->removeRow(row-1);
        ui->tableWidget_Device->removeRow(row-1);

        PALSDeviceTableWidgetItemCollector* collector_sigma = m_deviceWidgetCollection.takeAt(row-1);
        DDELETE_SAFETY(collector_sigma);

        PALSDeviceTableWidgetItemCollector* collector_mu = m_deviceWidgetCollection.takeAt(row-1);
        DDELETE_SAFETY(collector_mu);

        PALSDeviceTableWidgetItemCollector* collector_I = m_deviceWidgetCollection.takeAt(row-1);
        DDELETE_SAFETY(collector_I);
    }
    else //FWHM-row?
    {
        PALSFitParameter *sigmaParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row);
        PALSFitParameter *muParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row+1);
        PALSFitParameter *IParam = m_fitSet->getDeviceResolutionParamPtr()->getParameterAt(row+2);

        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(sigmaParam);
        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(muParam);
        m_fitSet->getDeviceResolutionParamPtr()->removeParameter(IParam);

        ui->tableWidget_Device->removeRow(row);
        ui->tableWidget_Device->removeRow(row);
        ui->tableWidget_Device->removeRow(row);

        PALSDeviceTableWidgetItemCollector* collector_sigma = m_deviceWidgetCollection.takeAt(row);
        DDELETE_SAFETY(collector_sigma);

        PALSDeviceTableWidgetItemCollector* collector_mu = m_deviceWidgetCollection.takeAt(row);
        DDELETE_SAFETY(collector_mu);

        PALSDeviceTableWidgetItemCollector* collector_I = m_deviceWidgetCollection.takeAt(row);
        DDELETE_SAFETY(collector_I);
    }

    updateDeviceComponentNames();

    ui->tableWidget_Device->selectRow(ui->tableWidget_Device->rowCount()-1);

    connect(ui->tableWidget_Device, SIGNAL(cellChanged(int,int)), this, SLOT(sendToInstantPreview(int, int)));

    emit dataChanged();
}

void ParameterListView::updateSourceComponentNames()
{
    int cnt = 0;
    for ( int i = 0 ; i < m_sourceWidgetCollection.size() ; i += 2 )
    {
        PALSSourceTableWidgetItemCollector *sourceTableRowWidget_tau = m_sourceWidgetCollection.at(i);
        PALSSourceTableWidgetItemCollector *sourceTableRowWidget_I = m_sourceWidgetCollection.at(i+1);

        const DString tauText(QString("<b>&#964;<sub>") % QVariant(cnt+1).toString() % QString("</sub></b> [ps]"));
        const DString iText(QString("<b>I<sub>") % QVariant(cnt+1).toString() % QString("</sub></b>"));

         sourceTableRowWidget_tau->setAlias(tauText);
         sourceTableRowWidget_I->setAlias(iText);

         cnt ++;

         if ( cnt % 2 )
         {
             sourceTableRowWidget_tau->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
             sourceTableRowWidget_I->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
         }
         else
         {
             sourceTableRowWidget_tau->setBackgroundColor(/*Qt::gray*/QColor(102, 134, 145));
             sourceTableRowWidget_I->setBackgroundColor(/*Qt::gray*/QColor(102, 134, 145));
         }
    }
}

void ParameterListView::updateSampleComponentNames()
{
    int cnt = 0;
    for ( int i = 0 ; i < m_sampleWidgetCollection.size() ; i += 2 )
    {
        PALSSampleTableWidgetItemCollector *sampleTableRowWidget_tau = m_sampleWidgetCollection.at(i);
        PALSSampleTableWidgetItemCollector *sampleTableRowWidget_I = m_sampleWidgetCollection.at(i+1);

        const DString tauText(QString("<b>&#964;<sub>") % QVariant(cnt+1).toString() % QString("</sub></b> [ps]"));
        const DString iText(QString("<b>I<sub>") % QVariant(cnt+1).toString() % QString("</sub></b>"));

         sampleTableRowWidget_tau->setAlias(tauText);
         sampleTableRowWidget_I->setAlias(iText);

         cnt ++;

         if ( cnt % 2 )
         {
             sampleTableRowWidget_tau->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
             sampleTableRowWidget_I->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
         }
         else
         {
             sampleTableRowWidget_tau->setBackgroundColor(/*Qt::gray*/QColor(102, 134, 145));
             sampleTableRowWidget_I->setBackgroundColor(/*Qt::gray*/QColor(102, 134, 145));
         }
    }
}

void ParameterListView::updateDeviceComponentNames()
{
    int cnt = 0;
    for ( int i = 0 ; i < m_deviceWidgetCollection.size() ; i += 3 )
    {
        PALSDeviceTableWidgetItemCollector *deviceTableRowWidget_sigma = m_deviceWidgetCollection.at(i);
        PALSDeviceTableWidgetItemCollector *deviceTableRowWidget_mu = m_deviceWidgetCollection.at(i+1);
        PALSDeviceTableWidgetItemCollector *deviceTableRowWidget_I = m_deviceWidgetCollection.at(i+2);

        const DString sigmaText(QString("<b>FWHM<sub>") % QVariant(cnt+1).toString() % QString("</sub></b> [ps]"));
        const DString muText(QString("<b>t0<sub>") % QVariant(cnt+1).toString() % QString("</sub></b> [ps]"));
        const DString IText(QString("<b>I<sub>") % QVariant(cnt+1).toString() % QString("</sub></b> [ps]"));

         deviceTableRowWidget_sigma->setAlias(sigmaText);
         deviceTableRowWidget_mu->setAlias(muText);
         deviceTableRowWidget_I->setAlias(IText);

         cnt ++;

         if ( cnt % 2 )
         {
             deviceTableRowWidget_sigma->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
             deviceTableRowWidget_mu->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
             deviceTableRowWidget_I->setBackgroundColor(/*Qt::lightGray*/QColor(145, 185, 199));
         }
         else
         {
             deviceTableRowWidget_sigma->setBackgroundColor(/*Qt::gray*/QColor(102, 134, 145));
             deviceTableRowWidget_mu->setBackgroundColor(/*Qt::gray*/QColor(102, 134, 145));
             deviceTableRowWidget_I->setBackgroundColor(/*Qt::gray*/QColor(102, 134, 145));
         }
    }
}

void ParameterListView::refreshBackgroundValue(double background)
{
    ui->doubleSpinBox_background->setValue(background);
}

void ParameterListView::setFitRangeLimits(int lower, int upper)
{
    disconnect(ui->widget, SIGNAL(rangeChanged(double, double)), this, SLOT(updateChannelRange(double, double)));

    ui->widget->setLimits(lower, upper);
    ui->widget->setLowerLevel(lower);
    ui->widget->setUpperLevel(upper);

    connect(ui->widget, SIGNAL(rangeChanged(double, double)), this, SLOT(updateChannelRange(double, double)));

    emit fitRangeChanged(lower, upper);
}

void ParameterListView::setFitRange(int lower, int upper)
{
    disconnect(ui->widget, SIGNAL(rangeChanged(double, double)), this, SLOT(updateChannelRange(double, double)));

    ui->widget->setLowerLevel(lower);
    ui->widget->setUpperLevel(upper);

    connect(ui->widget, SIGNAL(rangeChanged(double, double)), this, SLOT(updateChannelRange(double, double)));

    emit fitRangeChanged(lower, upper);
}

void ParameterListView::setBackgroundChannelRange(int range)
{
    ui->spinBox_backgroundChannel->setValue(range);
}

void ParameterListView::setBackgroundCalculationUsingFirstChannels(bool first)
{
    ui->checkBox_FirstChannelBkgrd->setChecked(first);
}

void PALSSourceTableWidgetItemCollector::updateValue(int row, int column)
{
    DUNUSED_PARAM(row);

    switch ( column ){
    /*case 0:
        m_param->setActive(m_activeItem->isChecked());
        break;*/

    case 0:
        m_param->setAlias(m_aliasItem->text());
        break;

    case 1:
        m_param->setName(m_nameItem->text());
        break;

    case 2:
    {
        bool ok = false;
        const double value = QVariant(m_startValueItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setStartValue(value);
        else
            m_startValueItem->setText(QString::number(m_param->getStartValue(), 'f', 5));
    }
        break;

    case 3:
    {
        bool ok = false;
        const double value = QVariant(m_lowerLimitItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setLowerBoundingValue(value);
        else
            m_lowerLimitItem->setText(QString::number(m_param->getLowerBoundingValue(), 'f', 5));

        m_param->setLowerBoundingEnabled(m_lowerLimitItem->checkState()==0?false:true);
    }
        break;

    case 4:
    {
        bool ok = false;
        const double value = QVariant(m_upperLimitItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setUpperBoundingValue(value);
        else
            m_upperLimitItem->setText(QString::number(m_param->getUpperBoundingValue(), 'f', 5));

        m_param->setUpperBoundingEnabled(m_upperLimitItem->checkState()==0?false:true);
    }
        break;

    case 5:
        m_param->setAsFixed(m_fixedItem->isChecked());
        break;

    default:
        break;
    }
}

void PALSSourceTableWidgetItemCollector::checkBoxStateChanged(CheckBoxTableWidgetItem *widget)
{
    if ( !widget )
        return;

    /*if ( widget == m_activeItem )
        m_param->setActive(m_activeItem->isChecked());
    else if ( widget == m_fixedItem )*/
        m_param->setAsFixed(m_fixedItem->isChecked());
}

void PALSSourceTableWidgetItemCollector::setAlias(const QString &name)
{
    if ( !m_param )
        return;

    m_aliasItem->setText(name);
    m_param->setAlias(name);
}

void PALSSourceTableWidgetItemCollector::setBackgroundColor(const QColor &color)
{
    //m_fixedItem->setBackgroundColor(color);
    m_nameItem->setBackgroundColor(color);
    //m_aliasItem->setBackgroundColor(color);
    m_lowerLimitItem->setBackgroundColor(color);
    m_upperLimitItem->setBackgroundColor(color);
    m_startValueItem->setBackgroundColor(color);
}

PALSSourceTableWidgetItemCollector::PALSSourceTableWidgetItemCollector(PALSFitParameter *fitParam, QTableWidget* tableWidget, int row) :
    m_param(fitParam)
{
    tableWidget->insertRow(row);

    m_fixedItem = new CheckBoxTableWidgetItem;
    m_aliasItem = new LabelTableWidgetItem; m_aliasItem->setColor(Qt::blue);
    m_nameItem = new QTableWidgetItem; m_nameItem->setTextAlignment(Qt::AlignCenter);
    m_startValueItem = new QTableWidgetItem;
    m_lowerLimitItem = new QTableWidgetItem; m_lowerLimitItem->setCheckState(Qt::Checked);
    m_upperLimitItem = new QTableWidgetItem; m_upperLimitItem->setCheckState(Qt::Checked);

    tableWidget->setCellWidget(row, 0, m_aliasItem);
    tableWidget->setItem(row, 1, m_nameItem);
    tableWidget->setItem(row, 2, m_startValueItem);
    tableWidget->setItem(row, 3, m_lowerLimitItem);
    tableWidget->setItem(row, 4, m_upperLimitItem);
    tableWidget->setCellWidget(row, 5, m_fixedItem);

    m_aliasItem->setText(m_param->getAlias());
    m_nameItem->setText(m_param->getName());
    m_startValueItem->setText(QString::number(m_param->getStartValue(), 'f', 5));
    m_lowerLimitItem->setText(QString::number(m_param->getLowerBoundingValue(), 'f', 5));
    m_lowerLimitItem->setCheckState(!m_param->isLowerBoundingEnabled()?Qt::Unchecked:Qt::Checked);
    m_upperLimitItem->setText(QString::number(m_param->getUpperBoundingValue(), 'f', 5));
    m_upperLimitItem->setCheckState(!m_param->isUpperBoundingEnabled()?Qt::Unchecked:Qt::Checked);
    m_fixedItem->setCheckedState(m_param->isFixed());

    connect(tableWidget,  SIGNAL(cellChanged(int, int)), this, SLOT(updateValue(int, int)));
    connect(m_fixedItem,  SIGNAL(stateChanged(CheckBoxTableWidgetItem*)), this, SLOT(checkBoxStateChanged(CheckBoxTableWidgetItem*)));

#if defined(Q_OS_WIN)
    tableWidget->setFont(WINDOWS_FONT(10));
    QFont font = tableWidget->font();
    font.setBold(true);
    tableWidget->setFont(font);
#endif
}

PALSSourceTableWidgetItemCollector::~PALSSourceTableWidgetItemCollector() {}


void PALSSampleTableWidgetItemCollector::updateValue(int row, int column)
{
    DUNUSED_PARAM(row);

    switch ( column ) {
    case 0:
        m_param->setAlias(m_aliasItem->text());
        break;

    case 1:
        m_param->setName(m_nameItem->text());
        break;

    case 2:
    {
        bool ok = false;
        const double value = QVariant(m_startValueItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setStartValue(value);
        else
            m_startValueItem->setText(QString::number(m_param->getStartValue(), 'f', 5));
    }
        break;

    case 3:
    {
        bool ok = false;
        const double value = QVariant(m_lowerLimitItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setLowerBoundingValue(value);
        else
            m_lowerLimitItem->setText(QString::number(m_param->getLowerBoundingValue(), 'f', 5));

        m_param->setLowerBoundingEnabled(m_lowerLimitItem->checkState()==0?false:true);
    }
        break;

    case 4:
    {
        bool ok = false;
        const double value = QVariant(m_upperLimitItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setUpperBoundingValue(value);
        else
            m_upperLimitItem->setText(QString::number(m_param->getUpperBoundingValue(), 'f', 5));

        m_param->setUpperBoundingEnabled(m_upperLimitItem->checkState()==0?false:true);
    }
        break;

    case 5:
        m_param->setAsFixed(m_fixedItem->isChecked());
        break;

    default:
        break;
    }
}

void PALSSampleTableWidgetItemCollector::checkBoxStateChanged(CheckBoxTableWidgetItem *widget)
{
    if ( !widget )
        return;

    m_param->setAsFixed(m_fixedItem->isChecked());
}

void PALSSampleTableWidgetItemCollector::setAlias(const QString &name)
{
     if ( !m_param )
         return;

     m_aliasItem->setText(name);
     m_param->setAlias(name);
}

void PALSSampleTableWidgetItemCollector::setBackgroundColor(const QColor &color)
{
    m_nameItem->setBackgroundColor(color);
    m_lowerLimitItem->setBackgroundColor(color);
    m_upperLimitItem->setBackgroundColor(color);
    m_startValueItem->setBackgroundColor(color);
}

PALSSampleTableWidgetItemCollector::PALSSampleTableWidgetItemCollector(PALSFitParameter *fitParam, QTableWidget* tableWidget, int row) :
    m_param(fitParam)
{
    tableWidget->insertRow(row); //insert the new row!

    m_fixedItem = new CheckBoxTableWidgetItem;
    m_nameItem = new QTableWidgetItem; m_nameItem->setTextAlignment(Qt::AlignCenter);
    m_aliasItem = new LabelTableWidgetItem; m_aliasItem->setColor(Qt::blue);
    m_startValueItem = new QTableWidgetItem;
    m_lowerLimitItem = new QTableWidgetItem; m_lowerLimitItem->setCheckState(Qt::Checked);
    m_upperLimitItem = new QTableWidgetItem; m_upperLimitItem->setCheckState(Qt::Checked);

    tableWidget->setCellWidget(row, 0, m_aliasItem);
    tableWidget->setItem(row, 1, m_nameItem);
    tableWidget->setItem(row, 2, m_startValueItem);
    tableWidget->setItem(row, 3, m_lowerLimitItem);
    tableWidget->setItem(row, 4, m_upperLimitItem);
    tableWidget->setCellWidget(row, 5, m_fixedItem);

    m_nameItem->setText(m_param->getName());
    m_aliasItem->setText(m_param->getAlias());
    m_startValueItem->setText(QString::number(m_param->getStartValue(), 'f', 5));
    m_lowerLimitItem->setText(QString::number(m_param->getLowerBoundingValue(), 'f', 5));
    m_lowerLimitItem->setCheckState(!m_param->isLowerBoundingEnabled()?Qt::Unchecked:Qt::Checked);
    m_upperLimitItem->setText(QString::number(m_param->getUpperBoundingValue(), 'f', 5));
    m_upperLimitItem->setCheckState(!m_param->isUpperBoundingEnabled()?Qt::Unchecked:Qt::Checked);
    m_fixedItem->setCheckedState(m_param->isFixed());

    connect(tableWidget,  SIGNAL(cellChanged(int, int)), this, SLOT(updateValue(int, int)));
    connect(m_fixedItem,  SIGNAL(stateChanged(CheckBoxTableWidgetItem*)), this, SLOT(checkBoxStateChanged(CheckBoxTableWidgetItem*)));

#if defined(Q_OS_WIN)
    tableWidget->setFont(WINDOWS_FONT(10));
    QFont font = tableWidget->font();
    font.setBold(true);
    tableWidget->setFont(font);
#endif
}

PALSSampleTableWidgetItemCollector::~PALSSampleTableWidgetItemCollector() {}


PALSBackgroundTableWidgetItemCollector::PALSBackgroundTableWidgetItemCollector(PALSFitParameter *fitParam, QTableWidget *tableWidget, int row) :
    m_param(fitParam)
{
    tableWidget->insertRow(row); //insert the new row!

    m_activeItem = new CheckBoxTableWidgetItem;
    m_fixedItem = new CheckBoxTableWidgetItem;

    m_startValueItem = new QTableWidgetItem;
    m_lowerLimitItem = new QTableWidgetItem; m_lowerLimitItem->setCheckState(Qt::Checked);
    m_upperLimitItem = new QTableWidgetItem; m_upperLimitItem->setCheckState(Qt::Checked);

    tableWidget->setCellWidget(row, 0, m_activeItem);
    tableWidget->setItem(row, 1, m_startValueItem);
    tableWidget->setItem(row, 2, m_lowerLimitItem);
    tableWidget->setItem(row, 3, m_upperLimitItem);
    tableWidget->setCellWidget(row, 4, m_fixedItem);

    m_activeItem->setCheckedState(m_param->isActive());
    m_startValueItem->setText(QString::number(m_param->getStartValue(), 'f', 5));
    m_lowerLimitItem->setText(QString::number(m_param->getLowerBoundingValue(), 'f', 5));
    m_lowerLimitItem->setCheckState(!m_param->isLowerBoundingEnabled()?Qt::Unchecked:Qt::Checked);
    m_upperLimitItem->setText(QString::number(m_param->getUpperBoundingValue(), 'f', 5));
    m_upperLimitItem->setCheckState(!m_param->isUpperBoundingEnabled()?Qt::Unchecked:Qt::Checked);
    m_fixedItem->setCheckedState(m_param->isFixed());

    connect(tableWidget,  SIGNAL(cellChanged(int, int)), this, SLOT(updateValue(int, int)));
    connect(m_activeItem, SIGNAL(stateChanged(CheckBoxTableWidgetItem*)), this, SLOT(checkBoxStateChanged(CheckBoxTableWidgetItem*)));
    connect(m_fixedItem,  SIGNAL(stateChanged(CheckBoxTableWidgetItem*)), this, SLOT(checkBoxStateChanged(CheckBoxTableWidgetItem*)));
}

PALSBackgroundTableWidgetItemCollector::~PALSBackgroundTableWidgetItemCollector() {}

void PALSBackgroundTableWidgetItemCollector::updateValue(int row, int column)
{
    DUNUSED_PARAM(row);

    switch ( column ){
    case 0:
        m_param->setActive(m_activeItem->isChecked());
        break;

    case 1:
        m_param->setStartValue(QVariant(m_startValueItem->text()).toDouble());
        break;

    case 2:
        m_param->setLowerBoundingValue(QVariant(m_lowerLimitItem->text()).toDouble());
        m_param->setLowerBoundingEnabled(m_lowerLimitItem->checkState()==0?false:true);
        break;

    case 3:
        m_param->setUpperBoundingValue(QVariant(m_upperLimitItem->text()).toDouble());
        m_param->setUpperBoundingEnabled(m_upperLimitItem->checkState()==0?false:true);
        break;

    case 4:
        m_param->setAsFixed(m_fixedItem->isChecked());
        break;

    default:
        break;
    }
}

void PALSBackgroundTableWidgetItemCollector::checkBoxStateChanged(CheckBoxTableWidgetItem *widget)
{
    if ( !widget )
        return;

    if ( widget == m_activeItem )
        m_param->setActive(m_activeItem->isChecked());
    else if ( widget == m_fixedItem )
        m_param->setAsFixed(m_fixedItem->isChecked());
}

void PALSBackgroundTableWidgetItemCollector::setBackground(double val)
{
    m_startValueItem->setText(QVariant(val).toString());
    m_param->setStartValue(val);
}


CheckBoxTableWidgetItem::CheckBoxTableWidgetItem(QWidget *parent) :
    QWidget(parent)
{
    m_pCheckBox = new QCheckBox;
    m_pLayout = new QHBoxLayout(this);

    m_pLayout->addWidget(m_pCheckBox);
    m_pLayout->setAlignment(Qt::AlignCenter);
    m_pLayout->setContentsMargins(0, 0, 0, 0);

    setLayout(m_pLayout);

    connect(m_pCheckBox, SIGNAL(clicked(bool)), this, SLOT(checked(bool)));
}

CheckBoxTableWidgetItem::~CheckBoxTableWidgetItem()
{
    DDELETE_SAFETY(m_pCheckBox);
    DDELETE_SAFETY(m_pLayout);
}

bool CheckBoxTableWidgetItem::isChecked() const
{
    return m_pCheckBox->isChecked();
}

void CheckBoxTableWidgetItem::setCheckedState(bool checked)
{
    m_pCheckBox->setCheckState(checked?Qt::Checked:Qt::Unchecked);
}

void CheckBoxTableWidgetItem::checked(bool checked)
{
    emit clicked(checked);
    emit stateChanged(this);
}

void CheckBoxTableWidgetItem::setBackgroundColor(const QColor &color)
{
    setStyleSheet("background: rgb(" % QVariant(color.red()).toString() % ", " %  QVariant(color.green()).toString() % ", " % QVariant(color.blue()).toString() % ");");
    m_pCheckBox->setStyleSheet("background: rgb(" % QVariant(color.red()).toString() % ", " %  QVariant(color.green()).toString() % ", " % QVariant(color.blue()).toString() % ");");
}


PALSDeviceTableWidgetItemCollector::PALSDeviceTableWidgetItemCollector(PALSFitParameter *fitParam, QTableWidget *tableWidget, int row) :
    m_param(fitParam)
{
    tableWidget->insertRow(row); //insert the new row!

    m_fixedItem = new CheckBoxTableWidgetItem;
    m_nameItem = new QTableWidgetItem; m_nameItem->setTextAlignment(Qt::AlignCenter);
    m_aliasItem = new LabelTableWidgetItem; m_aliasItem->setColor(Qt::blue);
    m_startValueItem = new QTableWidgetItem;
    m_lowerLimitItem = new QTableWidgetItem; m_lowerLimitItem->setCheckState(Qt::Checked);
    m_upperLimitItem = new QTableWidgetItem; m_upperLimitItem->setCheckState(Qt::Checked);

    tableWidget->setCellWidget(row, 0, m_aliasItem);
    tableWidget->setItem(row, 1, m_nameItem);
    tableWidget->setItem(row, 2, m_startValueItem);
    tableWidget->setItem(row, 3, m_lowerLimitItem);
    tableWidget->setItem(row, 4, m_upperLimitItem);
    tableWidget->setCellWidget(row, 5, m_fixedItem);

    m_nameItem->setText(m_param->getName());
    m_aliasItem->setText(m_param->getAlias());
    m_startValueItem->setText(QString::number(m_param->getStartValue(), 'f', 5));
    m_lowerLimitItem->setText(QString::number(m_param->getLowerBoundingValue(), 'f', 5));
    m_lowerLimitItem->setCheckState(!m_param->isLowerBoundingEnabled()?Qt::Unchecked:Qt::Checked);
    m_upperLimitItem->setText(QString::number(m_param->getUpperBoundingValue(), 'f', 5));
    m_upperLimitItem->setCheckState(!m_param->isUpperBoundingEnabled()?Qt::Unchecked:Qt::Checked);
    m_fixedItem->setCheckedState(m_param->isFixed());

    connect(tableWidget,  SIGNAL(cellChanged(int, int)), this, SLOT(updateValue(int, int)));
    connect(m_fixedItem,  SIGNAL(stateChanged(CheckBoxTableWidgetItem*)), this, SLOT(checkBoxStateChanged(CheckBoxTableWidgetItem*)));

#if defined(Q_OS_WIN)
    tableWidget->setFont(WINDOWS_FONT(10));
    QFont font = tableWidget->font();
    font.setBold(true);
    tableWidget->setFont(font);
#endif
}

PALSDeviceTableWidgetItemCollector::~PALSDeviceTableWidgetItemCollector() {}

void PALSDeviceTableWidgetItemCollector::updateValue(int row, int column)
{
    DUNUSED_PARAM(row);

    switch ( column ){
    case 0:
        m_param->setAlias(m_aliasItem->text());
        break;

    case 1:
        m_param->setName(m_nameItem->text());
        break;

    case 2:
    {
        bool ok = false;
        const double value = QVariant(m_startValueItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setStartValue(value);
        else
        {
            m_startValueItem->setText(QString::number(m_param->getStartValue(), 'f', 5));
        }
    }
        break;

    case 3:
    {
        bool ok = false;
        const double value = QVariant(m_lowerLimitItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setLowerBoundingValue(value);
        else
        {
            m_lowerLimitItem->setText(QString::number(m_param->getLowerBoundingValue(), 'f', 5));
        }

        m_param->setLowerBoundingEnabled(m_lowerLimitItem->checkState()==0?false:true);
    }
        break;

    case 4:
    {
        bool ok = false;
        const double value = QVariant(m_upperLimitItem->text()).toDouble(&ok);

        if ( ok )
            m_param->setUpperBoundingValue(value);
        else
        {
            m_upperLimitItem->setText(QString::number(m_param->getUpperBoundingValue(), 'f', 5));
        }

        m_param->setUpperBoundingEnabled(m_upperLimitItem->checkState()==0?false:true);
    }
        break;

    case 5:
        m_param->setAsFixed(m_fixedItem->isChecked());
        break;

    default:
        break;
    }
}

void PALSDeviceTableWidgetItemCollector::checkBoxStateChanged(CheckBoxTableWidgetItem *widget)
{
    if ( !widget )
        return;

    if ( widget == m_fixedItem )
        m_param->setAsFixed(m_fixedItem->isChecked());
}

void PALSDeviceTableWidgetItemCollector::setAlias(const QString &name)
{
    if ( !m_param )
        return;

    m_aliasItem->setText(name);
    m_param->setAlias(name);
}

void PALSDeviceTableWidgetItemCollector::setBackgroundColor(const QColor &color)
{
    m_nameItem->setBackgroundColor(color);
    m_lowerLimitItem->setBackgroundColor(color);
    m_upperLimitItem->setBackgroundColor(color);
    m_startValueItem->setBackgroundColor(color);
}


LabelTableWidgetItem::LabelTableWidgetItem(QWidget *parent) :
    QWidget(parent)
{
    m_label = new QLabel;
    m_pLayout = new QHBoxLayout(this);

    m_label->setTextFormat(Qt::RichText);

    m_pLayout->addWidget(m_label);
    m_pLayout->setAlignment(Qt::AlignCenter);
    m_pLayout->setContentsMargins(0, 0, 0, 0);

    setLayout(m_pLayout);
}

LabelTableWidgetItem::~LabelTableWidgetItem()
{
    DDELETE_SAFETY(m_label);
    DDELETE_SAFETY(m_pLayout);
}

DString LabelTableWidgetItem::text() const
{
    return (DString)m_label->text();
}

void LabelTableWidgetItem::setText(const DString &text)
{
    m_label->setText(text);
}

void LabelTableWidgetItem::setBackgroundColor(const QColor &color)
{
    setStyleSheet("background-color: rgb(" % QVariant(color.red()).toString() % ", " %  QVariant(color.green()).toString() % ", " % QVariant(color.blue()).toString() % ");");
    m_label->setStyleSheet("background: rgb(" % QVariant(color.red()).toString() % ", " %  QVariant(color.green()).toString() % ", " % QVariant(color.blue()).toString() % ");");
}

void LabelTableWidgetItem::setColor(const QColor &color)
{
    m_label->setStyleSheet("color: rgb(" % QVariant(color.red()).toString() % ", " %  QVariant(color.green()).toString() % ", " % QVariant(color.blue()).toString() % ");");
}
