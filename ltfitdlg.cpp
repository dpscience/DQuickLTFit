/****************************************************************************
**
**  DQuickLTFit, a software for the analysis of Positron-Lifetime Spectra
**  based on the Least-Square Optimization using the Levenberg-Marquardt
**  Algorithm.
**
**  Copyright (C) 2016-2019 Danny Petschke
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

#include "ltfitdlg.h"
#include "ui_ltfitdlg.h"

#ifndef WINDOWS_FONT
#define WINDOWS_FONT(__pointSize__)  QFont("Arial", __pointSize__)
#endif

DFastLTFitDlg::DFastLTFitDlg(const QString projectPath, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::DFastLTFitDlg),
    m_onStart(false)
{
    ui->setupUi(this);

    QApplication::setWindowIcon(QIcon(":/localImages/Images/IconPNGRounded.png"));

    if ( PALSProjectSettingsManager::sharedInstance()->load() )
    {
        if ( !PALSProjectSettingsManager::sharedInstance()->getLastProjectPathList().isEmpty() )
        {
            m_lastProjectsMenu = new QMenu(ui->menuLoad_file);
            m_lastProjectsMenu->setTitle("Recent Projects...");

            for ( int i = 0 ; i < PALSProjectSettingsManager::sharedInstance()->getLastProjectPathList().size() ; ++ i )
            {
                if ( PALSProjectSettingsManager::sharedInstance()->getLastProjectPathList().at(i).isEmpty() )
                    continue;

                QAction *action = new QAction(PALSProjectSettingsManager::sharedInstance()->getLastProjectPathList().at(i), this);
                action->setIcon(QIcon(":/localImages/Images/IconPNGRounded.png"));

                connect(action, SIGNAL(triggered()), this, SLOT(openProjectFromLastPath()));
                m_lastProjectActionList.append(action);
            }

            m_lastProjectsMenu->addActions(m_lastProjectActionList);
            ui->menuLoad_file->addAction(m_lastProjectsMenu->menuAction());
        }
        else
            m_lastProjectsMenu = nullptr;
    }
    else
        m_lastProjectsMenu = nullptr;


    m_plotWindow = new DFastPlotDlg;
    m_resultWindow = new DFastResultDlg;
    m_calculatorWindow = new DFastCalculatorDlg;

    m_fitEngineThread = new QThread;
    m_fitEngine = new LifeTimeDecayFitEngine;
    m_fitEngine->moveToThread(m_fitEngineThread);

    connect(m_fitEngineThread, SIGNAL(started()), m_fitEngine, SLOT(fit()));
    connect(ui->pushButtonRunFit, SIGNAL(clicked()), this, SLOT(runFit()));
    connect(m_fitEngine, SIGNAL(finished()), this, SLOT(fitHasFinished()));

    m_chiSquareLabel = new QLabel;
    m_integralCountInROI = new QLabel;

    ui->statusBar->setStyleSheet("background-color: lightgray");

    ui->statusBar->addPermanentWidget(m_chiSquareLabel);
    ui->statusBar->addPermanentWidget(m_integralCountInROI);

    connect(ui->actionLoad, SIGNAL(triggered()), this, SLOT(openProject()));
    connect(ui->actionSave, SIGNAL(triggered()), this, SLOT(saveProject()));
    connect(ui->actionNew, SIGNAL(triggered()), this, SLOT(newProject()));
    connect(ui->actionSaveAs, SIGNAL(triggered()), this, SLOT(saveProjectAs()));
    connect(ui->actionImport, SIGNAL(triggered()), this, SLOT(importASCII()));

    connect(ui->widget, SIGNAL(dataChanged()), this, SLOT(instantPreview()));

    connect(ui->widget, SIGNAL(fitRangeChanged(int,int)), m_plotWindow, SLOT(setFitRange(int,int)));
    connect(ui->widget, SIGNAL(fitRangeChanged(int,int)), SLOT(instantPreview()));

    connect(ui->widget->fixedBackgroundCheckBox(), SIGNAL(clicked(bool)), this, SLOT(changeFixedBackground(bool)));

    connect(ui->actionPlot_Window, SIGNAL(triggered(bool)), this, SLOT(changePlotWindowVisibility(bool)));
    connect(m_plotWindow, SIGNAL(visibilityChanged(bool)), this, SLOT(changePlotWindowVisibilityFromOutside(bool)));

    connect(ui->actionResult_Window, SIGNAL(triggered(bool)), this, SLOT(changeResultWindowVisibility(bool)));
    connect(m_resultWindow, SIGNAL(visibilityChanged(bool)), this, SLOT(changeResultWindowVisibilityFromOutside(bool)));

    connect(ui->actionOpen_Calculator, SIGNAL(triggered(bool)), this, SLOT(changeCalculatorWindowVisibility(bool)));
    connect(m_calculatorWindow, SIGNAL(visibilityChanged(bool)), this, SLOT(changeCalculatorWindowVisibilityFromOutside(bool)));

    connect(ui->actionRaw_Data_Trace_2, SIGNAL(triggered(bool)), this, SLOT(changeRawDataTraceVisibility(bool)));
    connect(ui->actionStart_Value_Trace_2, SIGNAL(triggered(bool)), this, SLOT(changeStartValueTraceVisibility(bool)));
    connect(ui->actionFit_Trace_2, SIGNAL(triggered(bool)), this, SLOT(changeFitTraceVisibility(bool)));

    connect(m_resultWindow, SIGNAL(resultListIsEmpty()), this, SLOT(disablePDFExport()));
    connect(m_resultWindow, SIGNAL(resultListHasResults()), this, SLOT(enablePDFExport()));

    connect(ui->actionExport_Current_Result_as_PDF, SIGNAL(triggered()), m_resultWindow, SLOT(printToPDF()));
    connect(ui->actionExport_Current_Result_as_HTML, SIGNAL(triggered()), m_resultWindow, SLOT(printToHTML()));
    connect(ui->actionSave_Plot_as_Image, SIGNAL(triggered()), m_plotWindow, SLOT(savePlotAsImage()));
    connect(ui->actionAbout, SIGNAL(triggered()), this, SLOT(showAbout()));

    ui->pushButtonRunFit->setLiteralSVG(":/localImages/Images/arrowRight");
    ui->pushButtonRunFit->setStatusTip("Fit Lifetime-Data...");

#if defined(Q_OS_WIN)
    ui->label->setFont(WINDOWS_FONT(10));
    ui->label_2->setFont(WINDOWS_FONT(10));
#endif

    ui->actionExport_Current_Result_as_PDF->setIcon(QIcon(":/localImages/Images/pdfExport.svg"));
    ui->actionExport_Current_Result_as_HTML->setIcon(QIcon(":/localImages/Images/htmlExport.svg"));
    ui->actionSave_Plot_as_Image->setIcon(QIcon(":/localImages/Images/pngExport.svg"));

    ui->actionLoad->setIcon(QIcon(":/localImages/Images/open.svg"));
    ui->actionSave->setIcon(QIcon(":/localImages/Images/save.svg"));
    ui->actionSaveAs->setIcon(QIcon(":/localImages/Images/save.svg"));
    ui->actionNew->setIcon(QIcon(":/localImages/Images/new.svg"));
    ui->actionImport->setIcon(QIcon(":/localImages/Images/plot.svg"));
    ui->actionAbout->setIcon(QIcon(":/localImages/Images/IconPNGRounded.png"));
    ui->actionOpen_Calculator->setIcon(QIcon(":/localImages/Images/calculator73.svg"));

    QPixmap redPixmap(20, 20), greenPixmap(20, 20), bluePixmap(20, 20);
    redPixmap.fill(Qt::red);
    greenPixmap.fill(Qt::green);
    bluePixmap.fill(Qt::blue);

    ui->actionFit_Trace_2->setIcon(greenPixmap);
    ui->actionRaw_Data_Trace_2->setIcon(redPixmap);
    ui->actionStart_Value_Trace_2->setIcon(bluePixmap);

    ui->pushButtonRunFit->setToolTip("Fit Lifetime-Data...");

#if defined(Q_OS_WIN)
    ui->actionLoad->setShortcut(QKeySequence("Ctrl+L"));
    ui->actionNew->setShortcut(QKeySequence("Ctrl+N"));
    ui->actionSave->setShortcut(QKeySequence("Ctrl+S"));
    ui->actionImport->setShortcut(QKeySequence("Ctrl+I"));

    ui->actionPlot_Window->setShortcut(QKeySequence("Ctrl+P"));
    ui->actionResult_Window->setShortcut(QKeySequence("Ctrl+R"));

    ui->actionRaw_Data_Trace_2->setShortcut(QKeySequence("Ctrl+D"));
    ui->actionStart_Value_Trace_2->setShortcut(QKeySequence("Ctrl+T"));
    ui->actionFit_Trace_2->setShortcut(QKeySequence("Ctrl+F"));
    ui->actionOpen_Calculator->setShortcut(QKeySequence("Ctrl+Alt+C"));
#endif

#if defined(Q_OS_WIN)
    m_calculatorWindow->setTextFont(WINDOWS_FONT(10));
#else
    m_calculatorWindow->setTextFont(QFont("Helvetica", 12));
#endif

    //if ( projectPath.isEmpty() )
    newProject();

    m_onStart = true;

    if ( !PALSProjectSettingsManager::sharedInstance()->isLinearLastScaling() )
        m_plotWindow->setLogarithmicScaling();
    else
        m_plotWindow->setLinearScaling();

    ui->widget->setBackgroundChannelRange(PALSProjectSettingsManager::sharedInstance()->getLastBackgroundChannelRange());
    ui->widget->setBackgroundCalculationUsingFirstChannels(PALSProjectSettingsManager::sharedInstance()->getBackgroundCalculationFromFirstChannels());

    m_plotWindow->setYRangeData(1, 10000);

    if ( !projectPath.isEmpty() )
        openProjectFromPath(projectPath);

    if ( PALSProjectSettingsManager::sharedInstance()->getPlotWindowWasShownOnExit() )
        m_plotWindow->showMaximized();
    else
    {
        m_plotWindow->show();
        m_plotWindow->hide();
    }

    if ( PALSProjectSettingsManager::sharedInstance()->getResultWindowWasShownOnExit() )
        m_resultWindow->show();
    else
    {
        m_resultWindow->show();
        m_resultWindow->hide();
    }

    show();
}

DFastLTFitDlg::~DFastLTFitDlg()
{
    PALSProjectSettingsManager::sharedInstance()->setLinearAsLastScaling(m_plotWindow->isLinearScalingEnabled());
    PALSProjectSettingsManager::sharedInstance()->save();

    if ( m_lastProjectsMenu )
    {
        while ( m_lastProjectActionList.size() > 0 )
        {
            m_lastProjectsMenu->removeAction(m_lastProjectActionList.first());

            QAction *action = m_lastProjectActionList.takeFirst();
            DDELETE_SAFETY(action);
        }

        DDELETE_SAFETY(m_lastProjectsMenu);
    }

    DDELETE_SAFETY(m_resultWindow);
    DDELETE_SAFETY(m_plotWindow);
    DDELETE_SAFETY(m_calculatorWindow);

    DDELETE_SAFETY(m_chiSquareLabel);
    DDELETE_SAFETY(m_integralCountInROI);

    DDELETE_SAFETY(m_fitEngine);
    DDELETE_SAFETY(m_fitEngineThread);

    DDELETE_SAFETY(ui);
}

void DFastLTFitDlg::closeEvent(QCloseEvent *event)
{
    event->ignore();

    const QMessageBox::StandardButton replyBtn = QMessageBox::question(this, "Closing DQuickLTFit?",
                                                                       "<nobr>Did you save the project?</nobr>",
                                                                       QMessageBox::Yes|QMessageBox::No);

    if ( replyBtn == QMessageBox::StandardButton::No )
    {
        event->ignore();
        return;
    }

    event->accept();

    if ( m_resultWindow->isVisible() )
        PALSProjectSettingsManager::sharedInstance()->setResultWindowWasShownOnExit(true);
    else
        PALSProjectSettingsManager::sharedInstance()->setResultWindowWasShownOnExit(false);

    if ( m_plotWindow->isVisible() )
        PALSProjectSettingsManager::sharedInstance()->setPlotWindowWasShownOnExit(true);
    else
        PALSProjectSettingsManager::sharedInstance()->setPlotWindowWasShownOnExit(false);

    PALSProjectSettingsManager::sharedInstance()->save();

    m_plotWindow->close();
    m_resultWindow->close();
    m_calculatorWindow->close();

    QMainWindow::closeEvent(event);
}

void DFastLTFitDlg::changeFixedBackground(bool fixed)
{
    PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getBackgroundParamPtr()->getParameter()->setAsFixed(fixed);
}

void DFastLTFitDlg::openProject()
{
    const QString fileName = QFileDialog::getOpenFileName(this, tr("Open a project"),
                                                    PALSProjectSettingsManager::sharedInstance()->getLastChosenPath(),
                                                    QString("DQuickLTFit Project File (*" % PROJECT_EXTENSION % ")"));

    openProjectFromPath(fileName);
}

void DFastLTFitDlg::openProjectFromLastPath()
{
    openProjectFromPath(((QAction*)sender())->text());
}

void DFastLTFitDlg::openProjectFromPath(const QString& fileName)
{
    if ( fileName.isEmpty() )
        return;

    PALSProjectSettingsManager::sharedInstance()->setLastChosenPath(QFileInfo(fileName).absoluteDir().absolutePath());

    if ( PALSProjectManager::sharedInstance()->getFileName() == fileName )
    {
        DMSGBOX("<nobr>This project is already open!</nobr>");
        return;
    }

    if ( !PALSProjectManager::sharedInstance()->load(fileName) )
    {
        DMSGBOX("<nobr>Sorry, an error occurred while loading this project!</nobr>");
        return;
    }
    else
    {
        if ( !PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getLifeTimeData().isEmpty() )
        {
            int minChn = INT_MAX;
            int maxChn = -INT_MAX;
            int minCnts = INT_MAX;
            int maxCnts = -INT_MAX;

            for ( int i = 0 ; i < PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getLifeTimeData().size() ; ++ i )
            {
                const int channel = (const int)PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getLifeTimeData().at(i).x();
                const int counts = (const int)PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getLifeTimeData().at(i).y();

                maxChn = qMax(channel, maxChn);
                minChn = qMin(channel, minChn);
                maxCnts = qMax(counts, maxCnts);
                minCnts = qMin(counts, minCnts);
            }

            PALSProjectManager::sharedInstance()->setChannelRanges(minChn, maxChn);

            m_plotWindow->clearAll();
            m_plotWindow->setXRange(minChn, maxChn);

            m_plotWindow->addRawData(PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getLifeTimeData());
            m_plotWindow->addFitData(PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getFitData());
            m_plotWindow->addResidualData(PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getResiduals());

            m_plotWindow->setXRange(minChn, maxChn);

            ui->widget->setFitRangeLimits(minChn, maxChn);
            ui->widget->setFitRange(minChn, maxChn);

            m_plotWindow->setYRangeData(1, 1.3*((double)maxCnts));
        }
        else
        {
            PALSProjectManager::sharedInstance()->setChannelRanges(PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel(), PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel());

            m_plotWindow->clearAll();
            m_plotWindow->setXRange(PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel(), PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel());
            m_plotWindow->setYRangeData(1, 10000);

            ui->widget->setFitRangeLimits(0, 10000);
            ui->widget->setFitRange(0, 10000);
        }

        PALSProjectManager::sharedInstance()->setFileName(fileName);
        PALSProjectSettingsManager::sharedInstance()->addLastProjectPathToList(fileName);

        ui->widget->updateParamterList();
        m_resultWindow->clearTabs();

        if ( !PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getFitData().isEmpty() )
            m_resultWindow->addResultTabsFromHistory();

        disconnect(ui->widget->fixedBackgroundCheckBox(), SIGNAL(clicked(bool)), this, SLOT(changeFixedBackground(bool)));
            ui->widget->fixedBackgroundCheckBox()->setChecked(PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getBackgroundParamPtr()->getParameter()->isFixed());
        connect(ui->widget->fixedBackgroundCheckBox(), SIGNAL(clicked(bool)), this, SLOT(changeFixedBackground(bool)));
    }

    updateLastProjectActionList();
    updateWindowTitle();    
}

void DFastLTFitDlg::saveProject()
{
    QString filename = "";
    if ( PALSProjectManager::sharedInstance()->getFileName().isEmpty() )
    {
        filename = QFileDialog::getSaveFileName(this, tr("Select or type a filename..."),
                                                PALSProjectSettingsManager::sharedInstance()->getLastChosenPath(),
                                                QString("DQuickLTFit Project File (*" % PROJECT_EXTENSION % ")"));

        if ( filename.isEmpty() )
            return;
        else
        {
            PALSProjectSettingsManager::sharedInstance()->setLastChosenPath(QFileInfo(filename).absoluteDir().absolutePath());

            PALSProjectManager::sharedInstance()->setFileName(filename);
            PALSProjectSettingsManager::sharedInstance()->addLastProjectPathToList(filename);

            updateLastProjectActionList();
        }
    }


    if ( PALSProjectManager::sharedInstance()->save(PALSProjectManager::sharedInstance()->getFileName()) ){
        //DMSGBOX("The project was saved successfully.");
    }else{
        DMSGBOX("Sorry, an error occurred while saving this project.");
    }

    updateWindowTitle();
}

void DFastLTFitDlg::saveProjectAs()
{
    const QString filename = QFileDialog::getSaveFileName(this, tr("Select or type a filename..."),
                                                          PALSProjectSettingsManager::sharedInstance()->getLastChosenPath(),
                                                          QString("DQuickLTFit Project File (*" % PROJECT_EXTENSION % ")"));

    if ( filename.isEmpty() )
        return;
    else
    {
        PALSProjectSettingsManager::sharedInstance()->setLastChosenPath(QFileInfo(filename).absoluteDir().absolutePath());

        PALSProjectManager::sharedInstance()->setFileName(filename);
        PALSProjectSettingsManager::sharedInstance()->addLastProjectPathToList(filename);

        updateLastProjectActionList();
    }


    if ( PALSProjectManager::sharedInstance()->save(PALSProjectManager::sharedInstance()->getFileName()) ){
        //DMSGBOX("The project was saved successfully.");
    }else{
        DMSGBOX("Sorry, an error occurred while saving this project.");
    }

    updateWindowTitle();
}

void DFastLTFitDlg::newProject()
{
    PALSProjectManager::sharedInstance()->createEmptyProject();
    PALSProjectManager::sharedInstance()->setFileName("");

    PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getBackgroundParamPtr()->getParameter()->setAsFixed(ui->widget->fixedBackgroundCheckBox()->isChecked());

    m_plotWindow->clearAll();
    m_plotWindow->setXRange(PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel(), PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel());
    m_plotWindow->setYRangeData(1, 10000);

    PALSProjectManager::sharedInstance()->setChannelRanges(PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel(), PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel());

    ui->widget->setFitRangeLimits(PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel(), PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel());
    ui->widget->setFitRange(PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel(), PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel());

    ui->widget->updateParamterList();

    m_resultWindow->clearTabs();

    updateWindowTitle();
}

QStringList DFastLTFitDlg::autoDetectDelimiter(const QString& row)
{
    if ( row.split(";").size() != 2 )
    {
        if ( row.split("|").size() != 2 )
        {
            if ( row.split("\t").size() < 2 )
            {
                if ( row.split(" ").size() >= 2 )
                {
                    QStringList returnList;
                    const QStringList list = row.split(" ");
                    const int splitSize = list.size();

                    int cnt = 0;
                    for ( int i = 0 ; i < splitSize ; ++ i )
                    {
                        if ( list.at(i).isEmpty() )
                            continue;
                        else
                        {
                            bool ok = false;
                            const int value = (int)QVariant(list.at(i).trimmed()).toDouble(&ok);
                            DUNUSED_PARAM(value);

                            if ( ok )
                            {
                                returnList.append(list.at(i).trimmed());
                                cnt ++;
                            }

                            if ( cnt == 2 )
                                break;
                        }
                    }


                    return returnList;
                }
                else
                {
                    QStringList returnList;

                    bool ok = false;
                    const int value = (int)QVariant(row.trimmed()).toDouble(&ok);
                    DUNUSED_PARAM(value);

                    if ( ok )
                    {
                        returnList.append(row.trimmed());
                    }


                    return returnList;
                }
            }
            else
            {
                QStringList returnList;
                const QStringList list = row.split("\t");
                const int splitSize = list.size();

                int cnt = 0;
                for ( int i = 0 ; i < splitSize ; ++ i )
                {
                    if ( list.at(i).isEmpty() )
                        continue;
                    else
                    {
                        bool ok = false;
                        const int value = (int)QVariant(list.at(i).trimmed()).toDouble(&ok);
                        DUNUSED_PARAM(value);

                        if ( ok )
                        {
                            returnList.append(list.at(i).trimmed());
                            cnt ++;
                        }

                        if ( cnt == 2 )
                            break;
                    }
                }


                return returnList;
            }
        }
        else
        {
            const QStringList list = row.split("|");
            const QString value1 = list.at(0).trimmed();
            const QString value2 = list.at(1).trimmed();

            QStringList returnList;
            returnList.append(value1);
            returnList.append(value2);


            return returnList;
        }
    }
    else
    {
        const QStringList list = row.split(";");
        const QString value1 = list.at(0).trimmed();
        const QString value2 = list.at(1).trimmed();

        QStringList returnList;
        returnList.append(value1);
        returnList.append(value2);


        return returnList;
    }

    return QStringList();
}

void DFastLTFitDlg::importASCII(const AccessType& type, const QString& fileNameFromSeq)
{
    QString fileName = "";
    int binFac = 1;

    if ( type == AccessType::FromOneFile ) {
        QFileDialog *fd = new QFileDialog;
        fd->setWindowTitle(tr("Import data from ASCII File..."));
        fd->setAcceptMode(QFileDialog::AcceptOpen);
        fd->setDirectory(PALSProjectSettingsManager::sharedInstance()->getLastChosenPath());
        fd->setFileMode(QFileDialog::ExistingFile);
        fd->setViewMode(QFileDialog::Detail);
        fd->setNameFilter(tr("Lifetime Data (*.dat *.txt *.log)"));

        fd->setOption(QFileDialog::DontUseNativeDialog, true);

        QGridLayout *layout = static_cast<QGridLayout*>(fd->layout());

        QHBoxLayout *hboxLayout = new QHBoxLayout;

        QLabel *label = new QLabel("Bin-Factor?");
        QSpinBox *binFacBox = new QSpinBox;
        binFacBox->setRange(1, 100);
        binFacBox->setSingleStep(1);
        binFacBox->setButtonSymbols(QAbstractSpinBox::NoButtons);

        hboxLayout->addWidget(label);
        hboxLayout->addWidget(binFacBox);

        layout->addLayout(hboxLayout, layout->rowCount()-1, layout->columnCount());

        if ( !fd->exec() )
            return;

        fileName = !fd->selectedFiles().isEmpty()?fd->selectedFiles().first():QString("");
        binFac = binFacBox->value();

        DDELETE_SAFETY(label);
        DDELETE_SAFETY(binFacBox);
        DDELETE_SAFETY(hboxLayout);
        DDELETE_SAFETY(fd);

        /*fileName = QFileDialog::getOpenFileName(widget, tr("Import data from ASCII File..."),
                                                PALSProjectSettingsManager::sharedInstance()->getLastChosenPath(),
                                                tr("Lifetime Data (*.dat *.txt *.log)"));*/


        if ( fileName.isEmpty() )
            return;


        PALSProjectSettingsManager::sharedInstance()->setLastChosenPath(QFileInfo(fileName).absoluteDir().absolutePath());
    }
    else
        fileName = fileNameFromSeq;


    QFile file(fileName);

    if ( file.open(QIODevice::ReadOnly) ) {
        QList<QPointF> dataSet;
        int minChn = INT_MAX;
        int  maxChn = -INT_MAX;
        int minCnts = INT_MAX;
        int maxCnts = -INT_MAX;

        int channelCounter = 0;
        int channel = 0;
        int counts = 0;

        while ( !file.atEnd() ) {
            QString dataRow = file.readLine();

            const QStringList dataSetString = autoDetectDelimiter(dataRow);

            if ( dataSetString.size() != 2
                 && dataSetString.size() != 1 )
                continue;

            bool ok_1 = true, ok_2 = false;

            channelCounter ++;

            if ( dataSetString.size() == 2 ) {
                int tchannel = (int)QVariant(dataSetString.at(0)).toInt(&ok_1);
                DUNUSED_PARAM(tchannel);
            }

            if ( dataSetString.size() == 2 )
                counts += (int)QVariant(dataSetString.at(1)).toInt(&ok_2);
            else if ( dataSetString.size() == 1 )
                counts += (int)QVariant(dataSetString.at(0)).toInt(&ok_2);

            if ( !ok_1 || !ok_2 )
                continue;

            maxChn = qMax(channel, maxChn);
            minChn = qMin(channel, minChn);
            maxCnts = qMax(counts, maxCnts);
            minCnts = qMin(counts, minCnts);

            if ( counts < 0 ) {
                if ( type == AccessType::FromOneFile ) {
                    DMSGBOX("Please correct the content of this file. Values lower than 0 detected.")
                }

                file.close();
                return;
            }

            if ( !(channelCounter%binFac) ) {
                dataSet.append(QPointF(channel, counts));
                counts = 0;
                channel ++;
            }
        }

        file.close();

        if ( dataSet.size() <= 2 ) {
            DMSGBOX("Either the Number of Data-Points was too low or the Bin-Factor is too high!");
            return;
        }


        PALSProjectManager::sharedInstance()->setChannelRanges(minChn, maxChn);

        PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->clearFitData();
        PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->clearResidualData();

        PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->setLifeTimeData(dataSet);
        PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->setBinFactor(binFac);

        m_plotWindow->clearAll();

        m_plotWindow->setXRange(minChn, maxChn);
        m_plotWindow->addRawData(dataSet);
        m_plotWindow->setXRange(minChn, maxChn);

        int newStartChannel = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel();
        if (  newStartChannel < minChn )
            newStartChannel = minChn;

        int newStopChannel = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel();
        if (  newStopChannel > maxChn )
            newStopChannel = maxChn;

        ui->widget->setFitRangeLimits(minChn, maxChn);
        ui->widget->setFitRange(newStartChannel, newStopChannel);

        m_plotWindow->setYRangeData(1, 1.3*((double)maxCnts));

        instantPreview();
    }
    else
    {
        if ( type == AccessType::FromOneFile )
        {
            DMSGBOX("Sorry, an error occurred while importing lifetime-data.")
        }
    }

    PALSProjectManager::sharedInstance()->setASCIIDataName(fileName);

    updateWindowTitle();
}

void DFastLTFitDlg::runFit()
{
    if ( PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getLifeTimeData().isEmpty() )
    {
        DMSGBOX("<nobr>No data for fitting: Please import lifetime data before.</nobr>");
        return;
    }

    QList<QString> sourceConflictList;
    QList<QString> sampleConflictList;
    QList<QString> deviceConflictList;

    for ( int a = 0 ; a < PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getSourceParamPtr()->getSize() ; ++ a )
    {
        const PALSFitParameter *param = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getSourceParamPtr()->getParameterAt(a);

        if ( param->isFixed() && (param->isLowerBoundingEnabled() || param->isUpperBoundingEnabled()) )
            sourceConflictList.append(param->getAlias());
    }

    for ( int a = 0 ; a < PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getLifeTimeParamPtr()->getSize() ; ++ a )
    {
        const PALSFitParameter *param = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getLifeTimeParamPtr()->getParameterAt(a);

        if ( param->isFixed() && (param->isLowerBoundingEnabled() || param->isUpperBoundingEnabled()) )
            sampleConflictList.append(param->getAlias());
    }

    for ( int a = 0 ; a < PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize() ; ++ a )
    {
        const PALSFitParameter *param = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(a);

        if ( param->isFixed() && (param->isLowerBoundingEnabled() || param->isUpperBoundingEnabled()) )
            deviceConflictList.append(param->getAlias());
    }

    if ( !sourceConflictList.isEmpty() || !sampleConflictList.isEmpty() || !deviceConflictList.isEmpty() )
    {
        QString conflictText = "There are parameter conflicts: <br><br><b>SOURCE:</b> ";

        for ( int i = 0 ; i <  sourceConflictList.size() ; ++ i)
        {
            conflictText.append(sourceConflictList.at(i));

            if ( i != sourceConflictList.size()-1 )
                conflictText.append(", ");
        }

        conflictText.append("<br><br><b>SAMPLE:</b> ");

        for ( int i = 0 ; i <  sampleConflictList.size() ; ++ i)
        {
            conflictText.append(sampleConflictList.at(i));

            if ( i != sampleConflictList.size()-1 )
                conflictText.append(", ");
        }

        conflictText.append("<br><br><b>IRF:</b> ");

        for ( int i = 0 ; i <  deviceConflictList.size() ; ++ i)
        {
            conflictText.append(deviceConflictList.at(i));

            if ( i != deviceConflictList.size()-1 )
                conflictText.append(", ");
        }

        conflictText.append("<br><br>The parameter can either be <b>fixed</b> or <b>has limits</b>.");

        DMSGBOX(conflictText);
        return;
    }

    enableGUI(false);

    m_fitEngine->init(PALSProjectManager::sharedInstance()->getDataStructure());
    m_fitEngineThread->start();
}

void DFastLTFitDlg::fitHasFinished()
{
    m_fitEngineThread->exit(0);

    instantPreview();

    m_plotWindow->clearFitData();
    m_plotWindow->addFitData(m_fitEngine->getFitPlotPoints());

    m_plotWindow->clearResidualData();
    m_plotWindow->addResidualData(PALSProjectManager::sharedInstance()->getDataStructure()->getDataSetPtr()->getResiduals());

    m_resultWindow->addResultTabFromLastFit();

    enableGUI(true);
}

void DFastLTFitDlg::updateWindowTitle()
{
    if ( !PALSProjectManager::sharedInstance()->getFileName().isEmpty() )
    {
        this->setWindowTitle("Scope - " % VERSION_STRING_AND_PROGRAM_NAME % " - " % PALSProjectManager::sharedInstance()->getFileName());
        m_plotWindow->setWindowTitle("Plot - " % VERSION_STRING_AND_PROGRAM_NAME % " - " % PALSProjectManager::sharedInstance()->getFileName());
        m_resultWindow->setWindowTitle("Results - " % VERSION_STRING_AND_PROGRAM_NAME % " - " % PALSProjectManager::sharedInstance()->getFileName());
    }
    else
    {
        this->setWindowTitle("Scope - " % VERSION_STRING_AND_PROGRAM_NAME % " - <empty project>");
        m_plotWindow->setWindowTitle("Plot - " % VERSION_STRING_AND_PROGRAM_NAME % " - <empty project>");
        m_resultWindow->setWindowTitle("Results - " % VERSION_STRING_AND_PROGRAM_NAME % " - <empty project>");
    }

    m_calculatorWindow->setWindowTitle("Calculator - " % VERSION_STRING_AND_PROGRAM_NAME);
}

void DFastLTFitDlg::updateLastProjectActionList()
{
    if ( m_lastProjectsMenu )
    {
        const int size = m_lastProjectActionList.size();
        int cnt = 0;

        if ( size != 0 ) {
            while ( cnt < size )
            {
                m_lastProjectsMenu->removeAction(m_lastProjectActionList.first());

                QAction *action = m_lastProjectActionList.takeFirst();

                if (action) {
                    disconnect(action, SIGNAL(triggered()), this, SLOT(openProjectFromLastPath()));
                    DDELETE_SAFETY(action);
                }

                cnt ++;
            }
        }

        ui->menuLoad_file->removeAction(m_lastProjectsMenu->menuAction());

        DDELETE_SAFETY(m_lastProjectsMenu);
        m_lastProjectActionList.clear();

        m_lastProjectsMenu = new QMenu(ui->menuLoad_file);
        m_lastProjectsMenu->setTitle("Recent Projects...");

        for ( int i = 0 ; i < PALSProjectSettingsManager::sharedInstance()->getLastProjectPathList().size() ; ++ i )
        {
            QAction *action = new QAction(PALSProjectSettingsManager::sharedInstance()->getLastProjectPathList().at(i), this);
            action->setIcon(QIcon(":/localImages/Images/IconPNGRounded.png"));

            connect(action, SIGNAL(triggered()), this, SLOT(openProjectFromLastPath()));
            m_lastProjectActionList.append(action);
        }

        m_lastProjectsMenu->addActions(m_lastProjectActionList);
        ui->menuLoad_file->addAction(m_lastProjectsMenu->menuAction());
    }
}

void DFastLTFitDlg::enablePDFExport()
{
    ui->actionExport_Current_Result_as_PDF->setEnabled(true);
    ui->actionExport_Current_Result_as_HTML->setEnabled(true);
}

void DFastLTFitDlg::disablePDFExport()
{
    ui->actionExport_Current_Result_as_PDF->setEnabled(false);
    ui->actionExport_Current_Result_as_HTML->setEnabled(false);
}

void DFastLTFitDlg::showAbout()
{
    const QString text = VERSION_STRING_AND_PROGRAM_NAME % " (" % VERSION_RELEASE_DATE % ") <br><br>" % COPYRIGHT_NOTICE % "<br><br>";
    const QString contact = "contact: <a href=\"danny.petschke@uni-wuerzburg.de\">danny.petschke@uni-wuerzburg.de</a><br><br>";

    const QString license = "<nobr>Fit-algorithm provided by: <br>MPFIT: A MINPACK-1 Least Squares Fitting Library in C</nobr><br><br>";
    const QString license2 = "<nobr>Icons provided by: <br>https://www.flaticon.com (flaticon)</nobr><br><br>";
    const QString license3 = "<nobr>Logo designed by Hannah Heil</nobr>";

    QMessageBox::about(this, "DQuickLTFit", text % contact % license % license2 % license3);
}

void DFastLTFitDlg::instantPreview()
{
    PALSDataStructure *dataStructure = PALSProjectManager::sharedInstance()->getDataStructure();

    if ( dataStructure->getDataSetPtr()->getLifeTimeData().isEmpty() )
        return;

    const int paramCnt = dataStructure->getFitSetPtr()->getComponentsCount()+dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize() + 1;

    double channelResolution = dataStructure->getFitSetPtr()->getChannelResolution();

    double startChannel = dataStructure->getFitSetPtr()->getStartChannel();
    double stopChannel = dataStructure->getFitSetPtr()->getStopChannel();
    double peakChannel = 0;
    double countsInPeak = -(double)INT_MAX;

     int startChannelIndex = 0;
     int stopChannelIndex = 0;
     int peakChannelIndex = 0;

     int channelCnt = 0;
     int dataCntInRange = (stopChannel-startChannel+1);

     double *x = new double[dataCntInRange];
     double *y = new double[dataCntInRange];
     double *ey = new double[dataCntInRange];

     int inRangeCnt = 0;
     int integralCounts = 0;
     int tZero = 0;

     for ( QPointF p : dataStructure->getDataSetPtr()->getLifeTimeData() ) {
         if ( ((int)p.x()) >= ((int)startChannel) && ((int)p.x()) <= ((int)stopChannel) ) {
             x[inRangeCnt] = p.x();
             y[inRangeCnt] = p.y();

             /* calculate error (weighting) (statistical weighting) */
             ey[inRangeCnt] = 1.0/sqrt(p.y() + 1.0); //prevent zero division

             integralCounts += (int)p.y();

             if ( ((int)p.x()) == ((int)startChannel) )
                 startChannelIndex = channelCnt;

             if ( ((int)p.x()) == ((int)stopChannel) )
                 stopChannelIndex = channelCnt;

             if ( ((int)p.y()) > ((int)countsInPeak) ) {
                 countsInPeak = p.y();
                 peakChannel = p.x();
                 peakChannelIndex = channelCnt;
                 tZero = peakChannel;
             }

             inRangeCnt ++;
         }

         channelCnt ++;
     }

     tZero -= startChannel;

    double *params = new double[paramCnt]; /* following order: source => sample => gaussian => bkgrd */

    int i = 0;
    for ( i = 0 ; i < dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i+=2 ) {
        const PALSFitParameter *fitParam = dataStructure->getFitSetPtr()->getSourceParamPtr()->getParameterAt(i);
        params[i] = fitParam->getStartValue()/channelResolution; //tau

        const PALSFitParameter *fitParam2 = dataStructure->getFitSetPtr()->getSourceParamPtr()->getParameterAt(i+1);
        params[i+1] = fitParam2->getStartValue(); //I
    }

    int cnt = 0;
    for ( i = dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i < dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i+=2 ) {
        const PALSFitParameter *fitParam = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getParameterAt(cnt);
        params[i] = fitParam->getStartValue()/channelResolution; //tau

        cnt ++;

        const PALSFitParameter *fitParam2 = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getParameterAt(cnt);
        params[i+1] = fitParam2->getStartValue(); //I

        cnt ++;
    }

    int cntGaussian= 0;
    for ( i = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() ; i < dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() + dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize() ; i+=3 ) {
        const PALSFitParameter *fitParam = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        params[i] = fitParam->getStartValue()/channelResolution; //sigma

        cntGaussian ++;

        const PALSFitParameter *fitParam2 = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        params[i+1] = fitParam2->getStartValue()/channelResolution; //mu

        cntGaussian ++;

        const PALSFitParameter *fitParam3 = dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getParameterAt(cntGaussian);
        params[i+2] = fitParam3->getStartValue(); //I

        cntGaussian ++;
    }

    const PALSFitParameter *bkgrd = dataStructure->getFitSetPtr()->getBackgroundParamPtr()->getParameter();

    const int bkgrdIndex = dataStructure->getFitSetPtr()->getLifeTimeParamPtr()->getSize()+dataStructure->getFitSetPtr()->getSourceParamPtr()->getSize() + dataStructure->getFitSetPtr()->getDeviceResolutionParamPtr()->getSize();

    params[bkgrdIndex] = bkgrd->getStartValue();

    countsInPeak -= bkgrd->getStartValue();

    const double bkgrdVal = params[bkgrdIndex];


    QList<QPointF> fitPlotSet;

    double residuals = 0.0;
    const int reducedCntInRange = (dataCntInRange - 1);
    const int reducedDevCount = (paramCnt - cntGaussian - 1);
    const int reducedParamCount = (paramCnt - 1);
    const double integralCountsWithoutBkgrd = (double)integralCounts-(double)dataCntInRange*bkgrdVal;

    for ( int i = 0 ; i < reducedCntInRange ; ++ i ) {
        double f = 0.0;

        x[i] -= startChannel;
        const double x_plus_1 = x[i+1]-startChannel;

        for ( int device = reducedDevCount ; device < reducedParamCount ; device += 3 ) {
            const double gaussianSigmaVal = params[device]/(2*sqrt(log(2))); /* transform FWHM to 1-sigma uncertainty */
            const double gaussianMuVal = params[device+1];

            const double gaussianIntensity = params[device+2]; /* IRF contribution */

            double valF = 0.0;

            /* Kirkegaard and Eldrup (1972) */
            for ( int param = 0 ; param <  reducedDevCount ; param += 2 ) { /* 1st param[0] = tau; 2nd param[1] = Intensity */
                const double yji = exp(-(x[i]-gaussianMuVal-(gaussianSigmaVal*gaussianSigmaVal)/(4*params[param]))/params[param])*(1-erf((0.5*gaussianSigmaVal/params[param])-(x[i]-gaussianMuVal)/gaussianSigmaVal));
                const double yji_plus_1 = exp(-(x_plus_1-gaussianMuVal-(gaussianSigmaVal*gaussianSigmaVal)/(4*params[param]))/params[param])*(1-erf((0.5*gaussianSigmaVal/params[param])-(x_plus_1-gaussianMuVal)/gaussianSigmaVal));

                valF += 0.5*params[param+1]*(yji-yji_plus_1-erf((x[i]-gaussianMuVal)/gaussianSigmaVal)+erf((x_plus_1-gaussianMuVal)/gaussianSigmaVal));
            }

            valF *= gaussianIntensity;
            f += valF;
        }

        x[i] += startChannel;

        f *= integralCountsWithoutBkgrd;
        f += (double)bkgrdVal;

        residuals += (y[i]-f)*(y[i]-f)*ey[i]*ey[i];

        fitPlotSet.append(QPointF(x[i], f));
    }

    /*approximated reduced chi-square (number of free parameters is not taken into account) */
    if ( dataCntInRange == 0 )
        residuals = -1;
    else
        residuals /= ((double)dataCntInRange);

    m_plotWindow->clearPreviewData();
    m_plotWindow->addPreviewData(fitPlotSet);
    m_plotWindow->updateBkgrdData();

    m_plotWindow->setFitRange(startChannel, stopChannel);

    if ( residuals == -1 ) {
        m_integralCountInROI->setText("");
        m_chiSquareLabel->setText("");
    }
    else {
        m_integralCountInROI->setText("estimated t<sub>0</sub>: <b>" % QVariant(tZero*channelResolution).toString() % "ps</b> Integral Cnts. ROI [" % QVariant(startChannel).toString() % ":" % QVariant(stopChannel).toString() % "]: <b>" % QVariant(integralCounts).toString() % "</b>");
        m_chiSquareLabel->setText("approx. &#967;<sub>&#957;</sub><sup>2</sup> ( @ start ): <b>" % QString::number(residuals, 'g', 3) % "</b>"); //appr. reduced chi-square
    }

    delete [] x;
    delete [] y;
    delete [] ey;
    delete [] params;
}

void DFastLTFitDlg::changePlotWindowVisibility(bool visible)
{
    if ( !visible )
        m_plotWindow->hide();
    else
        m_plotWindow->show();
}

void DFastLTFitDlg::changePlotWindowVisibilityFromOutside(bool visible)
{
    disconnect(ui->actionPlot_Window, SIGNAL(triggered(bool)), this, SLOT(changePlotWindowVisibility(bool)));

    ui->actionPlot_Window->setChecked(visible);

    connect(ui->actionPlot_Window, SIGNAL(triggered(bool)), this, SLOT(changePlotWindowVisibility(bool)));
}

void DFastLTFitDlg::changeResultWindowVisibility(bool visible)
{
    if ( !visible )
        m_resultWindow->hide();
    else
        m_resultWindow->show();
}

void DFastLTFitDlg::changeResultWindowVisibilityFromOutside(bool visible)
{
    disconnect(ui->actionResult_Window, SIGNAL(triggered(bool)), this, SLOT(changeResultWindowVisibility(bool)));

    ui->actionResult_Window->setChecked(visible);

    connect(ui->actionResult_Window, SIGNAL(triggered(bool)), this, SLOT(changeResultWindowVisibility(bool)));
}

void DFastLTFitDlg::changeCalculatorWindowVisibility(bool visible)
{
    if ( !visible )
        m_calculatorWindow->hide();
    else
        m_calculatorWindow->show();
}

void DFastLTFitDlg::changeCalculatorWindowVisibilityFromOutside(bool visible)
{
    disconnect(ui->actionOpen_Calculator, SIGNAL(triggered(bool)), this, SLOT(changeCalculatorWindowVisibility(bool)));

    ui->actionOpen_Calculator->setChecked(visible);

    connect(ui->actionOpen_Calculator, SIGNAL(triggered(bool)), this, SLOT(changeCalculatorWindowVisibility(bool)));
}

void DFastLTFitDlg::changeRawDataTraceVisibility(bool visible)
{
    m_plotWindow->setRawDataVisible(visible);
}

void DFastLTFitDlg::changeStartValueTraceVisibility(bool visible)
{
    m_plotWindow->setStartValueDataVisible(visible);
}

void DFastLTFitDlg::changeFitTraceVisibility(bool visible)
{
    m_plotWindow->setFitDataVisible(visible);
}

void DFastLTFitDlg::calculateBackground()
{
    ((ParameterListView*)ui->widget)->updateBackgroundValue();
}

void DFastLTFitDlg::enableGUI(bool enable)
{
    setEnabled(enable);
    ui->widget->setEnabled(enable);
    ui->pushButtonRunFit->enableWidget(enable);

    if (enable) {
        ui->pushButtonRunFit->setStatusTip("Fit Lifetime-Data...");
        ui->label->setText("Fit");
        ui->label_2->setText("Data");
        ui->pushButtonRunFit->setLiteralSVG(":/localImages/Images/arrowRight");
    }
    else {
        ui->pushButtonRunFit->setStatusTip("Fit is Running..");
        ui->label->setText("Fit is Running");
        ui->label_2->setText("!");
        ui->pushButtonRunFit->setLiteralSVG(":/localImages/Images/fit");
    }
}

void DFastLTFitDlg::printToFile(const QString &fileName, const QList<QPointF> &vec)
{
    QFile file(fileName + ".in");

    if (file.open(QIODevice::ReadWrite | QIODevice::Append)) {
        QTextStream stream(&file);
        stream << fileName + "\n\r";
        for (int i = 0 ; i < vec.size() ; i += 8) {
            QString text("       " + QVariant(i).toString() + "        \n\r");

            for (int u = 0 ; u < 8 ; ++ u) {
                text.append(QVariant(vec.at((i) + u).y()).toString() + "        \n\r");
            }

            stream << text << endl;
        }

        file.close();
    }
}
