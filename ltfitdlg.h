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

#ifndef DFASTLTFITDLG_H
#define DFASTLTFITDLG_H

#include <QMainWindow>
#include <QFileDialog>
#include <QInputDialog>
#include <QLabel>
#include <QSpinBox>
#include <QAction>
#include <QDebug>

#include "Settings/settings.h"
#include "Settings/projectmanager.h"
#include "Settings/projectsettingsmanager.h"

#include "ltplotdlg.h"
#include "ltresultdlg.h"
#include "ltcalculatordlg.h"

#include "Fit/lifetimedecayfit.h"

#include "ltdefines.h"

typedef enum {
    FromOneFile,
    FromSequence
} AccessType;

class DFastFileSequencerDlg;

namespace Ui {
class DFastLTFitDlg;
}

class DFastLTFitDlg : public QMainWindow
{
    Q_OBJECT
public:
    explicit DFastLTFitDlg(const QString projectPath = "", QWidget *parent = nullptr);
    virtual ~DFastLTFitDlg();

protected:
    virtual void closeEvent(QCloseEvent *event);

public slots:
    void changeFixedBackground(bool fixed);

    void openProject();
    void saveProject();
    void saveProjectAs();
    void newProject();

    void importASCII(const AccessType& type = AccessType::FromOneFile, const QString &fileNameFromSeq = "");

    void runFit();
    void instantPreview();

    void changePlotWindowVisibility(bool visible);
    void changePlotWindowVisibilityFromOutside(bool visible);

    void changeResultWindowVisibility(bool visible);
    void changeResultWindowVisibilityFromOutside(bool visible);

    void changeCalculatorWindowVisibility(bool visible);
    void changeCalculatorWindowVisibilityFromOutside(bool visible);

    void changeRawDataTraceVisibility(bool visible);
    void changeStartValueTraceVisibility(bool visible);
    void changeFitTraceVisibility(bool visible);

    void calculateBackground();

    void enableGUI(bool enable);

    void printToFile(const QString& fileName, const QList<QPointF>& vec);

private slots:
    void fitHasFinished();
    void updateWindowTitle();

    void openProjectFromPath(const QString& fileName);
    void openProjectFromLastPath();
    void updateLastProjectActionList();

    void enablePDFExport();
    void disablePDFExport();

    void showAbout();

public:
    static QStringList autoDetectDelimiter(const QString& row);

private:
    Ui::DFastLTFitDlg *ui;

    DFastPlotDlg *m_plotWindow;
    DFastResultDlg *m_resultWindow;
    DFastCalculatorDlg *m_calculatorWindow;

    LifeTimeDecayFitEngine *m_fitEngine;
    QThread *m_fitEngineThread;

    QLabel *m_chiSquareLabel;
    QLabel *m_integralCountInROI;

    QMenu *m_lastProjectsMenu;
    QList<QAction*> m_lastProjectActionList;

    bool m_onStart;
};

#endif // DFASTLTFITDLG_H
