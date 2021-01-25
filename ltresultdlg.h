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

#ifndef DFASTRESULTDLG_H
#define DFASTRESULTDLG_H

#include <QWidget>
#include <QTextEdit>
#include <QLayout>
#include <QFileDialog>
#include <QFileInfo>
#include <QPdfWriter>

#include "Settings/settings.h"
#include "Settings/projectmanager.h"
#include "Settings/projectsettingsmanager.h"

class ResultTab;

namespace Ui {
class DFastResultDlg;
}

class DFastResultDlg : public QWidget
{
    Q_OBJECT
public:
    explicit DFastResultDlg(QWidget *parent = nullptr);
    virtual ~DFastResultDlg();

protected:
    virtual void closeEvent(QCloseEvent *event);
    virtual void hideEvent(QHideEvent *event);
    virtual void showEvent(QShowEvent *event);

public slots:
    /* call this to add a new PALSFitResult automatically from the last results (directly after fitting) */
    void addResultTabFromLastFit();
    /* call this to add all results from the history (e.g. on loading a project) */
    void addResultTabsFromHistory();
    void clearTabs(bool fromButtonClick = false);
    void clearTabsFromButtonClick();

    void printToPDF();
    void printToHTML();

private slots:
    void closeTab(int);
    void renameTabs();

    void enablePDFExport();
    void disablePDFExport();

signals:
    void visibilityChanged(bool);
    void resultListIsEmpty();
    void resultListHasResults();

private:
    Ui::DFastResultDlg *ui;

    QList<ResultTab *> m_tabList;
};

class ResultTab : public QWidget
{
    Q_OBJECT
public:
    explicit ResultTab(QWidget *parent = nullptr);
    virtual ~ResultTab();

    QTextEdit *textEdit() const;

public slots:
    void addResultsFromLastFit();
    void addResult(PALSResult *result);

private:
    QHBoxLayout *m_pLayout;
    QTextEdit *m_textEdit;
};

#endif // DFASTRESULTDLG_H
