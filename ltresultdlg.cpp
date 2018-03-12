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

#include "ltresultdlg.h"
#include "ui_ltresultdlg.h"

#ifndef WINDOWS_FONT
#define WINDOWS_FONT(__pointSize__)  QFont("Arial", __pointSize__)
#endif

DFastResultDlg::DFastResultDlg(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DFastResultDlg)
{
    ui->setupUi(this);

    ui->tabWidget->setTabsClosable(true);
    ui->tabWidget->setUsesScrollButtons(true);

    connect(ui->tabWidget, SIGNAL(tabCloseRequested(int)), this, SLOT(closeTab(int)));
    connect(ui->pushButton_exportAsPDF, SIGNAL(clicked()), this, SLOT(printToPDF()));
    connect(ui->pushButton_exportAsHTML, SIGNAL(clicked()), this, SLOT(printToHTML()));
    connect(ui->pushButton_removeAllResults, SIGNAL(clicked()), this, SLOT(clearTabsFromButtonClick()));

    connect(this, SIGNAL(resultListIsEmpty()), this, SLOT(disablePDFExport()));
    connect(this, SIGNAL(resultListHasResults()), this, SLOT(enablePDFExport()));

    ui->pushButton_exportAsPDF->setLiteralSVG(":/localImages/Images/pdfExport");
    ui->pushButton_exportAsHTML->setLiteralSVG(":/localImages/Images/htmlExport");
    ui->pushButton_removeAllResults->setLiteralSVG(":/localImages/Images/remove");

    ui->pushButton_exportAsPDF->setToolTip("Export selected Results as PDF");
    ui->pushButton_exportAsHTML->setToolTip("Export selected Results as HTML");
    ui->pushButton_removeAllResults->setToolTip("Clear all Results");
}

DFastResultDlg::~DFastResultDlg()
{
    while ( ui->tabWidget->count() > 0 )
    {
        PALSProjectManager::sharedInstance()->getResultHistorie()->removeResult(0);
        ui->tabWidget->removeTab(0);

        ResultTab *resultTab = m_tabList.takeAt(0);
        DDELETE_SAFETY(resultTab);
    }

    DDELETE_SAFETY(ui);
}

void DFastResultDlg::closeEvent(QCloseEvent *event)
{
    event->ignore();

    QWidget::closeEvent(event);
}

void DFastResultDlg::hideEvent(QHideEvent *event)
{
    emit visibilityChanged(false);

    QWidget::hideEvent(event);
}

void DFastResultDlg::showEvent(QShowEvent *event)
{
    emit visibilityChanged(true);

    QWidget::showEvent(event);
}

void DFastResultDlg::addResultTabFromLastFit()
{
    ResultTab *tab = new ResultTab;
    tab->addResultsFromLastFit();

    m_tabList.append(tab);

    ui->tabWidget->addTab(tab, "");
    ui->tabWidget->setTabText(ui->tabWidget->count()-1, "Fit-Results " % QVariant(ui->tabWidget->count()).toString());
    ui->tabWidget->setTabToolTip(ui->tabWidget->count()-1, "Fit-Results " % QVariant(ui->tabWidget->count()).toString());

    ui->tabWidget->setCurrentWidget(tab);

    emit resultListHasResults();
}

void DFastResultDlg::addResultTabsFromHistory()
{
    for ( int i = 0 ; i < PALSProjectManager::sharedInstance()->getResultHistorie()->getSize() ; ++ i )
    {
        PALSResult *result = PALSProjectManager::sharedInstance()->getResultHistorie()->getResultAt(i);

        if ( !result )
            continue;

        ResultTab *tab = new ResultTab;
        tab->addResult(result);

        m_tabList.append(tab);

        ui->tabWidget->addTab(tab, "");
        ui->tabWidget->setTabText(ui->tabWidget->count()-1, "Fit-Results " % QVariant(ui->tabWidget->count()).toString());
        ui->tabWidget->setTabToolTip(ui->tabWidget->count()-1, "Fit-Results " % QVariant(ui->tabWidget->count()).toString());
    }

    ui->tabWidget->setCurrentIndex(ui->tabWidget->count()-1);

    if ( PALSProjectManager::sharedInstance()->getResultHistorie()->getSize() == 0 )
        emit resultListIsEmpty();
    else
        emit resultListHasResults();
}

void DFastResultDlg::clearTabs()
{
    while ( ui->tabWidget->count() > 0 )
    {
        ui->tabWidget->removeTab(0);

        ResultTab *resultTab = m_tabList.takeAt(0);
        DDELETE_SAFETY(resultTab);

        PALSProjectManager::sharedInstance()->getResultHistorie()->removeResult(0);
    }

    m_tabList.clear();

    emit resultListIsEmpty();
}

void DFastResultDlg::clearTabsFromButtonClick()
{
    if ( ui->tabWidget->count() ==  0 )
    {
        DMSGBOX("No results available!");
        return;
    }


    const QMessageBox::StandardButton replyBtn = QMessageBox::question(this, "Delete history?",
                                                                       "<nobr>Clearing the history cannot be undone. Are you sure?</nobr>",
                                                                       QMessageBox::Yes|QMessageBox::No);

    if ( replyBtn == QMessageBox::StandardButton::No )
        return;


    clearTabs();
}

#include <QDesktopWidget>

void DFastResultDlg::printToPDF()
{
    if ( PALSProjectManager::sharedInstance()->getResultHistorie()->getSize() == 0 )
    {
        DMSGBOX("<nobr>Sorry, no results available.</nobr>");
        return;
    }

    showMaximized();

    const QString filename = QFileDialog::getSaveFileName(this, tr("Select or type a filename..."),
                                                          PALSProjectSettingsManager::sharedInstance()->getLastChosenPath(),
                                                          tr("PDF (*.pdf)"));

    if ( filename.isEmpty() )
        return;

    PALSProjectSettingsManager::sharedInstance()->setLastChosenPath(QFileInfo(filename).absoluteDir().absolutePath());

    QPdfWriter writer(filename);

    writer.setCreator("Automatically generated by DQuickLTFit software using QPdfWriter-Plugin.");
    writer.setTitle("Fit-Results of DQuickLTFit software.");

    const double dpi = 300;

    writer.setResolution((int)dpi);
    writer.setPageOrientation(QPageLayout::Landscape);
    writer.setPageSize(QPagedPaintDevice::A4);

    const double pageSizeA4Width = writer.pageSizeMM().width();
    const double textEditSizeWidth = ((ResultTab*)ui->tabWidget->currentWidget())->textEdit()->widthMM();

    const double scalingWidth = textEditSizeWidth/pageSizeA4Width;

    QPainter painter(&writer);
    painter.scale(scalingWidth, scalingWidth);

    ((ResultTab*)ui->tabWidget->currentWidget())->textEdit()->document()->drawContents(&painter,  QRectF(QPointF(0, 0), ((ResultTab*)ui->tabWidget->currentWidget())->textEdit()->document()->size()));
    painter.end();
}

void DFastResultDlg::printToHTML()
{
    if ( PALSProjectManager::sharedInstance()->getResultHistorie()->getSize() == 0 )
    {
        DMSGBOX("<nobr>Sorry, no results available.</nobr>");
        return;
    }

    const QString filename = QFileDialog::getSaveFileName(this, tr("Select or type a filename..."),
                                                          PALSProjectSettingsManager::sharedInstance()->getLastChosenPath(),
                                                          tr("HTML (*.html)"));

    if ( filename.isEmpty() )
        return;

    PALSProjectSettingsManager::sharedInstance()->setLastChosenPath(QFileInfo(filename).absoluteDir().absolutePath());

    QFile file(filename);

    if (file.open(QIODevice::ReadWrite)) {
        QTextStream stream(&file);
        stream.setCodec("utf-8");

        stream << ((ResultTab*)ui->tabWidget->currentWidget())->textEdit()->document()->toHtml("utf-8");

        file.close();
    }
}

void DFastResultDlg::closeTab(int index)
{
    PALSProjectManager::sharedInstance()->getResultHistorie()->removeResult(index);
    ui->tabWidget->removeTab(index);

    ResultTab *resultTab = m_tabList.takeAt(index);
    DDELETE_SAFETY(resultTab);

    renameTabs();

    if ( PALSProjectManager::sharedInstance()->getResultHistorie()->getSize() == 0 )
        emit resultListIsEmpty();
    else
        emit resultListHasResults();
}

void DFastResultDlg::renameTabs()
{
    for ( int i = 0 ; i < ui->tabWidget->count() ; ++ i )
    {
        ui->tabWidget->setTabText(i, "Fit-Results " % QVariant(i+1).toString());
        ui->tabWidget->setTabToolTip(i, "Fit-Results " % QVariant(i+1).toString());
    }
}

void DFastResultDlg::enablePDFExport()
{
    ui->pushButton_exportAsPDF->setEnabled(true);
    ui->pushButton_exportAsHTML->setEnabled(true);
}

void DFastResultDlg::disablePDFExport()
{
    ui->pushButton_exportAsPDF->setEnabled(false);
    ui->pushButton_exportAsHTML->setEnabled(false);
}

ResultTab::ResultTab(QWidget *parent) :
    QWidget(parent)
{
    m_pLayout = new QHBoxLayout(this);
    m_textEdit = new QTextEdit;

    m_textEdit->setFrameStyle(QFrame::NoFrame);

#if defined(Q_OS_WIN)
    m_textEdit->setFont(WINDOWS_FONT(10));
#elif defined(Q_OS_MAC)
    m_textEdit->setFont(QFont("Arial", 12));
#endif

    m_pLayout->addWidget(m_textEdit);
    m_pLayout->setAlignment(Qt::AlignCenter);
    m_pLayout->setContentsMargins(0, 0, 0, 0);

    setLayout(m_pLayout);
}

ResultTab::~ResultTab()
{
    DDELETE_SAFETY(m_textEdit);
    DDELETE_SAFETY(m_pLayout);
}

QTextEdit *ResultTab::textEdit() const
{
    return m_textEdit;
}

void ResultTab::addResultsFromLastFit()
{
    PALSResult *result = PALSProjectManager::sharedInstance()->getResultHistorie()->getResultAt(PALSProjectManager::sharedInstance()->getResultHistorie()->getSize()-1);

    if ( !result )
        return;

    m_textEdit->append(result->getResultText());
}

void ResultTab::addResult(PALSResult *result)
{
    if ( !result )
        return;

    m_textEdit->append(result->getResultText());
}
