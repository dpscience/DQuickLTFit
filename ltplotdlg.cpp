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

#include "ltplotdlg.h"
#include "ui_ltplotdlg.h"

DFastPlotDlg::DFastPlotDlg(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DFastPlotDlg)
{
    ui->setupUi(this);

    ///raw-data:
    ui->widget->dataPlotView_1()->curve().at(0)->setCurveColor(Qt::red);
    ui->widget->dataPlotView_1()->curve().at(0)->setCurveStyle(plot2DXCurve::rect);
    ui->widget->dataPlotView_1()->curve().at(0)->setCurveWidth(2);

    ///preview-data:
    ui->widget->dataPlotView_1()->curve().at(1)->setCurveColor(Qt::blue);
    ui->widget->dataPlotView_1()->curve().at(1)->setCurveStyle(plot2DXCurve::line);
    ui->widget->dataPlotView_1()->curve().at(1)->setCurveWidth(2);

    ///fit-data:
    ui->widget->dataPlotView_1()->curve().at(2)->setCurveColor(Qt::green);
    ui->widget->dataPlotView_1()->curve().at(2)->setCurveStyle(plot2DXCurve::line);
    ui->widget->dataPlotView_1()->curve().at(2)->setCurveWidth(2);

    ///background-data:
    ui->widget->dataPlotView_1()->curve().at(5)->setCurveColor(Qt::black);
    ui->widget->dataPlotView_1()->curve().at(5)->setCurveStyle(plot2DXCurve::line);
    ui->widget->dataPlotView_1()->curve().at(5)->setCurveWidth(1);

    ///residual-data:
    ui->widget->dataPlotView_2()->curve().at(0)->setCurveColor(Qt::red);
    ui->widget->dataPlotView_2()->curve().at(0)->setCurveStyle(plot2DXCurve::line);
    ui->widget->dataPlotView_2()->curve().at(0)->setCurveWidth(2);

    ///fitrange-limits:
    ui->widget->dataPlotView_1()->curve().at(3)->setCurveColor(Qt::black);
    ui->widget->dataPlotView_1()->curve().at(3)->setCurveStyle(plot2DXCurve::line);
    ui->widget->dataPlotView_1()->curve().at(3)->setCurveWidth(1);

    ui->widget->dataPlotView_1()->curve().at(4)->setCurveColor(Qt::black);
    ui->widget->dataPlotView_1()->curve().at(4)->setCurveStyle(plot2DXCurve::line);
    ui->widget->dataPlotView_1()->curve().at(4)->setCurveWidth(1);

    setStyleSheet("background: white");

    connect(ui->widget->imageExportButton(), SIGNAL(clicked()), this, SLOT(savePlotAsImage()));
    connect(ui->widget->exportDataButton(), SIGNAL(clicked()), this, SLOT(exportResidualsFitAndRawData()));
    connect(ui->widget->dataPlotView_1()->yLeft(), SIGNAL(scalingPropertyChanged()), this, SLOT(updateROI()));
}

DFastPlotDlg::~DFastPlotDlg()
{
    DDELETE_SAFETY(ui);
}

bool DFastPlotDlg::isLinearScalingEnabled() const
{
    return ui->widget->isLinearScalingEnabled();
}

void DFastPlotDlg::closeEvent(QCloseEvent *event)
{
    event->ignore();

    QWidget::closeEvent(event);
}

void DFastPlotDlg::hideEvent(QHideEvent *event)
{
    emit visibilityChanged(false);

    QWidget::hideEvent(event);
}

void DFastPlotDlg::showEvent(QShowEvent *event)
{
    emit visibilityChanged(true);

    QWidget::showEvent(event);
}

void DFastPlotDlg::addRawData(const QList<QPointF> &datas)
{
    ui->widget->dataPlotView_1()->curve().at(0)->addData(datas);
    ui->widget->dataPlotView_1()->replot();

    ui->widget->dataPlotView_1()->autoscale();
}

void DFastPlotDlg::addPreviewData(const QList<QPointF> &datas)
{
    ui->widget->dataPlotView_1()->curve().at(1)->addData(datas);
    ui->widget->dataPlotView_1()->replot();
}

void DFastPlotDlg::addFitData(const QList<QPointF> &datas)
{
    ui->widget->dataPlotView_1()->curve().at(2)->addData(datas);
    ui->widget->dataPlotView_1()->replot();
}

void DFastPlotDlg::addResidualData(const QList<QPointF> &datas)
{
    ui->widget->dataPlotView_2()->curve().at(0)->addData(datas);
    ui->widget->dataPlotView_2()->replot();

    ui->widget->dataPlotView_2()->autoscale();
    ui->widget->dataPlotView_2()->yLeft()->setAxisRange(-4, 4);

    ui->widget->dataPlotView_2()->xBottom()->setAxisRange(ui->widget->dataPlotView_1()->xBottom()->getAxisMinValue(), ui->widget->dataPlotView_1()->xBottom()->getAxisMaxValue());
}

void DFastPlotDlg::setFitRange(int lower, int upper)
{
    ui->widget->dataPlotView_1()->curve().at(3)->clearCurveContent();
    ui->widget->dataPlotView_1()->curve().at(4)->clearCurveContent();

    double yAxisMin = 0;
    if ( ui->widget->dataPlotView_1()->yLeft()->getAxisMinValue() < 1 )
        yAxisMin = 1;
    else
        yAxisMin = ui->widget->dataPlotView_1()->yLeft()->getAxisMinValue();

    const QPointF pLower1(lower, yAxisMin);
    const QPointF pLower2(lower, ui->widget->dataPlotView_1()->yLeft()->getAxisMaxValue());

    const QPointF pUpper1(upper, yAxisMin);
    const QPointF pUpper2(upper, ui->widget->dataPlotView_1()->yLeft()->getAxisMaxValue());

    ui->widget->dataPlotView_1()->curve().at(3)->addData(pLower1.x(), pLower1.y());
    ui->widget->dataPlotView_1()->curve().at(3)->addData(pLower2.x(), pLower2.y());
    ui->widget->dataPlotView_1()->curve().at(3)->addData(pLower1.x(), pLower1.y());

    ui->widget->dataPlotView_1()->curve().at(4)->addData(pUpper1.x(), pUpper1.y());
    ui->widget->dataPlotView_1()->curve().at(4)->addData(pUpper2.x(), pUpper2.y());
    ui->widget->dataPlotView_1()->curve().at(4)->addData(pUpper1.x(), pUpper1.y());

    ui->widget->dataPlotView_1()->replot();
}

void DFastPlotDlg::updateBkgrdData()
{
    ui->widget->dataPlotView_1()->curve().at(5)->clearCurveContent();

    const QPointF pLeft(ui->widget->dataPlotView_1()->xBottom()->getAxisMinValue(), PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getBackgroundParamPtr()->getParameter()->getStartValue());
    const QPointF pRight(ui->widget->dataPlotView_1()->xBottom()->getAxisMaxValue(), PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getBackgroundParamPtr()->getParameter()->getStartValue());

    ui->widget->dataPlotView_1()->curve().at(5)->addData(pLeft.x(), pLeft.y());
    ui->widget->dataPlotView_1()->curve().at(5)->addData(pRight.x(), pRight.y());

    ui->widget->dataPlotView_1()->replot();
}

void DFastPlotDlg::setXRange(int min, int max)
{
    ui->widget->setXLimits(min, max);
}

void DFastPlotDlg::setYRangeData(int min, int max)
{
    ui->widget->setYLimits(min, max);
}

void DFastPlotDlg::setYRangeConvidenceLevel(double min, double max)
{
    //TODO!
}

void DFastPlotDlg::setRawDataVisible(bool visible)
{
    ui->widget->dataPlotView_1()->curve().at(0)->showCurve(visible);
}

void DFastPlotDlg::setStartValueDataVisible(bool visible)
{
    ui->widget->dataPlotView_1()->curve().at(1)->showCurve(visible);
}

void DFastPlotDlg::setFitDataVisible(bool visible)
{
    ui->widget->dataPlotView_1()->curve().at(2)->showCurve(visible);
}

void DFastPlotDlg::setLinearScaling()
{
    if ( !ui->widget->isLinearScalingEnabled() )
        ui->widget->changeYAxisScaling();
}

void DFastPlotDlg::setLogarithmicScaling()
{
    if ( ui->widget->isLinearScalingEnabled() )
        ui->widget->changeYAxisScaling();
}

void DFastPlotDlg::savePlotAsImage()
{
    showMaximized();

    const QString filename = QFileDialog::getSaveFileName(this, tr("Select or type a filename..."),
                                                          PALSProjectSettingsManager::sharedInstance()->getLastChosenPath(),
                                                          tr("PNG (*.png);;JPG (*.jpg);;JPEG (*.jpeg);; BMP (*.bmp);; PPM (*.ppm);; XBM (*.xbm);; XPM (*.xpm)"));

    if ( filename.isEmpty() )
        return;

    PALSProjectSettingsManager::sharedInstance()->setLastChosenPath(QFileInfo(filename).absoluteDir().absolutePath());

    ui->widget->setButtonsVisible(false);

    QPixmap map = ui->widget->grab();
    QPainter painter(&map);
    painter.setRenderHint(QPainter::Antialiasing);

    QTextOption o;
    o.setWrapMode(QTextOption::NoWrap);

    const QString str = PALSProjectManager::sharedInstance()->getFileName() % " [Saved: " % QDateTime::currentDateTime().toString() % "]";

    QFontMetrics metrics(painter.font());
    const QRect rect = metrics.boundingRect(str);

    painter.drawText(QRectF(20, 20, rect.width(), rect.height()), str, o);

    map.save(filename, 0, 100);
    ui->widget->setButtonsVisible(true);
}

void DFastPlotDlg::exportResidualsFitAndRawData()
{
    plot2DXWidget *rawAndFitData = ui->widget->dataPlotView_1();
    plot2DXWidget *residualsData = ui->widget->dataPlotView_2();

    if (!rawAndFitData || !residualsData) {
        DMSGBOX("Sorry, an unknown error occurred! Please save this project and restart this DQuickLTFit.");
        return;
    }

    plot2DXCurve *rawData = ui->widget->dataPlotView_1()->curve().at(0);

    if (rawData->getData().isEmpty()) {
        DMSGBOX("No data available!");
        return;
    }

    plot2DXCurve *fitData = ui->widget->dataPlotView_1()->curve().at(2);
    plot2DXCurve *residualData = ui->widget->dataPlotView_2()->curve().at(0);

    const QString filename = QFileDialog::getSaveFileName(this, tr("Select or type a filename..."),
                                                          PALSProjectSettingsManager::sharedInstance()->getLastChosenPath(),
                                                          tr("txt (*.txt)"));

    if ( filename.isEmpty() )
        return;

    PALSProjectSettingsManager::sharedInstance()->setLastChosenPath(QFileInfo(filename).absoluteDir().absolutePath());

    QString base_1 = filename.split(".txt").at(0);
    QString base_2 = filename.split(".txt").at(0);

    const QString fitTraceFileName = base_1.append("_fitData.txt");
    const QString residualsFileName = base_2.append("_residuals.txt");

    //save the rawData:
    QFile fileRawData(filename);

    if (fileRawData.open(QIODevice::ReadWrite)) {
        QTextStream stream(&fileRawData);

        stream << "channel [#]\tcounts[#]\r\n";

        for ( QPointF p : rawData->getData() ) {
            stream << QVariant((int)p.x()).toString() << "\t" << QVariant((int)p.y()).toString() << "\r\n";
        }

        fileRawData.close();
    }

    if (fitData->getData().isEmpty() || residualData->getData().isEmpty()) {
        DMSGBOX("Note: Residuals and Fitdata were not saved!");
        return;
    }

    //save the fitData:
    QFile fileFitData(fitTraceFileName);

    if (fileFitData.open(QIODevice::ReadWrite)) {
        QTextStream stream(&fileFitData);

        stream << "fraction of channel [#]\tfraction of counts[#]\r\n";

        for ( QPointF p : fitData->getData() ) {
            stream << QVariant(p.x()).toString() << "\t" << QVariant(p.y()).toString() << "\r\n";
        }

        fileFitData.close();
    }

    //save the residualData:
    QFile fileResidualData(residualsFileName);

    if (fileResidualData.open(QIODevice::ReadWrite)) {
        QTextStream stream(&fileResidualData);

        stream << "channel [#]\tresiduals [sigma]\r\n";

        for ( QPointF p : residualData->getData() ) {
            stream << QVariant((int)p.x()).toString() << "\t" << QVariant(p.y()).toString() << "\r\n";
        }

        fileResidualData.close();
    }
}

void DFastPlotDlg::updateROI()
{
    const int lower = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStartChannel();
    const int upper = PALSProjectManager::sharedInstance()->getDataStructure()->getFitSetPtr()->getStopChannel();

    ui->widget->dataPlotView_1()->curve().at(3)->clearCurveContent();
    ui->widget->dataPlotView_1()->curve().at(4)->clearCurveContent();

    double yAxisMin = 0;
    if ( ui->widget->dataPlotView_1()->yLeft()->getAxisMinValue() < 1 )
        yAxisMin = 1;
    else
        yAxisMin = ui->widget->dataPlotView_1()->yLeft()->getAxisMinValue();

    const QPointF pLower1(lower, yAxisMin);
    const QPointF pLower2(lower, ui->widget->dataPlotView_1()->yLeft()->getAxisMaxValue());

    const QPointF pUpper1(upper, yAxisMin);
    const QPointF pUpper2(upper, ui->widget->dataPlotView_1()->yLeft()->getAxisMaxValue());

    ui->widget->dataPlotView_1()->curve().at(3)->addData(pLower1.x(), pLower1.y());
    ui->widget->dataPlotView_1()->curve().at(3)->addData(pLower2.x(), pLower2.y());
    ui->widget->dataPlotView_1()->curve().at(3)->addData(pLower1.x(), pLower1.y());

    ui->widget->dataPlotView_1()->curve().at(4)->addData(pUpper1.x(), pUpper1.y());
    ui->widget->dataPlotView_1()->curve().at(4)->addData(pUpper2.x(), pUpper2.y());
    ui->widget->dataPlotView_1()->curve().at(4)->addData(pUpper1.x(), pUpper1.y());

    ui->widget->dataPlotView_1()->replot();
}

void DFastPlotDlg::clearAll()
{
    ui->widget->dataPlotView_1()->curve().at(0)->clearCurveContent();
    ui->widget->dataPlotView_1()->curve().at(1)->clearCurveContent();
    ui->widget->dataPlotView_1()->curve().at(2)->clearCurveContent();

    ui->widget->dataPlotView_2()->curve().at(0)->clearCurveContent();

    ui->widget->dataPlotView_1()->replot();
    ui->widget->dataPlotView_2()->replot();
}

void DFastPlotDlg::clearFitData()
{
    ui->widget->dataPlotView_1()->curve().at(2)->clearCurveContent();
}

void DFastPlotDlg::clearRawData()
{
    ui->widget->dataPlotView_1()->curve().at(0)->clearCurveContent();
}

void DFastPlotDlg::clearPreviewData()
{
    ui->widget->dataPlotView_1()->curve().at(1)->clearCurveContent();
}

void DFastPlotDlg::clearResidualData()
{
    ui->widget->dataPlotView_2()->curve().at(0)->clearCurveContent();
}
