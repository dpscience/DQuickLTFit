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

#include "ltfitplotresidualview.h"
#include "ui_ltfitplotresidualview.h"

#define WINDOWS_FONT(__pointSize__)  QFont("Arial", __pointSize__)

DSynchronizedDblPlotWindow::DSynchronizedDblPlotWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DSynchronizedDblPlotWindow)
{
    ui->setupUi(this);

    //buttons:
    ui->yAxisRangeButton->setLiteralSVG(":/localImages/Images/arrowUp");
    ui->xAxisRangeButton->setLiteralSVG(":/localImages/Images/arrowRight");
    ui->linLogButton->setLiteralSVG(":/localImages/Images/scaling");

    ui->saveAsPNGButton->setLiteralSVG(":/localImages/Images/pngExport");
    ui->saveFitAndResidualData->setLiteralSVG(":/localImages/Images/save");

    ui->yAxisRangeButton->setToolTip("Change the Vertical Axis Scaling");
    ui->xAxisRangeButton->setToolTip("Change the Horizontal Axis Scaling");
    ui->linLogButton->setToolTip("Switch between linear/logarithmic scaling");
    ui->saveAsPNGButton->setToolTip("Export Plot Window as PNG");
    ui->saveFitAndResidualData->setToolTip("Export Residuals, Fit- and Raw-Data");

    //plot-views:
    ui->plotWidget_1->yRight()->setVisible(false);
    ui->plotWidget_1->xTop()->setVisible(false);

    ui->plotWidget_2->yRight()->setVisible(false);
    ui->plotWidget_2->xTop()->setVisible(false);

    ui->plotWidget_1->showXBottomGrid(false);
    ui->plotWidget_1->showXTopGrid(false);
    ui->plotWidget_1->showYLeftGrid(false);
    ui->plotWidget_1->showYRightGrid(false);

    ui->plotWidget_2->showXBottomGrid(false);
    ui->plotWidget_2->showXTopGrid(false);
    ui->plotWidget_2->showYLeftGrid(true);
    ui->plotWidget_2->showYRightGrid(true);

    ui->plotWidget_1->yLeft()->setAxisDistribution(2);
    ui->plotWidget_2->yLeft()->setAxisDistribution(4);

    ui->plotWidget_1->yLeft()->setAxisLabelText("[#]");
    ui->plotWidget_2->xBottom()->setAxisLabelText("Channel [#]");
    ui->plotWidget_2->yLeft()->setAxisLabelText("Sigma");
    ui->plotWidget_1->xBottom()->setAxisLabelText("");

    ui->plotWidget_1->xBottom()->setAxisLabelPosition(plot2DXAxis::middle);
    ui->plotWidget_1->yLeft()->setAxisLabelPosition(plot2DXAxis::middle);

    ui->plotWidget_2->xBottom()->setAxisLabelPosition(plot2DXAxis::middle);
    ui->plotWidget_2->yLeft()->setAxisLabelPosition(plot2DXAxis::valueStart);

    ui->plotWidget_1->yLeft()->setAxisRange(1, 10000);
    ui->plotWidget_1->yLeft()->setNumberPrecision(0);
    ui->plotWidget_2->yLeft()->setAxisRange(-4, 4);
    ui->plotWidget_2->yLeft()->setNumberPrecision(0);

    ui->plotWidget_1->xBottom()->setAxisRange(0, 1024);
    ui->plotWidget_1->xBottom()->setNumberPrecision(0);
    ui->plotWidget_2->xBottom()->setAxisRange(0, 1024);
    ui->plotWidget_2->xBottom()->setNumberPrecision(0);

#if defined(Q_OS_WIN)
    ui->plotWidget_1->yLeft()->setFont(WINDOWS_FONT(8));
    ui->plotWidget_2->yLeft()->setFont(WINDOWS_FONT(8));
    ui->plotWidget_1->xBottom()->setFont(WINDOWS_FONT(8));
    ui->plotWidget_2->xBottom()->setFont(WINDOWS_FONT(8));
#endif

    m_dVRangeDblSlider = new DVerticalRangeDoubleSlider;
    m_dHRangeDblSlider = new DHorizontalRangeDblSlider;

#if defined(Q_OS_WIN)
    m_dVRangeDblSlider->setWindowFlags(Qt::Tool|Qt::WindowTitleHint|Qt::WindowCloseButtonHint|Qt::CustomizeWindowHint);
    m_dHRangeDblSlider->setWindowFlags(Qt::Tool|Qt::WindowTitleHint|Qt::WindowCloseButtonHint|Qt::CustomizeWindowHint);
#endif

#if defined(Q_OS_MAC) || defined(Q_OS_OSX)
    m_dVRangeDblSlider->setWindowFlags(Qt::Window | Qt::WindowTitleHint | Qt::CustomizeWindowHint | Qt::WindowCloseButtonHint | Qt::WindowMinimizeButtonHint);
    m_dHRangeDblSlider->setWindowFlags(Qt::Window | Qt::WindowTitleHint | Qt::CustomizeWindowHint | Qt::WindowCloseButtonHint | Qt::WindowMinimizeButtonHint);
#endif

    m_dHRangeDblSlider->setMaximumHeight(70);
    m_dHRangeDblSlider->setMinimumHeight(70);
    m_dHRangeDblSlider->setMinimumWidth(500);
    m_dHRangeDblSlider->setMaximumWidth(500);

    m_dVRangeDblSlider->setMaximumHeight(420);
    m_dVRangeDblSlider->setMinimumHeight(420);
    m_dVRangeDblSlider->setMinimumWidth(100);
    m_dVRangeDblSlider->setMaximumWidth(100);

#if defined(Q_OS_WIN)
    m_dVRangeDblSlider->setWindowTitle("Counts");
    m_dHRangeDblSlider->setWindowTitle("Channels");
#endif

#if defined(Q_OS_MAC) || defined(Q_OS_OSX)
    m_dVRangeDblSlider->setWindowTitle("");
    m_dHRangeDblSlider->setWindowTitle("Channels");
#endif

    connect(ui->yAxisRangeButton, SIGNAL(clicked()), this, SLOT(changeYRangeVisibility()));
    connect(ui->xAxisRangeButton, SIGNAL(clicked()), this, SLOT(changeXRangeVisibility()));
    connect(ui->linLogButton, SIGNAL(clicked()), this, SLOT(changeYAxisScaling()));

    connect(m_dHRangeDblSlider, SIGNAL(rangeChanged(double,double)), ui->plotWidget_1->xBottom(), SLOT(setAxisRange(double,double)));
    connect(m_dHRangeDblSlider, SIGNAL(rangeChanged(double,double)), ui->plotWidget_2->xBottom(), SLOT(setAxisRange(double,double)));
    connect(m_dVRangeDblSlider, SIGNAL(rangeChanged(double,double)), ui->plotWidget_1->yLeft(), SLOT(setAxisRange(double,double)));

    m_dHRangeDblSlider->setLimits(1, 1024);
    m_dVRangeDblSlider->setLimits(1, 10000);
}

DSynchronizedDblPlotWindow::~DSynchronizedDblPlotWindow()
{
    DDELETE_SAFETY(m_dVRangeDblSlider);
    DDELETE_SAFETY(m_dHRangeDblSlider);
    DDELETE_SAFETY(ui);
}

bool DSynchronizedDblPlotWindow::isLinearScalingEnabled() const
{
    if ( ui->plotWidget_1->yLeft()->getAxisScaling() == plot2DXAxis::linear )
        return true;
    else
        return false;
}

DSVGButton *DSynchronizedDblPlotWindow::imageExportButton() const
{
    return ui->saveAsPNGButton;
}

DSVGButton *DSynchronizedDblPlotWindow::exportDataButton() const
{
    return ui->saveFitAndResidualData;
}

plot2DXWidget *DSynchronizedDblPlotWindow::dataPlotView_1() const
{
    return ui->plotWidget_1;
}

plot2DXWidget *DSynchronizedDblPlotWindow::dataPlotView_2() const
{
    return ui->plotWidget_2;
}

void DSynchronizedDblPlotWindow::changeYAxisScaling()
{
    if ( ui->plotWidget_1->yLeft()->getAxisScaling() == plot2DXAxis::linear ){
        m_dVRangeDblSlider->setLimits(1.0, m_dVRangeDblSlider->upperLimit());
        ui->plotWidget_1->yLeft()->setAxisScaling(plot2DXAxis::logarithmic);
    }else{
        m_dVRangeDblSlider->setLimits(0.0, m_dVRangeDblSlider->upperLimit());
        ui->plotWidget_1->yLeft()->setAxisScaling(plot2DXAxis::linear);
    }
}

void DSynchronizedDblPlotWindow::setYLimits(double lower, double upper)
{
    if ( ui->plotWidget_1->yLeft()->getAxisScaling() == plot2DXAxis::linear )
    {
        m_dVRangeDblSlider->setLimits(lower, upper);
        m_dVRangeDblSlider->setLowerLevel(lower);
        m_dVRangeDblSlider->setUpperLevel(upper);
    }
    else if ( ui->plotWidget_1->yLeft()->getAxisScaling() == plot2DXAxis::logarithmic )
    {
        if ( lower <= 1.0f )
            lower = 1.0f;

        m_dVRangeDblSlider->setLimits(lower, upper);
        m_dVRangeDblSlider->setLowerLevel(lower);
        m_dVRangeDblSlider->setUpperLevel(upper);
    }
}

void DSynchronizedDblPlotWindow::setXLimits(double lower, double upper)
{
    m_dHRangeDblSlider->setLimits(lower, upper);
    m_dHRangeDblSlider->setUpperLevel(upper);
    m_dHRangeDblSlider->setLowerLevel(lower);
}

void DSynchronizedDblPlotWindow::changeYRangeVisibility()
{
    if ( !m_dVRangeDblSlider->isVisible() )
        m_dVRangeDblSlider->show();
    else
        m_dVRangeDblSlider->hide();
}

void DSynchronizedDblPlotWindow::setButtonsVisible(bool visible)
{
    ui->linLogButton->setVisible(visible);
    ui->yAxisRangeButton->setVisible(visible);
    ui->xAxisRangeButton->setVisible(visible);
    ui->saveAsPNGButton->setVisible(visible);
}

void DSynchronizedDblPlotWindow::autoscale()
{
    //TODO!!!
}

void DSynchronizedDblPlotWindow::changeXRangeVisibility()
{
    if ( !m_dHRangeDblSlider->isVisible() )
        m_dHRangeDblSlider->show();
    else
        m_dHRangeDblSlider->hide();
}
