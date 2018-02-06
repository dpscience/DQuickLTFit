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

#ifndef FITPLOTRESIDUALWINDOW_H
#define FITPLOTRESIDUALWINDOW_H

#include "DLib/DTypes/types.h"
#include "DLib/DPlot/plot2DXWidget.h"
#include "DLib/DGUI/svgbutton.h"
#include "DLib/DGUI/horizontalrangedoubleslider.h"
#include "DLib/DGUI/verticalrangedoubleslider.h"

namespace Ui {
class DSynchronizedDblPlotWindow;
}

class DSynchronizedDblPlotWindow : public QWidget
{
    Q_OBJECT
public:
    explicit DSynchronizedDblPlotWindow(QWidget *parent = 0);
    virtual ~DSynchronizedDblPlotWindow();

    bool isLinearScalingEnabled() const;

    DSVGButton *imageExportButton() const;
    DSVGButton *exportDataButton() const;

public slots:
    plot2DXWidget* dataPlotView_1() const;  //upper
    plot2DXWidget* dataPlotView_2() const;  //lower

    void changeYAxisScaling(); //switch between log and lin scaling (upper)

    void setYLimits(double lower, double upper); //only upper plot!
    void setXLimits(double lower, double upper); //both plots

    void setButtonsVisible(bool visible);

    void autoscale();

private slots:
    void changeXRangeVisibility();
    void changeYRangeVisibility();

private:
    Ui::DSynchronizedDblPlotWindow *ui;
    DVerticalRangeDoubleSlider *m_dVRangeDblSlider;
    DHorizontalRangeDblSlider *m_dHRangeDblSlider;
};

#endif // FITPLOTRESIDUALWINDOW_H
