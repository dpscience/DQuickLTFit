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

#ifndef PLOT2DXCURVE_H
#define PLOT2DXCURVE_H

#include <QObject>

#include <QList>
#include <QPair>
#include <QColor>
#include <QPointF>
#include <QPoint>
#include <QRectF>

class plot2DXCurve : QObject
{
    friend class plot2DXWidget;

    Q_OBJECT
public:
    typedef enum{
        yLeft_xBottom = 0,
        yLeft_xTop = 1,
        yRight_xBottom = 2,
        yRight_xTop = 3
    }scaleAxis;

    typedef enum{
        line = 0,
        point = 1,
        cross = 2,
        rect = 3,
        circle = 4
    }curveStyle;

    plot2DXCurve();
    virtual ~plot2DXCurve();

public slots:
    void setAxis(scaleAxis axis);
    void setCurveWidth(int width);
    void setCurveColor(const QColor& color);
    void setCurveStyle(curveStyle style);
    void showCurve(bool show);

    void addData(double x_value, double y_value);
    void addData(const QList<QPointF> &dataset);

    void clearCurveContent();
    void clearCurveContent(int from, int to);
    void clearCurveCache();
    void clearCurveCache(int from, int to);
    void setMaxContainerSize(int size);

signals:
    void curvePropertyChanged();
    void maxValueChanged(double minx, double miny, double maxx, double maxy, plot2DXCurve::scaleAxis axis);

public:
    scaleAxis getAxis() const;

    int getCurveWidth() const;
    QColor getCurveColor() const;
    curveStyle getCurveStyle() const;

    bool isCurveShown() const;

    QList<QPointF> getData() const;
    QList<QPointF> getCache() const;

    int setMaxContainerSize() const;

private slots:
    void reset();
    void setLastValueBeforeReplot(const QPointF& lastValue);
    void swapToContainer();

private:
    QPointF getLastValueBeforeReplot() const;

private:
    QList<QPointF> m_cache;
    QList<QPointF> m_dataContainer;

    scaleAxis m_scale;

    int m_lineWidth;
    QColor m_lineColor;

    curveStyle m_curveStyle;

    bool m_shown;

    QPointF m_lastValueBeforeReplot;

    int m_maxCount;

    double x_min, x_max;
    double y_min, y_max;
};

#endif // PLOT2DXCURVE_H
