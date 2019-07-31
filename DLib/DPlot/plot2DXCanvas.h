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

#ifndef PLOT2DXCANVAS_H
#define PLOT2DXCANVAS_H

#include <QWidget>

#include <QColor>
#include <QPainter>
#include <QPen>
#include <QList>
#include <QPoint>
#include <QPixmap>

#include "plot2DXCurve.h"

class plot2DXCanvas : public QWidget
{
    friend class plot2DXWidget;

    Q_OBJECT
public:
    plot2DXCanvas(QWidget *parent = 0);

public slots:
    void setBackgroundColor(const QColor& color);
    void clear();

signals:
    void canvasPropertyChanged();

public:
    QColor getBackgroundColor() const;
    void drawCurve(int curveWidth, const QColor &curveColor, plot2DXCurve::curveStyle style, const QList<QPoint> &pixelList);

protected:
    virtual void paintEvent(QPaintEvent *event);

private:
    QPixmap* pixmap();

private slots:
    void shiftPixmap(int shift);

    void drawPoints(const QPoint& pixel, QPainter* painter);
    void drawCross(const QPoint& pixel, QPainter* painter);
    void drawRect(const QPoint& pixel, QPainter* painter);
    void drawLine(const QList<QPoint>& pixelList, QPainter* painter);
    void drawCircle(const QPoint& pixel, QPainter* painter);

    void drawYLeftGrid(const QList<double>& yPxList, const QPen& pen);
    void drawYRightGrid(const QList<double>& yPxList, const QPen& pen);
    void drawXBottomGrid(const QList<double>& xPxList, const QPen &pen);
    void drawXTopGrid(const QList<double>& xPxList, const QPen &pen);

private:
    QPixmap m_canvasPixmap;
    QColor m_bgdColor;
};

#endif // PLOT2DXCANVAS_H
