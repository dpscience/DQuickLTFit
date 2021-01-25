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

#ifndef DFASTPLOTDLG_H
#define DFASTPLOTDLG_H

#include <QWidget>
#include <QLabel>
#include <QFileDialog>
#include <QKeyEvent>
#include <QTextItem>
#include <QPainter>

#include "Settings/projectmanager.h"
#include "Settings/projectsettingsmanager.h"

namespace Ui {
class DFastPlotDlg;
}

class DFastPlotDlg : public QWidget
{
    Q_OBJECT
public:
    explicit DFastPlotDlg(QWidget *parent = 0);
    virtual ~DFastPlotDlg();

    bool isLinearScalingEnabled() const;

protected:
    virtual void closeEvent(QCloseEvent *event);
    virtual void hideEvent(QHideEvent *event);
    virtual void showEvent(QShowEvent *event);

public slots:
    void addRawData(const QList<QPointF>& datas);
    void addPreviewData(const QList<QPointF>& datas);
    void addFitData(const QList<QPointF>& datas);
    void addResidualData(const QList<QPointF>& datas);

    void clearAll();

    void clearFitData();
    void clearRawData();
    void clearPreviewData();
    void clearResidualData();

    void setFitRange(int lower, int upper);
    void updateBkgrdData();

    void setXRange(int min, int max);
    void setYRangeData(int min, int max);
    void setYRangeConvidenceLevel(double min, double max);

    void setRawDataVisible(bool visible);
    void setStartValueDataVisible(bool visible);
    void setFitDataVisible(bool visible);

    void setLinearScaling();
    void setLogarithmicScaling();

    void savePlotAsImage();
    void exportResidualsFitAndRawData();

private slots:
    void updateROI();

signals:
    void visibilityChanged(bool);

private:
    Ui::DFastPlotDlg *ui;
};

#endif // DFASTPLOTDLG_H
