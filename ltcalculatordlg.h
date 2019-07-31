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

#ifndef DFASTCALCULATORDLG_H
#define DFASTCALCULATORDLG_H

#include <QWidget>

#include "DLib/DLib.h"

namespace Ui {
class DFastCalculatorDlg;
}

class DFastCalculatorDlg : public QWidget
{
    Q_OBJECT
public:
    explicit DFastCalculatorDlg(QWidget *parent = 0);
    virtual ~DFastCalculatorDlg();

public slots:
    void setTextFont(const QFont& font);

protected:
    virtual void closeEvent(QCloseEvent *event);
    virtual void hideEvent(QHideEvent *event);
    virtual void showEvent(QShowEvent *event);

signals:
    void visibilityChanged(bool);

private:
    Ui::DFastCalculatorDlg *ui;
};

#endif // DFASTCALCULATORDLG_H
