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

#include "ltcalculatordlg.h"
#include "ui_ltcalculatordlg.h"

#ifndef WINDOWS_FONT
#define WINDOWS_FONT(__pointSize__)  QFont("Arial", __pointSize__)
#endif

DFastCalculatorDlg::DFastCalculatorDlg(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DFastCalculatorDlg)
{
    ui->setupUi(this);
}

DFastCalculatorDlg::~DFastCalculatorDlg()
{
    DDELETE_SAFETY(ui);
}

void DFastCalculatorDlg::setTextFont(const QFont &font)
{
    ui->textEdit->setTextFont(font);
}

void DFastCalculatorDlg::closeEvent(QCloseEvent *event)
{
    event->ignore();

    QWidget::closeEvent(event);
}

void DFastCalculatorDlg::hideEvent(QHideEvent *event)
{
    emit visibilityChanged(false);

    QWidget::hideEvent(event);
}

void DFastCalculatorDlg::showEvent(QShowEvent *event)
{
    emit visibilityChanged(true);

    QWidget::showEvent(event);
}
