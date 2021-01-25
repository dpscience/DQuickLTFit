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

#ifndef DFASTLICENSETEXTBOX_H
#define DFASTLICENSETEXTBOX_H

#include <QWidget>
#include <QIODevice>
#include <QFile>

#include "DLib/DTypes/defines.h"

namespace Ui {
class DFastLicenseTextBox;
}

class DFastLicenseTextBox : public QWidget {
    Q_OBJECT
public:
    explicit DFastLicenseTextBox(QWidget *parent = nullptr);
    virtual ~DFastLicenseTextBox();

    void addLicense(const QString& license, const QString& header);

private:
    Ui::DFastLicenseTextBox *ui;
};

#endif // DFASTLICENSETEXTBOX_H
