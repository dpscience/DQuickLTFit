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

#include "ltlicensetextbox.h"
#include "ui_ltlicensetextbox.h"

DFastLicenseTextBox::DFastLicenseTextBox(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DFastLicenseTextBox) {
    ui->setupUi(this);

    connect(ui->pushButton, SIGNAL(clicked()), this, SLOT(close()));    
}

DFastLicenseTextBox::~DFastLicenseTextBox() {
    DDELETE_SAFETY(ui)
}

void DFastLicenseTextBox::addLicense(const QString &license, const QString &header) {
    QFile file(license);

    if (file.open(QIODevice::ReadOnly)) {
        ui->plainTextEdit->appendPlainText(QString(file.readAll()));

        file.close();
    }

    setWindowTitle(header);

    ui->plainTextEdit->moveCursor(QTextCursor::Start);
    ui->plainTextEdit->setReadOnly(true);
}
