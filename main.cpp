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

#include "ltfitdlg.h"
#include "ltdefines.h"

#include <QApplication>
#include <QSplashScreen>
#include <QDesktopWidget>
#include <QSharedMemory>
#include <QStringBuilder>

int main(int argc, char *argv[])
{
    QSharedMemory sharedMemory;
    sharedMemory.setKey("DQuickLTFit0123456789qwetzuioasdfghjklerfgbnpokjn,.-234567890weuhcq8934cn43q8DQuickLTFit");

    if ( !sharedMemory.create(1) ) {
        DMSGBOX("An instance of DQuickLTFit is already running!");

        exit(0);
        return 0;
    }

    QApplication a(argc, argv);
    a.setApplicationName("DQuickLTFit");

    QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);

    QSplashScreen splash;
    splash.setPixmap(QPixmap::fromImage(QImage(":/localImages/Images/PALS.JPG").scaledToWidth(QApplication::desktop()->availableGeometry().width()*0.5)));
    splash.show();

#if !defined(Q_OS_WIN)
    splash.setFont(QFont("Arial", 12));
#endif

    splash.showMessage((QString(QString("<b>") % VERSION_STRING_AND_PROGRAM_NAME % QString(" (") % VERSION_RELEASE_DATE % QString(") </b><br>(C) Copyright 2016-2018 by Danny Petschke. All rights reserved."))), Qt::AlignLeft | Qt::AlignTop, Qt::darkGray);

    const QTime dieTime= QTime::currentTime().addSecs(3);
    while ( QTime::currentTime() < dieTime )
       QCoreApplication::processEvents();

    DFastLTFitDlg w;
    w.show();
    splash.finish(&w);

    return a.exec();
}
