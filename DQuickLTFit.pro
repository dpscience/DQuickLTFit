# ****************************************************************************
#
#  DQuickLTFit, a software for the analysis of Positron-Lifetime Spectra
#  based on the Least-Square Optimization using the Levenberg-Marquardt
#  Algorithm.
#
#  Copyright (C) 2016-2021 Danny Petschke
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see http://www.gnu.org/licenses/.
#
# *****************************************************************************
#
#  @author: Danny Petschke
#  @contact: danny.petschke@uni-wuerzburg.de
#
# *****************************************************************************

QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

macx {
QMAKE_MAC_SDK = macosx10.11
}

TARGET = DQuickLTFit_4_2
TEMPLATE = app

win32{
RC_FILE = myapp.rc
}

macx{
ICON = myappicon.icns
}

SOURCES += main.cpp\
        Settings/projectmanager.cpp \
        Settings/projectsettingsmanager.cpp \
        Settings/settings.cpp \
        Fit/mpfit.c\
        Fit/lifetimedecayfit.cpp \
        ltfitdlg.cpp \
        ltresultdlg.cpp \
        ltplotdlg.cpp \
        ltparameterlistview.cpp \
        ltfitplotresidualview.cpp \
        ltcalculatordlg.cpp \
        ltlicensetextbox.cpp \

HEADERS  += \
                    Settings/projectmanager.h \
                    Settings/projectsettingsmanager.h \
                    Settings/settings.h \
                    Fit/mpfit.h \
                    Fit/mpfit_DISCLAIMER \
                    Fit/lifetimedecayfit.h \
                    ltfitdlg.h \
                    ltresultdlg.h \
                    ltplotdlg.h \
                    ltparameterlistview.h \
                    ltfitplotresidualview.h \
                    ltcalculatordlg.h \
                    ltlicensetextbox.h \
                    ltdefines.h

FORMS    += \
                   ltfitdlg.ui \
                   ltresultdlg.ui \
                   ltplotdlg.ui \
                   ltparameterlistview.ui \
                   ltfitplotresidualview.ui \
                   ltcalculatordlg.ui \
                   ltlicensetextbox.ui \


#DLib-import <START>:
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
QT       += xml svg

macx {
QMAKE_MAC_SDK = macosx10.11
}

CONFIG   += c++11

HEADERS  += DLib/DLib.h\
            DLib/DTypes/defines.h\
            DLib/DTypes/types.h\
            DLib/DXML/simplexml.h\
            DLib/DGUI/svgbutton.h\
            DLib/DGUI/slider.h\
            DLib/DGUI/verticalrangedoubleslider.h\
            DLib/DGUI/horizontalrangedoubleslider.h\
            DLib/DGUI/constantexplanations.h \
            DLib/DGUI/mathconsoletextbox.h \
            DLib/DPlot/plot2DXWidget.h\
            DLib/DPlot/plot2DXCurve.h\
            DLib/DPlot/plot2DXAxis.h\
            DLib/DPlot/plot2DXCanvas.h\
            DLib/DCompression/compressionwrapper.h

SOURCES  += DLib/DTypes/defines.cpp\
            DLib/DTypes/types.cpp\
            DLib/DXML/simplexml.cpp\
            DLib/DGUI/svgbutton.cpp\
            DLib/DGUI/slider.cpp\
            DLib/DGUI/verticalrangedoubleslider.cpp\
            DLib/DGUI/horizontalrangedoubleslider.cpp\
            DLib/DGUI/constantexplanations.cpp \
            DLib/DGUI/mathconsoletextbox.cpp \
            DLib/DPlot/plot2DXWidget.cpp\
            DLib/DPlot/plot2DXCurve.cpp\
            DLib/DPlot/plot2DXAxis.cpp\
            DLib/DPlot/plot2DXCanvas.cpp\
            DLib/DCompression/miniz.c\
            DLib/DCompression/compressionwrapper.cpp

FORMS    += DLib/DGUI/horizontalrangedoubleslider.ui\
            DLib/DGUI/verticalrangedoubleslider.ui

RESOURCES += \
            DLib/DResources/res.qrc \
                            Images.qrc \
    license.qrc
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#DLib-import <END>:
