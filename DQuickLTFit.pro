#-------------------------------------------------
#
# Project created by QtCreator 2015-12-06T15:24:47
#
#-------------------------------------------------

QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

macx {
QMAKE_MAC_SDK = macosx10.11
}

TARGET = DQuickLTFit_4_1
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
        ltcalculatordlg.cpp

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
                    ltdefines.h

FORMS    += \
                   ltfitdlg.ui \
                   ltresultdlg.ui \
                   ltplotdlg.ui \
                   ltparameterlistview.ui \
                   ltfitplotresidualview.ui \
                   ltcalculatordlg.ui


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
                            Images.qrc
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#DLib-import <END>:
