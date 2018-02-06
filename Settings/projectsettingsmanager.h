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

#ifndef PALSPROJECTSETTINGSMANAGER_H
#define PALSPROJECTSETTINGSMANAGER_H

#define MAX_PROJECT_PATH_CNT 15

#include "settings.h"

class PALSProjectSettingsManager
{
    DSimpleXMLNode *m_rootNode;
    DSimpleXMLNode *m_lastProjectNode;
    DSimpleXMLNode *m_linLogOnExitNode;
    DSimpleXMLNode *m_lastPathNode;
    DSimpleXMLNode *m_lastBackgroundChannelRangeNode;
    DSimpleXMLNode *m_resultWindowWasShown;
    DSimpleXMLNode *m_backgroundCalculationWithLastChannels;

    QStringList m_projectPathList;

public:
    static PALSProjectSettingsManager *sharedInstance();

    bool load();
    bool save();

    void addLastProjectPathToList(const QString& path);
    void setLinearAsLastScaling(bool lin);
    void setLastChosenPath(const QString& path);
    void setLastBackgroundChannelRange(int range);
    void setBackgroundCalculationFromFirstChannels(bool first);
    void setResultWindowWasShownOnExit(bool on);

    QStringList getLastProjectPathList() const;
    bool isLinearLastScaling() const;
    QString getLastChosenPath() const;
    int getLastBackgroundChannelRange() const;
    bool getBackgroundCalculationFromFirstChannels() const;
    bool getResultWindowWasShownOnExit() const;

private:
    PALSProjectSettingsManager();
    virtual ~PALSProjectSettingsManager();
};

#endif // PALSPROJECTSETTINGSMANAGER_H
