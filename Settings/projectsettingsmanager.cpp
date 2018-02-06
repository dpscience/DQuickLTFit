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

#include "projectsettingsmanager.h"

static PALSProjectSettingsManager *pSettingsManagerInstance = nullptr;

PALSProjectSettingsManager *PALSProjectSettingsManager::sharedInstance()
{
    if ( !pSettingsManagerInstance )
        pSettingsManagerInstance = new PALSProjectSettingsManager();

    return pSettingsManagerInstance;
}

bool PALSProjectSettingsManager::load()
{
    DSimpleXMLReader reader("dquickltfitsetup.dqltfsettings");

    DSimpleXMLTag contentTag;
    if ( reader.readFromFile(&contentTag) )
    {
        bool ok = false;
        DSimpleXMLTag tag = contentTag.getTag(m_rootNode, &ok);

        if ( ok )
        {
            DSimpleXMLTag valueTag = tag.getTag("last-project-path", &ok);

            if ( ok ) m_lastProjectNode->setValue(valueTag.getValue());
            else m_lastProjectNode->setValue("");

            valueTag = tag.getTag("linear-as-last-scaling", &ok);

            if ( ok ) m_linLogOnExitNode->setValue(valueTag.getValue());
            else m_linLogOnExitNode->setValue(true);

            valueTag = tag.getTag("last-chosen-path");

            if ( ok ) m_lastPathNode->setValue(valueTag.getValue());
            else m_lastPathNode->setValue("/home");

            valueTag = tag.getTag("last-background-channel-range");

            if ( ok ) m_lastBackgroundChannelRangeNode->setValue(valueTag.getValue());
            else m_lastBackgroundChannelRangeNode->setValue(1000);

            valueTag = tag.getTag("background-calculation-using-first-channels");

            if ( ok ) m_backgroundCalculationWithLastChannels->setValue(valueTag.getValue());
            else m_backgroundCalculationWithLastChannels->setValue(false);

            valueTag = tag.getTag("result-window-was-shown");

            if ( ok ) m_resultWindowWasShown->setValue(valueTag.getValue());
            else m_resultWindowWasShown->setValue(true);

            const QStringList pathList = ((DString)m_lastProjectNode->getValue().toString()).parseBetween2("{", "}");

            m_projectPathList.clear();

            for ( QString item : pathList )
                m_projectPathList.append(item);

            return true;
        }
        else
        {
            m_lastProjectNode->setValue("");
            m_linLogOnExitNode->setValue(true);
            m_lastPathNode->setValue("/home");
            m_lastBackgroundChannelRangeNode->setValue(1000);
            m_resultWindowWasShown->setValue(true);
            m_projectPathList.clear();

            return false;
        }
    }
    else
    {
        m_lastProjectNode->setValue("");
        m_linLogOnExitNode->setValue(true);
        m_lastPathNode->setValue("/home");
        m_lastBackgroundChannelRangeNode->setValue(1000);
        m_resultWindowWasShown->setValue(true);
        m_projectPathList.clear();

        return false;
    }
}

bool PALSProjectSettingsManager::save()
{
    QString list = "";
    for ( QString item : m_projectPathList )
        list = list % "{" % item % "}";

    m_lastProjectNode->setValue(list);

    DSimpleXMLWriter writer("dquickltfitsetup.dqltfsettings");

    return writer.writeToFile(m_rootNode);
}

void PALSProjectSettingsManager::addLastProjectPathToList(const QString &path)
{
    if ( path.isEmpty() )
        return;


    if ( m_projectPathList.contains(path) )
    {
        QString existingPath = m_projectPathList.takeAt(m_projectPathList.indexOf(path));
        DUNUSED_PARAM(existingPath);
    }

    if ( m_projectPathList.size() == MAX_PROJECT_PATH_CNT )
    {
        QString lastPath = m_projectPathList.takeLast();
        DUNUSED_PARAM(lastPath);

        m_projectPathList.prepend(path);
    }
    else
        m_projectPathList.prepend(path);
}

void PALSProjectSettingsManager::setLinearAsLastScaling(bool lin)
{
    m_linLogOnExitNode->setValue(lin);
}

void PALSProjectSettingsManager::setLastChosenPath(const QString &path)
{
    m_lastPathNode->setValue(path);
}

void PALSProjectSettingsManager::setLastBackgroundChannelRange(int range)
{
    m_lastBackgroundChannelRangeNode->setValue(range);
}

void PALSProjectSettingsManager::setBackgroundCalculationFromFirstChannels(bool first)
{
    m_backgroundCalculationWithLastChannels->setValue(first);
}

void PALSProjectSettingsManager::setResultWindowWasShownOnExit(bool on)
{
    m_resultWindowWasShown->setValue(on);
}

QStringList PALSProjectSettingsManager::getLastProjectPathList() const
{
    return m_projectPathList;
}

bool PALSProjectSettingsManager::isLinearLastScaling() const
{
    return m_linLogOnExitNode->getValue().toBool();
}

QString PALSProjectSettingsManager::getLastChosenPath() const
{
    return m_lastPathNode->getValue().toString();
}

int PALSProjectSettingsManager::getLastBackgroundChannelRange() const
{
    return m_lastBackgroundChannelRangeNode->getValue().toInt();
}

bool PALSProjectSettingsManager::getBackgroundCalculationFromFirstChannels() const
{
    return m_backgroundCalculationWithLastChannels->getValue().toBool();
}

bool PALSProjectSettingsManager::getResultWindowWasShownOnExit() const
{
    return m_resultWindowWasShown->getValue().toBool();
}

PALSProjectSettingsManager::PALSProjectSettingsManager()
{
    m_rootNode = new DSimpleXMLNode("project-settings");
    m_lastProjectNode = new DSimpleXMLNode("last-project-path");
    m_linLogOnExitNode = new DSimpleXMLNode("linear-as-last-scaling");
    m_lastPathNode = new DSimpleXMLNode("last-chosen-path");
    m_lastBackgroundChannelRangeNode = new DSimpleXMLNode("last-background-channel-range");
    m_backgroundCalculationWithLastChannels = new DSimpleXMLNode("background-calculation-using-first-channels");
    m_resultWindowWasShown = new DSimpleXMLNode("result-window-was-shown");

    (*m_rootNode) << m_lastProjectNode << m_linLogOnExitNode << m_lastPathNode << m_lastBackgroundChannelRangeNode << m_backgroundCalculationWithLastChannels << m_resultWindowWasShown;
}

PALSProjectSettingsManager::~PALSProjectSettingsManager()
{
    save();

    DDELETE_SAFETY(m_resultWindowWasShown);
    DDELETE_SAFETY(m_lastProjectNode);
    DDELETE_SAFETY(m_linLogOnExitNode);
    DDELETE_SAFETY(m_lastPathNode);
    DDELETE_SAFETY(m_lastBackgroundChannelRangeNode);
    DDELETE_SAFETY(m_backgroundCalculationWithLastChannels);
    DDELETE_SAFETY(m_rootNode);
    DDELETE_SAFETY(pSettingsManagerInstance);
}
