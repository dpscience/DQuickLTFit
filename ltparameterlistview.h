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

#ifndef PARAMETERLISTVIEW_H
#define PARAMETERLISTVIEW_H

#include <QWidget>
#include <QTableWidget>
#include <QCheckBox>
#include <QHBoxLayout>

#include "Settings/settings.h"
#include "Settings/projectmanager.h"
#include "Settings/projectsettingsmanager.h"

class PALSSourceTableWidgetItemCollector;
class PALSSampleTableWidgetItemCollector;
class PALSBackgroundTableWidgetItemCollector;
class PALSDeviceTableWidgetItemCollector;

class CheckBoxTableWidgetItem;
class LabelTableWidgetItem;

namespace Ui {
class ParameterListView;
}

class ParameterListView : public QWidget
{
    Q_OBJECT

    friend class DFastLTFitDlg;
public:
    ParameterListView(QWidget *parent = 0);
    virtual ~ParameterListView();

    void updateParamterList();

    void setEnabled(bool enable);

private:
    void initializeSourceTableWidget();
    void initializeSampleTableWidget();
    void initializeDeviceTableWidget();

    QCheckBox *fixedBackgroundCheckBox() const;

private slots:
    void updateChannelResolution(double value);
    void updateIterations(int iter);
    void updateBackground(double value);
    void updateChannelRange(double lower, double upper);
    void sendToInstantPreview(int row, int col);
    void updateBackgroundValue();
    void saveBackgroundChannelRanges(int channels);
    void setUsingFirstChannelsForBkgrdCalc();

public slots:
    void addSourceComponent();
    void addSampleComponent();
    void addDeviceResolutionComponent();

    void removeSourceComponent();
    void removeSampleComponent();
    void removeDeviceResolutionComponent();

    void updateSourceComponentNames();
    void updateSampleComponentNames();
    void updateDeviceComponentNames();

    void refreshBackgroundValue(double background);

    void setFitRangeLimits(int lower, int upper);
    void setFitRange(int lower, int upper);

    void setBackgroundChannelRange(int range);
    void setBackgroundCalculationUsingFirstChannels(bool first);

signals:
    void dataChanged();
    void calculateBackground(int);
    void fitRangeChanged(int, int);

private:
    Ui::ParameterListView *ui;

    PALSFitSet *m_fitSet;

    QList<PALSSourceTableWidgetItemCollector*> m_sourceWidgetCollection;
    QList<PALSSampleTableWidgetItemCollector*> m_sampleWidgetCollection;
    QList<PALSDeviceTableWidgetItemCollector*> m_deviceWidgetCollection;
};


class PALSSourceTableWidgetItemCollector : public QObject
{
    PALSFitParameter *m_param;

    CheckBoxTableWidgetItem *m_fixedItem;
    //CheckBoxTableWidgetItem *m_activeItem;
    QTableWidgetItem *m_nameItem;
    LabelTableWidgetItem *m_aliasItem;
    QTableWidgetItem *m_lowerLimitItem;
    QTableWidgetItem *m_upperLimitItem;
    QTableWidgetItem *m_startValueItem;

    Q_OBJECT
public:
    PALSSourceTableWidgetItemCollector(PALSFitParameter *fitParam, QTableWidget *tableWidget, int row);
    virtual ~PALSSourceTableWidgetItemCollector();

public slots:
    void updateValue(int row, int column);
    void checkBoxStateChanged(CheckBoxTableWidgetItem *widget);
    void setAlias(const QString& name);
    void setBackgroundColor(const QColor& color);
};


class PALSSampleTableWidgetItemCollector : public QObject
{
    PALSFitParameter *m_param;

    CheckBoxTableWidgetItem *m_fixedItem;
    //CheckBoxTableWidgetItem *m_activeItem;
    QTableWidgetItem *m_nameItem;
    LabelTableWidgetItem *m_aliasItem;
    QTableWidgetItem *m_lowerLimitItem;
    QTableWidgetItem *m_upperLimitItem;
    QTableWidgetItem *m_startValueItem;

    Q_OBJECT
public:
    PALSSampleTableWidgetItemCollector(PALSFitParameter *fitParam, QTableWidget *tableWidget, int row);
    virtual ~PALSSampleTableWidgetItemCollector();

public slots:
    void updateValue(int row, int column);
    void checkBoxStateChanged(CheckBoxTableWidgetItem *widget);
    void setAlias(const QString& name);
    void setBackgroundColor(const QColor& color);
};


class PALSDeviceTableWidgetItemCollector : public QObject
{
    PALSFitParameter *m_param;

    CheckBoxTableWidgetItem *m_fixedItem;
    QTableWidgetItem *m_nameItem;
    LabelTableWidgetItem *m_aliasItem;
    QTableWidgetItem *m_lowerLimitItem;
    QTableWidgetItem *m_upperLimitItem;
    QTableWidgetItem *m_startValueItem;

    Q_OBJECT
public:
    PALSDeviceTableWidgetItemCollector(PALSFitParameter *fitParam, QTableWidget *tableWidget, int row);
    virtual ~PALSDeviceTableWidgetItemCollector();

public slots:
    void updateValue(int row, int column);
    void checkBoxStateChanged(CheckBoxTableWidgetItem *widget);
    void setAlias(const QString& name);
    void setBackgroundColor(const QColor& color);
};


class PALSBackgroundTableWidgetItemCollector : public QObject
{
    PALSFitParameter *m_param;

    CheckBoxTableWidgetItem *m_activeItem;
    CheckBoxTableWidgetItem *m_fixedItem;
    QTableWidgetItem *m_lowerLimitItem;
    QTableWidgetItem *m_upperLimitItem;
    QTableWidgetItem *m_startValueItem;

    Q_OBJECT
public:
    PALSBackgroundTableWidgetItemCollector(PALSFitParameter *fitParam, QTableWidget *tableWidget, int row);
    virtual ~PALSBackgroundTableWidgetItemCollector();

public slots:
    void updateValue(int row, int column);
    void checkBoxStateChanged(CheckBoxTableWidgetItem *widget);

    void setBackground(double val);
};


class LabelTableWidgetItem : public QWidget
{
    QLabel *m_label;
    QHBoxLayout *m_pLayout;

    Q_OBJECT
public:
    LabelTableWidgetItem(QWidget *parent = 0);
    virtual ~LabelTableWidgetItem();

    DString text() const;

public slots:
    void setText(const DString& text);
    void setBackgroundColor(const QColor& color);
    void setColor(const QColor& color);
};

class CheckBoxTableWidgetItem : public QWidget
{
    QCheckBox   *m_pCheckBox;
    QHBoxLayout *m_pLayout;

    Q_OBJECT
public:
    CheckBoxTableWidgetItem(QWidget *parent = 0);
    virtual ~CheckBoxTableWidgetItem();

    bool isChecked() const;
    void setCheckedState(bool checked);

private slots:
    void checked(bool checked);

public slots:
    void setBackgroundColor(const QColor& color);

signals:
    void clicked(bool checked);
    void stateChanged(CheckBoxTableWidgetItem* widget);
};

#endif // PARAMETERLISTVIEW_H
