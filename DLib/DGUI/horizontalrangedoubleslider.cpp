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

#include "horizontalrangedoubleslider.h"
#include "ui_horizontalrangedoubleslider.h"

DHorizontalRangeDblSlider::DHorizontalRangeDblSlider(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DHorizontalRangeDblSlider),
    m_type(INT)
{
    ui->setupUi(this);

    if ( m_type == INT )
        m_numberValidator = new QIntValidator;
    else
        m_numberValidator = new QDoubleValidator;

    ui->lineEdit->setValidator(m_numberValidator);
    ui->lineEdit_2->setValidator(m_numberValidator);

    connect(ui->widget, SIGNAL(lowerValueChanged(double)), this, SLOT(updateLowerLevel(double)));
    connect(ui->widget, SIGNAL(upperValueChanged(double)), this, SLOT(updateUpperLevel(double)));

    connect(ui->lineEdit_2, SIGNAL(returnPressed()), this, SLOT(updateUpperLevel()));
    connect(ui->lineEdit, SIGNAL(returnPressed()), this, SLOT(updateLowerLevel()));
    connect(ui->lineEdit_2, SIGNAL(editingFinished()), this, SLOT(updateUpperLevel()));
    connect(ui->lineEdit, SIGNAL(editingFinished()), this, SLOT(updateLowerLevel()));

    connect(ui->widget, SIGNAL(valuesChanged(double,double)), this, SLOT(rangeHasChanged(double,double)));

    setLimits(DHorizontalRangeDoubleSlider_WIDGET_MIN, DHorizontalRangeDoubleSlider_WIDGET_MAX);

    updateUpperLevel(DHorizontalRangeDoubleSlider_WIDGET_MAX);
    updateLowerLevel(DHorizontalRangeDoubleSlider_WIDGET_MIN);
}

DHorizontalRangeDblSlider::~DHorizontalRangeDblSlider()
{
    DDELETE_SAFETY(m_numberValidator);
    DDELETE_SAFETY(ui);
}

double DHorizontalRangeDblSlider::upperLimit() const
{
    if ( m_type == INT )
        return (int)m_upperLimit;
    else
        return m_upperLimit;
}

double DHorizontalRangeDblSlider::lowerLimit() const
{
    if ( m_type == INT )
        return (int)m_lowerLimit;
    else
        return m_lowerLimit;
}

double DHorizontalRangeDblSlider::upperValue() const
{
    if ( m_type == INT )
        return (int)ui->widget->valueUpperLimit();
    else
        return ui->widget->valueUpperLimit();
}

double DHorizontalRangeDblSlider::lowerValue() const
{
    if ( m_type == INT )
        return (int)ui->widget->valueLowerLimit();
    else
        return ui->widget->valueLowerLimit();
}

DHorizontalRangeDblSlider::VALUE_TYPE DHorizontalRangeDblSlider::valueType() const
{
    return m_type;
}

void DHorizontalRangeDblSlider::setLimits(double min, double max)
{
    if ( max <= min )
        return;

    if ( m_type == INT )
    {
        min = (int)min;
        max = (int)max;
    }

    if ( m_type == INT )
    {
        ui->label->setText("(" % QVariant((int)min).toString() % ")"); //min
        ui->label_2->setText("(" % QVariant((int)max).toString() % ")"); //max
    }
    else
    {
        ui->label->setText("(" % QVariant(min).toString() % ")"); //min
        ui->label_2->setText("(" % QVariant(max).toString() % ")"); //max
    }

    ui->widget->setRange(min, max);

    ui->widget->setValueLowerSlider(min);
    ui->widget->setValueUpperSlider(max);

    update();
}

void DHorizontalRangeDblSlider::updateLowerLevel(double value)
{
    if ( m_type == INT )
        ui->lineEdit->setText(QVariant((int)value).toString());
    else
        ui->lineEdit->setText(QString::number(value, 'f', 2));

    update();
}

void DHorizontalRangeDblSlider::updateUpperLevel(double value)
{
    if ( m_type == INT )
        ui->lineEdit_2->setText(QVariant((int)value).toString());
    else
        ui->lineEdit_2->setText(QString::number(value, 'f', 2));

    update();
}

void DHorizontalRangeDblSlider::updateLowerLevel()
{
    bool ok = false;

    QString value = ui->lineEdit->text();
    value.replace(",", ".");

    const double dVal = QVariant(value).toDouble(&ok);

    if ( ok )
    {
        if ( m_type == DOUBLE )
        {
            if ( !ui->widget->setValueLowerSlider(dVal) )
                ui->lineEdit->setText(QString::number(ui->widget->valueLowerLimit(), 'f', 2));
            else
                ui->lineEdit->setText(QString::number(dVal, 'f', 2));
        }
        else
        {
            if ( !ui->widget->setValueLowerSlider((int)dVal) )
                ui->lineEdit->setText(QVariant((int)ui->widget->valueLowerLimit()).toString());
            else
                ui->lineEdit->setText(QVariant((int)dVal).toString());
        }
    }
    else
    {
        if ( m_type == DOUBLE )
            ui->lineEdit->setText(QString::number(ui->widget->valueLowerLimit(), 'f', 2));
        else
            ui->lineEdit->setText(QVariant(ui->widget->valueLowerLimit()).toString());
    }

    update();
}

void DHorizontalRangeDblSlider::updateUpperLevel()
{
    bool ok = false;

    QString value = ui->lineEdit_2->text();
    value.replace(",", ".");

    const double dVal = QVariant(value).toDouble(&ok);

    if ( ok )
    {
        if ( m_type == DOUBLE )
        {
            if ( !ui->widget->setValueUpperSlider(dVal) )
                ui->lineEdit_2->setText(QString::number(ui->widget->valueUpperLimit(), 'f', 2));
            else
                ui->lineEdit_2->setText(QString::number(dVal, 'f', 2));
        }
        else
        {
            if ( !ui->widget->setValueUpperSlider((int)dVal) )
                ui->lineEdit_2->setText(QVariant((int)ui->widget->valueUpperLimit()).toString());
            else
                ui->lineEdit_2->setText(QVariant((int)dVal).toString());
        }
    }
    else
    {
        if ( m_type == DOUBLE )
            ui->lineEdit_2->setText(QString::number(ui->widget->valueUpperLimit(), 'f', 2));
        else
            ui->lineEdit_2->setText(QVariant(ui->widget->valueUpperLimit()).toString());
    }

    update();
}

void DHorizontalRangeDblSlider::setLowerLevel(double value)
{
    if ( m_type == INT )
        value = (int)value;

    emit ui->widget->lowerValueChanged(value);

    updateLowerLevel(value);
    updateLowerLevel();

    update();
}

void DHorizontalRangeDblSlider::setUpperLevel(double value)
{
    if ( m_type == INT )
        value = (int)value;

    emit ui->widget->upperValueChanged(value);

    updateUpperLevel(value);
    updateUpperLevel();

    update();
}

void DHorizontalRangeDblSlider::setValueType(DHorizontalRangeDblSlider::VALUE_TYPE type)
{
    m_type = type;

    DDELETE_SAFETY(m_numberValidator);

    if ( m_type == INT )
        m_numberValidator = new QIntValidator;
    else
        m_numberValidator = new QDoubleValidator;
}

void DHorizontalRangeDblSlider::rangeHasChanged(double lowerLimit, double upperLimit)
{
    if ( m_type == INT )
    {
        lowerLimit = (int)lowerLimit;
        upperLimit = (int)upperLimit;
    }

    emit rangeChanged(lowerLimit, upperLimit);
}
