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

#include "svgbutton.h"

DSVGButton::DSVGButton(QWidget *parent) :
    QWidget(parent),
    m_enabled(true),
    m_state(undefined),
    m_bgColor("transparent")
{
    setStyleSheet("QWidget{background-color: " + m_bgColor + "}");
}

DSVGButton::DSVGButton(const QString &pathLiteral, QWidget *parent) :
    QWidget(parent),
    m_defSVGPath(pathLiteral % "_default.svg"),
    m_hoverSVGPath(pathLiteral % "_hover.svg"),
    m_clickedSVGPath(pathLiteral % "_click.svg"),
    m_enabled(true),
    m_state(undefined),
    m_bgColor("transparent")
{
    setStyleSheet("QWidget{background-color: " + m_bgColor + "}");
}

DSVGButton::DSVGButton(const QString &defaultStateSVGPath, const QString &hoverStateSVGPath, const QString &clickStateSVGPath, QWidget *parent) :
    QWidget(parent),
    m_defSVGPath(defaultStateSVGPath),
    m_hoverSVGPath(hoverStateSVGPath),
    m_clickedSVGPath(clickStateSVGPath),
    m_enabled(true),
    m_state(undefined),
    m_bgColor("transparent")
{
    setStyleSheet("QWidget{background-color: " % m_bgColor % "}");
}

DSVGButton::~DSVGButton() {}

QString DSVGButton::customStatusTip() const
{
    return m_statusTip;
}

void DSVGButton::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    painter.setRenderHints(QPainter::Antialiasing|QPainter::HighQualityAntialiasing);

    QImage img;

    if ( m_state == hover || m_state == release || !m_enabled )
        img = DSVGImage::getImage(m_hoverSVGPath, geometry().height()-2*DSVGButton_SVG_OFFSET, geometry().height()-2*DSVGButton_SVG_OFFSET);
    else if ( m_state == click && m_enabled )
        img = DSVGImage::getImage(m_clickedSVGPath, geometry().height()-2*DSVGButton_SVG_OFFSET, geometry().height()-2*DSVGButton_SVG_OFFSET);
    else if ( (m_state == leave || m_state == undefined) && m_enabled )
        img = DSVGImage::getImage(m_defSVGPath, geometry().height()-2*DSVGButton_SVG_OFFSET, geometry().height()-2*DSVGButton_SVG_OFFSET);

    if ( !img.isNull() )
        painter.drawImage(DSVGButton_SVG_OFFSET, DSVGButton_SVG_OFFSET, img);

    QWidget::paintEvent(event);
}

void DSVGButton::resizeEvent(QResizeEvent *event){ QWidget::resizeEvent(event); }

bool DSVGButton::event(QEvent *event)
{
    if ( event->type() == QEvent::Enter )
    {
        m_state = hover;
        update();

        if ( !m_statusTip.isEmpty() )
            emit statusChanged(m_statusTip);
    }
    else if ( event->type() == QEvent::Leave )
    {
        m_state = leave;
        update();

        emit statusChanged("");
    }
    else if ( event->type() == QEvent::MouseButtonPress )
    {
        m_state = click;
        update();

        if ( !m_statusTip.isEmpty() )
            emit statusChanged(m_statusTip);
    }
    else if ( event->type() == QEvent::MouseButtonRelease )
    {
        m_state = release;
        update();

        if ( !m_statusTip.isEmpty() )
            emit statusChanged(m_statusTip);

        emit clicked();
    }

    return QWidget::event(event);
}

void DSVGButton::mousePressEvent(QMouseEvent *event)
{
    event->ignore();
}

void DSVGButton::mouseReleaseEvent(QMouseEvent *event)
{
    event->ignore();
}

void DSVGButton::mouseMoveEvent(QMouseEvent *event)
{
    event->ignore();
}

void DSVGButton::setLiteralSVG(const QString &pathLiteral)
{
    m_defSVGPath     = pathLiteral % "_default.svg";
    m_hoverSVGPath   = pathLiteral % "_hover.svg";
    m_clickedSVGPath = pathLiteral % "_click.svg";
}

void DSVGButton::setDefaultStateSVG(const QString &path)
{
    m_defSVGPath = path;
    update();
}

void DSVGButton::setHoverStateSVG(const QString &path)
{
    m_hoverSVGPath = path;
    update();
}

void DSVGButton::setClickedStateSVG(const QString &path)
{
    m_clickedSVGPath = path;
    update();
}

void DSVGButton::setBackgroundColor(const QString &cssName)
{
    m_bgColor = cssName;
    setStyleSheet("QWidget{background-color: " + m_bgColor + "}");
}

void DSVGButton::setBackgroundColor(const QColor &color)
{
    m_bgColor = QString("rgb(" % QVariant(color.red()).toString() % ", " + QVariant(color.green()).toString() % ", " % QVariant(color.blue()).toString() % ")");
    setStyleSheet("QWidget{background-color: " % m_bgColor % "}");
}

void DSVGButton::setCustomStatusTip(const QString& statusTip)
{
    m_statusTip = statusTip;
    update();
}

void DSVGButton::enableWidget(bool enabled)
{
    m_enabled = enabled;
    QWidget::setEnabled(enabled);
    update();
}


DSVGToolButton::DSVGToolButton(QWidget *parent) :
    QToolButton(parent),
    m_enabled(true),
    m_state(undefined)
{
    setStyleSheet("QWidget{background-color: transparent}");
}

DSVGToolButton::DSVGToolButton(const QString &pathLiteral, QWidget *parent) :
    QToolButton(parent),
    m_defSVGPath(pathLiteral % "_default.svg"),
    m_hoverSVGPath(pathLiteral % "_hover.svg"),
    m_clickedSVGPath(pathLiteral % "_click.svg"),
    m_enabled(true),
    m_state(undefined)
{
    setStyleSheet("QWidget{background-color: transparent}");
}

DSVGToolButton::DSVGToolButton(const QString &defaultStateSVGPath, const QString &hoverStateSVGPath, const QString &clickStateSVGPath, QWidget *parent) :
    QToolButton(parent),
    m_defSVGPath(defaultStateSVGPath),
    m_hoverSVGPath(hoverStateSVGPath),
    m_clickedSVGPath(clickStateSVGPath),
    m_enabled(true),
    m_state(undefined)
{
    setStyleSheet("QWidget{background-color: transparent}");
}

DSVGToolButton::~DSVGToolButton() {}

QString DSVGToolButton::customStatusTip() const
{
    return m_statusTip;
}

void DSVGToolButton::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    painter.setRenderHints(QPainter::Antialiasing|QPainter::HighQualityAntialiasing);

    QImage img;

    if ( m_state == hover || m_state == release || !m_enabled )
        img = DSVGImage::getImage(m_hoverSVGPath, geometry().height()-2*DSVGButton_SVG_OFFSET, geometry().height()-2*DSVGButton_SVG_OFFSET);
    else if ( m_state == click && m_enabled )
        img = DSVGImage::getImage(m_clickedSVGPath, geometry().height()-2*DSVGButton_SVG_OFFSET, geometry().height()-2*DSVGButton_SVG_OFFSET);
    else if ( (m_state == leave || m_state == undefined) && m_enabled )
        img = DSVGImage::getImage(m_defSVGPath, geometry().height()-2*DSVGButton_SVG_OFFSET, geometry().height()-2*DSVGButton_SVG_OFFSET);

    if ( !img.isNull() )
        painter.drawImage(DSVGButton_SVG_OFFSET, DSVGButton_SVG_OFFSET, img);

    QToolButton::paintEvent(event);
}

bool DSVGToolButton::event(QEvent *event)
{
    if ( event->type() == QEvent::Enter )
    {
        m_state = hover;
        update();

        if ( !m_statusTip.isEmpty() )
            emit statusChanged(m_statusTip);
    }
    else if ( event->type() == QEvent::Leave )
    {
        m_state = leave;
        update();

        emit statusChanged("");
    }
    else if ( event->type() == QEvent::MouseButtonPress )
    {
        m_state = click;
        update();

        if ( !m_statusTip.isEmpty() )
            emit statusChanged(m_statusTip);
    }
    else if ( event->type() == QEvent::MouseButtonRelease )
    {
        m_state = release;
        update();

        if ( !m_statusTip.isEmpty() )
            emit statusChanged(m_statusTip);

        emit clicked();
    }

    return QToolButton::event(event);
}

void DSVGToolButton::mousePressEvent(QMouseEvent *event)
{
    event->ignore();
}

void DSVGToolButton::mouseReleaseEvent(QMouseEvent *event)
{
    event->ignore();
}

void DSVGToolButton::mouseMoveEvent(QMouseEvent *event)
{
    event->ignore();
}

void DSVGToolButton::setDefaultStateSVG(const QString &path)
{
    m_defSVGPath = path;
    update();
}

void DSVGToolButton::setHoverStateSVG(const QString &path)
{
    m_hoverSVGPath = path;
    update();
}

void DSVGToolButton::setClickedStateSVG(const QString &path)
{
    m_clickedSVGPath = path;
    update();
}

void DSVGToolButton::setCustomStatusTip(const QString& statusTip)
{
    m_statusTip = statusTip;
    update();
}

void DSVGToolButton::enableWidget(bool enabled)
{
    m_enabled = enabled;
    QToolButton::setEnabled(enabled);
    update();
}
