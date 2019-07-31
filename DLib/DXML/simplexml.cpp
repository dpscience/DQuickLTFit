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

#include "simplexml.h"

static int dnodeCnt = 0;

QString getTabJumpString(int cnt)
{
    QString tabString = "";
    for ( int i = 0 ; i < cnt ; ++ i )
        tabString += "    ";

    return tabString;
}

DSimpleXMLWriter::DSimpleXMLWriter() :
    QFile(),
    m_tabCnt(0)
{
    dnodeCnt = 0;
}

DSimpleXMLWriter::DSimpleXMLWriter(const QString &fileName) :
    QFile(fileName),
    m_tabCnt(0)
{
    dnodeCnt = 0;
}

DSimpleXMLWriter::DSimpleXMLWriter(const DString &fileName) :
    QFile((QString)fileName),
    m_tabCnt(0)
{
    dnodeCnt = 0;
}

DSimpleXMLWriter::~DSimpleXMLWriter() {}

void DSimpleXMLWriter::writeNode(const QList<DSimpleXMLNode*>& nodeList, QTextStream *stream)
{
    m_tabCnt ++;
    (*stream) << QString("\r\n");

    for ( int index = 0 ; index < nodeList.size() ; ++ index )
    {
        if ( !nodeList.at(index) )
            continue;


        if ( nodeList.at(index)->hasValue() )
        {
            (*stream) << DSIMPLEXML_TABJUMP(m_tabCnt) % DSIMPLEXML_STARTNODE(nodeList.at(index)->nodeName());
                (*stream) << nodeList.at(index)->getValue().toString();
            (*stream) << DSIMPLEXML_ENDNODE(nodeList.at(index)->nodeName());
        }
        else if ( nodeList.at(index)->hasChilds() )
        {
            (*stream) << DSIMPLEXML_TABJUMP(m_tabCnt) % DSIMPLEXML_STARTNODE(nodeList.at(index)->nodeName());
                writeNode(nodeList.at(index)->getChilds(), stream);
            (*stream) << DSIMPLEXML_TABJUMP(m_tabCnt) % DSIMPLEXML_ENDNODE(nodeList.at(index)->nodeName());
        }

#ifdef QT_DEBUG
            qDebug() << nodeList.at(index)->nodeName();
#endif
        (*stream) << QString("\r\n");
    }

    m_tabCnt --;
}

bool DSimpleXMLWriter::writeToFile(DSimpleXMLNode *rootNode)
{
    if ( !rootNode )
    {
        DERRORLOG("DSimpleXMLNode: !Root-Node is nullptr.");
        return false;
    }


    QTextStream xmlStream(this);

    if ( open(QIODevice::ReadWrite|QIODevice::Truncate) )
    {
        xmlStream << DSIMPLEXML_STARTNODE(rootNode->nodeName());

        if ( rootNode->hasValue() )
            xmlStream << rootNode->getValue().toString();
        else if ( rootNode->hasChilds() )
            writeNode(rootNode->getChilds(), &xmlStream);

        xmlStream << DSIMPLEXML_ENDNODE(rootNode->nodeName());

        close();
    }
    else
        return false;


    return true;
}

bool DSimpleXMLWriter::writeToFile(const QList<DSimpleXMLNode *> &rootNodeList)
{
    QTextStream xmlStream(this);

    if ( open(QIODevice::ReadWrite|QIODevice::Truncate) )
    {
        for ( int index = 0 ; index < rootNodeList.size() ; ++ index )
        {
            if ( !rootNodeList.at(index) )
                continue;

            xmlStream << DSIMPLEXML_STARTNODE(rootNodeList.at(index)->nodeName());

            if ( rootNodeList.at(index)->hasValue() )
                xmlStream << rootNodeList.at(index)->getValue().toString();
            else if ( rootNodeList.at(index)->hasChilds() )
                writeNode(rootNodeList.at(index)->getChilds(), &xmlStream);

            xmlStream << DSIMPLEXML_ENDNODE(rootNodeList.at(index)->nodeName());

            xmlStream << QString("\r\n");
        }

        close();
    }
    else
        return false;


    return true;
}


DSimpleXMLNode::DSimpleXMLNode() :
       m_nodeName(""),
       m_value(""),
       m_parent(nullptr)
{
    if ( m_nodeName.isEmpty() )
    {
        m_nodeName = "node_" % QVariant(dnodeCnt).toString();
        dnodeCnt ++;

        DERRORLOG("DSimpleXMLNode-Name (temporarely) set to '%s'.\n Use function 'setNodeName(...)'.", m_nodeName.toStdString().c_str());
    }
}

DSimpleXMLNode::DSimpleXMLNode(const QString &nodeName) :
    m_nodeName(nodeName),
    m_value(""),
    m_parent(nullptr)
{
    if ( m_nodeName.isEmpty() )
    {
        m_nodeName = "node_" % QVariant(dnodeCnt).toString();
        dnodeCnt ++;

        DERRORLOG("DSimpleXMLNode-Name (temporarely) set to '%s'.\n Use function 'setNodeName(...)'.", m_nodeName.toStdString().c_str());
    }
}

DSimpleXMLNode::DSimpleXMLNode(const DString &nodeName) :
    m_nodeName((QString)nodeName),
    m_value(""),
    m_parent(nullptr)
{
    if ( m_nodeName.isEmpty() )
    {
        m_nodeName = "node_" % QVariant(dnodeCnt).toString();
        dnodeCnt ++;

        DERRORLOG("DSimpleXMLNode-Name (temporarely) set to '%s'.\n Use function 'setNodeName(...)'.", m_nodeName.toStdString().c_str());
    }
}

DSimpleXMLNode& DSimpleXMLNode::operator<<(DSimpleXMLNode *childNode)
{
    if ( !childNode )
    {
        DERRORLOG("DSimpleXMLNode: !Child-Pointer is nullptr.");
        return *this;
    }
    else
    {
        addChild(childNode);
        childNode->setParent(this);
    }

    return *this;
}

DSimpleXMLNode::~DSimpleXMLNode()
{
    const int iNum = m_childs.count();
    for ( int i = 0 ; i < iNum ; ++ i )
    {
        DSimpleXMLNode* child = m_childs.takeFirst();
        DDELETE_SAFETY(child);
    }

    m_childs.clear();

    if ( m_parent )
        m_parent->removeFromParent(this);
}

QString DSimpleXMLNode::nodeName() const
{
    return m_nodeName;
}

bool DSimpleXMLNode::hasChilds() const
{
    if ( m_childs.isEmpty() )
        return false;
    else
        return true;
}

bool DSimpleXMLNode::hasParent() const
{
    return m_parent ? true : false;
}

bool DSimpleXMLNode::hasValue() const
{
    if ( !hasChilds() )
        return true;
    else
        return false;
}

void DSimpleXMLNode::addChild(DSimpleXMLNode *childNode)
{
    if ( childNode )
    {
        m_childs.append(childNode);
        childNode->setParent(this);
    }
    else
    {
        DERRORLOG("DSimpleXMLNode: !Child-Pointer is nullptr.");
    }
}

QList<DSimpleXMLNode*> DSimpleXMLNode::getChilds() const
{
     return m_childs;
}

QVariant DSimpleXMLNode::getValue() const
{
     return m_value;
}

void DSimpleXMLNode::setValue(const QVariant &value)
{
    m_value = value;
}

DSimpleXMLNode *DSimpleXMLNode::getParent() const
{
    return m_parent;
}

void DSimpleXMLNode::setParent(DSimpleXMLNode *parentNode)
{
    m_parent = parentNode;
}

void DSimpleXMLNode::removeFromParent(DSimpleXMLNode *childNode)
{
    if ( !childNode )
        return;


    const int index = m_childs.indexOf(childNode);
    m_childs.takeAt(index);
}

void DSimpleXMLNode::XMLMessageBox()
{
     DMSGBOX((QString)DSimpleXMLString(this));
}

void DSimpleXMLNode::setNodeName(const QString &nodeName)
{
     m_nodeName = nodeName;
}

bool DSimpleXMLNode::isValid() const
{
     if ( hasChilds() || hasValue() )
         return true;
     else
         return false;
}


DSimpleXMLTag::DSimpleXMLTag() :
    DString("") {}

DSimpleXMLTag::DSimpleXMLTag(const QString &value) :
    DString(value) {}

DSimpleXMLTag::DSimpleXMLTag(const DString &value) :
    DString(value) {}

DSimpleXMLTag::DSimpleXMLTag(const DSimpleXMLString &value) :
    DString(value) {}

DSimpleXMLTag::~DSimpleXMLTag() {}

DSimpleXMLTag DSimpleXMLTag::getTag(const QString &tagName, bool *ok) const
{
    const DString pValue = parseBetween(DSIMPLEXML_STARTNODE(tagName), DSIMPLEXML_ENDNODE(tagName));

    if ( pValue.isEmpty() )
    {
        if ( ok ) *ok = false;
        return DSimpleXMLTag("");
    }
    else
    {
        if ( ok ) *ok = true;
        return DSimpleXMLTag(pValue);
    }
}

DSimpleXMLTag DSimpleXMLTag::getTag(const DString &tagName, bool *ok) const
{
    return getTag((QString)tagName, ok);
}

DSimpleXMLTag DSimpleXMLTag::getTag(const DSimpleXMLString &tagName, bool *ok) const
{
    return getTag((QString)tagName, ok);
}

DSimpleXMLTag DSimpleXMLTag::getTag(const DSimpleXMLNode &node, bool *ok) const
{
    return getTag(node.nodeName(), ok);
}

DSimpleXMLTag DSimpleXMLTag::getTag(DSimpleXMLNode *node, bool *ok) const
{
    if ( !node )
        return DSimpleXMLTag("");
    else
        return getTag(node->nodeName(), ok);
}

QVariant DSimpleXMLTag::getValue() const
{
    return QVariant(*this);
}

QVariant DSimpleXMLTag::getValueAt(const QString &tagName, bool *ok) const
{
    return getTag(tagName, ok).getValue();
}

QVariant DSimpleXMLTag::getValueAt(const DString &tagName, bool *ok) const
{
    return getTag(tagName, ok).getValue();
}

QVariant DSimpleXMLTag::getValueAt(const DSimpleXMLNode &node, bool *ok) const
{
    return getTag(node, ok).getValue();
}

QVariant DSimpleXMLTag::getValueAt(DSimpleXMLNode *node, bool *ok) const
{
    if ( node )
        return getTag(node, ok).getValue();
    else
        return QVariant("");
}

void DSimpleXMLTag::XMLMessageBox()
{
    DMSGBOX((QString)*this);
}

QVariant DSimpleXMLTag::getValueAt(const DSimpleXMLString &tagName, bool *ok) const
{
    return getTag(tagName, ok).getValue();
}


DSimpleXMLReader::DSimpleXMLReader() :
    QFile("") {}

DSimpleXMLReader::DSimpleXMLReader(const QString &fileName) :
    QFile(fileName) {}

DSimpleXMLReader::DSimpleXMLReader(const DString &fileName) :
    QFile((QString)fileName) {}

DSimpleXMLReader::~DSimpleXMLReader() {}

bool DSimpleXMLReader::readFromFile(DSimpleXMLTag *xmlContent)
{
    if ( !xmlContent )
        return false;


    QTextStream xmlStream(this);

    if ( open(QIODevice::ReadOnly|QIODevice::Text) )
    {
        *xmlContent = DSimpleXMLTag(xmlStream.readAll());

        close();
    }
    else
        return false;


    return true;
}


DSimpleXMLString::DSimpleXMLString() :
    DString(""),
    m_tabCnt(0)
{
    dnodeCnt = 0;
}

DSimpleXMLString::DSimpleXMLString(DSimpleXMLNode *rootNode) :
    DString(""),
    m_tabCnt(0)
{
    dnodeCnt = 0;
    setXMLNode(rootNode);
}

DSimpleXMLString::DSimpleXMLString(const DSimpleXMLNode *rootNode) :
    DString(""),
    m_tabCnt(0)
{
    dnodeCnt = 0;
    setXMLNode(rootNode);
}

DSimpleXMLString::DSimpleXMLString(const QList<DSimpleXMLNode *> &rootNodeList) :
    DString(""),
    m_tabCnt(0)
{
    dnodeCnt = 0;
    setXMLNode(rootNodeList);
}

DSimpleXMLString::~DSimpleXMLString() {}

void DSimpleXMLString::setXMLNode(DSimpleXMLNode *rootNode)
{
    if ( !rootNode )
    {
        DERRORLOG("DSimpleXMLString: !Root-Node is nullptr.");
        return;
    }


    QTextStream xmlStream(this);

    xmlStream << DSIMPLEXML_STARTNODE(rootNode->nodeName());

    if ( rootNode->hasValue() )
        xmlStream << rootNode->getValue().toString();
    else if ( rootNode->hasChilds() )
        writeNode(rootNode->getChilds(), &xmlStream);

    xmlStream << DSIMPLEXML_ENDNODE(rootNode->nodeName());
}

void DSimpleXMLString::setXMLNode(const DSimpleXMLNode *rootNode)
{
    if ( !rootNode )
    {
        DERRORLOG("DSimpleXMLString: !Root-Node is nullptr.");
        return;
    }


    QTextStream xmlStream(this);

    xmlStream << DSIMPLEXML_STARTNODE(rootNode->nodeName());

    if ( rootNode->hasValue() )
        xmlStream << rootNode->getValue().toString();
    else if ( rootNode->hasChilds() )
        writeNode(rootNode->getChilds(), &xmlStream);

    xmlStream << DSIMPLEXML_ENDNODE(rootNode->nodeName());
}

void DSimpleXMLString::setXMLNode(const QList<DSimpleXMLNode *> &rootNodeList)
{
    QTextStream xmlStream(this);

    for ( int index = 0 ; index < rootNodeList.size() ; ++ index )
    {
        if ( !rootNodeList.at(index) )
            continue;

        xmlStream << DSIMPLEXML_STARTNODE(rootNodeList.at(index)->nodeName());

        if ( rootNodeList.at(index)->hasValue() )
            xmlStream << rootNodeList.at(index)->getValue().toString();
        else if ( rootNodeList.at(index)->hasChilds() )
            writeNode(rootNodeList.at(index)->getChilds(), &xmlStream);

        xmlStream << DSIMPLEXML_ENDNODE(rootNodeList.at(index)->nodeName());

        xmlStream << QString("\r\n");
    }
}

void DSimpleXMLString::writeNode(const QList<DSimpleXMLNode *> &nodeList, QTextStream *stream)
{
    m_tabCnt ++;
    (*stream) << QString("\r\n");

    for ( int index = 0 ; index < nodeList.size() ; ++ index )
    {
        if ( !nodeList.at(index) )
            continue;


        if ( nodeList.at(index)->hasValue() )
        {
            (*stream) << DSIMPLEXML_TABJUMP(m_tabCnt) % DSIMPLEXML_STARTNODE(nodeList.at(index)->nodeName());
                (*stream) << nodeList.at(index)->getValue().toString();
            (*stream) << DSIMPLEXML_ENDNODE(nodeList.at(index)->nodeName());
        }
        else if ( nodeList.at(index)->hasChilds() )
        {
            (*stream) << DSIMPLEXML_TABJUMP(m_tabCnt) % DSIMPLEXML_STARTNODE(nodeList.at(index)->nodeName());
                writeNode(nodeList.at(index)->getChilds(), stream);
            (*stream) << DSIMPLEXML_TABJUMP(m_tabCnt) % DSIMPLEXML_ENDNODE(nodeList.at(index)->nodeName());
        }

        (*stream) << QString("\r\n");
    }

    m_tabCnt --;
}
