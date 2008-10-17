#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "MOOT config reporting base classes"
__author__   = "J. Panetta <panetta@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/01/25 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"

import logging
import os
import getpass
import time

from Ft.Xml.Domlette                 import NonvalidatingReader, implementation, PrettyPrint
from Ft.Xml                          import EMPTY_NAMESPACE, InputSource
from Ft.Lib                          import Uri
from Ft.Xml.Xvif                     import RelaxNgValidator
from Ft.Xml.Xslt         import Transform
from Ft.Lib.Uri          import OsPathToUri

from ConfigReport import *

_log = logging.getLogger('offline.configData')

CONFIG_NAMESPACE = EMPTY_NAMESPACE
CONFIG_DEST_PRODUCTION = "/nfs/farm/g/glast/u03/"

CONFIG_XSL_TRANSFORM = os.path.expandvars("$CONFIGDATAROOT/python/configTransform.xsl")

TAG_CFGRPT     = "ConfigReport"
TAG_CINFO      = "ConfigInfo"
TAG_BINFO      = "BaselineInfo"
TAG_PRECINCTS  = "Precincts"
TAG_PRECINCT   = "PrecinctReportLink"
TAG_COMMENT    = "Comment"
TAG_INTENT     = "Intent"
TAG_LINKTO     = "LinkTo"
TAG_PRCRPT     = "PrecinctReport"
TAG_PRCINFO    = "PrecinctInfo"
TAG_ANCINFO    = "AncillaryInfo"
TAG_PARINFO    = "ParamInfo"
TAG_SECTION    = "ReportSection"
TAG_IMG        = "Image"
TAG_TABLE      = "Table"
TAG_INCTEXT    = "TextInclude"
TAG_TEXTLINE   = "TextLine"

ATTR_ACTIVE    = "Active"
ATTR_ALIAS     = "Alias"
ATTR_BKEY      = "BaselineKey"
ATTR_CKEY      = "ConfigKey"
ATTR_CNAME     = "ConfigName"
ATTR_DATE      = "Date"
ATTR_DESC      = "Description"
ATTR_FILELINK  = "FileLink"
ATTR_FNAME     = "FileName"
ATTR_KEY       = "Key"
ATTR_NAME      = "Name"
ATTR_NLINES    = "NLines"
ATTR_NOCONV    = "NoConvert"
ATTR_PKEY      = "Key"
ATTR_PNAME     = "PrecinctName"
ATTR_STATUS    = "Status"
ATTR_TEXT      = "Text"
ATTR_TITLE     = "Title"
ATTR_USER      = "User"
ATTR_VKEY      = "VoteKey"

PRECINCT_HANDLERS = {}

class XmlReport(object):
    ROOT_TAG_TYPE = ""
    def __init__(self):
        self.__doc = implementation.createDocument(CONFIG_NAMESPACE, None, None)
        self.__root = None

    @property
    def rootNode(self):
        return self.__root

    def createHeader(self):
        self.__root = makeChildNode(self.__doc, self.ROOT_TAG_TYPE)
        setAttribute(self.rootNode, ATTR_CKEY, self.data.config)
        setAttribute(self.rootNode, ATTR_BKEY, self.data.baseline)
        setAttribute(self.rootNode, ATTR_USER, getpass.getuser())
        setAttribute(self.rootNode, ATTR_DATE, time.asctime())

    def writeReport(self, force=False):
        """!@brief write this report to a file.

        @param force   Force the creation of a file
        """
        if force or not self.isBuilt:
            f = open(self.fileName, 'w')
            PrettyPrint(self.__doc, f)
            f.close()
        return self.fileName
        
    def addConfigInfo(self):
        cInfoNode = makeChildNode(self.rootNode, TAG_CINFO)
        self.__addCInfoTags(cInfoNode, self.data.configInfo)

    def addBaselineInfo(self):
        bInfoNode = makeChildNode(self.rootNode, TAG_BINFO)
        self.__addCInfoTags(bInfoNode, self.data.baselineInfo)

    def __addCInfoTags(self, node, info):
        if not info: return
        setAttribute(node, ATTR_CNAME, info.getName())
        setAttribute(node, ATTR_KEY, info.getKey())
        setAttribute(node, ATTR_DATE, info.getCreationTime())
        setAttribute(node, ATTR_STATUS, info.getStatus())
        setAttribute(node, ATTR_ACTIVE, info.getActive())
        setAttribute(node, ATTR_DESC, info.getDescrip())
        if hasattr(info, 'getVoteKey'):
            setAttribute(node, ATTR_VKEY, int(info.getVoteKey()))
            vInfo = self.data.db.getVoteInfo(int(info.getVoteKey()))
            setAttribute(node, ATTR_FNAME, vInfo.getSrc())
        
    def addSection(self, sectionName, parent=None):
        if not parent:
            parent = self.rootNode
        sec = makeChildNode(parent, TAG_SECTION)
        setAttribute(sec, ATTR_NAME, sectionName)
        return sec

    def addComment(self, parent, commentText):
        comment = makeChildNode(parent, TAG_COMMENT)
        setAttribute(comment, ATTR_USER, getpass.getuser())
        setAttribute(comment, ATTR_DATE, time.asctime())
        makeTextChildNode(comment, commentText)

    def addLink(self, parent, linkTarget, linkText, noConvert=False):
        "!@brief add a link to another file.  Links will be made relative to the report directory"
        # noConvert:  Don't convert 'xml' to 'html' in links
        relTarget = relpath(linkTarget, self.path)
        link = makeChildNode(parent, TAG_LINKTO)
        setAttribute(link, ATTR_FNAME, str(relTarget))
        setAttribute(link, ATTR_NOCONV, str(noConvert))
        setAttribute(link, ATTR_DESC, str(linkText))

    def addImage(self, parent, fileName, title, caption=""):
        imgNode = makeChildNode(parent, TAG_IMG)
        setAttribute(imgNode, ATTR_FNAME, relpath(fileName, self.path))
        setAttribute(imgNode, ATTR_TITLE, title)
        if caption:
            makeTextChildNode(imgNode, caption)

    def addTable(self, parentSection,  border='1', width='', caption=''):
        """!@brief  Add a table under a section

        Border and width are optional.

        \param parentSection  section node that this image will be inserted in.
        \param border    Specifies whether the table will have a border or not. Default value is 1.
        \param width     Specifies the width of the border in percentage relative to the width of the page.

        \return          created node
        """
        tabNode = makeChildNode(parentSection, TAG_TABLE)
        if border != '':
            setAttribute(tabNode, 'border', border)
        if width != '':
            setAttribute(tabNode, 'width', width)
        if caption:
            cap = makeChildNode(tabNode, 'CAPTION')
            setAttribute(cap, 'align', 'bottom')
            makeTextChildNode(cap, str(caption))
        return tabNode

    def addTableHeader(self, tableNode, text, align=''):
        """!@brief Adds one column header to the table.
    
        \param tableNode tableNode where column headers are being added.
        \param text    Column header text.
        \param align   Optional text alignment.
        """
        firstRow = getChildNode(tableNode, 'TR')
        if not firstRow:
            firstRow = makeChildNode(tableNode, 'TR')
        header = makeChildNode(firstRow, 'TH')
        if align:
            setAttribute(header, 'align', align)
        makeTextChildNode(header, text)

    def addTableRow(self, tableNode):
        """!@brief start a table row
        \param tableNode   tableNode to add a row to

        \return the created row
        """
        return makeChildNode(tableNode, 'TR')

    def addTableData(self, row, text, align='', width='', bold=False, underline=False):
        """!@brief Adds a data column to a table row

        \param row   the table row
        \param text  Text of the data column.
        \param align Optional text alignment.
        \param width Optional column width in percentage.
        """
        element = makeChildNode(row, 'TD')
        if align:
            setAttribute(element, 'align', align)
        if width:
            setAttribute(element, 'width', width)
        if bold:
            element = makeChildNode(element, 'B')
        if underline:
            element = makeChildNode(element, 'U')
        makeTextChildNode(element, text)
        

    @property
    def isBuilt(self):
        """!@brief   Is this report already built?  Return t/f
        """
        # Simply implement as check if self.fileName exists
        try:
            return os.path.exists(self.fileName)
        except OSError, e:
            _log.error("Saw exception in isBuilt: %s"%e)
            return False


class ConfigXmlReport(XmlReport, ConfigReport):
    def __init__(self, configData):
        ConfigReport.__init__(self, configData)
        XmlReport.__init__(self)
        self.ROOT_TAG_TYPE = TAG_CFGRPT
        self.__precRpts = []

    def createHeader(self):
        super(ConfigXmlReport, self).createHeader()
        setAttribute(self.rootNode, ATTR_CNAME, self.data.configInfo.getName())

    def addPrecincts(self, rebuild):
        self.__getAliases()
        precsNode = makeChildNode(self.rootNode, TAG_PRECINCTS)
        # The heritage report was requested to be first
        heritage = PRECINCT_HANDLERS['Heritage'](self.data)
        heritage.createReport()
        rptName = heritage.writeReport(rebuild)
        self.addLink(precsNode, rptName, "Heritage Report")
        self.__precRpts.append(rptName)
        # Loop over all precincts, alphabetically
        pSet = self.data.precinctInfo
        pKeys = pSet.keys()
        pKeys.sort()
        for pName in pKeys:
            pInfo = pSet[pName]
            if pInfo.getPrecinct() in PRECINCT_HANDLERS:
                handler = PRECINCT_HANDLERS[pInfo.getPrecinct()]
            else:
                handler = PRECINCT_HANDLERS['DEFAULT']
            pRpt = handler(pInfo, self.data)
            pRpt.createReport(rebuild)
            rptName = pRpt.writeReport(rebuild)
            self.addLink(precsNode, rptName, "Precinct: %s;  Alias: %s"%(pInfo.getPrecinct(), pInfo.alias))
            self.__precRpts.append(rptName)

    def __getAliases(self):
        # determine the set of precinct info aliases.
        votePath = self.data.db.getVoteInfo(int(self.data.configInfo.getVoteKey())).getSrc()
        confPath = os.path.join(os.environ['MOOT_ARCHIVE'], votePath)
        tmp,confDoc = openXmlFileByName(confPath)
        for pInfo in self.data.precinctInfo.values():
            precinctNode = confDoc.xpath('/descendant::%s-vote'%pInfo.getPrecinct())[0]
            precinctAlias = str(precinctNode.xpath('child::text()')[0].data)
            pInfo.alias = precinctAlias
        

    def addDataFiles(self, rebuild=False):
        if self.data.configRootFileName(rebuild):
            self.addLink(self.rootNode, self.data.configRootFileName(), "Configuration ROOT data file")
        if self.data.baselineRootFileName(rebuild):
            self.addLink(self.rootNode, self.data.baselineRootFileName(), "Baseline ROOT data file")
        if self.data.compareRootFileName(rebuild):
            self.addLink(self.rootNode, self.data.compareRootFileName(), "Comparison ROOT data file")

    def precinctXml(self):
        return self.__precRpts

    def writeReport(self, rebuild=False):
        return super(ConfigXmlReport, self).writeReport(rebuild)


class PrecinctXmlReport(XmlReport, PrecinctReport):
    def __init__(self, precinctInfo, configData):
        PrecinctReport.__init__(self, precinctInfo, configData)
        XmlReport.__init__(self)
        self.ROOT_TAG_TYPE = TAG_PRCRPT

    def createHeader(self):
        super(PrecinctXmlReport, self).createHeader()
        setAttribute(self.rootNode, ATTR_NAME, self.info.getPrecinct())
        setAttribute(self.rootNode, ATTR_ALIAS, self.info.alias)
        setAttribute(self.rootNode, ATTR_VKEY, self.info.getKey())
        # precinct report info
        pInfo = makeChildNode(self.rootNode, TAG_PRCINFO)
        setAttribute(pInfo, ATTR_PNAME, self.info.getPrecinct())
        setAttribute(pInfo, ATTR_PKEY, int(self.info.getKey()))
        setAttribute(pInfo, ATTR_USER, self.info.getCreator())
        setAttribute(pInfo, ATTR_DATE, self.info.getTs())
        setAttribute(pInfo, ATTR_STATUS, self.info.getStatus())
        setAttribute(pInfo, ATTR_DESC, self.info.getDescrip())
        setAttribute(pInfo, ATTR_FNAME, self.info.getSrc())
        setAttribute(pInfo, ATTR_ALIAS, self.info.alias)
        for anc in self.info.ancillaries:
            ancInfo = makeChildNode(pInfo, TAG_ANCINFO)
            setAttribute(ancInfo, ATTR_NAME, anc.getClass())
            setAttribute(ancInfo, ATTR_KEY, anc.getKey())
            setAttribute(ancInfo, ATTR_USER, anc.getCreator())
            setAttribute(ancInfo, ATTR_DATE, anc.getRowCreation())
            setAttribute(ancInfo, ATTR_STATUS, anc.getStatus())
            setAttribute(ancInfo, ATTR_DESC, anc.getDescrip())
            setAttribute(ancInfo, ATTR_FNAME, anc.getSrc())
        for par in self.info.parameters:
            parInfo = makeChildNode(pInfo, TAG_PARINFO)
            setAttribute(parInfo, ATTR_NAME, par.getClass())
            setAttribute(parInfo, ATTR_KEY, par.getKey())
            setAttribute(parInfo, ATTR_USER, par.getCreator())
            setAttribute(parInfo, ATTR_DATE, par.getRowCreation())
            setAttribute(parInfo, ATTR_STATUS, par.getStatus())
            setAttribute(parInfo, ATTR_DESC, par.getDescrip())
            setAttribute(parInfo, ATTR_FNAME, par.getSrc())        
        
    def addIntent(self, parent, intentText=""):
        intent = makeChildNode(parent, TAG_INTENT)
        makeTextChildNode(intent, intentText)

    def includeText(self, parent, includeFile, nLines=-1, isHtml=False):
        tInc = makeChildNode(parent, TAG_INCTEXT)
        setAttribute(tInc, ATTR_FNAME, relpath(includeFile, self.path))
        setAttribute(tInc, ATTR_NOCONV, str(not isHtml))
        count = 0
        for line in open(includeFile, 'r').readlines():
            if count < nLines or nLines == -1:
                n = makeChildNode(tInc, TAG_TEXTLINE)
                setAttribute(n, ATTR_TEXT, line)
                setAttribute(n, ATTR_NOCONV, str(not isHtml))
            count += 1
        if nLines>=0 and count > nLines:
            setAttribute(tInc, ATTR_NLINES, nLines)

    def writeReport(self, force=False):
        return super(PrecinctXmlReport, self).writeReport(force)



def transformToFile(xslFileName, xmlFileName, htmlFileName):
    """!\brief Based on the input template and output filename parameters
    creates an HTML report.
    
    \param xslFileName  Input XSLT template file name.
    \param htmlFileName Output HTML filename.
    """
    htmlFile = file(htmlFileName, 'w+')
    xslAsUri = OsPathToUri(xslFileName)
    xmlStream = open(xmlFileName, 'r').read()
    Transform(xmlStream, xslAsUri, output=htmlFile)
    htmlFile.close()



def makePrecinctHandlers():
    PRECINCT_HANDLERS["DEFAULT"] = PrecinctXmlReport
    # ACD reports
    from AcdXmlReport import AcdBiasXmlReport, AcdHldXmlReport, AcdModeXmlReport, AcdPHAXmlReport, AcdVetoXmlReport ## These are now trivial implementations and the base class would do below
    PRECINCT_HANDLERS["ACD_Bias"] = AcdBiasXmlReport
    PRECINCT_HANDLERS["ACD_Hld"] = AcdHldXmlReport
    PRECINCT_HANDLERS["ACD_Mode"] = AcdModeXmlReport
    PRECINCT_HANDLERS["ACD_PHA"] = AcdPHAXmlReport
    PRECINCT_HANDLERS["ACD_Veto"] = AcdVetoXmlReport
    # Cal reports
    from CalXmlReport import CalLacXmlReport, CalFleXmlReport, CalFheXmlReport, CalUldXmlReport, CalModeXmlReport
    PRECINCT_HANDLERS["CAL_LAC"]    = CalLacXmlReport
    PRECINCT_HANDLERS["CAL_FLE"]    = CalFleXmlReport
    PRECINCT_HANDLERS["CAL_FHE"]    = CalFheXmlReport
    PRECINCT_HANDLERS["CAL_ULD"]    = CalUldXmlReport
    PRECINCT_HANDLERS["CAL_Mode"]   = CalModeXmlReport
    # Timing Reports
    from TimingXmlReport import AcdTimingXmlReport, GnlModeXmlReport, GnlTimingXmlReport, CalTimingXmlReport, TkrTimingXmlReport
    PRECINCT_HANDLERS["ACD_Timing"] = AcdTimingXmlReport
    PRECINCT_HANDLERS["GNL_Mode"] = GnlModeXmlReport
    PRECINCT_HANDLERS["GNL_Timing"] = GnlTimingXmlReport
    PRECINCT_HANDLERS["CAL_Timing"] = CalTimingXmlReport
    PRECINCT_HANDLERS["TKR_Timing"] = TkrTimingXmlReport
    # Trigger reports
    from TrgGemXmlReport import TrgGemXmlReport
    PRECINCT_HANDLERS["TRG_GEM"] = TrgGemXmlReport
    from TrgRoiXmlReport import TrgRoiXmlReport
    PRECINCT_HANDLERS["TRG_ROI"] = TrgRoiXmlReport
    # Tracker reports
    from TkrXmlReport import TkrModeXmlReport, TkrStripsXmlReport, TkrThreshXmlReport
    PRECINCT_HANDLERS["TKR_Mode"] = TkrModeXmlReport
    PRECINCT_HANDLERS["TKR_Strips"] = TkrStripsXmlReport
    PRECINCT_HANDLERS["TKR_Thresh"] = TkrThreshXmlReport
    # Non LATC reports
    from LciXmlReport import LciXmlReport
    PRECINCT_HANDLERS["ACD_LCI"] = LciXmlReport
    PRECINCT_HANDLERS["CAL_LCI"] = LciXmlReport
    PRECINCT_HANDLERS["TKR_LCI"] = LciXmlReport
    from IgnoreXmlReport import IgnoreXmlReport
    PRECINCT_HANDLERS["latc_ignore"] = IgnoreXmlReport
    from AssociateXmlReport import AssociateXmlReport
    PRECINCT_HANDLERS["Filter_Associate"] = AssociateXmlReport
    from HeritageXmlReport import HeritageXmlReport
    PRECINCT_HANDLERS["Heritage"] = HeritageXmlReport
    
    


def makeChildNode(parentNode, childName):
    """ makes a dom node and attaches it as a child to another node

    \param parentNode attachment point
    \param childName name of the child node
    """
    childNode = parentNode.rootNode.createElementNS(CONFIG_NAMESPACE, childName)
    parentNode.appendChild(childNode)
    return childNode

def getChildNode(parentNode, childName):
    "!@brief Get one child node with name childName"
    n = getChildNodes(parentNode, childName)
    if n:
        return n[0]
    
def getChildNodes(parentNode, childName):
    """ Finds the child nodes of the parent node with name childName
    \param childName

    return None if childName not found, else the child
    """
    cNodes = parentNode.xpath('child::%s'%childName)
    if not cNodes:
        return None
    return cNodes

def getAttribute(node, attrName):
    """ Return the value of attribute 'attrName' as a string

    return None if attrName is not an attribute of node
    """
    if node.hasAttributeNS(CONFIG_NAMESPACE, attrName):
        attr =  node.getAttributeNS(CONFIG_NAMESPACE, attrName)
        if attr is not None:
            return str(attr)
    return None

def setAttribute(node, attr, value):
    """ Wrap the setAttribute of a node.  Different XML engines use different interfaces
    """
    return node.setAttributeNS(CONFIG_NAMESPACE, attr, str(value))

def makeTextChildNode(parentNode, value):
    """ makes a dom text node

    \param parentNode attachment point
    \param value text to put in child node 
    """
    childNode = parentNode.rootNode.createTextNode(str(value))
    parentNode.appendChild(childNode)
    return childNode

def openXmlFileByName(fileName):    
    """ open an xml file for reading

    \param fileName
    \return (file,doc) the opened file object and the associated xml doc object

    after the user is done with the file object the should remember to clean up by
    calling file.close()
    """
    theFile = open(fileName,'r')
    input = Uri.OsPathToUri(fileName, attemptAbsolute=1)
    try:
        doc = NonvalidatingReader.parseUri(input)
    except:
        theFile.close()
        theFile = None
        doc = None
    return (theFile,doc)    

def setLogging(verbosity=None):
    """
    """
    import logging
    if verbosity is None:
        logging.basicConfig()
    elif verbosity == 1:
        logging.basicConfig(level=logging.WARNING)
    elif verbosity == 2:
        logging.basicConfig(level=logging.INFO)
    elif verbosity >= 3:
        logging.basicConfig(level=logging.DEBUG)
    return

def optparse():
    from optparse import OptionParser, OptionGroup
    version = __version__[__version__.find(': ')+2:-2]
    parser = OptionParser(usage = "%prog [options]",
                          version = "%prog" + " %s" % version )
    required = OptionGroup(parser, "Required", "Required flags must be specified")
    parser.add_option_group(required)
    optional = OptionGroup(parser, "Optional")
    parser.add_option_group(optional)
    expert = OptionGroup(parser, "Expert", "Expert-only option flags")
    parser.add_option_group(expert)
    required.add_option("-k", "--mootKey", dest="configKey", action="store",
                      type="int", help="MOOT configuration key")
    optional.add_option("-b", "--baselineKey", dest="baselineKey", action="store",
                      type="int", help="Baseline comparison key")
    optional.add_option("-v", "--verbose", dest="verbose", action="count",
                        help="verbosity level")
    expert.add_option("--configDir", dest="configDir", action="store",
                      type="string", help="Destination directory for config reports")
    expert.add_option("--xslTransform", dest="xslTransform", action="store",
                      type="string", help="Use a different transform file")
    expert.add_option("--rebuild", dest="rebuild", action="store_true",
                      help="do rebuild")

    parser.set_defaults(rebuild=False)
    parser.set_defaults(configDir=CONFIG_DEST_PRODUCTION)
    parser.set_defaults(xslTransform=CONFIG_XSL_TRANSFORM)
    options,args = parser.parse_args(sys.argv[1:])
    if not options.configKey:
        logging.error("Option --mootKey (-k) must be specified")
        print
        parser.print_help()
        sys.exit(-1)
    return options, args

if __name__ == '__main__':
    import sys, os

    options, args = optparse()
    setLogging(options.verbose)

    try:
        makePrecinctHandlers()
        holder = ConfigDataHolder(options.configKey,
                                  options.baselineKey,
                                  configDirBase=options.configDir)
        cr = ConfigXmlReport(holder)
        cr.createReport(options.rebuild)
        xmlName = cr.writeReport(options.rebuild)

        if options.rebuild or not os.path.exists(xmlName[:-3]+"html"):
            transformToFile(options.xslTransform,
                            xmlName,
                            xmlName[:-3]+"html")
        for pR in cr.precinctXml():
            if options.rebuild or not os.path.exists(pR[:-3]+"html"):
                transformToFile(options.xslTransform,
                                pR,
                                pR[:-3]+"html")
    except ConfigReportError, e:
        _log.error(e.message)
        sys.exit(-1)
    except Exception, e:
        _log.error("Unexpected exception recieved:")
        _log.error(e)
        raise
