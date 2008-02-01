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
CONFIG_DEST_PRODUCTION = "/nfs/farm/g/glast/u03/ConfigReports/"
CONFIG_DEST_JIMTEST = "/u/ec/panetta/public_html/glast/configReports/"

CONFIG_XSL_TRANSFORM = os.path.expandvars("$CONFIGDATAROOT/python/configTransform.xsl")

TAG_CFGRPT     = "ConfigReport"
TAG_CINFO      = "ConfigInfo"
TAG_BINFO      = "BaselineInfo"
TAG_PRECINCTS  = "Precincts"
TAG_PRECINCT   = "PrecinctReportLink"
TAG_COMMENT    = "Comment"
TAG_INTENT     = "Intent"
             

ATTR_CKEY      = "ConfigKey"
ATTR_BKEY      = "BaselineKey"
ATTR_USER      = "User"
ATTR_DATE      = "Date"
ATTR_KEY       = "Key"
ATTR_CNAME     = "ConfigName"
ATTR_STATUS    = "Status"
ATTR_ACTIVE    = "Active"
ATTR_DESC      = "Description"
ATTR_NAME      = "Name"
ATTR_FILELINK  = "FileLink"
ATTR_VKEY      = "VoteKey"
ATTR_FNAME     = "FileName"

PRECINCT_HANDLERS = {}

class ConfigXmlReport(ConfigReport):
    def __init__(self, configData):
        ConfigReport.__init__(self, configData)

        self.__doc = implementation.createDocument(CONFIG_NAMESPACE, None, None)
        self.__root = None
        self.__precRpts = []

    def createHeader(self):
        self.__root = makeChildNode(self.__doc, TAG_CFGRPT)
        setAttribute(self.__root, ATTR_CKEY, self.data.config)
        setAttribute(self.__root, ATTR_BKEY, self.data.baseline)
        setAttribute(self.__root, ATTR_USER, getpass.getuser())
        setAttribute(self.__root, ATTR_DATE, time.asctime())
        self.data.configRootFileName()
        self.data.baselineRootFileName()
        self.data.compareRootFileName()

    def addConfigInfo(self):
        cInfoNode = makeChildNode(self.__root, TAG_CINFO)
        self.__addCInfoTags(cInfoNode, self.data.configInfo)

    def addBaselineInfo(self):
        bInfoNode = makeChildNode(self.__root, TAG_BINFO)
        self.__addCInfoTags(bInfoNode, self.data.baselineInfo)

    def __addCInfoTags(self, node, info):
        setAttribute(node, ATTR_CNAME, info.getName())
        setAttribute(node, ATTR_KEY, info.getKey())
        setAttribute(node, ATTR_DATE, info.getCreationTime())
        setAttribute(node, ATTR_STATUS, info.getStatus())
        setAttribute(node, ATTR_ACTIVE, info.getActive())
        setAttribute(node, ATTR_DESC, info.getDescrip())
        if hasattr(info, 'getVoteKey'):
          setAttribute(node, ATTR_VKEY, info.getVoteKey())
          vInfo = self.data.db.getVoteInfo(info.getVoteKey())
          setAttribute(node, ATTR_FNAME, vInfo.getSrc())
        
    def addPrecincts(self):
        precsNode = makeChildNode(self.__root, TAG_PRECINCTS)
        for pInfo in self.data.precinctInfo:
            if pInfo.getPrecinct() in PRECINCT_HANDLERS:
                handler = PRECINCT_HANDLERS[pInfo.getPrecinct()]
            else:
                handler = PRECINCT_HANDLERS['DEFAULT']
            pRpt = handler(pInfo, self.data)
            pRpt.createReport()
            rName = self.data.makeRelative(pRpt.writeReport())
            pNode = makeChildNode(precsNode, TAG_PRECINCT)
            setAttribute(pNode, ATTR_NAME, pInfo.getPrecinct())
            setAttribute(pNode, ATTR_FILELINK, rName[:-4]) # chop '.xml' off end
            self.__precRpts.append(rName)

    def addComment(self, parent, commentText):
        comment = makeChildNode(parent, TAG_COMMENT)
        setAttribute(comment, ATTR_USER, getpass.getuser())
        setAttribute(comment, ATTR_DATE, time.asctime())
        makeTextChildNode(comment, commentText)

    def addDataFiles(self):
        pass

    def precinctXml(self):
        return self.__precRpts

    def writeReport(self):
        fileName = os.path.join(self.data.configDir, 'ConfigReport.xml')
        f = open(fileName, 'w')
        PrettyPrint(self.__doc,f)
        f.close()
        return self.data.makeRelative(fileName)

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
    from PrecinctReport import PrecinctXmlReport
    PRECINCT_HANDLERS["DEFAULT"] = PrecinctXmlReport
    # Trigger reports
    from TrgGemXmlReport import TrgGemXmlReport
    PRECINCT_HANDLERS["TRG_GEM"] = TrgGemXmlReport
    # Tracker reports
    from TkrXmlReport import TkrModeXmlReport, TkrStripsXmlReport, TkrThreshXmlReport, TkrTimingXmlReport
    PRECINCT_HANDLERS["TKR_Mode"] = TkrModeXmlReport
    PRECINCT_HANDLERS["TKR_Strips"] = TkrStripsXmlReport
    PRECINCT_HANDLERS["TKR_Thresh"] = TkrThreshXmlReport
    PRECINCT_HANDLERS["TKR_Timing"] = TkrTimingXmlReport
    from AcdXmlReport import AcdBiasXmlReport, AcdHldXmlReport, AcdModeXmlReport, AcdPHAXmlReport, AcdTimingXmlReport, AcdVetoXmlReport ## These are now trivial implementations and the base class would do below
    PRECINCT_HANDLERS["ACD_Bias"] = AcdBiasXmlReport
    PRECINCT_HANDLERS["ACD_Hld"] = AcdHldXmlReport
    PRECINCT_HANDLERS["ACD_Mode"] = AcdModeXmlReport
    PRECINCT_HANDLERS["ACD_PHA"] = AcdPHAXmlReport
    PRECINCT_HANDLERS["ACD_Veto"] = AcdVetoXmlReport
    PRECINCT_HANDLERS["ACD_Timing"] = AcdTimingXmlReport



def makeChildNode(parentNode, childName):
    """ makes a dom node and attaches it as a child to another node

    \param parentNode attachment point
    \param childName name of the child node
    """
    childNode = parentNode.rootNode.createElementNS(CONFIG_NAMESPACE, childName)
    parentNode.appendChild(childNode)
    return childNode

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
    parser.set_defaults(configDir=CONFIG_DEST_PRODUCTION)
    #parser.set_defaults(configDir=CONFIG_DEST_JIMTEST)
    parser.set_defaults(xslTransform=CONFIG_XSL_TRANSFORM)
    options,args = parser.parse_args(sys.argv[1:])
    if not options.configKey:
        logging.error("Option --mootKey (-k) must be specified")
        print
        parser.print_help()
        sys.exit(-1)
    return options, args

if __name__ == '__main__':
    import sys

    options, args = optparse()
    setLogging(options.verbose)    
    makePrecinctHandlers()
    holder = ConfigDataHolder(options.configKey,
                              options.baselineKey,
                              configDirBase=options.configDir)
    cr = ConfigXmlReport(holder)
    cr.createReport()
    xmlName = cr.writeReport()
    print xmlName
    transformToFile(options.xslTransform,
                    os.path.join(holder.configDir, xmlName),
                    os.path.join(holder.configDir, xmlName[:-3]+"html"))
    for pR in cr.precinctXml():
        transformToFile(options.xslTransform,
                        os.path.join(holder.configDir, pR),
                        os.path.join(holder.configDir, pR[:-3]+"html"))
