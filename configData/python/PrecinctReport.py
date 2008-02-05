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

from Ft.Xml.Domlette                 import NonvalidatingReader, implementation, PrettyPrint
from Ft.Xml                          import EMPTY_NAMESPACE, InputSource
from Ft.Lib                          import Uri
from Ft.Xml.Xvif                     import RelaxNgValidator

from ConfigXmlReport                 import *

_log = logging.getLogger('offline.configData')

# PrecinctReport is a class that attempts to abstract the
# creation of configuration difference reports on a per precinct
# basis.

TAG_PRCRPT     = "PrecinctReport"
TAG_PRCINFO    = "PrecinctInfo"
TAG_SECTION    = "ReportSection"
TAG_IMG        = "Image"
TAG_INCTEXT    = "TextInclude"
TAG_TEXTLINE   = "TextLine"


ATTR_TITLE     = "Title"
ATTR_TEXT      = "Text"
ATTR_PNAME     = "PrecinctName"
ATTR_PKEY      = "Key"
ATTR_NLINES    = "NLines"




class PrecinctReport(object):
    def __init__(self, precinctInfo, configData):
        self.__pInfo = precinctInfo
        self.__confData = configData
    @property
    def info(self):
        return self.__pInfo
    @property
    def data(self):
        return self.__confData
    
    def createReport(self):
        self.createHeader()

    def createHeader(self):
        pass

    def writeReport(self):
        raise NotImplementedError




class PrecinctXmlReport(PrecinctReport):
    def __init__(self, precinctInfo, configData):
        PrecinctReport.__init__(self, precinctInfo, configData)
        self.__doc = implementation.createDocument(CONFIG_NAMESPACE, None, None)
        self.__root = None

    def createHeader(self):
        self.__root = makeChildNode(self.__doc, TAG_PRCRPT)
        setAttribute(self.__root, ATTR_CKEY, self.data.config)
        setAttribute(self.__root, ATTR_BKEY, self.data.baseline)
        setAttribute(self.__root, ATTR_USER, getpass.getuser())
        setAttribute(self.__root, ATTR_DATE, time.asctime())
        setAttribute(self.__root, ATTR_NAME, self.info.getPrecinct())
        # precinct report info
        pInfo = makeChildNode(self.__root, TAG_PRCINFO)
        setAttribute(pInfo, ATTR_PNAME, self.info.getPrecinct())
        setAttribute(pInfo, ATTR_PKEY, int(self.info.getKey()))
        setAttribute(pInfo, ATTR_USER, self.info.getCreator())
        setAttribute(pInfo, ATTR_DATE, self.info.getTs())
        setAttribute(pInfo, ATTR_STATUS, self.info.getStatus())
        setAttribute(pInfo, ATTR_DESC, self.info.getDescrip())
        setAttribute(pInfo, ATTR_FNAME, self.info.getSrc())
        self.addLink(self.__root, "ConfigReport.xml", "Parent Report")
        
    def addSection(self, sectionName):
        sec = makeChildNode(self.__root, TAG_SECTION)
        setAttribute(sec, ATTR_NAME, sectionName)
        return sec

    def addIntent(self, parent, intentText=""):
        intent = makeChildNode(parent, TAG_INTENT)
        makeTextChildNode(intent, intentText)

    def addImage(self, parent, fileName, title, caption=""):
        imgNode = makeChildNode(parent, TAG_IMG)
        setAttribute(imgNode, ATTR_FNAME, self.data.makeRelative(fileName))
        setAttribute(imgNode, ATTR_TITLE, title)
        if caption:
            makeTextChildNode(imgNode, caption)

    def addComment(self, parent, commentText):
        comment = makeChildNode(parent, TAG_COMMENT)
        setAttribute(comment, ATTR_USER, getpass.getuser())
        setAttribute(comment, ATTR_DATE, time.asctime())
        makeTextChildNode(comment, commentText)

    def addLink(self, parent, linkTarget, linkText, noConvert=False):
        "!@brief add a link to another file.  Links are all relative to $(CWD)"
        # noConvert:  Don't convert 'xml' to 'html' in links
        link = makeChildNode(parent, TAG_LINKTO)
        setAttribute(link, ATTR_FNAME, str(linkTarget))
        setAttribute(link, ATTR_NOCONV, str(noConvert))
        setAttribute(link, ATTR_DESC, str(linkText))

    def includeText(self, parent, includeFile, nLines=-1, isHtml=False):
        tInc = makeChildNode(parent, TAG_INCTEXT)
        setAttribute(tInc, ATTR_FNAME, self.data.makeRelative(includeFile))
        setAttribute(tInc, ATTR_NOCONV, str(not isHtml))
        if nLines>=0:
            setAttribute(tInc, ATTR_NLINES, nLines)
        count = 0
        for line in open(includeFile, 'r').readlines():
            if count < nLines or nLines == -1:
                n = makeChildNode(tInc, TAG_TEXTLINE)
                setAttribute(n, ATTR_TEXT, line)
                count += 1

    def writeReport(self):
        fileName = '%s_report.xml'%self.info.getPrecinct()
        f = open(os.path.join(self.data.configDir, fileName), 'w')
        PrettyPrint(self.__doc,f)
        f.close()
        return fileName


if __name__ == '__main__':
    pass


