#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "Tkr config reporting classes"
__author__   = "P.A.Hart <philiph@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/01/25 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"

import logging

from ConfigXmlReport import *
from RootRptGenerator import SystemCommand
from TkrRegisterChecker import *

FN_SHORTSUM = "%sTKR_%s_shortSum.txt"

class TkrXmlReport(PrecinctXmlReport):
    def __init__(self, precinctInfo, configData, type):
        PrecinctXmlReport.__init__(self, precinctInfo, configData)
        self.__type = type
        self.__confRootFile = configData.configRootFileName()
        self.__baseRootFile = configData.baselineRootFileName()
        self.__pngFileInfos = []
        
    def createReport(self):
        self.createHeader()
        summary = self.addSection("TKR_%s_Summary" %(self.__type))

        self.addIntent(summary)  # blank intent node for later?

        self.shortSummary(self.data.configDir, rebuild=True) ## rebuild for now
        self.includeText(summary, self.__builtName, nLines=100)
        for info in self.__pngFileInfos:
            file = info.fileName
            caption = info.caption
            title = info.title
            #print 'add image: file, title, caption', file, title, caption
            self.addImage(summary, file, title, caption)
        #print 'add comment'
        self.addComment(summary, "empty comment")

    def shortSummary(self, outputStub="", rebuild=False):
        self.__builtName = FN_SHORTSUM % (outputStub, self.__type)
        if not os.path.exists(self.__builtName) or rebuild:
            tkrRegisterChecker = TkrRegisterChecker(self.__type, self.__confRootFile, self.__baseRootFile, self.__builtName,self.info.alias)
            tkrRegisterChecker.doChecks()
            self.__pngFileInfos = tkrRegisterChecker.makePngs(outputStub)
            tkrRegisterChecker.cleanupRoot()
            
class TkrModeXmlReport(TkrXmlReport):
    def __init__(self, precinctInfo, configData):
        TkrXmlReport.__init__(self, precinctInfo, configData, 'Mode')

class TkrTimingXmlReport(TkrXmlReport):
    def __init__(self, precinctInfo, configData):
        TkrXmlReport.__init__(self, precinctInfo, configData, 'Timing')

class TkrStripsXmlReport(TkrXmlReport):
    def __init__(self, precinctInfo, configData):
        TkrXmlReport.__init__(self, precinctInfo, configData, 'Strips')

class TkrThreshXmlReport(TkrXmlReport):
    def __init__(self, precinctInfo, configData):
        TkrXmlReport.__init__(self, precinctInfo, configData, 'Thresh')

