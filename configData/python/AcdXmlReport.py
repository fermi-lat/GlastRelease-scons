#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "ACD config reporting classes"
__author__   = "P.A.Hart <philiph@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/01/25 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"

import logging, os

from ConfigXmlReport import *
from RootRptGenerator import SystemCommand
from AcdPlotChecker import *

FN_SHORTSUM = "%s_shortSum.txt"

class AcdXmlReport(PrecinctXmlReport):
    def __init__(self, precinctInfo, configData):
        PrecinctXmlReport.__init__(self, precinctInfo, configData)
        self.__precinctName = precinctInfo.getPrecinct()
        self.__confRootFile = configData.configRootFileName()
        self.__baseRootFile = configData.baselineRootFileName()
        self.__compRootFile = configData.compareRootFileName()
        self.__pngFileInfos = []
        
    def createReport(self, rebuild=False):
        self.createHeader()
        summary = self.addSection("%s_Summary" %(self.__precinctName))
        self.shortSummary(self.path, rebuild)
        self.includeText(summary, self.__builtName, nLines=10)
        for info in self.__pngFileInfos:
            file = info.fileName
            caption = info.caption
            title = info.title
            self.addImage(summary, file, title, caption)
        self.longReport(rebuild)

    def shortSummary(self, path, rebuild=False):
        self.__builtName = os.path.join(path, FN_SHORTSUM % (self.__precinctName))
        if not os.path.exists(self.__builtName) or rebuild:
            fp = open(self.__builtName, 'w')
            fp.write("no comparison information available currently")
            fp.close()
            acdPlotChecker = AcdPlotChecker(self.__precinctName, self.__compRootFile)
            self.__pngFileInfos = acdPlotChecker.makePngs(path)

    def longReport(self, rebuild=False):
        #longSec = self.addSection("Details")
        pass
            
class AcdBiasXmlReport(AcdXmlReport):
    def __init__(self, precinctInfo, configData):
        AcdXmlReport.__init__(self, precinctInfo, configData)

class AcdHldXmlReport(AcdXmlReport):
    def __init__(self, precinctInfo, configData):
        AcdXmlReport.__init__(self, precinctInfo, configData)
        
class AcdModeXmlReport(AcdXmlReport):
    def __init__(self, precinctInfo, configData):
        AcdXmlReport.__init__(self, precinctInfo, configData)

class AcdPHAXmlReport(AcdXmlReport):
    def __init__(self, precinctInfo, configData):
        AcdXmlReport.__init__(self, precinctInfo, configData)

class AcdTimingXmlReport(AcdXmlReport):
    def __init__(self, precinctInfo, configData):
        AcdXmlReport.__init__(self, precinctInfo, configData)
        ## Martin's domain - just #s
        
class AcdVetoXmlReport(AcdXmlReport):
    def __init__(self, precinctInfo, configData):
        AcdXmlReport.__init__(self, precinctInfo, configData)


if __name__ == '__main__':
    from ConfigReport import *
    configKey = 100
    baseKey   = 136

    cr = ConfigDataHolder(configKey, baseKey, './')
    pInfo = cr.precinctInfo['ACD_Bias']
    pInfo.alias = 'nomGSFC'

    tr = AcdBiasXmlReport(pInfo, cr)
    tr.createReport()
    pR = tr.writeReport(force=True)

    transformToFile(CONFIG_XSL_TRANSFORM, pR, pR[:-3]+'html')
