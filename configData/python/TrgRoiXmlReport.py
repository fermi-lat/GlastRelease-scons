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

import logging, os
from difflib import context_diff

from ConfigXmlReport import *
from RootRptGenerator import SystemCommand

CMD_SHORTSUM = "dumpGemConfiguration -r -m -M -o %s %d" # (outputStub, configKey)
FN_SHORTSUM = "TRG_ROI_shortSum_%d.txt"
FN_SHORTDIFF = "TRG_ROI_shortDiff.txt"


class TrgRoiXmlReport(PrecinctXmlReport):
    def __init__(self, precinctInfo, configData):
        PrecinctXmlReport.__init__(self, precinctInfo, configData)

    def createReport(self, rebuild=False):
        #logging.debug("Creating report for %s" % self.info.getPrecinct())
        self.createHeader()
        #self.addIntent(self.rootNode, 'intent is empty...')  # blank intent node for later?

        shortFileCfg = shortSummary(self.data.config, self.path, rebuild)
        self.includeText(self.rootNode, shortFileCfg)
        shortFileBsl = shortSummary(self.data.baseline, self.path, rebuild)
        self.addLink(self.rootNode, shortFileBsl, "Baseline summary")

        cfgLines = file(shortFileCfg, 'r').readlines()
        baseLines= file(shortFileBsl, 'r').readlines()
        cdShort = context_diff(cfgLines, baseLines)
        if cdShort:
            diffName = os.path.join(self.path,FN_SHORTDIFF)
            if rebuild or not os.path.exists(diffName):
                diffOut = file(diffName, 'w')
                diffOut.writelines(cdShort)
                diffOut.close()
            self.includeText(self.rootNode, diffName)

        self.addComment(self.rootNode, "empty comment")



def shortSummary(configKey, path, rebuild=False):
    builtName = os.path.join(path, FN_SHORTSUM % (configKey))
    if not os.path.exists(builtName) or rebuild:
        cmd = SystemCommand(CMD_SHORTSUM % (builtName, configKey))
        cmd.handle()
    return builtName


if __name__ == '__main__':
    from ConfigReport import *
    configKey = 8
    baseKey   = 9
    voteKey = 111

    cr = ConfigDataHolder(configKey, baseKey, './')

    pInfo = cr.db.getVoteInfo(voteKey)
    pInfo.alias = 'foo'

    tr = TrgGemXmlReport(pInfo, cr)

    tr.createReport()
    pR = tr.writeReport()

    transformToFile(CONFIG_XSL_TRANSFORM, pR,
                            pR[:-3]+"html")

    
