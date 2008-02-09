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
from difflib import context_diff

from ConfigXmlReport import *
from RootRptGenerator import SystemCommand

CMD_SHORTSUM = "dumpGemConfiguration.exe -m -M -o %s %d" # (outputStub, configKey)
FN_SHORTSUM = "%sTRG_GEM_shortSum_%d.txt"
FN_SHORTDIFF = "%sTRG_GEM_shortDiff.txt"

CMD_LONGSUM = "dumpGemConfiguration.exe -g -f -m -M -o %s %d" # (outputStub, configKey)
FN_LONGSUM = "%sTRG_GEM_longSum_%d.txt"
FN_LONGDIFF = "%sTRG_GEM_longDiff.txt"

class TrgGemXmlReport(PrecinctXmlReport):
    def __init__(self, precinctInfo, configData):
        PrecinctXmlReport.__init__(self, precinctInfo, configData)

    def createReport(self):
        #logging.debug("Creating report for %s" % self.info.getPrecinct())
        self.createHeader()
        summary = self.addSection("TRG_GEM_Summary")
        #self.addIntent(summary, 'intent is empty...')  # blank intent node for later?

        shortFileCfg = shortSummary(self.data.config, self.data.configDir)
        self.includeText(summary, shortFileCfg)
        shortFileBsl = shortSummary(self.data.baseline, self.data.configDir)
        self.addLink(summary, shortFileBsl, "Baseline summary")

        cfgLines = file(shortFileCfg, 'r').readlines()
        baseLines= file(shortFileBsl, 'r').readlines()
        cdShort = context_diff(cfgLines, baseLines)
        if cdShort:
            diffName = FN_SHORTDIFF%self.data.configDir
            diffOut = file(diffName, 'w')
            diffOut.writelines(cdShort)
            diffOut.close()
            self.includeText(summary, diffName)

        longSect = self.addSection("TRG_GEM Long Comparison")
        longFileCfg = longSummary(self.data.config, self.data.configDir)
        longFileBsl = longSummary(self.data.baseline, self.data.configDir)
        self.addLink(longSect, longFileCfg, "Full Configuration data")
        self.addLink(longSect, longFileBsl, "Full Baseline data")
        cfgLines = file(longFileCfg, 'r').readlines()
        baseLines= file(longFileBsl, 'r').readlines()
        cdLong = context_diff(cfgLines, baseLines)
        if cdLong:
            diffName = FN_LONGDIFF%self.data.configDir
            diffOut = file(diffName, 'w')
            diffOut.writelines(cdLong)
            diffOut.close()
            self.includeText(longSect, diffName)
        

        self.addComment(self.rootNode, "empty comment")



def shortSummary(configKey, outputStub="", rebuild=False):
    builtName = FN_SHORTSUM % (outputStub, configKey)
    if not os.path.exists(builtName) or rebuild:
        cmd = SystemCommand(CMD_SHORTSUM % (builtName, configKey))
        cmd.handle()
    return builtName

def longSummary(configKey, outputStub="", rebuild=False):
    builtName = FN_LONGSUM % (outputStub, configKey)
    if not os.path.exists(builtName) or rebuild:
        cmd = SystemCommand(CMD_LONGSUM % (builtName, configKey))
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

    
