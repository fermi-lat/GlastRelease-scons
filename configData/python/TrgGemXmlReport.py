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

from PrecinctReport import *
from RootRptGenerator import SystemCommand

CMD_SHORTSUM = "dumpGemConfiguration.exe -m -M -o %s %d" # (outputStub, configKey)
FN_SHORTSUM = "%sTRG_GEM_shortSum.txt"

class TrgGemXmlReport(PrecinctXmlReport):
    def __init__(self, precinctInfo, configData):
        PrecinctXmlReport.__init__(self, precinctInfo, configData)

    def createReport(self):
        self.createHeader()
        summary = self.addSection("TRG_GEM_Summary")

        self.addIntent(summary, 'intent is empty...')  # blank intent node for later?

        shortFile = shortSummary(self.data.config, self.data.configDir)
        self.includeText(summary, shortFile, nLines=10)
        #self.includeText(summary, shortFile)

        self.addComment(summary, "empty comment")




def shortSummary(configKey, outputStub="", rebuild=False):
    builtName = FN_SHORTSUM % (outputStub)
    if not os.path.exists(builtName) or rebuild:
        cmd = SystemCommand(CMD_SHORTSUM % (builtName, configKey))
        cmd.handle()
    return builtName
