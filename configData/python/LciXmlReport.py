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

from ConfigXmlReport import *
from RootRptGenerator import SystemCommand


class LciXmlReport(PrecinctXmlReport):
    def __init__(self, precinctInfo, configData):
        PrecinctXmlReport.__init__(self, precinctInfo, configData)

    def createReport(self, rebuild=False):
        logging.debug("Creating report for %s" % self.info.getPrecinct())
        self.createHeader()
        summary = self.addSection("LCI_File")
        lciFile = os.path.join(os.environ['MOOT_ARCHIVE'],self.info.getSrc())
        self.includeText(summary, lciFile)
        self.addComment(summary, "empty comment")

