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
        voteSec = self.addSection("LCI_Vote")
        lciFile = os.path.join(os.environ['MOOT_ARCHIVE'],self.info.getSrc())
        self.includeText(voteSec, lciFile)

        parmSec = self.addSection("LCI Parameter file")
        parm = [p for p in self.info.parameters if 'lci' in p.getClass()][0] 
        parmFile = os.path.join(os.environ['MOOT_ARCHIVE'], parm.getSrc())
        self.includeText(parmSec, parmFile)


if __name__ == '__main__':
    from ConfigReport import *
    configKey = 2527
    lciVoteKey = 713
    cr = ConfigDataHolder(configKey, None, './')
    pInfo = cr.precinctInfo['CAL_LCI']
    pInfo.alias = "cal_calibGen_flt"

    report = LciXmlReport(pInfo, cr)
    report.createReport()
    xml = report.writeReport()

    transformToFile(CONFIG_XSL_TRANSFORM, xml, xml[:-3]+'html')
