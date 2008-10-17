#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "Just grab some subsystem config plots and make png files"
__author__   = "P.A.Hart <philiph@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/02/01 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"

#
#

from BasicPlotChecker import *

class AcdPlotChecker(BasicPlotChecker):
    def __init__(self, precinctName, compareRootFile, baseRootFile=None, configurationRootFile=None, outputFile=None):
        BasicPlotChecker.__init__(self, precinctName, compareRootFile, baseRootFile, configurationRootFile, outputFile)

    def setupPlotSelection(self):
        self.plotSelection = [('ACD_Bias', 'bias_dac', 'bias dac plot', 'bias dac caption', 'comp'),
                              ('ACD_Veto', 'veto_dac', 'veto dac', 'config veto dac settings', 'comp'),
                              ('ACD_Hld', 'hld_dac', 'hld dac', 'veto dac dist', 'comp'),
                              ('ACD_PHA', 'pha_threshold', 'pha threshold', 'pha threshold' , 'comp'),
                              ('ACD_PHA', 'pha_threshold_v', 'pha threshold v', 'pha threshold v' , 'comp'),
                              ('ACD_Bias', 'delta_bias_dac', 'difference in bias dac', 'change in bias dac setting, new - old', 'comp'),
                              ('ACD_Veto', 'delta_veto_dac', 'veto dac difference', 'change in veto dac settings, new - old', 'comp'),
                              ('ACD_Hld', 'delta_hld_dac', 'hld dac difference', 'change in hold dac settings, new - old', 'comp'),
                              ('ACD_PHA', 'delta_pha_threshold', 'pha threshold difference', 'change in pha threshold setting, new - old' , 'comp'),
                              ('ACD_PHA', 'delta_pha_threshold_v', 'pha threshold v difference', 'change in pha threshold v, new - old' , 'comp')
                              ]
