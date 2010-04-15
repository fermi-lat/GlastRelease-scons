#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "MOOT config reporting tool"
__author__   = "J. Panetta <panetta@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/01/25 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"

import logging
import subprocess
import os

_log = logging.getLogger()

# system commands we are using inside SystemCommand
CMD_XML2ROOT = "ConfigXml2Root -M -o %s %d"  # (outputStub, configKey)
CMD_COMPLATC = "RootConfigCompare -f -r -o %s %s %s"  # (outputStub, configName, baselineName)

# file name formats for above commands
FN_XML2ROOT  = "%sMootConfig_%d.root"  # (outputStub, configKey)
FN_COMPLATC  = "%sreport.root" # (outputStub)

# RootRptGenerator uses offline tools (through popen) to
# create the ROOT files needed for reports

class SystemCommand(subprocess.Popen):
    """ Subclass subprocess.Popen to change the defaults
    """
    def __init__(self, args):
        subprocess.Popen.__init__(self, args,
                                shell     = True,
                                stdin     = subprocess.PIPE,
                                stdout    = subprocess.PIPE,
                                stderr    = subprocess.PIPE,
                                close_fds = True)
        self.args = args
        _log.debug( "Executing system command: [%s] " % args )
    
    def communicate(self, input=None):
        # Return array of lines of stdout, stderr rather than strings containing \n
        out, err = subprocess.Popen.communicate(self, input)
        out = out.splitlines()
        err = err.splitlines()
        return out, err
    
    def handle(self):
        # Handle the output and errors.  Return stdout, stderr.
        # Can raise RuntimeError
        out, err = self.communicate()
        if self.returncode != 0:
            raise RuntimeError('Error seen in executing system command %s.  stdout: %s, stderr: %s' %\
                               (self.args, out, err))
        return out, err
    

def genLatcRoot(configKey, path, rebuild=False):
    "!@brief given a configuration key, generate the root repr for later processing"
    builtName = FN_XML2ROOT % (path, configKey)
    if not os.path.exists(builtName) or rebuild:
        cmd = SystemCommand(CMD_XML2ROOT % (path, configKey))
        cmd.handle()
    return builtName

def compareLatcRoot(configRootName, baselineName, path, rebuild=False):
    "!@brief given two file names, run through offline comparator"
    builtName = FN_COMPLATC % path
    if not os.path.exists(builtName) or rebuild:
        cmd = SystemCommand(CMD_COMPLATC % (path, configRootName, baselineName))
        cmd.handle()
    return builtName
    



if __name__ == '__main__':

    print "Regression for root wrapper code"

    n1 = genLatcRoot(52, 'nomSciOps_')
    n2 = genLatcRoot(42, 'lacCalib_HiNeg_')


    compareLatcRoot(n1, n2, '%dvs%d_' % (52,42))
