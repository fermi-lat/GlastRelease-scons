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
import os, sys

from py_mootCore import MootQuery

_log = logging.getLogger()


class ConfigReportError(Exception):
    pass

class ConfigReport(object):
    def __init__(self, configData):
        self.__confData = configData

    @property
    def data(self):
        return self.__confData

    def createReport(self):
        self.createHeader()
        self.addConfigInfo()
        self.addBaselineInfo()
        self.addPrecincts()

    def createHeader(self):
        pass

    def addConfigInfo(self):
        pass

    def addBaselineInfo(self):
        pass

    def addPrecincts(self):
        pass

    def writeReport(self):
        raise NotImplementedError





class ConfigDataHolder(object):
    """ Holder for configuration data 
    """
    def __init__(self, configKey, baselineKey=None, configDirBase="", mootQuery=None):
        self.__mq = mootQuery
        if not self.__mq:
            self.__mq = MootQuery()
        # base and delta configInfo for comparison
        self.__base = None
        self.__conf = None
        if configKey is not None:
            self.setConfig(configKey)
        if baselineKey is not None:
            self.setBaseline(baselineKey)
        self.__confRoot = ""
        self.__baseRoot = ""
        self.__compRoot = ""
        self.__confBase = configDirBase
        self.__confDir  = configDirBase + "%s/" % str(configKey)
        if not os.path.isdir(self.__confDir):
            try:
                os.mkdir(self.__confDir)
            except OSError, e:
                _log.error("ConfigReport failed in creating the configuration directory %s." % self.__confDir)
            except Exception, e:
                _log.error("An unknown error: '%s' occurred ConfigReport.__init__" % e)
            sys.exit(-1)

    @property
    def db(self):
        return self.__mq
    
    @property
    def baselineInfo(self):
        return self.__base
    def setBaseline(self, baselineKey):
        self.__base = self.__mq.getConfigInfo(baselineKey)
    def getBaseline(self):
        return int(self.baselineInfo.getKey()) if self.baselineInfo else -1
    baseline = property(getBaseline, setBaseline, doc='the baseline configuration key to compare to')

    @property
    def configInfo(self):
        return self.__conf
    def setConfig(self, configKey):
        self.__conf = self.__mq.getConfigInfo(configKey)
    def getConfig(self):
        return int(self.configInfo.getKey()) if self.configInfo else -1
    config = property(getConfig, setConfig, doc='the config key of interest')

    @property
    def configDir(self):
        return self.__confDir

    @property
    def precinctInfo(self):
        params = map(self.__mq.getParmInfo, self.__mq.configParmsUsed(self.config))
        vKeys = set([int(i.getVoteFk()) for i in params])
        return [self.__mq.getVoteInfo(int(i)) for i in vKeys]

    def configRootFileName(self, rebuild=False):
        "!@brief return the config ROOT file name.  rebuild if necessary."
        if not self.configInfo: return None
        if self.__confRoot and not rebuild: return self.__confRoot
        from RootRptGenerator import genLatcRoot
        self.__confRoot = genLatcRoot(self.config, self.configDir, rebuild)
        return self.__confRoot

    def baselineRootFileName(self, rebuild=False):
        "!@brief return the baseline ROOT file name.  rebuild if necessary."
        if not self.baselineInfo: return None
        if self.__baseRoot and not rebuild: return self.__baseRoot
        from RootRptGenerator import genLatcRoot
        self.__baseRoot = genLatcRoot(self.baseline, self.configDir, rebuild)
        return self.__baseRoot

    def compareRootFileName(self, rebuild=False):
        "!@brief return the config/baseline comparison ROOT file name.  rebuild if necessary."
        if self.__compRoot and not rebuild: return None
        if not self.configRootFileName(rebuild): return None
        if not self.baselineRootFileName(rebuild): return None
        from RootRptGenerator import compareLatcRoot
        self.__compRoot = compareLatcRoot(self.configRootFileName(),
                                          self.baselineRootFileName(),
                                          os.path.join(self.configDir, "%dvs%d_" % (self.config, self.baseline)),
                                          rebuild)
        return self.__compRoot

    def makeRelative(self, aPath):
        # strip the configDir from a path, making it relative to configDir
        if aPath.startswith(self.configDir):
            return aPath[len(self.configDir):].strip('/')
        return aPath


if __name__ == '__main__':

    baselineKey = 52
    testKey = 42
    rebuild = False

    cr = ConfigDataHolder(testKey, baselineKey, './')
    print cr.configRootFileName(srebuild)
    print cr.baselineRootFileName(rebuild)
    print cr.compareRootFileName(rebuild)

    print cr.precinctInfo()
