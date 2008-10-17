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

from py_mootCore import MootQuery, VoteInfo

_log = logging.getLogger()


class ConfigReportError(Exception):
    pass


class ReportBase(object):
    def __init__(self, data):
        self.__data = data
    @property
    def data(self):
        return self.__data

    @property
    def fileName(self):
        """!@brief  Generate a filename for this report object.

        Note:  The filename will be an absolute path.
        Note:  This function is abstract.
        """
        raise NotImplementedError

    @property
    def path(self):
        """!@brief Return a fully-qualified path for storing report data

        Note:  This function is abstract
        """
        raise NotImplementedError

    def createReport(self, rebuild=False):
        raise NotImplementedError

    def createHeader(self):
        raise NotImplementedError

    def addConfigInfo(self):
        raise NotImplementedError

    def addBaselineInfo(self):
        raise NotImplementedError

    def addSection(self, sectionName, parent):
        raise NotImplementedError

    def writeReport(self, force):
        raise NotImplementedError


class ConfigReport(ReportBase):
    CFGDIRSTR = 'ConfigReports'
    def __init__(self, configData):
        ReportBase.__init__(self, configData)
        self.__path  = None

    def createReport(self, rebuild=False):
        self.createHeader()
        self.addConfigInfo()
        self.addBaselineInfo()
        self.addDataFiles(rebuild)
        self.addPrecincts(rebuild)

    def addPrecincts(self, rebuild=False):
        raise NotImplementedError

    def addDataFiles(self, rebuild=False):
        raise NotImplementedError

    @property
    def info(self):
        return self.data.configInfo

    @property
    def path(self):
        # return and (if necessary) create directory path where files are stored
        # path is of format /TOP/ConfigReports/<configKey>
        if not self.__path:
            self.__path = os.path.join(self.data.configBase, self.CFGDIRSTR, self.info.getKey())
            if not os.path.exists(self.__path):
                os.makedirs(self.__path)
        return self.__path

    @property
    def fileName(self):
        relFile = "ConfigReport.xml"
        return os.path.join(self.path, relFile)


class PrecinctReport(ReportBase):
    PDIRSTR = 'ConfigReports/Precincts'
    def __init__(self, precinctInfo, configData):
        ReportBase.__init__(self, configData)
        self.__pInfo = precinctInfo
        self.__path  = None
    @property
    def info(self):
        return self.__pInfo

    def createReport(self, rebuild=False):
        # implemented to create a base report object
        self.createHeader()

    @property
    def path(self):
        # return and (if necessary) create directory path where files are stored
        # path is of format: /TOP/Precinct/<pName>/<pName>-<voteKey>{-<ancType>-<ancKey>}*/
        if not self.__path:
            baseDir = os.path.join(self.data.configBase, self.PDIRSTR, self.info.getPrecinct())
            relDir = "%s-%s" % (self.info.getPrecinct(), self.info.getKey())
            # for each ancillary, add a bit to the path
            #   Doing this because if ancillary changes and vote doesn't, report still is different
            for anc in self.info.ancillaries:
                relDir = "".join([relDir, "+%s-%s" % (anc.getClass(), anc.getKey())])
            self.__path = os.path.join(baseDir,relDir)
        if not os.path.exists(self.__path):
            os.makedirs(self.__path)
        return self.__path
        
    @property
    def fileName(self):
        relFile = "%s.xml" % self.info.getPrecinct()
        return os.path.join(self.path, relFile)
        
    

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
        self.setConfig(configKey)
        if baselineKey is not None:
            self.setBaseline(baselineKey)
        self.__confRoot = ""
        self.__baseRoot = ""
        self.__compRoot = ""
        self.__confBase = configDirBase
        self.__confDir  = os.path.join(configDirBase, "ConfigReports", "%s/" % str(configKey))
        if not os.path.isdir(self.__confDir):
            try:
                # create the output directory and chmod it g+rwx, o+rx
                import stat
                os.makedirs(self.__confDir)
                mode = os.stat(self.__confDir).st_mode
                os.chmod(self.__confDir, mode | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)
            except OSError, e:
                raise ConfigReportError("ConfigReport failed in creating the configuration directory %s." % self.__confDir)
        self.__pInfo = {}
        self.__fswInfo = []
        self.__vParMap = {}
        self.__vAncMap = {}

    @property
    def db(self):
        return self.__mq
    
    @property
    def baselineInfo(self):
        return self.__base
    def setBaseline(self, baselineKey):
        self.__base = self.__retrieveConfigInfo(baselineKey)
    def getBaseline(self):
        return int(self.baselineInfo.getKey()) if self.baselineInfo else -1
    baseline = property(getBaseline, setBaseline, doc='the baseline configuration key to compare to')

    @property
    def configInfo(self):
        return self.__conf
    def setConfig(self, configKey):
        self.__conf = self.__retrieveConfigInfo(configKey)
    def getConfig(self):
        return int(self.configInfo.getKey()) if self.configInfo else -1
    config = property(getConfig, setConfig, doc='the config key of interest')

    def __retrieveConfigInfo(self, key):
        # wrapper function to get all known info about a specific config
        if not isinstance(key, int):
            try:
                key = int(key)
            except TypeError:
                raise ConfigReportError("Non-integer (or coercable to integer) config key %s is disallowed"%configKey)
        cInfo = self.db.getConfigInfo(key)
        # attach FSW constituent information
        #cInfo.fswConstituents = set(self.db.constituentInfoBySbs(int(cInfo.getSbsKey())))
        cInfo.fswConstituents = set(self.db.constituentInfoBySbs(6906))
        # attach FSW uploadables (LATC/LCI files)
        cInfo.fswInputs = set(self.db.configFswInfo(key))
        # attach parameter file info (LATC fragments)
        cInfo.parameters = set(map(self.db.getParmInfo, self.db.configParmsRequest(key)))
        # attach vote file info
        votes = dict([(p.getVoteFk(), self.db.getVoteInfo(int(p.getVoteFk()))) for p in cInfo.parameters])
        cInfo.votes = set(votes.values())
        # attach ancillary file info
        from py_mootCore import vectorOfUnsigned
        pVec = vectorOfUnsigned()
        aVec = vectorOfUnsigned()
        [pVec.append(int(p.getKey())) for p in cInfo.parameters]
        self.db.getAncsFromParms(pVec, aVec)
        ancs = list(set(aVec))
        ancs.sort()
        cInfo.ancillaries = set([self.db.getAncInfo(int(a)) for a in ancs])
        return cInfo

    @property
    def configDir(self):
        return self.__confDir

    @property
    def configBase(self):
        return self.__confBase

    @property
    def precinctInfo(self):
        # retrieve or fill precinct info dict
        if not self.__pInfo:
            # base creation of precinct info.
            db = self.db
            # find out what ancils were used to build a specific param
            # and which params hook to which votes
            par2Anc = {}
            vote2Par = {}
            from py_mootCore import vectorOfUnsigned
            for p in self.configInfo.parameters:
                pVec = vectorOfUnsigned()
                aVec = vectorOfUnsigned()
                pVec.append(int(p.getKey()))
                db.getAncsFromParms(pVec, aVec) # aVec is output...
                par2Anc[p.getKey()] = [db.getAncInfo(int(i)) for i in aVec]
                if p.getVoteFk() not in vote2Par:
                    vote2Par[p.getVoteFk()] = []
                vote2Par[p.getVoteFk()].append(p)
            # Get the vote info, and hook the ancillary info to the votes
            vKeys = set(vote2Par.keys())
            pList = [self.__mq.getVoteInfo(int(i)) for i in vKeys]
            pList.sort(key=VoteInfo.getPrecinct) # sort in precinct alphabetical
            self.__pInfo = dict([(p.getPrecinct(), p) for p in pList])
            for v in self.__pInfo.values():
                v.ancillaries = []
                v.parameters  = vote2Par[v.getKey()]
                for p in v.parameters:
                    v.ancillaries.extend(par2Anc[p.getKey()])
        return self.__pInfo

    @property
    def fswInfo(self):
        if not self.__fswInfo:
            self.__fswInfo = self.__mq.configFswInfo(self.config)
        return self.__fswInfo

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
        if self.__compRoot and not rebuild: return self.__compRoot
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
        return relpath(aPath, self.configDir)



def relpath(target, base=os.curdir):
    """
    Return a relative path to the target from either the current dir or an optional base dir.
    Base can be a directory specified either as absolute or relative to current dir.
    """

    if not os.path.exists(target):
        raise OSError, 'Target does not exist: '+target

    if not os.path.isdir(base):
        raise OSError, 'Base is not a directory or does not exist: '+base

    base_list = (os.path.abspath(base)).split(os.sep)
    target_list = (os.path.abspath(target)).split(os.sep)

    # On the windows platform the target may be on a completely different drive from the base.
    if os.name in ['nt','dos','os2'] and base_list[0] <> target_list[0]:
        raise OSError, 'Target is on a different drive to base. Target: '+target_list[0].upper()+', base: '+base_list[0].upper()

    # Starting from the filepath root, work out how much of the filepath is
    # shared by base and target.
    for i in range(min(len(base_list), len(target_list))):
        if base_list[i] <> target_list[i]: break
    else:
        # If we broke out of the loop, i is pointing to the first differing path elements.
        # If we didn't break out of the loop, i is pointing to identical path elements.
        # Increment i so that in all cases it points to the first differing path elements.
        i+=1

    rel_list = [os.pardir] * (len(base_list)-i) + target_list[i:]
    return os.path.join(*rel_list)


if __name__ == '__main__':

    baselineKey = 52
    testKey = 42
    rebuild = False

    cr = ConfigDataHolder(testKey, baselineKey, './')
    print cr.configRootFileName(rebuild)
    print cr.baselineRootFileName(rebuild)
    print cr.compareRootFileName(rebuild)

    print cr.precinctInfo()
