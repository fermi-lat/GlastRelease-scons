#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "Timing config reporting classes"
__author__   = "P.A.Hart <philiph@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/01/25 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"

import logging, os, operator

from ConfigXmlReport import *

LPA_MODENAMES = ("NORM", "TOO", "GRB0", "GRB1", "GRB2",
                 "SOLAR", "CALIB", "DIAG")

class AssociateXmlReport(PrecinctXmlReport):
    def __init__(self, precinctInfo, configData):
        PrecinctXmlReport.__init__(self, precinctInfo, configData)
        if not self._PrecinctReport__path:
            baseDir = os.path.join(self.data.configBase, self.PDIRSTR, self.info.getPrecinct())
            relDir = ("%s-%s" % (self.info.getPrecinct(), self.info.getKey()) +
                      "-sbsKey-%d" % int(self.data.configInfo.getSbsKey()))
            self._PrecinctReport__path = os.path.join(baseDir, relDir)
        
    def createReport(self, rebuild=False):
        self.createHeader()
        self.addFilterSection(rebuild)
        self.addConstituentSection(rebuild)

    def addFilterSection(self, rebuild=False):
        assocSec = self.addSection("FSW Filter Associations")
        self.addIntent(assocSec)
        # Create a table showing the intent of this filter association
        # To do this, we must read/parse the param file
        if not self.info.parameters:
            raise RuntimeError('FilterAssociate has no parameter files!  This should be impossible.')
        f,dom = openXmlFileByName(os.path.join(os.environ['MOOT_ARCHIVE'],
                                               self.info.parameters[0].getSrc()))
        assoc = getChildNode(dom, 'Associations')
        lpa_db = getChildNodes(assoc, 'LPA_DB')
        if len(lpa_db)>1:
            logging.warning("At this time, only one LPA_DB is supported for associate reporting")
        lpa_db = lpa_db[0]
        # determine which handlers are active for which modes
        # modeTrix is a dict of dicts, keyed on modename, then id
        #   inner dict content is (handler, instance, isActive).  
        #   Ex: {'NORM': {0: ('pfc_data', 'NO_INSTANCE', False), ...}}
        modeTrix = dict([(mode, {}) for mode in LPA_MODENAMES])
        maxID = -1
        modeNodes = getChildNodes(lpa_db, 'LPA_Mode')
        for modeNode in modeNodes:
            modeSet = getAttribute(modeNode, 'name').split()
            actNode = getChildNode(modeNode, 'Filters')
            activeFilters = getAttribute(actNode, 'active').split()
            for n in getChildNodes(modeNode, 'Handler'):
                hName = getAttribute(n,'name')
                hInst = getAttribute(n,'instance')
                hID   = int(getAttribute(n, 'id'))
                for aMode in modeSet:
                    modeTrix[aMode][hID] = (hName, hInst, hName in activeFilters)
                if hID > maxID: maxID = hID
        # and finally, create and populate the table
        table = self.addTable(assocSec, caption='Note: Bold items are active filters for that mode')
        self.addTableHeader(table, '')  # empty header for ID column
        for modeName in LPA_MODENAMES:  # one column per LPA mode
            self.addTableHeader(table, modeName)
        # now, the data rows
        for anID in range(maxID+1):
            theRow = self.addTableRow(table)
            self.addTableData(theRow, anID)
            for modeName in LPA_MODENAMES:
                (hName, hInst, isActive) = modeTrix[modeName][anID] if anID in modeTrix[modeName] else ('','',False)
                if hName.endswith('_data'):
                    rowElem = hName[:-5].upper()
                elif hName.endswith('_master'):
                    rowElem = '%s:%s' % (hName[:-7].upper(), hInst)
                else:
                    rowElem = '%s:%s' % (hName.upper(), hInst)
                self.addTableData(theRow, rowElem, bold=isActive)


    def addConstituentSection(self, rebuild=False):
        # create a table containing FSW constituents
        relSec = self.addSection('FSW Constituents (relocatable)')
        cdmSec = self.addSection('FSW Constituents (CDM)')
        constInfoList = self.data.configInfo.fswConstituents
        # break constituents into 'rel' and 'cdm'
        relInfoList = [(r.getPkg(), r.getName(), r) for r in constInfoList if r.getDir() == 'rel']
        relInfoList.sort()
        relInfoList = [r for _, _, r in relInfoList]
        cdmInfoList = [(c.getPkg(), c.getSchemaId(), c.getSchemaVersionId(), c.getInstanceId(), c) for c in constInfoList if c.getDir() == 'cdm']
        cdmInfoList.sort()
        cdmInfoList = [c for _, _, _, _, c in cdmInfoList]
        
        relTable = self.addTable(relSec)
        cdmTable = self.addTable(cdmSec)
        # relocatables are pretty much dir/pkg/ver/name/logical
        # actual paths aren't that interesting
        self.addTableHeader(relTable, 'dir')
        self.addTableHeader(relTable, 'pkg')
        self.addTableHeader(relTable, 'version')
        self.addTableHeader(relTable, 'name')
        self.addTableHeader(relTable, 'logical key')
        for rel in relInfoList:
            row = self.addTableRow(relTable)
            self.addTableData(row, rel.getDir())
            self.addTableData(row, rel.getPkg())
            self.addTableData(row, rel.getPkgVersion())
            self.addTableData(row, rel.getName())
            self.addTableData(row, rel.getFswId())
        # CDMs have dir/pkg/ver/name/logical + schema_id, version_id, instance_id
        # actual paths aren't that interesting
        self.addTableHeader(cdmTable, 'dir')
        self.addTableHeader(cdmTable, 'pkg')
        self.addTableHeader(cdmTable, 'version')
        self.addTableHeader(cdmTable, 'name')
        self.addTableHeader(cdmTable, 'logical key')
        self.addTableHeader(cdmTable, 'schema id')
        self.addTableHeader(cdmTable, 'version id')
        self.addTableHeader(cdmTable, 'instance id')
        for cdm in cdmInfoList:
            row = self.addTableRow(cdmTable)
            self.addTableData(row, cdm.getDir())
            self.addTableData(row, cdm.getPkg())
            self.addTableData(row, cdm.getPkgVersion())
            self.addTableData(row, cdm.getName())
            self.addTableData(row, cdm.getFswId())
            self.addTableData(row, cdm.getSchemaId())
            self.addTableData(row, cdm.getSchemaVersionId())
            self.addTableData(row, cdm.getInstanceId())


     

            
if __name__ == '__main__':
    from ConfigReport import *
    import sys
    configKey = int(sys.argv[1])
    baseKey   = int(sys.argv[2])

    cr = ConfigDataHolder(configKey, baseKey, './')
    pInfo = cr.precinctInfo['Filter_Associate']
    pInfo.alias = 'passthru'

    tr = AssociateXmlReport(pInfo, cr)
    tr.createReport()
    pR = tr.writeReport(force=True)

    transformToFile(CONFIG_XSL_TRANSFORM, pR, pR[:-3]+'html')

