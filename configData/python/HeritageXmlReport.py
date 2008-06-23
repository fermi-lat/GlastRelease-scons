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

import logging, os

from ConfigXmlReport import *

from py_mootCore import vectorOfUnsigned, ParmInfo, VoteInfo, AncInfo

class HeritageXmlReport(XmlReport, ReportBase):
    CFGDIRSTR = 'ConfigReports'
    ROOT_TAG_TYPE = 'HeritageReport'
    def __init__(self, configData):
        XmlReport.__init__(self)
        ReportBase.__init__(self, configData)
        self.__path = None

    @property
    def path(self):
        # return and (if necessary) create directory path where files are stored
        # path is of format /TOP/ConfigReports/<configKey>
        if not self.__path:
            self.__path = os.path.join(self.data.configBase, self.CFGDIRSTR, str(self.data.config))
            if not os.path.exists(self.__path):
                os.makedirs(self.__path)
        return self.__path

    @property
    def fileName(self):
        relFile = "HeritageReport.xml"
        return os.path.join(self.path, relFile)

    def createReport(self, rebuild=False):
        self.createHeader()
        self.uploadSection(rebuild)
        self.paramSection(rebuild)
        self.voteSection(rebuild)
        self.ancilSection(rebuild)


    def uploadSection(self, rebuild=False):
        section = self.addSection('FSW Uploads')
        cInfo = self.data.configInfo
        table = self.addTable(section)
        self.addTableHeader(table, 'FMX Logical Key')
        self.addTableHeader(table, 'FMX path')
        self.addTableHeader(table, 'MOOT src path')
        for upload in cInfo.fswInputs:
            theRow = self.addTableRow(table)
            self.addTableData(theRow, upload.getFmxKey())
            self.addTableData(theRow, upload.getFmxPath())
            self.addTableData(theRow, upload.getSrc())


    def paramSection(self, rebuild=False):
        section = self.addSection('Parameter Files')
        cInfo = self.data.configInfo
        table = self.addTable(section)
        self.addTableHeader(table, 'MOOT key')
        self.addTableHeader(table, 'MOOT Path')
        self.addTableHeader(table, 'Status')
        self.addTableHeader(table, 'Type')
        self.addTableHeader(table, 'Creation')
        self.addTableHeader(table, 'Vote Key')
        pList = list(cInfo.parameters)
        pList.sort(key=_intKey)
        for param in pList:
            theRow = self.addTableRow(table)
            self.addTableData(theRow, param.getKey())
            self.addTableData(theRow, param.getSrc())
            self.addTableData(theRow, param.getStatus())
            self.addTableData(theRow, param.getClass())
            self.addTableData(theRow, "%s @ %s" %(param.getCreator(), param.getRowCreation()))
            self.addTableData(theRow, param.getVoteFk())

    def voteSection(self, rebuild=False):
        section = self.addSection('Vote Files')
        cInfo = self.data.configInfo
        table = self.addTable(section)
        self.addTableHeader(table, 'MOOT Key')
        self.addTableHeader(table, 'MOOT Path')
        self.addTableHeader(table, 'Status')
        self.addTableHeader(table, 'Precinct')
        self.addTableHeader(table, 'Creation')
        vList = list(cInfo.votes)
        vList.sort(key=_intKey)
        for vote in vList:
            theRow = self.addTableRow(table)
            self.addTableData(theRow, vote.getKey())
            self.addTableData(theRow, vote.getSrc())
            self.addTableData(theRow, vote.getStatus())
            self.addTableData(theRow, vote.getPrecinct())
            self.addTableData(theRow, "%s @ %s" %(vote.getCreator(), vote.getTs()))

    def ancilSection(self, rebuild=False):
        section = self.addSection('Ancillary Inputs')
        cInfo = self.data.configInfo
        table = self.addTable(section)
        self.addTableHeader(table, 'MOOT Key')
        self.addTableHeader(table, 'MOOT Path')
        self.addTableHeader(table, 'Status')
        self.addTableHeader(table, 'Type')
        self.addTableHeader(table, 'Creation')
        aList = list(cInfo.ancillaries)
        aList.sort(key=_intKey)
        for ancil in aList:
            theRow = self.addTableRow(table)
            self.addTableData(theRow, ancil.getKey())
            self.addTableData(theRow, ancil.getSrc())
            self.addTableData(theRow, ancil.getStatus())
            self.addTableData(theRow, ancil.getClass())
            self.addTableData(theRow, "%s @ %s" %(ancil.getCreator(), ancil.getRowCreation()))



def _intKey(obj):
    return int(obj.getKey())

if __name__ == '__main__':
    from ConfigReport import *
    configKey = 2007
    baseKey   = 2006

    cr = ConfigDataHolder(configKey, baseKey, './')

    tr = HeritageXmlReport(cr)
    tr.createReport()
    pR = tr.writeReport(force=True)

    transformToFile(CONFIG_XSL_TRANSFORM, pR, pR[:-3]+'html')


