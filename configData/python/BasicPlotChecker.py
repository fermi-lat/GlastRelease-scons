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

import ROOT, os

class BasicPlotChecker(object):
    def __init__(self, precinctName, compareRootFile, baseRootFile=None, configurationRootFile=None, outputFile=None):
        self.__precinctName = precinctName
        self.__compareRootFile = compareRootFile
        self.__configurationRootFile = configurationRootFile
        self.__baseRootFile = baseRootFile
        if outputFile == None:
            self.__outputFile = None
        else:
            self.setupOutputFile(outputFile)

        self.__resetLogs()
        
        self.__maxErrors = 10
        self.__maxWarns = 10
        self.__maxInfos = 10

        self.setupPlotSelection()
    
    def __resetLogs(self):
        self.errors = []
        self.warns = []
        self.infos = []

    def setupPlotSelection(self):
        if self.__precinctName == 'all':
            self.plotSelection = [('foo', 'bias_dac', 'bias blah blah', 'bias dac caption', 'comp'),
                                  ('foo', 'veto_dac', 'veto blah blah', 'veto dac caption', None),
                                  ('bar', 'fhe', 'fhe title blah', 'fhe blah blah caption', 'comp')
                                  ]
        else:
            raise RuntimeError, 'Please define plots'
    
    def doChecks(self, doAll=False):
        if self.__outputFile:
            self.closeOutputFile()
      
    def makePngs(self, path):
        ## for the moment put all ROOT objects here because of garbage collection issue
        if self.__compareRootFile is not None:
            compRootFile = ROOT.TFile(self.__compareRootFile) 
##        compRootFile.ls()
        if self.__configurationRootFile is not None:
            confRootFile = ROOT.TFile(self.__configurationRootFile)
##            confRootFile.ls()
        if self.__baseRootFile is not None:
            baseRootFile = ROOT.TFile(self.__baseRootFile)        
##            baseRootFile.ls()
        
        self.plotInfos = []
        self.__filePath = path
        for (precinct, plot, title, caption, option) in self.plotSelection:
            if self.__precinctName == precinct or self.__precinctName == 'all':
                if option == 'comp' and self.__compareRootFile is not None:
                    histo = compRootFile.Get(plot)
                    if histo is not None:
                        compRootFile.Get(plot).Draw()
                        self.savePng(plot+".png", title, caption)
                    else:
                        raise RuntimeError, 'No plot named %s in %s' %(plot, self.__compareRootFile)
                else:
                    if self.__configurationRootFile is None or self.__baseRootFile is None:
                        print 'unable to plot base/conf overlay plot - no root file defined'
                    else:
                        self.plotFromTree(baseRootFile, plot)
                        if self.__compareRootFile is not None:
                            self.plotFromTree(compRootFile, plot, same=True)
                        self.savePng(plot+".png", title, caption)
        return self.plotInfos

    def plotFromTree(self, file, leafName):
        raise RuntimeError, 'implement this in base class'
    
    def cleanupRoot(self):
##        import gc
##        gc.enable()
        pass
    
    def dump(self, msgType):
        if msgType == 'errors':
            list = self.errors
            max = self.__maxErrors
        elif msgType == 'warnings':
            list = self.warns
            max = self.__maxWarns
        elif msgType == 'informations':
            list = self.infos
            max = self.__maxInfos
        else:
            raise RuntimeError, 'no such message type %s' %(msgType)
        
        if len(list)>0:
            typeStr = "# of %s: %d;" %(msgType, len(list))
            if len(list) > max:
                self.write("%s: display first %d %s below." %(typeStr, max, msgType))
            else:
                self.write("%s: display all %s below" %(typeStr, msgType))
            for msg in list[:10]:
                self.write(msg)

    def write(self, msg):
        if self.__outputFile:
            self.__outputFile.write("%s\n" %(msg))
        else:
            print msg

    def savePng(self, name, title, caption):
        fileName = os.path.join(self.__filePath, self.__precinctName + "_" + name)
        plotInfo = PlotInfo()
        ROOT.gPad.SaveAs(fileName)
        plotInfo.save(fileName, title, caption)
        self.plotInfos.append(plotInfo)

    def setupOutputFile(self, fileName):
        self.__outputFile = open(fileName, 'w')

    def closeOutputFile(self):
        self.__outputFile.close()

class PlotInfo(object):
    def __init__(self):
        pass

    def save(self, name, title, caption):
        self.fileName = name
        self.title = title
        self.caption = caption

# **********************
# main
# **********************
if __name__ == '__main__':
    ROOT.gStyle.SetOptStat( 0 )
    ROOT.gStyle.SetPalette( 1 )

    precinctName = 'all'#'Mode' #'all'
    import sys
    if len(sys.argv) == 2:
        fname = str( sys.argv[1] )
        print "open difference root file:", fname
        bpc = BasicPlotChecker(precinctName, fname)
    elif len(sys.argv) == 4:
        compName = str( sys.argv[1] )
        baseName = str( sys.argv[2] )
        confName = str( sys.argv[3] )
        print "open difference root file:", compName
        bpc = BasicPlotChecker(precinctName, compName, baseName, confName)
    else:
        print "usage: checkRegisters.py [compare root file] [config root file] (reference root file)"
        sys.exit()

    bpc.makePngs('bar')
