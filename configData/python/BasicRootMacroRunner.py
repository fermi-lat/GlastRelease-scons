#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "Run root macro, get text"
__author__   = "P.A.Hart <philiph@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/02/01 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"

#
#

import ROOT
import os  ## drop

class BasicRootMacroRunner(object):
  def __init__(self, precinctName, configurationRootFile=None, outputFile=None):
    self.__precinctName = precinctName
    self.__configurationRootFile = configurationRootFile
    self.__outputFile = outputFile
    
  def doChecks(self):
    if self.__outputFile:
      cmd = "root  -l -b -q 'dumpPrecinct.C(\"%s\", \"%s\")' > %s" %(self.__configurationRootFile, self.__precinctName,self.__outputFile)
    else:
      cmd = "root  -l -b -q 'dumpPrecinct.C(\"%s\", \"%s\")'" %(self.__configurationRootFile, self.__precinctName)
    print cmd
    os.system(cmd) ## do something smarter
    
  def write(self, msg):
    if self.__outputFile:
      self.__outputFile.write("%s\n" %(msg))
    else:
      print msg

  def savePng(self, name, title, caption):
    fileName = self.__outputStub + self.__precinctName + "_" + name
    plotInfo = PlotInfo()
    ROOT.gPad.SaveAs(fileName)
    plotInfo.save(fileName, title, caption)
    self.plotInfos.append(plotInfo)

  def setupOutputFile(self, fileName):
    self.__outputFile = open(fileName, 'w')

  def closeOutputFile(self):
    self.__outputFile.close()

# **********************
# main
# **********************
if __name__ == '__main__':

  import sys
  if len(sys.argv) == 3:
    precinctName = str( sys.argv[1] )
    rootFile = str( sys.argv[2] )
    brmr = BasicRootMacroRunner(precinctName, rootFile)
  elif len(sys.argv) == 4:
    precinctName = str( sys.argv[1] )
    rootFile = str( sys.argv[2] )
    outFile = str( sys.argv[3] )
    print "open difference root file:", rootFile
    brmr = BasicRootMacroRunner(precinctName, rootFile, outFile)
  else:
    print "usage: BasicRootMacroRunner.py "
    sys.exit()

  brmr.doChecks()
