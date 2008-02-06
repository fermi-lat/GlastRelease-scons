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
import time, re
from RootRptGenerator import SystemCommand

class BasicRootMacroRunner(object):
  def __init__(self, precinctName, configurationRootFile=None, outputFile=None):
    self.__precinctName = precinctName
    self.__configurationRootFile = configurationRootFile
    self.__outputFile = outputFile
    self.__macro = os.path.join(os.environ['CONFIGDATAROOT'], 'python', 'Config_check.C')
    self.__exemacro = os.path.join(os.environ['CONFIGDATAROOT'], 'python', 'dumpPrecinct.C')
    
  def doChecks(self):
    tmpFile = "/tmp/rootOutput_%d.tmp" %(int(time.time()*1000000))
    if self.__outputFile:
       cmd = "root  -l -b '%s(\"%s\", \"%s\", \"%s\")'" % (self.__exemacro, self.__configurationRootFile, self.__macro, self.__precinctName)#, tmpFile)
    else:
      cmd = "root  -l -b -q 'dumpPrecinct.C(\"%s\", \"%s\")'" %(self.__configurationRootFile, self.__precinctName)
    newCmd = SystemCommand(cmd)
    out, err = newCmd.handle()
    for line in out+err:
      print 'doChecks: ', line
    if self.__outputFile:
      self.htmlize(out+err, self.__outputFile)

  def htmlize(self, lines, outFile):
    #fp = open(inFile, 'r')
    #lines = fp.readlines()
    #fp.close()
    fp = open(outFile, 'w')
    for line in lines:
      formattedLine = re.sub(" ", "&nbsp;", line)
      formattedLine = re.sub("font&nbsp;", "font ", formattedLine)
      formattedLine = re.sub("\n", "<br>", formattedLine)
      fp.write('%s\n' %(formattedLine))
    fp.close()
    
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
