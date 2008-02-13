#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "Cal precinct report code"
__author__   = "Z.Fewtrell, based on TkrRegisterChecker by P.A.Hart <philiph@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/01/25 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"


import ROOT
# set batch mode (avoid plotting)
ROOT.gROOT.SetBatch(True)


class CalPrecinctReport(object):
  def __init__(self,
               precinctName,
               rootDataPath,
               baselineRootPath):
    """
    rootDataPath - path to LATC Root data for 'current' config
    baselineRootPath - path to LATC Root data for 'baseline' comparison
    """
    object.__init__(self)
    
    self._precinctName = precinctName

    # read in root files
    # path to 'current' LATC ROOT File
    self.__rootDataPath = rootDataPath
    # current LATC register data dictionary
    import LATCRootData
    self._cfgRootData = LATCRootData.PrecinctData(self.__rootDataPath, self._precinctName)

    # optional baseline ROOT data
    if baselineRootPath == None:
      self.__baselineRootPath = None
      self._baselineRootData = None
    else:
      self.__baselineRootPath = baselineRootPath
      self._baselineRootData = LATCRootData.PrecinctData(self.__baselineRootPath, self._precinctName)

  def createTXTReport(self, filename):
    """
    Write TXT report data to output file.
    """
    txtFile = open(filename,'w')
    txtFile.close()

  def makeImgs(self, outputDir):
    """
    Generate image report files & return list describing each output img
    """
    return []


  def _genCCCDiff(self, regName, outputDir):
    """
    plot the difference (cfg - baseline) between two configs against GCCC id
    
    return:
    BasicPlotChecker.PlotInfo
    """

    import LATCRootData
    regInfo = LATCRootData.registerInfo[regName]
    title = "%s_%s_diff_vs_GCCC"%(self._precinctName, regName)
    caption = title
    
    regData = self._cfgRootData.getRegisterData(regName)
    baselineData = self._baselineRootData.getRegisterData(regName)
    diffData = [new - old for (new,old) in zip(regData,baselineData)]
    # extract GCCC index from global GCFE indeces
    import calConstant
    cccData = [calConstant.getCFEId(idx)[1] for idx in range(calConstant.GLOBAL_N_GCFE)]
    
    import ROOTPlotUtils
    return ROOTPlotUtils.make2DHist(cccData,
                                    diffData,
                                    title,
                                    caption,
                                    "GCCC",
                                    "%s diff (cfg - baseline)"%regName,
                                    outputDir,
                                    calConstant.NUM_GCCC,
                                    2*(regInfo.maxVal+1),
                                    calConstant.NUM_GCCC,
                                    regInfo.maxVal+1,
                                    0,
                                    -1*(regInfo.maxVal+1))

  def _genCCCPlot(self, regName, outputDir):
    """
    Plot CFE data against CCC id - write to img file

    return:
    BasicPlotChecker.PlotInfo
    """
    import LATCRootData
    regInfo = LATCRootData.registerInfo[regName]
    title = "%s_%s_vs_GCCC"%(self._precinctName, regName)
    caption = title

    regData = self._cfgRootData.getRegisterData(regName)
    # extract GCCC index from global GCFE indeces
    import calConstant
    cccData = [calConstant.getCFEId(idx)[1] for idx in range(calConstant.GLOBAL_N_GCFE)]

    import ROOTPlotUtils
    return ROOTPlotUtils.make2DHist(cccData,
                                    regData,
                                    title,
                                    caption,
                                    "GCCC",
                                    regName,
                                    outputDir,
                                    calConstant.NUM_GCCC,
                                    regInfo.maxVal+1,
                                    calConstant.NUM_GCCC,
                                    regInfo.maxVal+1)
                                    
  def _genDiffHist(self, regName, outputDir):
    """
    Generated 1D histogram of (new - old) difference between settings

    return:
    BasicPlotChecker.PlotInfo object
    """

    import LATCRootData
    regInfo = LATCRootData.registerInfo[regName]

    title = "%s_%s_diff"%(self._precinctName, regName)
    caption = title

    cfgData = self._cfgRootData.getRegisterData(regName)
    baselineData = self._baselineRootData.getRegisterData(regName)
    diffData = [cfg - baseline for (cfg, baseline) in zip(cfgData, baselineData)]

    import ROOTPlotUtils
    return ROOTPlotUtils.make1DHist(diffData,
                                    title,
                                    caption,
                                    "%s diff (cfg - baseline)"%regName,
                                    "",
                                    outputDir,
                                    (regInfo.maxVal+1)*2,
                                    regInfo.maxVal+1,
                                    -1*(regInfo.maxVal+1))

  def _genScatterPlot(self, regName, outputDir):
    """
    Generate scatter plot between cfg and baseline settings

    return:
    BasicPlotChecker.PlotInfo object
    """

    import LATCRootData
    regInfo = LATCRootData.registerInfo[regName]

    title = "%s_%s_scatter"%(self._precinctName, regName)
    caption = title

    cfgData = self._cfgRootData.getRegisterData(regName)
    baselineData = self._baselineRootData.getRegisterData(regName)

    import ROOTPlotUtils
    return ROOTPlotUtils.make2DHist(baselineData,
                                    cfgData,
                                    title,
                                    caption,
                                    "BASELINE",
                                    "CFG",
                                    outputDir,
                                    regInfo.maxVal+1,
                                    regInfo.maxVal+1,
                                    regInfo.maxVal+1,
                                    regInfo.maxVal+1)

class CalThreshReport(CalPrecinctReport):
  """
  Code shared by various Cal threhold precincts (LAC,FLE,FHE,ULD)
  """
  def __init__(self,
               precinctName,
               rootDataPath,
               baselineRootPath):
    CalPrecinctReport.__init__(self,
                               precinctName,
                               rootDataPath,
                               baselineRootPath)

  def createTXTReport(self, filename):
    """
    Write TXT rpeort to otuput file
    """

    txtFile = open(filename,'w')

    import calConstant
    regName = calConstant.THRESH_REG_NAME[self._precinctName]
    txtFile.write("Nominal DAC slope: %f (%s)\n"%
                  (calConstant.NOMINAL_THRESH_DAC_SLOPE[regName],
                   calConstant.THRESH_DAC_SLOPE_UNITS[regName]))

    bcast_regName = "%s_bcast"%regName
    bcast_value = self._cfgRootData.getRegisterData(bcast_regName)[0]
    txtFile.write("Broadcast value: %d\n"%bcast_value)

    cfgData = self._cfgRootData.getRegisterData(regName)
    nBcast = cfgData.count(bcast_value)
    txtFile.write("Num Broadcast channels: %d of %d\n"%(nBcast, len(cfgData)))

    txtFile.close()
  

class CalLACReport(CalThreshReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalThreshReport.__init__(self,
                             "CAL_LAC",
                             rootDataPath,
                             baselineRootPath)
    
  def makeImgs(self, outputDir):
    imgList = []
    import calConstant
    regName = calConstant.THRESH_REG_NAME[self._precinctName]

    import LATCRootData
    regInfo = LATCRootData.registerInfo[regName]
    import ROOTPlotUtils
    imgList.append(ROOTPlotUtils.make1DHist(self._cfgRootData.getRegisterData(regName),
                                            "%s_%s_cfg_summary"%(self._precinctName,regName),
                                            "%s_%s_cfg_summary"%(self._precinctName,regName),
                                            regName,
                                            "",
                                            outputDir,
                                            regInfo.maxVal+1,
                                            regInfo.maxVal+1))
    
                   
                                            

    if self._baselineRootData is not None:
      imgList.append(ROOTPlotUtils.make1DHist(self._baselineRootData.getRegisterData(regName),
                                              "%s_%s_baseline_summary"%(self._precinctName,regName),
                                              "%s_%s_baseline_summary"%(self._precinctName,regName),
                                              regName,
                                              "",
                                              outputDir,
                                              regInfo.maxVal+1,
                                              regInfo.maxVal+1))
      imgList.append(self._genScatterPlot(regName, outputDir))
      imgList.append(self._genDiffHist(regName, outputDir))

    imgList.append(self._genCCCPlot(regName, outputDir))
    if self._baselineRootData is not None:
      imgList.append(self._genCCCDiff(regName, outputDir))
    

      
    return imgList

class CalFLEReport(CalThreshReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalThreshReport.__init__(self,
                               "CAL_FLE",
                               rootDataPath,
                               baselineRootPath)

  def makeImgs(self, outputDir):
    imgList = []

    import calConstant
    regName = calConstant.THRESH_REG_NAME[self._precinctName]
    import LATCRootData
    regInfo = LATCRootData.registerInfo[regName]
    import ROOTPlotUtils
    imgList.append(ROOTPlotUtils.make1DHist(self._cfgRootData.getRegisterData(regName),
                                            "%s_%s_cfg_summary"%(self._precinctName,regName),
                                            "%s_%s_cfg_summary"%(self._precinctName,regName),
                                            regName,
                                            "",
                                            outputDir,
                                            regInfo.maxVal+1,
                                            regInfo.maxVal+1))

    if self._baselineRootData is not None:
      imgList.append(ROOTPlotUtils.make1DHist(self._baselineRootData.getRegisterData(regName),
                                              "%s_%s_baseline_summary"%(self._precinctName,regName),
                                              "%s_%s_baseline_summary"%(self._precinctName,regName),
                                              regName,
                                              "",
                                              outputDir,
                                              regInfo.maxVal+1,
                                              regInfo.maxVal+1))
      imgList.append(self._genScatterPlot(regName, outputDir))
      imgList.append(self._genDiffHist(regName, outputDir))

    imgList.append(self._genCCCPlot(regName, outputDir))
    if self._baselineRootData is not None:
      imgList.append(self._genCCCDiff(regName, outputDir))

    return imgList
  

class CalFHEReport(CalThreshReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalThreshReport.__init__(self,
                             "CAL_FHE",
                             rootDataPath,
                             baselineRootPath)
    
    
  def makeImgs(self, outputDir):
    imgList = []

    import calConstant
    regName = calConstant.THRESH_REG_NAME[self._precinctName]

    import LATCRootData
    regInfo = LATCRootData.registerInfo[regName]
    import ROOTPlotUtils
    imgList.append(ROOTPlotUtils.make1DHist(self._cfgRootData.getRegisterData(regName),
                                            "%s_%s_cfg_summary"%(self._precinctName,regName),
                                            "%s_%s_cfg_summary"%(self._precinctName,regName),
                                            regName,
                                            "",
                                            outputDir,
                                            regInfo.maxVal+1,
                                            regInfo.maxVal+1))

    if self._baselineRootData is not None:
      imgList.append(ROOTPlotUtils.make1DHist(self._baselineRootData.getRegisterData(regName),
                                              "%s_%s_baseline_summary"%(self._precinctName,regName),
                                              "%s_%s_baseline_summary"%(self._precinctName,regName),
                                              regName,
                                              "",
                                              outputDir,
                                              regInfo.maxVal+1,
                                              regInfo.maxVal+1))
      imgList.append(self._genScatterPlot(regName, outputDir))
      imgList.append(self._genDiffHist(regName, outputDir))

    return imgList


class CalULDReport(CalThreshReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalThreshReport.__init__(self,
                               "CAL_ULD",
                               rootDataPath,
                               baselineRootPath)


  def makeImgs(self, outputDir):
    imgList = []

    import calConstant
    regName = calConstant.THRESH_REG_NAME[self._precinctName]

    import LATCRootData
    regInfo = LATCRootData.registerInfo[regName]
    import ROOTPlotUtils
    imgList.append(ROOTPlotUtils.make1DHist(self._cfgRootData.getRegisterData(regName),
                                            "%s_%s_cfg_summary"%(self._precinctName,regName),
                                            "%s_%s_cfg_summary"%(self._precinctName,regName),
                                            regName,
                                            "",
                                            outputDir,
                                            regInfo.maxVal+1,
                                            regInfo.maxVal+1))

    if self._baselineRootData is not None:
      imgList.append(ROOTPlotUtils.make1DHist(self._baselineRootData.getRegisterData(regName),
                                              "%s_%s_baseline_summary"%(self._precinctName,regName),
                                              "%s_%s_baseline_summary"%(self._precinctName,regName),
                                              regName,
                                              "",
                                              outputDir,
                                              regInfo.maxVal+1,
                                              regInfo.maxVal+1))
      imgList.append(self._genScatterPlot(regName, outputDir))
      imgList.append(self._genDiffHist(regName, outputDir))

    return imgList


class CalModeReport(CalPrecinctReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalPrecinctReport.__init__(self,
                               "CAL_Mode",
                               rootDataPath,
                               baselineRootPath)
