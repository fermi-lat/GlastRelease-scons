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

class RegisterInfo:
  """
  store data type & shape information for given LATC Register

  args:
  nInstances - number of register instances in LAT
  typeCode - ROOT.TBranch data type code character
  maxVal - maximum integer value for register
  """
  def __init__(self,
               registerName,
               precinctName,
               nInstances,
               typeCode,
               maxVal=-1):
    self.name = registerName
    self.precinctName = precinctName
    self.nInstances = nInstances
    self.typeCode = typeCode
    self.maxVal = maxVal

# describe cal LATC registers
_registerInfo = {}

import calConstant
# CAL_LAC
_registerInfo["log_acpt"] = RegisterInfo("log_acpt", "CAL_LAC", calConstant.GLOBAL_N_GCFE, 'i',127)
_registerInfo["log_acpt_bcast"] = RegisterInfo("log_acpt_bcast", "CAL_LAC",1, 'i',127)

# CAL_FLE
_registerInfo["fle_dac"] = RegisterInfo("fle_dac", "CAL_FLE", calConstant.GLOBAL_N_GCFE, 'i',127)
_registerInfo["fle_dac_bcast"] = RegisterInfo("fle_dac_bcast", "CAL_FLE", 1, 'i',127)

# CAL_FHE
_registerInfo["fhe_dac"] = RegisterInfo("fhe_dac", "CAL_FHE", calConstant.GLOBAL_N_GCFE, 'i',127)
_registerInfo["fhe_dac_bcast"] = RegisterInfo("fhe_dac_bcast", "CAL_FHE", 1, 'i',127)

# CAL_ULD
_registerInfo["rng_uld_dac"] = RegisterInfo("rng_uld_dac", "CAL_ULD", calConstant.GLOBAL_N_GCFE, 'i',127)
_registerInfo["rng_uld_dac_bcast"] = RegisterInfo("rng_uld_dac_bcast", "CAL_ULD", 1, 'i',127)

# CAL_Mode
_registerInfo["layer_mask_0"] = RegisterInfo("layer_mask_0", "CAL_Mode", calConstant.GLOBAL_N_GCCC, 'i')
_registerInfo["layer_mask_1"] = RegisterInfo("layer_mask_1", "CAL_Mode", calConstant.GLOBAL_N_GCCC, 'i')
_registerInfo["ccc_configuration"] = RegisterInfo("ccc_configuration", "CAL_Mode", calConstant.GLOBAL_N_GCCC, 'i')
_registerInfo["crc_dac"] = RegisterInfo("crc_dac", "CAL_Mode", calConstant.GLOBAL_N_GCRC, 'i')
_registerInfo["config"] = RegisterInfo("config", "CAL_Mode", calConstant.GLOBAL_N_GCRC, 'i')
_registerInfo["config_0"] = RegisterInfo("config_0", "CAL_Mode", calConstant.GLOBAL_N_GCFE, 'i')
_registerInfo["config_1"] = RegisterInfo("config_1", "CAL_Mode", calConstant.GLOBAL_N_GCFE, 'i')
_registerInfo["ref_dac"] = RegisterInfo("ref_dac", "CAL_Mode", calConstant.GLOBAL_N_GCFE, 'i')
_registerInfo["layer_mask_0_bcast"] = RegisterInfo("layer_mask_0_bcast", "CAL_Mode", 1, 'i')
_registerInfo["layer_mask_1_bcast"] = RegisterInfo("layer_mask_1_bcast", "CAL_Mode", 1, 'i')
_registerInfo["ccc_configuration_bcast"] = RegisterInfo("ccc_configuration_bcast", "CAL_Mode", 1, 'i')
_registerInfo["crc_dac_bcast"] = RegisterInfo("crc_dac_bcast", "CAL_Mode", 1, 'i')
_registerInfo["config_bcast"] = RegisterInfo("config_bcast", "CAL_Mode", 1, 'i')
_registerInfo["config_0_bcast"] = RegisterInfo("config_0_bcast", "CAL_Mode", 1, 'i')
_registerInfo["config_1_bcast"] = RegisterInfo("config_1_bcast", "CAL_Mode", 1, 'i')
_registerInfo["ref_dac_bcast"] = RegisterInfo("ref_dac_bcast", "CAL_Mode", 1, 'i')


# associate ROOT typecodes to python array typcodes
ROOT_TYPECODE_2_ARRAY = {}
ROOT_TYPECODE_2_ARRAY['i'] =  'L'

class LATCRootFile:
  """
  Provide LATC ROOT file data in form of python arrays
  """

  def __init__(self, rootPath):
    """
    open ROOT File

    args:
    rootPath - path to input LATC ROOT file
    """
    import ROOT
    self.__rootFile = ROOT.TFile(rootPath)
    self.__rootTree = self.__rootFile.Get("Config")
    # all LATC root data is in first event
    self.__rootTree.GetEntry(0)


  def getRegisterData(self, registerName):
    """
    retrieve single dimensional python array with data for specified register
    """
    leaf = self.__rootTree.GetLeaf(registerName)

    import array
    regInfo = _registerInfo[registerName]
    regData = array.array(ROOT_TYPECODE_2_ARRAY[regInfo.typeCode],
                          [0]*regInfo.nInstances)
    # loop through each value
    for idx in range(regInfo.nInstances):
      regData[idx] = int(leaf.GetValue(idx))

    return regData

class PrecinctData:
  """
  register_name -> register_value dictionary for all registers associated with given precinct.
  """
  def __init__(self, filename, precinctName):
    # associate register names with value arrays
    self._registerData = {}

    # load rootFile
    rootFile = LATCRootFile(filename)
    
    # find all registers for given precinct
    for (name, regInfo) in _registerInfo.iteritems():
      if regInfo.precinctName == precinctName:
        self._registerData[name] = rootFile.getRegisterData(name)

  def getRegisterData(self, registerName):
    """
    retrieve single dimensional python array with data for specified register
    """
    return self._registerData[registerName]


def _array2hist(title, name, data, maxVal, minVal=0):
  """
  Return ROOT TH1S from python array data
  """
  
  import ROOT
  hist = ROOT.TH1S(title,title,maxVal-minVal+1, minVal, maxVal)

  import array
  hist.FillN(len(data),
             array.array('d',data),
             array.array('d',[1]*len(data)))
  
  return hist

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
    
    self.__precinctName = precinctName

    # read in root files
    import ROOT
    # path to 'current' LATC ROOT File
    self.__rootDataPath = rootDataPath
    # current LATC register data dictionary
    self._currentRootData = PrecinctData(self.__rootDataPath, precinctName)

    # optional baseline ROOT data
    if baselineRootPath == None:
      self.__baselineRootPath == None
    else:
      self.__baselineRootPath = baselineRootPath
      self._baselineRootData = PrecinctData(self.__baselineRootPath, precinctName)

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


  def __saveHistToFile(self, hist, title, outputDir):
    """
    Save histogram to image (png) file.
    filename is autogenerated from plot title, precinct options

    args:
    hist - ROOT histogram object
    outputDir - output directory for image file
    
    return:
    final img path
    """
    imgPath = outputDir + "_" + title + ".png"
    
    hist.Draw()
    import ROOT
    ROOT.gPad.SaveAs(imgPath)

    return imgPath
    

  def __genSummaryHist(self, registerName, rootData, outputDir, title, caption):
    """
    Generate summary histogram for alli nstances of given register in given rootData

    args:
    min, max = valid range of register values 
    return:
    BasicPlotChecker.PlotInfo
    """

    regData = rootData.getRegisterData(registerName)
    regInfo = _registerInfo[registerName]

    hist = _array2hist(title, title, regData, regInfo.maxVal, 0)

    imgPath = self.__saveHistToFile(hist, title, outputDir)

    # save record of image file
    import BasicPlotChecker
    plotInfo = BasicPlotChecker.PlotInfo()
    plotInfo.save(imgPath, title, caption)

    return plotInfo
    

  def _genSummaryHists(self, registerName, outputDir):
    """
    Generate summary histograms for all instances of given register type

    args:

    return:
    list of BasicPlotChecker.PlotInfo objects
    """

    regInfo = _registerInfo[registerName]

    imgList = []
    # summarize current settings
    title = "%s_%s_summary"%(self.__precinctName, registerName)
    caption = title
    imgList.append(self.__genSummaryHist(registerName,
                                         self._currentRootData,
                                         outputDir,
                                         title,
                                         caption))

    if self._baselineRootData is not None:
      title = "%s_%s_baseline"%(self.__precinctName, registerName)
      caption = title
      imgList.append(self.__genSummaryHist(registerName,
                                           self._baselineRootData,
                                           outputDir,
                                           title,
                                           caption))


    return imgList

  def _genCCCPlot(self, registerName, outputDir):
    """
    Plot CFE data against CCC id - write to img file

    return:
    BasicPlotChecker.PlotInfo
    """
    regInfo = _registerInfo[registerName]

    title = "%s_%s_vs_GCCC"%(self.__precinctName, registerName)
    caption = title

    import ROOT
    import calConstant
    hist = ROOT.TH2S(title, title,
                     calConstant.NUM_GCCC, 0, calConstant.NUM_GCCC,
                     regInfo.maxVal+1, 0, regInfo.maxVal+1)

    hist.SetXTitle("GCCC")
    hist.SetYTitle(registerName)
    hist.SetMarkerStyle(ROOT.kFullTriangleUp)
    hist.SetMarkerSize(2)

    regData = self._currentRootData.getRegisterData(registerName)
    import array
    for idx in range(calConstant.GLOBAL_N_GCFE):
      (tem,ccc,crc,cfe) = calConstant.getCFEId(idx)
      hist.Fill(ccc, regData[idx])

    imgPath = self.__saveHistToFile(hist, title, outputDir)

    # save record of image file
    import BasicPlotChecker
    plotInfo = BasicPlotChecker.PlotInfo()
    plotInfo.save(imgPath, title, caption)

    return plotInfo

  def _genScatterPlot(self, registerName, outputDir):
    """
    Generate scatter plot between current and

    return:
    BasicPlotChecker.PlotInfo object
    """

    regInfo = _registerInfo[registerName]

    title = "%s_%s_scatter"%(self.__precinctName, registerName)
    caption = title

    import ROOT
    hist = ROOT.TH2S(title, title,
                     regInfo.maxVal+1, 0, regInfo.maxVal,
                     regInfo.maxVal+1, 0, regInfo.maxVal)

    hist.SetXTitle("BASELINE")
    hist.SetYTitle("CURRENT_CFG")
    hist.SetMarkerStyle(ROOT.kFullTriangleUp)
    hist.SetMarkerSize(2)

    cfgData = self._currentRootData.getRegisterData(registerName)
    baselineData = self._currentRootData.getRegisterData(registerName)

    import array
    hist.FillN(len(cfgData),
               array.array('d',baselineData),
               array.array('d',cfgData),
               array.array('d',[1]*len(cfgData))
               )

    imgPath = self.__saveHistToFile(hist, title, outputDir)

    # save record of image file
    import BasicPlotChecker
    plotInfo = BasicPlotChecker.PlotInfo()
    plotInfo.save(imgPath, title, caption)

    return plotInfo
    
    

class CalLACReport(CalPrecinctReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalPrecinctReport.__init__(self,
                               "CAL_LAC",
                               rootDataPath,
                               baselineRootPath)

  def makeImgs(self, outputDir):
    imgList = []
    imgList += self._genSummaryHists("log_acpt", outputDir)

    if self._baselineRootData is not None:
      imgList.append(self._genScatterPlot("log_acpt", outputDir))

    imgList.append(self._genCCCPlot("log_acpt", outputDir))

      
    return imgList

class CalFLEReport(CalPrecinctReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalPrecinctReport.__init__(self,
                               "CAL_FLE",
                               rootDataPath,
                               baselineRootPath)

  def makeImgs(self, outputDir):
    imgList = []

    imgList += self._genSummaryHists("fle_dac", outputDir)

    if self._baselineRootData is not None:
      imgList.append(self._genScatterPlot("fle_dac", outputDir))

    imgList.append(self._genCCCPlot("fle_dac", outputDir))
    return imgList
  

class CalFHEReport(CalPrecinctReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalPrecinctReport.__init__(self,
                               "CAL_FHE",
                               rootDataPath,
                               baselineRootPath)
    
    
  def makeImgs(self, outputDir):
    imgList = []
    imgList += self._genSummaryHists("fhe_dac", outputDir)
    if self._baselineRootData is not None:
      imgList.append(self._genScatterPlot("fhe_dac", outputDir))
    return imgList


class CalULDReport(CalPrecinctReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalPrecinctReport.__init__(self,
                               "CAL_ULD",
                               rootDataPath,
                               baselineRootPath)


  def makeImgs(self, outputDir):
    imgList = []
    imgList += self._genSummaryHists("rng_uld_dac", outputDir)
    if self._baselineRootData is not None:
      imgList.append(self._genScatterPlot("rng_uld_dac", outputDir))
    return imgList


class CalModeReport(CalPrecinctReport):
  def __init__(self,
               rootDataPath,
               baselineRootPath):
    CalPrecinctReport.__init__(self,
                               "CAL_Mode",
                               rootDataPath,
                               baselineRootPath)
