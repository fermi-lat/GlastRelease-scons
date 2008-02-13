#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "Utilties for reading LATC ROOT files in python"
__author__   = "Z.Fewtrell"
__date__     = "2008/01/25 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"


import ROOT
# set batch mode (avoid plotting)
ROOT.gROOT.SetBatch(True)

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
registerInfo = {}

import calConstant
# CAL_LAC
registerInfo["log_acpt"] = RegisterInfo("log_acpt", "CAL_LAC", calConstant.GLOBAL_N_GCFE, 'i',127)
registerInfo["log_acpt_bcast"] = RegisterInfo("log_acpt_bcast", "CAL_LAC",1, 'i',127)

# CAL_FLE
registerInfo["fle_dac"] = RegisterInfo("fle_dac", "CAL_FLE", calConstant.GLOBAL_N_GCFE, 'i',127)
registerInfo["fle_dac_bcast"] = RegisterInfo("fle_dac_bcast", "CAL_FLE", 1, 'i',127)

# CAL_FHE
registerInfo["fhe_dac"] = RegisterInfo("fhe_dac", "CAL_FHE", calConstant.GLOBAL_N_GCFE, 'i',127)
registerInfo["fhe_dac_bcast"] = RegisterInfo("fhe_dac_bcast", "CAL_FHE", 1, 'i',127)

# CAL_ULD
registerInfo["rng_uld_dac"] = RegisterInfo("rng_uld_dac", "CAL_ULD", calConstant.GLOBAL_N_GCFE, 'i',127)
registerInfo["rng_uld_dac_bcast"] = RegisterInfo("rng_uld_dac_bcast", "CAL_ULD", 1, 'i',127)

# CAL_Mode
registerInfo["layer_mask_0"] = RegisterInfo("layer_mask_0", "CAL_Mode", calConstant.GLOBAL_N_GCCC, 'i')
registerInfo["layer_mask_1"] = RegisterInfo("layer_mask_1", "CAL_Mode", calConstant.GLOBAL_N_GCCC, 'i')
registerInfo["ccc_configuration"] = RegisterInfo("ccc_configuration", "CAL_Mode", calConstant.GLOBAL_N_GCCC, 'i')
registerInfo["crc_dac"] = RegisterInfo("crc_dac", "CAL_Mode", calConstant.GLOBAL_N_GCRC, 'i')
registerInfo["config"] = RegisterInfo("config", "CAL_Mode", calConstant.GLOBAL_N_GCRC, 'i')
registerInfo["config_0"] = RegisterInfo("config_0", "CAL_Mode", calConstant.GLOBAL_N_GCFE, 'i')
registerInfo["config_1"] = RegisterInfo("config_1", "CAL_Mode", calConstant.GLOBAL_N_GCFE, 'i')
registerInfo["ref_dac"] = RegisterInfo("ref_dac", "CAL_Mode", calConstant.GLOBAL_N_GCFE, 'i')
registerInfo["layer_mask_0_bcast"] = RegisterInfo("layer_mask_0_bcast", "CAL_Mode", 1, 'i')
registerInfo["layer_mask_1_bcast"] = RegisterInfo("layer_mask_1_bcast", "CAL_Mode", 1, 'i')
registerInfo["ccc_configuration_bcast"] = RegisterInfo("ccc_configuration_bcast", "CAL_Mode", 1, 'i')
registerInfo["crc_dac_bcast"] = RegisterInfo("crc_dac_bcast", "CAL_Mode", 1, 'i')
registerInfo["config_bcast"] = RegisterInfo("config_bcast", "CAL_Mode", 1, 'i')
registerInfo["config_0_bcast"] = RegisterInfo("config_0_bcast", "CAL_Mode", 1, 'i')
registerInfo["config_1_bcast"] = RegisterInfo("config_1_bcast", "CAL_Mode", 1, 'i')
registerInfo["ref_dac_bcast"] = RegisterInfo("ref_dac_bcast", "CAL_Mode", 1, 'i')


# associate ROOT typecodes to python standard library array module typecodes
ROOT_TYPECODE_2_ARRAY_TYPECODE = {}
ROOT_TYPECODE_2_ARRAY_TYPECODE['i'] =  'L'

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
    regInfo = registerInfo[registerName]
    regData = array.array(ROOT_TYPECODE_2_ARRAY_TYPECODE[regInfo.typeCode],
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
    for (name, regInfo) in registerInfo.iteritems():
      if regInfo.precinctName == precinctName:
        self._registerData[name] = rootFile.getRegisterData(name)

  def getRegisterData(self, registerName):
    """
    retrieve single dimensional python array with data for specified register
    """
    return self._registerData[registerName]

