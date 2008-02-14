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

class RegisterFieldInfo:
  """
  Describe meaning of a register bit-set
  """
  def __init__(self,
               minBit,
               nBits,
               desc):
    """
    fieldName - name of bitfield withing register
    minBit - LSB of field (numbered from 0)
    nBits - size of bit-field
    desc - txt description
    """
    self.minBit = minBit
    self.nBits = nBits
    self.desc = desc

  def extractVal(self, regVal):
    """
    Extract register field value (as integer) from full register value, shift LSB of field to bit 0.
    """
    reg_shift = regVal/2**self.minBit

    reg_mod = reg_shift % 2**self.nBits

    return reg_mod

class RegisterInfo:
  """
  store data type & shape information for given LATC Register

  args:
  nInstances - number of register instances in LAT
  typeCode - ROOT.TBranch data type code character
  maxVal - maximum integer value for register
  fieldDict - dictionary of field_name -> RegisterFieldInfo
  """
  def __init__(self,
               nInstances,
               typeCode,
               maxVal=-1,
               fieldDict=None):
    self.nInstances = nInstances
    self.typeCode = typeCode
    self.maxVal = maxVal

    if fieldDict is not None:
      self.fieldDict = fieldDict
    else:
      self.fieldDict = {}

# describe cal LATC precincts
PRECINCT_INFO = {}

import calConstant
PRECINCT_INFO["CAL_LAC"] = {}
PRECINCT_INFO["CAL_LAC"]["log_acpt"] = RegisterInfo(calConstant.GLOBAL_N_GCFE, 'i',127)
PRECINCT_INFO["CAL_LAC"]["log_acpt_bcast"] = RegisterInfo(1, 'i',127)

PRECINCT_INFO["CAL_FLE"] = {}
PRECINCT_INFO["CAL_FLE"]["fle_dac"] = RegisterInfo(calConstant.GLOBAL_N_GCFE, 'i',127)
PRECINCT_INFO["CAL_FLE"]["fle_dac_bcast"] = RegisterInfo(1, 'i',127)

PRECINCT_INFO["CAL_FHE"] = {}
PRECINCT_INFO["CAL_FHE"]["fhe_dac"] = RegisterInfo(calConstant.GLOBAL_N_GCFE, 'i',127)
PRECINCT_INFO["CAL_FHE"]["fhe_dac_bcast"] = RegisterInfo(1, 'i',127)

PRECINCT_INFO["CAL_ULD"] = {}
PRECINCT_INFO["CAL_ULD"]["rng_uld_dac"] = RegisterInfo(calConstant.GLOBAL_N_GCFE, 'i',127)
PRECINCT_INFO["CAL_ULD"]["rng_uld_dac_bcast"] = RegisterInfo(1, 'i',127)

PRECINCT_INFO["CAL_Mode"] = {}
PRECINCT_INFO["CAL_Mode"]["layer_mask_0"] = RegisterInfo(calConstant.GLOBAL_N_GCCC, 'i')
PRECINCT_INFO["CAL_Mode"]["layer_mask_1"] = RegisterInfo(calConstant.GLOBAL_N_GCCC, 'i')
PRECINCT_INFO["CAL_Mode"]["ccc_configuration"] = RegisterInfo(calConstant.GLOBAL_N_GCCC, 'i')
PRECINCT_INFO["CAL_Mode"]["crc_dac"] = RegisterInfo(calConstant.GLOBAL_N_GCRC, 'i')
PRECINCT_INFO["CAL_Mode"]["config"] = RegisterInfo(calConstant.GLOBAL_N_GCRC, 'i')
PRECINCT_INFO["CAL_Mode"]["config_0"] = RegisterInfo(calConstant.GLOBAL_N_GCFE, 'i')
PRECINCT_INFO["CAL_Mode"]["config_1"] = RegisterInfo(calConstant.GLOBAL_N_GCFE, 'i')
PRECINCT_INFO["CAL_Mode"]["ref_dac"] = RegisterInfo(calConstant.GLOBAL_N_GCFE, 'i')
PRECINCT_INFO["CAL_Mode"]["layer_mask_0_bcast"] = RegisterInfo(1, 'i')
PRECINCT_INFO["CAL_Mode"]["layer_mask_1_bcast"] = RegisterInfo(1, 'i')
PRECINCT_INFO["CAL_Mode"]["ccc_configuration_bcast"] = RegisterInfo(1, 'i')
PRECINCT_INFO["CAL_Mode"]["crc_dac_bcast"] = RegisterInfo(1, 'i')
PRECINCT_INFO["CAL_Mode"]["config_bcast"] = RegisterInfo(1, 'i')
PRECINCT_INFO["CAL_Mode"]["config_0_bcast"] = RegisterInfo(1, 'i')
PRECINCT_INFO["CAL_Mode"]["config_1_bcast"] = RegisterInfo(1, 'i')
PRECINCT_INFO["CAL_Mode"]["ref_dac_bcast"] = RegisterInfo(1, 'i')

# describe LATC register fields
PRECINCT_INFO["CAL_Mode"]["config_0"].fieldDict["LE_GN_SEL"] = RegisterFieldInfo(0, 3, "Gain_Select mode for charge amplifier low-energy range")
PRECINCT_INFO["CAL_Mode"]["config_0"].fieldDict["HE_GN_SEL"] = RegisterFieldInfo(3, 4, "Gain_Select mode for charge amplifier high-energy range")
PRECINCT_INFO["CAL_Mode"]["config_0"].fieldDict["LE_RNG_ENA"] = RegisterFieldInfo(7, 1, "Range_Enable bit 0 for auto-range circuit (low-energy range enable)")
PRECINCT_INFO["CAL_Mode"]["config_0"].fieldDict["HE_RNG_ENA"] = RegisterFieldInfo(8, 1, "Range_Enable bit 1 for auto-range circuit (high-energy range enable)")
PRECINCT_INFO["CAL_Mode"]["config_0"].fieldDict["USE_FRST_RNG"] = RegisterFieldInfo(9, 1, "Use first range for range selection circuit, given by First_Range bits")
PRECINCT_INFO["CAL_Mode"]["config_0"].fieldDict["FIRST_RNG"] = RegisterFieldInfo(10, 2, "First Range for range-selection circuit")
PRECINCT_INFO["CAL_Mode"]["config_0"].fieldDict["OVERWRITE"] = RegisterFieldInfo(12, 1, "Overwrite bit 1 for range-selection circuit")
PRECINCT_INFO["CAL_Mode"]["config_0_bcast"].fieldDict = PRECINCT_INFO["CAL_Mode"]["config_0"].fieldDict

PRECINCT_INFO["CAL_Mode"]["config_1"].fieldDict["PREAMP_AUTO_RESET_ENA"] = RegisterFieldInfo(0, 1, "Enable bit for Preamp to automatically reset for large input pulses. Bit set to 1 enables automatic preamp reset.")
PRECINCT_INFO["CAL_Mode"]["config_1"].fieldDict["LE_TRG_ENA"] = RegisterFieldInfo(1, 1, "Trigger_Enable for LE trigger circuits")
PRECINCT_INFO["CAL_Mode"]["config_1"].fieldDict["HE_TRG_ENA"] = RegisterFieldInfo(2, 1, "Trigger_Enable for HE trigger circuits")
PRECINCT_INFO["CAL_Mode"]["config_1"].fieldDict["CALIB_GAIN"] = RegisterFieldInfo(3, 1, "Calibration_Gain bit (1) for calibration circuit . Bit set to 0 for autoranging test.")
PRECINCT_INFO["CAL_Mode"]["config_1"].fieldDict["CALIB_LE_EN"] = RegisterFieldInfo(4, 1, "Calibration_Enable for low-energy calibration circuit")
PRECINCT_INFO["CAL_Mode"]["config_1"].fieldDict["CALIB_HE_EN"] = RegisterFieldInfo(5, 1, "Calibration_Enable for high-energy calibration circuit")
PRECINCT_INFO["CAL_Mode"]["config_1_bcast"].fieldDict = PRECINCT_INFO["CAL_Mode"]["config_1"].fieldDict

PRECINCT_INFO["CAL_Mode"]["ccc_configuration"].fieldDict["DATA_ALMOST_FULL"] = RegisterFieldInfo(0, 7, "FIFO full condition (data)")
PRECINCT_INFO["CAL_Mode"]["ccc_configuration"].fieldDict["ERROR_ALMOST_FULL"] = RegisterFieldInfo(8, 6, "FIFO full condition (error)")
PRECINCT_INFO["CAL_Mode"]["ccc_configuration"].fieldDict["DIAG_ALMOST_FULL"] = RegisterFieldInfo(16, 6, "FIFO full condition (summary/diagnostic)")
PRECINCT_INFO["CAL_Mode"]["ccc_configuration"].fieldDict["LISTEN_MSB"] = RegisterFieldInfo(27, 4, "listen MSB")
PRECINCT_INFO["CAL_Mode"]["ccc_configuration"].fieldDict["OUTPUT_ENA"] = RegisterFieldInfo(31, 1, "controller output enable")
PRECINCT_INFO["CAL_Mode"]["ccc_configuration_bcast"].fieldDict = PRECINCT_INFO["CAL_Mode"]["ccc_configuration"].fieldDict

PRECINCT_INFO["CAL_Mode"]["config"].fieldDict["GCRC_PARITY"] = RegisterFieldInfo(0,3, "GCRC Parity (3 bits)")
PRECINCT_INFO["CAL_Mode"]["config"].fieldDict["GCRC_WRITE_FLAG"] = RegisterFieldInfo(8,8, "GCRC write flag")
PRECINCT_INFO["CAL_Mode"]["config_bcast"].fieldDict = PRECINCT_INFO["CAL_Mode"]["config"].fieldDict







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


  def getRegisterData(self, precinctName, registerName):
    """
    retrieve single dimensional python array with data for specified register
    """
    leaf = self.__rootTree.GetLeaf(registerName)

    import array
    regInfo = PRECINCT_INFO[precinctName][registerName]
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
    for (regName, regInfo) in PRECINCT_INFO[precinctName].iteritems():
      self._registerData[regName] = rootFile.getRegisterData(precinctName, regName)

  def getRegisterData(self, registerName):
    """
    retrieve single dimensional python array with data for specified register
    """
    return self._registerData[registerName]



  
