#!/bin/env python
#
#                               Copyright 2008
#                                     by
#                        The Board of Trustees of the
#                     Leland Stanford Junior University.
#                            All rights reserved.
#

__facility__ = "Online"
__abstract__ = "Tkr register checking/comparing code, based on Hiro's work"
__author__   = "P.A.Hart <philiph@SLAC.Stanford.edu> SLAC - GLAST LAT I&T/Online"
__date__     = "2008/01/25 00:00:00"
__updated__  = "$Date$"
__version__  = "$Revision$"
__release__  = "$Name$"
__credits__  = "SLAC"

#
#
import os, sys, array, struct
import ROOT, numarray
import tkrUtils
import logging

##rndm = ROOT.TRandom()

g_nTowers = tkrUtils.g_nTowers
g_nUniPlanes = tkrUtils.g_nUniPlanes
g_nTCC = tkrUtils.g_nTCC
g_nTRC = tkrUtils.g_nTRC
g_nTFE = tkrUtils.g_nFE
g_nChannels = tkrUtils.g_nChannels

tem_tkr_trgseqArray = array.array( 'I', [0]*g_nTowers )
tcc_configurationArray = array.array( 'I', [0]*g_nTowers*g_nTCC )
tcc_trg_alignmentArray = array.array( 'I', [0]*g_nTowers*g_nTCC )
trc_csrArray = array.array( 'd', [0]*g_nTowers*g_nTCC*g_nTRC )
tfe_thrdacArray = array.array( 'I', [0]*g_nTowers*g_nUniPlanes*g_nTFE )
tfe_calibdacArray = array.array( 'I', [0]*g_nTowers*g_nUniPlanes*g_nTFE )
tfe_data_maskArray = array.array( 'd', [0L]*g_nTowers*g_nUniPlanes*g_nTFE )
tfe_trg_maskArray = array.array( 'd', [0L]*g_nTowers*g_nUniPlanes*g_nTFE )
tfe_calib_maskArray = array.array( 'd', [0L]*g_nTowers*g_nUniPlanes*g_nTFE )
spt_lowArray = array.array( 'I', [0]*g_nTowers*g_nUniPlanes )
spt_highArray = array.array( 'I' , [0]*g_nTowers*g_nUniPlanes )

altSPT = [0,24,0,24,0,24,0,24,12]
altSize = [26,26,26,26,26,26,26,26,19]
taperedSize = [19,19,19,10,10,10,10,10,10]

VERBOSE = False ## revisit to make smarter

#
# load register values
#
#
# get TFE index
#
def getTFEid( index ):
  tem = index / (g_nUniPlanes*g_nTFE)
  ispt = (index%(g_nUniPlanes*g_nTFE)) / g_nTFE
  spt = tkrUtils.getSPTFromISPT( ispt )
  lname = tkrUtils.getLayerNameFromSPT( spt )
  tfe = index % g_nTFE
  feid = "TEM %d, %s (%s), FE %d" % (tem, lname, spt, tfe)
  return feid

#
# get TFE index
#
def getTRCid( index ):
  tem = index / (g_nTCC*g_nTRC)
  tcc = (index%(g_nTCC*g_nTRC)) / g_nTRC
  trc = index % g_nTRC
  lname = tkrUtils.getLayerNameFromCCRC( tcc, trc )
  trcid = "TEM %d, %s (CC %d, RC %d)" % (tem, lname, tcc, trc)
  return trcid

#
# get SPT id
#
def getSPTid( index ):
  tem = index / g_nUniPlanes
  ispt = index % g_nUniPlanes
  spt = tkrUtils.getSPTFromISPT( ispt )
  lname = tkrUtils.getLayerNameFromSPT( spt )
  sptid = "TEM %d, %s (%s)" % (tem, lname, spt)
  return sptid

#
# decode TCC configuraton register
#
def decodeTCC_CONFIGURATION( value ):
  enable = (value & 0x80000000) >> 31
  nTRC = (value & 0x78000000) >> 27
  diagFIFO = (value & 0x003F0000) >> 16
  errFIFO = (value & 0x00007F00) >> 8
  dataFIFO = (value & 0x0000007F)
  return (enable, nTRC, diagFIFO, errFIFO, dataFIFO)

#
# decode TRC CSR register
#
def decodeTRC_CSR( value ):
  LD_FT        = (value & 0x200000000L) >> 33
  LD_delay     = (value & 0x100000000L) >> 32
  LD_stretch   = (value & 0x080000000L) >> 31
  LD_CNT       = (value & 0x040000000L) >> 30
  LD_size      = (value & 0x020000000L) >> 29
  TOT_EN       = (value & 0x000200000L) >> 21
  force_no_ERR = (value & 0x000100000L) >> 20
  read_delay   = (value & 0x0000e0000L) >> 17
  OR_stretch   = (value & 0x00001f000L) >> 12
  TFE_CNT      = (value & 0x000000f80L) >> 7
  size         = (value & 0x00000007fL)
  return (LD_FT, LD_delay, LD_stretch, LD_CNT, LD_size, TOT_EN, force_no_ERR, \
          read_delay, OR_stretch, TFE_CNT, size)


#
# check register values
#
class TkrRegisterChecker(object):
  def __init__(self, precinctName, compareRootFile, baselineRootFile=None, outputFile=None, alias=None):
##    import gc
##    gc.disable()
    self.__precinctName = precinctName
    if outputFile == None:
      self.__outputFile = None
    else:
      self.setupOutputFile(outputFile)
      
    self.names = [ ( 'Timing', "tkr_trgseq", tem_tkr_trgseqArray, self.tem_tkr_trgseq ),
                   ( 'Mode', "tcc_configuration", tcc_configurationArray, self.tcc_configuration ),
                   ( 'Timing', "tcc_trg_align", tcc_trg_alignmentArray, self.tcc_trg_alignment ),
                   ( 'Mode', "trc_csr", trc_csrArray, self.trc_csr ),
                   ( 'Thresh', "threshold", tfe_thrdacArray, self.tfe_thrdac ),
                   ( 'Thresh', "injection", tfe_calibdacArray, self.tfe_calibdac ),
                   ( 'Strips', "data_mask", tfe_data_maskArray, self.tfe_data_mask ),
                   ( 'Strips', "trig_enable", tfe_trg_maskArray, self.tfe_trg_mask ),
                   ( 'Strips', "calib_mask", tfe_calib_maskArray, self.tfe_calib_mask ),
                   ( 'Mode', "low", spt_lowArray, self.tfe_spt_low ),
                   ( 'Mode', "high", spt_highArray, self.tfe_spt_high )
                   ]

    self.rFiles = []
    self.ingestRoot(compareRootFile, baselineRootFile)
    self.__values = self.getValues(self.rFiles)
    
    self.__alternatingSplits = "sptAlt" in alias
    self.__taperedSize = "taper" in alias
    self.__maxErrors = 10
    self.__maxWarns = 10
    self.__maxInfos = 10
    ROOT.gROOT.SetBatch(True)
    
  def __resetLogs(self):
    self.errors = []
    self.warns = []
    self.infos = []

  def getValues(self, rFiles):
    values = {}
    for rFile in rFiles:
      tree = rFile.Get( "Config" )
      for ( precinct, key, varray, name ) in self.names:
        tree.SetBranchAddress( key, varray )

      tree.GetEntry( 0 )

      for ( precinct, key, varray, regMethod ) in self.names:
        name = regMethod.__name__
        branch = tree.GetBranch( key )
        branch.GetEntry( 0 )
        if name in ("tfe_data_mask", "tfe_trg_mask","tfe_calib_mask","trc_csr"):
          narray = numarray.array( varray, numarray.UInt64 )
          for i in range(len(varray)):
            st = struct.pack( 'd', varray[i] ) # convert double to unsigned int64
            narray[i] = struct.unpack( 'Q', st )[0]
        else: narray = numarray.array( varray )
        if VERBOSE:
          print name, len(narray), narray[0:4], narray[-4:]
        if values.has_key( name ):
          values[name].append( narray )
        else: values[name] = [ narray ]
    return values

  def setupHistos(self):
    self.hndmask = ROOT.TH1F( "ndmask", "# of masked channels (data)", \
                              g_nChannels, 0, g_nChannels )
    self.hnewdmask = ROOT.TH1F( "newdmask", "# of new masked channels (data)", \
                                g_nChannels, -g_nChannels, g_nChannels )
    self.hntmask = ROOT.TH1F( "ntmask", "# of masked channels (trig)", \
                              g_nChannels, 0, g_nChannels )
    self.hnewtmask = ROOT.TH1F( "newtmask", "# of new masked channels (trig)", \
                                g_nChannels, -g_nChannels, g_nChannels )
    self.hthrdac = ROOT.TH1F( "thrdac", "TFE threshold DAC", 64, 0, 64 )
    self.hthrdacdiff = ROOT.TH1F( "thrdacdiff", "TFE thresh DAC change", 64, -32, 32 )
    self.hcalibdac = ROOT.TH1F( "calibdac", "TFE calib DAC", 64, 0, 64 )
    self.hcalibdacdiff = ROOT.TH1F( "calibdacdiff", "TFE calib DAC change", 64, -32, 32 )
    
  def doChecks(self, doAll=False):
    self.setupHistos()
    for (precinct, key, varray, method ) in self.names:
      if self.__precinctName == "all" or doAll or precinct == self.__precinctName:
        self.checkValues(key, varray, method)
    if self.__outputFile:
      self.closeOutputFile()
      
  def ingestRoot(self, comp, base):
    if comp:
      self.rFiles.append(ROOT.TFile(comp))
    if base:
      self.rFiles.append(ROOT.TFile(base))
    
  def checkValues(self, key, varray, registerMethod):## revisit to divide into compValues and CompVsRef
    self.__resetLogs() ## revisit - may not be the right idea
    name = registerMethod.__name__
    self.write("\n***** register: %s *****" % name)
    varray = self.__values[name][0]
    for index in range(len(varray)):
      value = varray[index]
      registerMethod(name, value, index)
    for msgType in ['errors', 'warnings', 'informations']:
      self.dump(msgType)
      
    # compare with reference if available
    #
    if len(self.__values[name]) == 1: return  ## end of checkCompValues
    
    #if name == "tfe_dac":
    #  for i in range(len(values[name][-1])):
    #    values[name][-1][i] += rndm.Gaus()*3+1.3
    diffs = varray - self.__values[name][-1]
    elms = diffs.nonzero()[0]
    if len(elms) == 0:
      self.write('ending %s check because no differences were found' %(name))
      return # nothing is different
    if name == "tfe_thrdac" or name == "tfe_calibdac":
      if name == "tfe_thrdac": hdacdiff = self.hthrdacdiff
      else: hdacdiff = self.hcalibdacdiff
      for value in diffs:
        hdacdiff.Fill( value )
      mean = hdacdiff.GetMean()
      rms = hdacdiff.GetRMS()
      self.write("INFO: %s is changed in average %.1f with rms %.1f" \
            % (name, mean, rms))
      if rms == 0: rms = 1.0

    self.__resetLogs()
    total = 0
    for index in elms:
      if name[:3] == "tcc": ID = "TEM %d TCC %d" % (index/g_nTCC, index%g_nTCC)
      elif name[:3] == "tem": ID = "TEM %d" % (index)
      elif name[:3] == "trc": ID = getTRCid( index )
      elif name[:3] == "tfe": ID = getTFEid( index )
      elif name[:3] == "spt": ID = getSPTid( index )
      if name == "tfe_thrdac" or name == "tfe_calibdac":
        self.infos.append( "INFO: %s changed in %s: %d, %d (0x%X, 0x%X)" \
                           % (name, ID, varray[index], self.__values[name][-1][index], \
                              varray[index], self.__values[name][-1][index] ) )
        if abs(diffs[index]-mean)/rms > 4.0:
          self.warns.append( "WARN: %s change is >4sigma in %s: %d, %d (0x%X, 0x%X)"\
                             % (name, ID, varray[index], self.__values[name][-1][index], \
                                varray[index], self.__values[name][-1][index] ) )          
      elif name == "tfe_data_mask" or name == "tfe_trg_mask":
        nmask = len( tkrUtils.getStripListFromMask( str(varray[index]) ) )
        rmask = len( tkrUtils.getStripListFromMask( str(self.__values[name][-1][index]) ) )
        dmask = nmask - rmask
        total += dmask
        if name == "tfe_data_mask": self.hnewdmask.Fill( dmask )
        elif name == "tfe_trg_mask": self.hnewtmask.Fill( dmask )
        if dmask > 5: 
          self.warns.append( "WARN: increase of masked strips >5 in %s: %d %d" \
                             % (ID,nmask,rmask) + "\n      %X, %X"\
                             % (varray[index],self.__values[name][-1][index]) )
        elif dmask < 0: 
          feid = getTFEid( index )
          self.errors.append( "ERR: decrease of masked strips in %s: %d %d" \
                              % (ID,nmask,rmask) + "\n     %X, %X"\
                              % (varray[index],self.__values[name][-1][index]) )
        else:
          self.infos.append( "INFO: increase of masked strips in %s: %d %d" \
                             % (ID,nmask,rmask) + "\n      %X, %X"\
                             % (varray[index],self.__values[name][-1][index]) )
      else:
        self.infos.append( "INFO: %s changed in %s: 0x%X, 0x%X" \
                           % (name, ID, varray[index], self.__values[name][-1][index]) )
        
    if total != 0:
      self.write("INFO: change of %s is %d in LAT total" % (name,total))
        
    for msgType in ['errors', 'warnings', 'informations']:
      self.dump(msgType)

  def makePngs(self, outputPath):
    #print 'make imagefiles for ', self.__precinctName, outputPath
    self.plotInfos = []
    self.__outputPath = outputPath
    if self.__precinctName == 'Strips':
      self.hndmask.Draw()
##      ROOT.gPad.SetLogy( 1 )
      self.savePng("ndmask.png", "Masked data channel plot", "masked data channel plot")
      self.hntmask.Draw()
      self.savePng("ntmask.png", "Masked trigger channel plot", "masked trigger channel plot")
    if self.__precinctName == 'Thresh' and self.hthrdac.GetEntries()>0:
      self.hthrdac.Draw()
      self.savePng("thrdac.png", "Threshold dac plot", "threshold dac plot")
    if self.__precinctName == 'Strips':
      if self.hcalibdac.GetEntries()>0:
        self.hcalibdac.Draw()
        self.savePng("calibdac.png", "TFE calib dac plot", "tfe calib dac plot")
      if self.__precinctName == 'Thresh' and self.hthrdacdiff.GetEntries()>0:
        self.hthrdacdiff.Draw()
        self.savePng("thrdacdiff.png", "TFE threshold dac change", "tfe threshold dac change")
      if self.hcalibdacdiff.GetEntries()>0:
        self.hcalibdacdiff.Draw()
        self.savePng("calibdacdiff.png", "TFE calib DAC change", "tfe calib dac change")
      if self.hnewtmask.GetEntries()>0:
        self.hnewtmask.Draw()
        self.savePng("newtmask.png", "New masked trigger channel plot", "new masked trigger channel plot")
      if self.hnewdmask.GetEntries()>0:
        self.hnewdmask.Draw()
        self.savePng("newdmask.png", "New masked data channels plot", "new masked data channel plot")
    return self.plotInfos

  def cleanupRoot(self):
    del self.hnewdmask
    del self.hndmask
    del self.hnewtmask
    del self.hntmask
    del self.hcalibdac
    del self.hthrdacdiff
    del self.hthrdac
    del self.hcalibdacdiff
##    import gc
##    gc.enable()
    
  def tem_tkr_trgseq(self, name, value, index):
    if value != 0:
      self.infos.append( "INFO: %s is not 0 in TEM %d: %d" \
                    % (name,index,value) )
    if value<0 or (value>30 and value!=64):
      self.errors.append( "ERR: %s is out of range (0,30) or 64 in TEM %d: %d" \
                          % (name,index,value) )

  def tcc_configuration(self, name, value, index):
    (enable, nTRC, diagFIFO, errFIFO, dataFIFO) \
             = decodeTCC_CONFIGURATION( value )
    if enable != 1:
      self.warns.append( "WARN: TEM %d TCC %d is disabled: %d (%X)" \
                         % (index/g_nTCC, index%g_nTCC, enable, value) )
    if nTRC != 0:
      self.warns.append( "WARN: # of TRC in TEM %d TCC %d is not 9: %d (%X)" \
                         % (index/g_nTCC, index%g_nTCC, (9-nTRC), value) )
    if diagFIFO != 7:
      self.warns.append( "WARN: diag FIFO in TEM %d TCC %d is not 7: %d (%X)" \
                         % (index/g_nTCC, index%g_nTCC, diagFIFO, value) )
    if errFIFO != 70:
      self.warns.append( "WARN: err FIFO in TEM %d TCC %d is not 70: %d (%X)" \
                         % (index/g_nTCC, index%g_nTCC, errFIFO, value) )
    if dataFIFO != 124:
      self.warns.append( "WARN: dataFIFO in TEM %d TCC %d is not 124: %d (%X)"\
                         % (index/g_nTCC, index%g_nTCC, dataFIFO, value) )

  def tcc_trg_alignment(self, name, value, index):# only check TREQ
    treq = value & 0xF
    if treq < 4 or treq > 6:
      self.infos.append( "INFO: %s is out of (4,6) in TEM %d TCC %d: %d (%X)" \
                    % (name,index/g_nTCC,index%g_nTCC,treq,value) )

  def trc_csr(self, name, value, index):
    (LD_FT, LD_delay, LD_stretch, LD_CNT, LD_size, TOT_EN, force_no_ERR, \
     read_delay, OR_stretch, TFE_CNT, size) = decodeTRC_CSR( value )
    ID = getTRCid( index )
    if LD_FT != 1:
      self.errors.append( "ERR: LD_FT of CSR in %s is not 0: %d (%X)" \
                          % (ID, LD_FT, value) )
    if LD_delay != 0:
      self.errors.append( "ERR: LD_delay of CSR in %s is not 0: %d (%X)" \
                          % (ID, LD_delay, value) )
    if LD_stretch != 1:
      self.errors.append( "ERR: LD_stretch of CSR in %s is not 1: %d (%X)" \
                          % (ID, LD_stretch, value) )
    if LD_CNT != 0:
      self.errors.append( "ERR: LD_CNT of CSR in %s is not 0: %d (%X)" \
                          % (ID, LD_CNT, value) )
    if LD_size != 1:
      self.errors.append( "ERR: LD_size of CSR in %s is not 1: %d (%X)" \
                          % (ID, LD_size, value) )
    if TOT_EN != 1:
      self.errors.append( "ERR: TOT_EN of CSR in %s is not 1: %d (%X)" \
                          % (ID, TOT_EN, value) )
    if force_no_ERR != 0:
      self.errors.append( "ERR: force_no_ERR of CSR in %s is not 0: %d (%X)" \
                          % (ID, force_no_ERR, value) )
    if read_delay != 0:
      self.errors.append( "ERR: read_delay of CSR in %s is not 0: %d (%X)" \
                          % (ID, read_delay, value) )
    if OR_stretch != 31 and OR_stretch != 16:
      self.warns.append( "WARN: OR_stretch of CSR in %s is not 16 or 31: %d (%X)" \
                         % (ID, OR_stretch, value) )
    if TFE_CNT != 0:
      self.errors.append( "ERR: TFE_CNT of CSR in %s is not 0: %d (%X)" \
                          % (ID, TFE_CNT, value) )
    if size != 13:# and size != 20 and size != 27:
##      self.warns.append( "WARN: size of CSR in %s is not 13, 20 or 27: %d (%X)"\
      self.warns.append( "WARN: size of CSR in %s is not 13: %d (%X)"\
                         % (ID, size, value) )
    trc = index % g_nTRC

    if self.__alternatingSplits:
      if size != altSize[trc]:
        self.errors.append( "ERROR: TRC buffer is inconsistent with alternated splits: for trc %d, saw %d, expected %d" %(trc, size,altSize[trc]))
      elif VERBOSE:
        self.infos.append( "INFO: TRC buffer is consistent with alternated splits: for trc %d, saw %d, expected %d" %(trc, size,altSize[trc]))

    if self.__taperedSize:
      if size != taperedSize[trc]:
        self.errors.append( "ERROR: TRC buffer is inconsistent with tapered buffer size: for trc %d, saw %d, expected %d" %(trc, size,taperedSize[trc]))
      elif VERBOSE:
        self.infos.append( "INFO: TRC buffer is consistent with tapered buffer size: for trc %d, saw %d, expected %d" %(trc, size,taperedSize[trc]))
      
  def tfe_thrdac(self, name, value, index):
    self.hthrdac.Fill( value )
    if value<22 or value>44:
      feid = getTFEid( index )
      self.errors.append( "ERR: %s is out of range (22,44) in %s: %d" \
                          % (name,feid,value) )

  def tfe_calibdac(self, name, value, index):
    self.hcalibdac.Fill( value )
    if value==0: return
    feid = getTFEid( index )
    self.infos.append( "INFO: %s is not zero in %s: %d" \
                  % (name,feid,value) )
    if value<10 or value>26:
      self.warns.append( "WARN: %s is out of range (10,26) in %s: %d" \
                         % (name,feid,value) )

  def tfe_masks(self, nmask, value, index):
    if nmask > 10: 
      feid = getTFEid( index )
      self.infos.append( "INFO: # of masked strips >10 in %s: %d" \
                         % (feid,nmask) + ", (%X)"%value )

  def tfe_data_mask(self, name, value, index):
    nmask = len( tkrUtils.getStripListFromMask( str(value) ) )
    self.hndmask.Fill( nmask )
    self.tfe_masks(nmask, value, index)
    
  def tfe_trg_mask(self, name, value, index):
    nmask = len( tkrUtils.getStripListFromMask( str(value) ) )
    self.hntmask.Fill( nmask)
    
  def tfe_calib_mask(self, name, value, index):
    if value != 0:
      feid = getTFEid( index )
      self.warns.append( "WARN: %s is not 0 in %s: %d" % (name,value,feid) )

  def tfe_spt_low(self, name, value, index):
    self.__tfe_spt(name, value, index)
    
  def tfe_spt_high(self, name, value, index):
    self.__tfe_spt(name, value, index)
    
  def __tfe_spt(self, name, value, index):
    vmax = 24
    vmin = 0
    vnom = 12
    if name.endswith("low"): offset = -1
    elif name.endswith("high"): offset = 0
    else: raise RuntimeError, 'method %d does not exist' %(name)
    
    faltSPT = True
    st = struct.pack( 'i', value ) # convert int to unsigned short
    value = struct.unpack( 'hh', st )[0]
    spt = getSPTid(index)
    trc = index % g_nTRC
    if value>vmax+offset or value<vmin+offset:
      self.errors.append( "ERR: %s is invalid in %s: %d" % (name,spt,value) )
    elif value != vnom+offset:
      self.warns.append( "WARN: %s is not %d in %s: %d" \
                         % (name,vnom+offset,spt,value) )
    if self.__alternatingSplits and value != altSPT[trc]+offset:
      self.errors.append( "INFO: SPT pattern is not consistent with alternated splits: found %d, expected %d at index %d" %(value, altSPT[trc]+offset, index))

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
        self.write("%s: display first %d %s below." %(typeStr, max, msgType), msgType)
      else:
        self.write("%s: display all %s below" %(typeStr, msgType), msgType)
      for msg in list[:10]:
        self.write(msg, msgType)


  def write(self, msg, msgType='informations'):
    if self.__outputFile:
      self.__outputFile.write("%s\n" %(msg))
    else:
      if msgType == 'errors':
        logging.error(msg)
      elif msgType == 'warnings':
        logging.warning(msg)
      elif msgType == 'informations':
        logging.debug(msg)


  def savePng(self, name, title, caption):
    fileName = os.path.join(self.__outputPath, "TKR_" + self.__precinctName + "_" + name)
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
    print "open root file:", fname
    trc = TkrRegisterChecker(precinctName, fname, alias='sptAlt')
  elif len(sys.argv) == 3:
    fname = str( sys.argv[1] )
    print "open root file:", fname
    refname = str( sys.argv[2] )
    print "open root file:", refname
    trc = TkrRegisterChecker(precinctName, fname, refname)
  else:
    print "usage: checkRegisters.py [config root file] (reference root file)"


  trc.doChecks()
  trc.makePngs()
