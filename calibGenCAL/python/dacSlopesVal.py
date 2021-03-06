"""
Validate CAL DacSlopes calibration data in XML format.  The command
line is:

dacSlopesVal [-V] [-L <log_file>] [-r] [-R <root_file>] <xml_file>

where:
    -r              - generate ROOT output with default name
    -R <root_file>  - output validation diagnostics in ROOT file
    -L <log_file>   - save console output to log text file
    -V              - verbose; turn on debug output
    <xml_file>      - The CAL DacSlopes calibration XML file to validate.    
"""


__facility__  = "Offline"
__abstract__  = "Validate CAL DacSlopes calibration data in XML format"
__author__    = "D.L.Wood"
__date__      = "$Date$"
__version__   = "$Revision$, $Author$"
__release__   = "$Name$"
__credits__   = "NRL code 7650"



import sys, os
import getopt
import logging

import numarray
import numarray.mlab

import calCalibXML
import calConstant
                  
                  

# validation limits


lacFineSlopeWarnLow    = 0.21
lacFineSlopeWarnHigh   = 0.45

lacFineSlopeErrLow     = 0.17
lacFineSlopeErrHigh    = 0.49

lacCoarseSlopeWarnLow  = 0.44
lacCoarseSlopeWarnHigh = 0.78

lacCoarseSlopeErrLow   = 0.37
lacCoarseSlopeErrHigh  = 0.85

fleFineSlopeWarnLow    = 0.74
fleFineSlopeWarnHigh   = 1.26

fleFineSlopeErrLow     = 0.61
fleFineSlopeErrHigh    = 1.39

fleCoarseSlopeWarnLow  = 3.0
fleCoarseSlopeWarnHigh = 4.8 

fleCoarseSlopeErrLow   = 2.7
fleCoarseSlopeErrHigh  = 5.1

fheFineSlopeWarnLow    = 36.0
fheFineSlopeWarnHigh   = 72.0

fheFineSlopeErrLow     = 30.0
fheFineSlopeErrHigh    = 78.0

fheCoarseSlopeWarnLow  = 81.0
fheCoarseSlopeWarnHigh = 113.0

fheCoarseSlopeErrLow   = 55.0
fheCoarseSlopeErrHigh  = 135.0

uldLex8SlopeWarnLow    = 0.67
uldLex8SlopeWarnHigh   = 1.02

uldLex8SlopeErrLow     = 0.61
uldLex8SlopeErrHigh    = 1.07

uldLex1SlopeWarnLow    = 6.0
uldLex1SlopeWarnHigh   = 9.0

uldLex1SlopeErrLow     = 5.0
uldLex1SlopeErrHigh    = 10.0

uldHex8SlopeWarnLow    = 47.5
uldHex8SlopeWarnHigh   = 74.5

uldHex8SlopeErrLow     = 43.0
uldHex8SlopeErrHigh    = 79.0



def rootHists(dacData, uldData, rangeData, fileName):

    # create summary ULD slope histograms

    sumHists = [None, None, None]
    cs = ROOT.TCanvas('c_Slopes_ULD', 'Slopes_ULD', -1)
    cs.SetGrid()
    cs.SetLogx()
    cs.SetLogy()
    sumLeg = ROOT.TLegend(0.88, 0.88, 0.99, 0.99)

    for erng in range(3):

        hName = "h_Slopes_ULD_%s" % calConstant.CRNG[erng]       
        hs = ROOT.TH1F(hName, 'DacSlopes_Slope_ULD: %s' % fileName, 100, 0, 100)
        hs.SetLineColor(erng + 1)
        hs.SetStats(False)
        axis = hs.GetXaxis()
        axis.SetTitle('Slope (MeV/DAC)')
        axis.CenterTitle()
        axis = hs.GetYaxis()
        axis.SetTitle('Counts')
        axis.CenterTitle()
        sumHists[erng] = hs
        sumLeg.AddEntry(hs, calConstant.CRNG[erng], 'L')
        cs.Update()
        
     
    for tem in towers:          
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                    for erng in range(3):
                        hs = sumHists[erng]
                        uld = uldData[erng,tem,row,end,fe,0]                        
                        hs.Fill(uld)

   
    hMax = 0
    for erng in range(3):
        hs = sumHists[erng]
        if hs.GetMaximum() > hMax:
            hMax = hs.GetMaximum()    
        
    for erng in range(3):

        hs = sumHists[erng]
        if erng == 0:
            dopt = ''
        else:
            dopt = 'SAME'
        hs.SetMaximum(hMax)
        hs.Draw(dopt)
        cs.Update()

    sumLeg.Draw()
    cs.Update()
    cs.Write()

   
    # create summary ULD saturation histograms

    sumHists = [None, None, None]
    cs = ROOT.TCanvas('c_Saturation_ULD', 'Saturation_ULD', -1)
    cs.SetGrid()
    cs.SetLogx()
    cs.SetLogy()
    sumLeg = ROOT.TLegend(0.88, 0.88, 0.99, 0.99)

    for erng in range(3):

        hName = "h_Saturation_ULD_%s" % calConstant.CRNG[erng]       
        hs = ROOT.TH1F(hName, 'DacSlopes_Saturation_ULD: %s' % fileName, 100, 0, 10000)
        hs.SetLineColor(erng + 1)
        hs.SetStats(False)
        axis = hs.GetXaxis()
        axis.SetTitle('Saturation (MeV)')
        axis.CenterTitle()
        axis = hs.GetYaxis()
        axis.SetTitle('Counts')
        axis.CenterTitle()
        sumHists[erng] = hs
        sumLeg.AddEntry(hs, calConstant.CRNG[erng], 'L')
        cs.Update()
        
     
    for tem in towers:          
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                    for erng in range(3):
                        hs = sumHists[erng]
                        uld = uldData[erng,tem,row,end,fe,2]                        
                        hs.Fill(uld)

   
    hMax = 0
    for erng in range(3):
        hs = sumHists[erng]
        if hs.GetMaximum() > hMax:
            hMax = hs.GetMaximum()    
        
    for erng in range(3):

        hs = sumHists[erng]
        if erng == 0:
            dopt = ''
        else:
            dopt = 'SAME'
        hs.SetMaximum(hMax)
        hs.Draw(dopt)
        cs.Update()

    sumLeg.Draw()
    cs.Update()
    cs.Write()
    

    # create summary DAC slope histograms

    sumHists = [None, None, None]
    cs = ROOT.TCanvas('c_Slopes_DAC', 'Slopes_DAC', -1)
    cs.SetGrid()
    cs.SetLogx()
    cs.SetLogy()
    sumLeg = ROOT.TLegend(0.88, 0.88, 0.99, 0.99)
    
    dacNames = ('LAC', 'FLE', 'FHE')

    for dac in range(3):

        hName = "h_Slopes_DAC_%s" % dacNames[dac]       
        hs = ROOT.TH1F(hName, 'DacSlopes_Slopes_DAC: %s' % fileName, 100, 0, 100)
        hs.SetLineColor(dac + 1)
        hs.SetStats(False)
        axis = hs.GetXaxis()
        axis.SetTitle('Slope (MeV/DAC)')
        axis.CenterTitle()
        axis = hs.GetYaxis()
        axis.SetTitle('Counts')
        axis.CenterTitle()
        sumHists[dac] = hs
        sumLeg.AddEntry(hs, dacNames[dac], 'L')
        cs.Update()
        
     
    for tem in towers:          
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                    for dac in range(3):
                        hs = sumHists[dac]
                        d = dacData[tem,row,end,fe, (dac * 2)]                        
                        hs.Fill(d)

   
    hMax = 0
    for dac in range(3):
        hs = sumHists[dac]
        if hs.GetMaximum() > hMax:
            hMax = hs.GetMaximum()    
        
    for dac in range(3):

        hs = sumHists[dac]
        if dac == 0:
            dopt = ''
        else:
            dopt = 'SAME'
        hs.SetMaximum(hMax)
        hs.Draw(dopt)
        cs.Update()

    sumLeg.Draw()
    cs.Update()
    cs.Write()
    
    
    # create summary DAC offset histograms

    sumHists = [None, None, None]
    cs = ROOT.TCanvas('c_Offsets_DAC', 'Offsets_DAC', -1)
    cs.SetGrid()
    cs.SetLogx()
    cs.SetLogy()
    sumLeg = ROOT.TLegend(0.88, 0.88, 0.99, 0.99)
    
    dacNames = ('LAC', 'FLE', 'FHE')

    for dac in range(3):

        hName = "h_Offsets_DAC_%s" % dacNames[dac]       
        hs = ROOT.TH1F(hName, 'DacSlopes_Offsets_DAC: %s' % fileName, 100, 0, 1000)
        hs.SetLineColor(dac + 1)
        hs.SetStats(False)
        axis = hs.GetXaxis()
        axis.SetTitle('Offset (MeV)')
        axis.CenterTitle()
        axis = hs.GetYaxis()
        axis.SetTitle('Counts')
        axis.CenterTitle()
        sumHists[dac] = hs
        sumLeg.AddEntry(hs, dacNames[dac], 'L')
        cs.Update()
        
     
    for tem in towers:          
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                    for dac in range(3):
                        hs = sumHists[dac]
                        d = dacData[tem,row,end,fe, (dac * 2) + 1]                        
                        hs.Fill(d)

   
    hMax = 0
    for dac in range(3):
        hs = sumHists[dac]
        if hs.GetMaximum() > hMax:
            hMax = hs.GetMaximum()    
        
    for dac in range(3):

        hs = sumHists[dac]
        if dac == 0:
            dopt = ''
        else:
            dopt = 'SAME'
        hs.SetMaximum(hMax)
        hs.Draw(dopt)
        cs.Update()

    sumLeg.Draw()
    cs.Update()
    cs.Write()



def calcError(dacData, uldData, rangeData):

    status = 0
    
    # check LAC slopes and offsets
    
    for tem in towers:
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                
                    slope = dacData[tem, row, end, fe, 0]
                    offset = dacData[tem, row, end, fe, 1]
                    rng = int(rangeData[tem, row, end, fe, 0])
                    
                    if rng == calConstant.CDAC_FINE:
                        warnLow = lacFineSlopeWarnLow
                        warnHigh = lacFineSlopeWarnHigh
                        errLow = lacFineSlopeErrLow
                        errHigh = lacFineSlopeErrHigh
                    else:
                        log.warning('using COARSE range for T%d,%s%s,%d,LAC', tem, 
                            calConstant.CROW[row], calConstant.CPM[end], fe)
                        warnLow = lacCoarseSlopeWarnLow
                        warnHigh = lacCoarseSlopeWarnHigh
                        errLow = lacCoarseSlopeErrLow
                        errHigh = lacCoarseSlopeErrHigh 
                        
                    if slope > warnHigh or slope < warnLow:
                    
                        if slope > errHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,LAC,%s' % \
                                (slope, errHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.error(msg)     
                            status = 1
                        elif slope < errLow:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,LAC,%s' % \
                                (slope, errLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.error(msg)     
                            status = 1
                            
                        elif slope > warnHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,LAC,%s' % \
                                (slope, warnHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.warning(msg)
                        else:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,LAC,%s' % \
                                (slope, warnLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng]) 
                        log.warning(msg)
                        
                    if offset > 0:
                        log.warning("offset %0.3f > 0 for for T%d,%s%s,%d,LAC,%s" % \
                            (offset, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe, calConstant.CDAC[rng]))
                        
    # check FLE slopes and offsets
    
    for tem in towers:
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                
                    slope = dacData[tem, row, end, fe, 2]
                    offset = dacData[tem, row, end, fe, 3]
                    rng = int(rangeData[tem, row, end, fe, 1])
                    
                    if rng == calConstant.CDAC_FINE:
                        log.warning('using FINE range for T%d,%s%s,%d,FLE', tem, 
                            calConstant.CROW[row], calConstant.CPM[end], fe)
                        warnLow = fleFineSlopeWarnLow
                        warnHigh = fleFineSlopeWarnHigh
                        errLow = fleFineSlopeErrLow
                        errHigh = fleFineSlopeErrHigh  
                    else:
                        warnLow = fleCoarseSlopeWarnLow
                        warnHigh = fleCoarseSlopeWarnHigh
                        errLow = fleCoarseSlopeErrLow
                        errHigh = fleCoarseSlopeErrHigh 
                        
                    if slope > warnHigh or slope < warnLow:
                    
                        if slope > errHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,FLE,%s' % \
                                (slope, errHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.error(msg)     
                            status = 1
                        elif slope < errLow:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,FLE,%s' % \
                                (slope, errLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.error(msg)     
                            status = 1
                            
                        elif slope > warnHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,FLE,%s' % \
                                (slope, warnHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.warning(msg)
                        else:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,FLE,%s' % \
                                (slope, warnLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng]) 
                        log.warning(msg)
                        
                    if offset > 0:
                        log.warning("offset %0.3f > 0 for for T%d,%s%s,%d,FLE,%s" % \
                            (offset, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe, calConstant.CDAC[rng]))
                        
    # check FHE slopes and offsets
    
    for tem in towers:
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                
                    slope = dacData[tem, row, end, fe, 4]
                    offset = dacData[tem, row, end, fe, 5]
                    rng = int(rangeData[tem, row, end, fe, 2])
                    
                    if rng == calConstant.CDAC_COARSE:
                        log.warning('using COARSE range for T%d,%s%s,%d,FHE', tem, 
                            calConstant.CROW[row], calConstant.CPM[end], fe)
                        warnLow = fheCoarseSlopeWarnLow
                        warnHigh = fheCoarseSlopeWarnHigh
                        errLow = fheCoarseSlopeErrLow
                        errHigh = fheCoarseSlopeErrHigh    
                    else:
                        warnLow = fheFineSlopeWarnLow
                        warnHigh = fheFineSlopeWarnHigh
                        errLow = fheFineSlopeErrLow
                        errHigh = fheFineSlopeErrHigh 
                        
                    if slope > warnHigh or slope < warnLow:
                    
                        if slope > errHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,FHE,%s' % \
                                (slope, errHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.error(msg)     
                            status = 1
                        elif slope < errLow:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,FHE,%s' % \
                                (slope, errLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.error(msg)     
                            status = 1
                            
                        elif slope > warnHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,FHE,%s' % \
                                (slope, warnHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng])
                            log.warning(msg)
                        else:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,FHE,%s' % \
                                (slope, warnLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe, calConstant.CDAC[rng]) 
                        log.warning(msg)
                        
                    if offset > 0:
                        log.warning("offset %0.3f > 0 for for T%d,%s%s,%d,FHE,%s" % \
                            (offset, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe, calConstant.CDAC[rng]))
                             
                        
    # check ULD LEX8 slopes and offsets and saturation points
    
    for tem in towers:
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                
                    slope = uldData[calConstant.CRNG_LEX8, tem, row, end, fe, 0]
                    offset = uldData[calConstant.CRNG_LEX8, tem, row, end, fe, 1]
                    sat = uldData[calConstant.CRNG_LEX8, tem, row, end, fe, 2]
                    
                    if rangeData[tem, row, end, fe, 3] == calConstant.CDAC_FINE:
                        log.error('using FINE range for T%d,%s%s,%d,ULD,LEX8', tem, row, end, fe)
                        status = 1
                        continue   
                    else:
                        warnLow = uldLex8SlopeWarnLow
                        warnHigh = uldLex8SlopeWarnHigh
                        errLow = uldLex8SlopeErrLow
                        errHigh = uldLex8SlopeErrHigh 
                        
                    if slope > warnHigh or slope < warnLow:
                    
                        if slope > errHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,ULD,LEX8' % \
                                (slope, errHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.error(msg)     
                            status = 1
                        elif slope < errLow:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,ULD,LEX8' % \
                                (slope, errLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.error(msg)     
                            status = 1
                            
                        elif slope > warnHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,ULD,LEX8' % \
                                (slope, warnHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.warning(msg)
                        else:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,ULD,LEX8' % \
                                (slope, warnLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe) 
                        log.warning(msg)
                        
                    if offset < 0:
                        log.warning('offset %0.3f < 0 for T%d,%s%s,%d,ULD,LEX8' % \
                            (offset, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe))
                             
                    if sat < 0:
                        log.warning('saturation %0.3f < 0 for T%d,%s%s,%d,ULD,LEX8' % \
                            (sat, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe))
                 
                        
    # check ULD LEX1 slopes
    
    for tem in towers:
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                
                    slope = uldData[calConstant.CRNG_LEX1, tem, row, end, fe, 0]
                    offset = uldData[calConstant.CRNG_LEX1, tem, row, end, fe, 1]
                    sat = uldData[calConstant.CRNG_LEX1, tem, row, end, fe, 2]
                    
                    if rangeData[tem, row, end, fe, 4] == 0:
                        log.error('using FINE range for T%d,%s%s,%d,ULD,LEX1', tem, row, end, fe)
                        status = 1
                        continue   
                    else:
                        warnLow = uldLex1SlopeWarnLow
                        warnHigh = uldLex1SlopeWarnHigh
                        errLow = uldLex1SlopeErrLow
                        errHigh = uldLex1SlopeErrHigh 
                        
                    if slope > warnHigh or slope < warnLow:
                    
                        if slope > errHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,ULD,LEX1' % \
                                (slope, errHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.error(msg)     
                            status = 1
                        elif slope < errLow:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,ULD,LEX1' % \
                                (slope, errLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.error(msg)     
                            status = 1
                            
                        elif slope > warnHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,ULD,LEX1' % \
                                (slope, warnHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.warning(msg)
                        else:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,ULD,LEX1' % \
                                (slope, warnLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe) 
                            log.warning(msg)
                            
                    if offset < 0:
                        log.warning('offset %0.3f < 0 for T%d,%s%s,%d,ULD,LEX1' % \
                            (offset, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe))
                             
                    if sat < 0:
                        log.warning('saturation %0.3f < 0 for T%d,%s%s,%d,ULD,LEX1' % \
                            (sat, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe))
                             
                            
    # check ULD HEX8 slopes
    
    for tem in towers:
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):
                
                    slope = uldData[calConstant.CRNG_HEX8, tem, row, end, fe, 0]
                    offset = uldData[calConstant.CRNG_HEX8, tem, row, end, fe, 1]
                    sat = uldData[calConstant.CRNG_HEX8, tem, row, end, fe, 2]
                    
                    if rangeData[tem, row, end, fe, 5] == 0:
                        log.error('using FINE range for T%d,%s%s,%d,ULD,HEX8', tem, row, end, fe)
                        status = 1
                        continue   
                    else:
                        warnLow = uldHex8SlopeWarnLow
                        warnHigh = uldHex8SlopeWarnHigh
                        errLow = uldHex8SlopeErrLow
                        errHigh = uldHex8SlopeErrHigh 
                        
                    if slope > warnHigh or slope < warnLow:
                    
                        if slope > errHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,ULD,HEX8' % \
                                (slope, errHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.error(msg)     
                            status = 1
                        elif slope < errLow:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,ULD,HEX8' % \
                                (slope, errLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.error(msg)     
                            status = 1
                            
                        elif slope > warnHigh:
                            msg = 'slope %0.3f > %0.3f for T%d,%s%s,%d,ULD,HEX8' % \
                                (slope, warnHigh, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe)
                            log.warning(msg)
                        else:
                            msg = 'slope %0.3f < %0.3f for T%d,%s%s,%d,ULD,HEX8' % \
                                (slope, warnLow, tem, calConstant.CROW[row], 
                                 calConstant.CPM[end], fe) 
                        log.warning(msg) 
                        
                    if offset < 0:
                        log.warning('offset %0.3f < 0 for T%d,%s%s,%d,ULD,HEX8' % \
                            (offset, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe))
                             
                    if sat < 0:
                        log.warning('saturation %0.3f < 0 for T%d,%s%s,%d,ULD,HEX8' % \
                            (sat, tem, calConstant.CROW[row], 
                             calConstant.CPM[end], fe))                                                                                                                                              

    return status




if __name__ == '__main__':


    rootOutput = False
    logName = None
    
    # setup logger

    logging.basicConfig()
    log = logging.getLogger('dacSlopesVal')
    log.setLevel(logging.INFO)

    # check command line

    try:
        opts = getopt.getopt(sys.argv[1:], "-R:-L:-V-r")
    except getopt.GetoptError:
        log.error(__doc__)
        sys.exit(1)

    optList = opts[0]
    for o in optList:
        if o[0] == '-R':
            rootName = o[1]
            rootOutput = True
        elif o[0] == '-r':
            rootName = None
            rootOutput = True    
        elif o[0] == '-L':
            logName = o[1]
        elif o[0] == '-V':
            log.setLevel(logging.DEBUG)    
        
    args = opts[1]
    if len(args) != 1:
        log.error(__doc__)
        sys.exit(1)    

    xmlName = args[0]
    
    ext = os.path.splitext(xmlName)
    if rootOutput and rootName is None:
        rootName = "%s.val.root" % ext[0]
    if logName is None:
        logName = "%s.val.log" % ext[0]
        
    if os.path.exists(logName):
        log.debug('Removing old log file %s', logName)
        os.remove(logName)    

    hdl = logging.FileHandler(logName)
    fmt = logging.Formatter('%(levelname)s %(message)s')
    hdl.setFormatter(fmt)
    log.addHandler(hdl)

    # open and read XML DacSlopes file

    log.info('Reading file %s', xmlName) 
    xmlFile = calCalibXML.calDacSlopesCalibXML(xmlName)
    (dacData, uldData, rangeData) = xmlFile.read()
    towers = xmlFile.getTowers()
    xmlFile.close()

    # validate calibration data

    valStatus = calcError(dacData, uldData, rangeData)    

    # create ROOT output file
    
    if rootOutput:

        import ROOT

        log.info('Creating file %s' % rootName)
        ROOT.gROOT.Reset()
        rootFile = ROOT.TFile(rootName, "recreate")

        # write error histograms

        rootHists(dacData, uldData, rangeData, xmlName)        

        # clean up

        rootFile.Close()
        
        
    # do simple stats
    
    # LAC
    
    sf = (rangeData[...,0] == 0)
    sc = (rangeData[...,0] == 1)
    
    fineCount = numarray.sum(numarray.ravel(sf))
    coarseCount = numarray.sum(numarray.ravel(sc))
    log.info("LAC FINE count = %d, COARSE count = %d", fineCount, coarseCount)  
    
    data = numarray.compress(numarray.ravel(sf), numarray.ravel(dacData[...,0]))
    if len(data) == 0:
        av = 0
        sd = 0
        mn = 0
        mx = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
        mn = min(data)
        mx = max(data)
    log.info("LAC FINE Mev/DAC slope average=%f, stddev=%f", av, sd)
    log.info("LAC FINE Mev/DAC slope minumum=%f, maximum=%f", mn, mx) 
    
    data = numarray.compress(numarray.ravel(sf), numarray.ravel(dacData[...,1]))
    if len(data) == 0:
        av = 0
        sd = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
    log.info("LAC FINE Mev offset average = %f, stddev=%f", av, sd)
    
    data = numarray.compress(numarray.ravel(sc), numarray.ravel(dacData[...,0]))
    if len(data) == 0:
        av = 0
        sd = 0
        mn = 0
        mx = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
        mn = min(data)
        mx = max(data)
    log.info("LAC COARSE Mev/DAC slope average=%f, stddev=%f", av, sd)
    log.info("LAC COARSE Mev/DAC slope minumum=%f, maximum=%f", mn, mx) 
    
    data = numarray.compress(numarray.ravel(sc), numarray.ravel(dacData[...,1]))
    if len(data) == 0:
        av = 0
        sd = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
    log.info("LAC COARSE Mev offset average = %f, stddev=%f", av, sd)
    
    # FLE    
        
    sf = (rangeData[...,1] == 0)
    sc = (rangeData[...,1] == 1)
    
    fineCount = numarray.sum(numarray.ravel(sf))
    coarseCount = numarray.sum(numarray.ravel(sc))
    log.info("FLE FINE count = %d, COARSE COUNT = %d", fineCount, coarseCount)    
    
    data = numarray.compress(numarray.ravel(sf), numarray.ravel(dacData[...,2]))  
    if len(data) == 0:
        av = 0
        sd = 0
        mn = 0
        mx = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
        mn = min(data)
        mx = max(data)
    log.info("FLE FINE Mev/DAC slope average = %f, stddev=%f", av, sd)
    log.info("FLE FINE Mev/DAC slope minumum=%f, maximum=%f", mn, mx) 
    
    data = numarray.compress(numarray.ravel(sf), numarray.ravel(dacData[...,3]))
    if len(data) == 0:
        av = 0
        sd = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
    log.info("FLE FINE Mev offset average = %f, stddev=%f", av, sd)
    
    data = numarray.compress(numarray.ravel(sc), numarray.ravel(dacData[...,2])) 
    if len(data) == 0:
        av = 0
        sd = 0
        mn = 0
        mx = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
        mn = min(data)
        mx = max(data)
    log.info("FLE COARSE Mev/DAC slope average = %f, stddev=%f", av, sd)
    log.info("FLE COARSE Mev/DAC slope minumum=%f, maximum=%f", mn, mx) 
    
    data = numarray.compress(numarray.ravel(sc), numarray.ravel(dacData[...,3]))
    if len(data) == 0:
        av = 0
        sd = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
    log.info("FLE COARSE Mev offset average = %f, stddev=%f", av, sd)
    
    # FHE
    
    sf = (rangeData[...,2] == 0)
    sc = (rangeData[...,2] == 1)
    
    fineCount = numarray.sum(numarray.ravel(sf))
    coarseCount = numarray.sum(numarray.ravel(sc))
    log.info("FHE FINE count = %d, COARSE count = %d", fineCount, coarseCount)
    
    data = numarray.compress(numarray.ravel(sf), numarray.ravel(dacData[...,4]))
    if len(data) == 0:
        av = 0
        sd = 0
        mn = 0
        mx = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
        mn = min(data)
        mx = max(data)
    log.info("FHE FINE Mev/DAC slope average = %f, stddev=%f", av, sd)
    log.info("FHE FINE Mev/DAC slope minumum=%f, maximum=%f", mn, mx) 
    
    data = numarray.compress(numarray.ravel(sf), numarray.ravel(dacData[...,5]))
    if len(data) == 0:
        av = 0
        sd = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
    log.info("FHE FINE Mev offset average = %f, stddev=%f", av, sd)
    
    data = numarray.compress(numarray.ravel(sc), numarray.ravel(dacData[...,4]))
    if len(data) == 0:
        av = 0
        sd = 0
        mn = 0
        mx = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
        mn = min(data)
        mx = max(data)
    log.info("FHE COARSE Mev/DAC slope average = %f, stddev=%f", av, sd)
    log.info("FHE COARSE Mev/DAC slope minumum=%f, maximum=%f", mn, mx) 
    
    data = numarray.compress(numarray.ravel(sc), numarray.ravel(dacData[...,5]))
    if len(data) == 0:
        av = 0
        sd = 0
    else:    
        av = numarray.average(data)
        sd = numarray.mlab.std(data)
    log.info("FHE COARSE Mev offset average = %f, stddev=%f", av, sd)
    
    # ULD LEX8
    
    data = numarray.ravel(uldData[calConstant.CRNG_LEX8,...,0])
    av = numarray.average(data)
    sd = numarray.mlab.std(data)
    mn = min(data)
    mx = max(data)
    log.info("ULD LEX8 Mev/DAC slope average = %f, stddev=%f", av, sd)
    log.info("ULD LEX8 Mev/DAC slope minumum = %f, maximum = %f", mn, mx) 
    
    av = numarray.average(uldData[calConstant.CRNG_LEX8,...,1], axis = None)
    sd = numarray.mlab.std(numarray.ravel(uldData[calConstant.CRNG_LEX8,...,1])) 
    log.info("ULD LEX8 MeV offset average = %f, stddev=%f", av, sd)
    
    av = numarray.average(uldData[calConstant.CRNG_LEX8,...,2], axis = None)
    sd = numarray.mlab.std(numarray.ravel(uldData[calConstant.CRNG_LEX8,...,2]))
    log.info("ULD LEX8 MeV saturation average = %f, stddev=%f", av, sd)
    
    mn = min(numarray.ravel(uldData[calConstant.CRNG_LEX8,...,2]))
    log.info("ULD LEX8 MeV saturation minimum = %f", mn)
    
    # ULD LEX1
    
    data = numarray.ravel(uldData[calConstant.CRNG_LEX1,...,0])
    av = numarray.average(data)
    sd = numarray.mlab.std(data)
    mn = min(data)
    mx = max(data)
    log.info("ULD LEX1 Mev/DAC slope average = %f, stddev=%f", av, sd)
    log.info("ULD LEX1 Mev/DAC slope minumum = %f, maximum = %f", mn, mx)
     
    av = numarray.average(uldData[calConstant.CRNG_LEX1,...,1], axis = None) 
    sd = numarray.mlab.std(numarray.ravel(uldData[calConstant.CRNG_LEX1,...,1]))
    log.info("ULD LEX1 MeV offset average = %f, stddev=%f", av, sd)
    
    av = numarray.average(uldData[calConstant.CRNG_LEX1,...,2], axis = None)
    sd = numarray.mlab.std(numarray.ravel(uldData[calConstant.CRNG_LEX1,...,2]))
    log.info("ULD LEX1 MeV saturation average = %f, stddev=%f", av, sd)
    
    mn = min(numarray.ravel(uldData[calConstant.CRNG_LEX1,...,2]))
    log.info("ULD LEX1 MeV saturation minimum = %f", mn)
    
    # ULD HEX8
    
    data = numarray.ravel(uldData[calConstant.CRNG_HEX8,...,0])
    av = numarray.average(data)
    sd = numarray.mlab.std(data)
    mn = min(data)
    mx = max(data)
    log.info("ULD HEX8 Mev/DAC slope average = %f, stddev=%f", av, sd) 
    log.info("ULD HEX8 Mev/DAC slope minumum = %f, maximum = %f", mn, mx)
    
    av = numarray.average(uldData[calConstant.CRNG_HEX8,...,1], axis = None) 
    sd = numarray.mlab.std(numarray.ravel(uldData[calConstant.CRNG_HEX8,...,1]))
    log.info("ULD HEX8 MeV offset average = %f, stddev=%f", av, sd)
    
    av = numarray.average(uldData[calConstant.CRNG_HEX8,...,2], axis = None)
    sd = numarray.mlab.std(numarray.ravel(uldData[calConstant.CRNG_HEX8,...,2]))
    log.info("ULD HEX8 MeV saturation average = %f, stddev=%f", av, sd)
    
    mn = min(numarray.ravel(uldData[calConstant.CRNG_HEX8,...,2]))
    log.info("ULD HEX8 MeV saturation minimum = %f", mn)

    # report results

    if valStatus == 0:
        statusStr = 'PASSED'
    else:
        statusStr = 'FAILED'

    log.info('Validation %s for file %s', statusStr, xmlName)        
        
    sys.exit(valStatus)

    
