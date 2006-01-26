"""
Validate CAL DAC/ADC characterization calibration data in XML format.  The command
line is:

charVal [-V] [-E <err_limit>] [-W <warn_limit>] [-R <root_file>] [-L <log_file>] <xml_file>

where:

    -R <root_file> - output validation diagnostics in ROOT file
    -E <err_limit> - error limit
    -W <warn_limit> - warning limit
    -L <log_file>  - save console output to log text file
    -V             - verbose; turn on debug output
    <xml_file> The CAL characterization XML file to validate.    
"""


__facility__  = "Offline"
__abstract__  = "Validate CAL adc2nrg calibration data in XML format"
__author__    = "D.L.Wood"
__date__      = "$Date$"
__version__   = "$Revision$, $Author$"
__release__   = "$Name$"
__credits__   = "NRL code 7650"


import sys, os
import getopt
import logging

import Numeric
import mpfit

import calFitsXML
import calConstant




def residuals(p, y, x, fjac = None):
  (a, b) = p
  err = y - a*x -b
  return (0, err)



def rootHists(errData):

    # create summary histogram

    cs = ROOT.TCanvas('c_Summary', 'Summary', -1)
    cs.SetGrid()
    cs.SetLogy()

    hists = [None, None]
    leg = ROOT.TLegend(0.88, 0.88, 0.99, 0.99)

    for rng in range(2):

        hName = "h_Summary_%s" % calConstant.CDAC[rng]      
        hs = ROOT.TH1F(hName, '', 100, 0.0, (dnormErrLimit * 2))
        hs.SetLineColor(rng + 1)
        hs.SetStats(False)
        for e in errData[rng]:
            hs.Fill(e)
                        
        hists[rng] = hs
        leg.AddEntry(hs, calConstant.CDAC[rng], 'L')
        cs.Update()

    hMax = 0
    for rng in range(2):
        hs = hists[rng]
        if hs.GetMaximum() > hMax:
            hMax = hs.GetMaximum()

    for rng in range(2):
        if rng == 0:
            dopt = ''
        else:
            dopt = 'SAME'
        hs = hists[rng]
        hs.SetMaximum(hMax)
        hs.Draw(dopt)
        cs.Update()
    
    leg.Draw()
    cs.Update()

    cs.Write()


if __name__ == '__main__':

    usage = "charVal [-V] [-L <log_file>] [-E <err_limit>] [-W <warn_limit>] [-R <root_file>] <xml_file>"

    rootOutput = False
    valStatus = 0
    dnormWarnLimit = 400.0
    dnormErrLimit = 800.0

    # setup logger

    logging.basicConfig()
    log = logging.getLogger('charVal')
    log.setLevel(logging.INFO)

    # check command line

    try:
        opts = getopt.getopt(sys.argv[1:], "-R:-E:-W:-L:-V")
    except getopt.GetoptError:
        log.error(usage)
        sys.exit(1)

    optList = opts[0]
    for o in optList:
        if o[0] == '-R':
            rootName = o[1]
            rootOutput = True
        elif o[0] == '-E':
            dnormErrLimit = float(o[1])
        elif o[0] == '-W':
            dnormWarnLimit = float(o[1])
        elif o[0] == '-L':
            if os.path.exists(o[1]):
                log.warning('Deleting old log file %s', o[1])
                os.remove(o[1])
            hdl = logging.FileHandler(o[1])
            fmt = logging.Formatter('%(levelname)s %(message)s')
            hdl.setFormatter(fmt)
            log.addHandler(hdl)
        elif o[0] == '-V':
            log.setLevel(logging.DEBUG)
        
    args = opts[1]
    if len(args) != 1:
        log.error(usage)
        sys.exit(1)    

    xmlName = args[0]

    log.debug('Using input file %s', xmlName)

    # open and read XML adc2nrg file

    xmlFile = calFitsXML.calFitsXML(fileName = xmlName)
    data = xmlFile.read()
    towers = xmlFile.getTowers()
    info = xmlFile.info()
    type = info['TTYPE1']
    xmlFile.close()

    if type != 'fle_dac' and type != 'fhe_dac' and type != 'log_acpt':
        log.error('XML file type %s not supported', type)
        sys.exit(1)

    log.debug('Validating file type %s', type)

    pi = {'fixed' : 0, 'limited' : (0, 0), 'mpprint' : 0}
    pinfo = [pi, pi]
    x0 = Numeric.arange(0.0, 64.0, 1.0)

    errData = ([], [])    

    for tem in towers:
        for row in range(calConstant.NUM_ROW):
            for end in range(calConstant.NUM_END):
                for fe in range(calConstant.NUM_FE):

                    fineData = data[tem, row, end, fe, 0:64]
                    coarseData = data[tem, row, end, fe, 64:128]                    

                    # fit FINE range data                    
                    
                    z = Numeric.nonzero(fineData)
                    y = Numeric.take(fineData, z)
                    x = Numeric.take(x0, z)
                    p0 = (20.0, -200.0)
                    fkw = {'x': x, 'y' : y}

                    if len(x) < 3:
                        log.error('Too little data: %d,%s%s,%d,FINE', tem, calConstant.CROW[row],
                                      calConstant.CPM[end], fe)
                        valStatus = 1
                    else:

                        fit = mpfit.mpfit(residuals, p0, functkw = fkw, parinfo = pinfo, quiet = 1)
                        if fit.status <= 0:
                            log.error('mpfit error - %s', fit.errmsg)
                            sys.exit(1)
                        dnorm = (fit.fnorm / len(x))
                        errData[0].append(dnorm)
                        log.debug("%d,%s%s,%d,FINE: %0.1f %0.1f %0.2f", tem, calConstant.CROW[row],
                                calConstant.CPM[end], fe, fit.params[0], fit.params[1], dnorm)

                        if dnorm > dnormWarnLimit:
                            if dnorm > dnormErrLimit:
                                log.error('dnorm %0.2f > %0.2f for %d,%s%s,%d,FINE', dnorm, dnormErrLimit, tem,
                                          calConstant.CROW[row], calConstant.CPM[end], fe)
                                valStatus = 1
                            else:
                                log.warning('dnorm %0.2f > %0.2f for %d,%s%s,%d,FINE', dnorm, dnormWarnLimit, tem,
                                          calConstant.CROW[row], calConstant.CPM[end], fe)

                    # fit coarse range data

                    z = Numeric.nonzero(coarseData)
                    y = Numeric.take(coarseData, z)
                    x = Numeric.take(x0, z)
                    p0 = (40, -400)
                    fkw = {'x': x, 'y' : y}
                    
                    if len(x) < 3:
                        log.error('Too little data: %d,%s%s,%d,COARSE', tem, calConstant.CROW[row],
                                      calConstant.CPM[end], fe)
                    else:

                        fit = mpfit.mpfit(residuals, p0, functkw = fkw, parinfo = pinfo, quiet = 1)
                        if fit.status <= 0:
                            log.error('mpfit error - %s', fit.errmsg)
                            sys.exit(1)
                        dnorm = (fit.fnorm / len(x))
                        errData[1].append(dnorm)
                        log.debug("%d,%s%s,%d,COARSE: %0.1f %0.1f %0.2f", tem, calConstant.CROW[row],
                                calConstant.CPM[end], fe, fit.params[0], fit.params[1], dnorm)

                        if dnorm > dnormWarnLimit:
                            if dnorm > dnormErrLimit:
                                log.error('dnorm %0.2f > %0.2f for %d,%s%s,%d,COARSE', dnorm, dnormErrLimit, tem,
                                          calConstant.CROW[row], calConstant.CPM[end], fe)
                                valStatus = 1
                            else:
                                log.warning('dnorm %0.2f > %0.2f for %d,%s%s,%d,COARSE', dnorm, dnormWarnLimit, tem,
                                          calConstant.CROW[row], calConstant.CPM[end], fe)                        


    # create ROOT output file
    
    if rootOutput:

        import ROOT

        ROOT.gROOT.Reset()
        rootFile = ROOT.TFile(rootName, "recreate")

        # write error histograms

        rootHists(errData)        

        # clean up

        rootFile.Close()
        

    # report results

    if valStatus == 0:
        statusStr = 'PASSED'
    else:
        statusStr = 'FAILED'

    log.info('Validation %s for file %s', statusStr, xmlName)
    sys.exit(valStatus)
    
    