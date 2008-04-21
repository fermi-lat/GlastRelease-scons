"""
Create plots of Cal intNonlin curves

inlXML2TXT <input_xml_file> <output_root_file

where:
    <input_xml_file> = input intNonlin GLAST Cal offline calibration file
    <output_root_file> = output ROOT file with plots

"""


__facility__  = "Offline"
__abstract__  = "Create plots of Cal intNonlin curves"
__author__    = "Z. Fewtrell"
__date__      = "$Date$"
__version__   = "$Revision$, $Author$"
__release__   = "$Name$"
__credits__   = "NRL code 7650"


import getopt
import sys
import calCalibXML
import ROOT
import cgc_util
import array
import numarray
import calConstant

if __name__ == '__main__':
    # check commandline
    try:
        (opts,args) = getopt.getopt(sys.argv[1:], "")
    except getopt.GetoptError:
        log.error(__doc__)
        sys.exit(1)

    if len(args) != 2:
        # should just be the one input file.
        print __doc__
        sys.exit(1)

    # retrieve commandline parms
    inName  = args[0]
    outName = args[1]

    # open and read XML IntNonlin file
    xmlFile = calCalibXML.calIntNonlinCalibXML(inName)
    (lenData, dacData, adcData) = xmlFile.read()
    towers = xmlFile.getTowers()
    xmlFile.close()

    # create output file
    rootFile = ROOT.TFile(outName, "RECREATE")

    for twr in towers:
        for offline_lyr in range(8):
            # calCalibXML uses 'row' indexing, not layer
            online_row = calCalibXML.layerToRow(offline_lyr)
            for col in range(12):
                for offline_face in range(2):
                    online_face = calConstant.offline_face_to_online[offline_face]
                    for rng in range(4):
                        nPts = lenData[rng][twr,online_row,online_face, col,0]
                        adcs = array.array('d',adcData[rng][twr,online_row,online_face, col])
                        dacs = array.array('d',dacData[rng][twr,online_row,online_face, col])
                                           
                        
                        # plot spline method
                        channel_str = "T%dL%dC%dF%dR%d"%(twr,offline_lyr,col,offline_face,rng)
                        spline = ROOT.TSpline3(channel_str, dacs, adcs, nPts)
                        
                        c = ROOT.TCanvas(channel_str, channel_str,-1)

                        spline.Draw("C")

                        g = ROOT.TGraph(nPts,
                                        dacs,
                                        adcs)

                        g.Fit("pol1","Q")
                        g.Draw("*")
                                         
                        # save plot to file
                        c.Write()

                        

    rootFile.Close()
