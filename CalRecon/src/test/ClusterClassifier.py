import ROOT
import time
from array import array
from pXmlWriter import *

ROOT.gStyle.SetCanvasColor(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPaintTextFormat('.3f')
ROOT.gStyle.SetOptStat(111111)


LOG_E_MIN = 1
LOG_E_MAX = 6
NUM_E_BINS = 10
PRE_CUT   = 'CalNumXtals > 4'
INI_NUM_BINS = 4000

ENERGY_BINS = [(1.0,1.5),(1.5,2.0),(2.0,2.5),(2.5,3.0),(3.0,3.5),
              (3.5,4.0),(4.0,4.5),(4.5,5.0),(5.0,5.5),(5.5,6.0)]


class ClusterVariable:

    def __init__(self, expression, minValue, maxValue, numBins = 50,
                 label = None):
        self.Expression = expression
        self.MinValue = minValue
        self.MaxValue = maxValue
        self.NumBins = numBins
        self.Label = label or expression

    def __str__(self):
        return 'Variable "%s" (%.3f--%.3f in %d bins)' %\
               (self.Expression, self.MinValue, self.MaxValue, self.NumBins)
        
# There's a subtle thing in the AnalysisTuple package that need to be
# taken into account!
#
#float logRLn = 0; 
#    if ((CAL_LAT_RLn - CAL_Cntr_RLn) > 0.0) {
#        logRLn = log(CAL_LAT_RLn - CAL_Cntr_RLn);
#    }
#    if (logRLn > 0.0) {
#        CAL_Long_Rms  = sqrt(calCluster->getRmsLong()/CAL_EnergyRaw) / logRLn;
#    }
#
# where:
#
# addItem("CalLATRLn",     &CAL_LAT_RLn);
# addItem("CalCntRLn",     &CAL_Cntr_RLn);

#Remove info from tkr to get the calongrms only depending on cal. Now no longer
#needed, variable Cal1LongRms has the right info.
#CAL_LONG_EXP = 'CalLongRms*log(CalLATRLn - CalCntRLn)*((CalLATRLn - CalCntRLn) != 0)'


VARIABLE_LIST  = [ClusterVariable('Cal1TransRms', 0, 100),
                  ClusterVariable('Cal1LongRmsAsym', 0, 0.25),
                  ClusterVariable('log10(Cal1LongRms/Cal1TransRms)',
                                  0, 2.5, label = 'MomentRatio'),
                  ClusterVariable('Cal1CoreEneFrac', 0, 1),
                  ClusterVariable('Cal1XtalEneRms/CalEnergyRaw',0, 0.4,
                                  label = 'XtalEneRms'),
                  ClusterVariable('Cal1XtalEneSkewness',-2, 10),
                  ClusterVariable('Cal1dEdxAve/Cal1FullLength', 0, 100,
                                  label = 'dEdxperLength')
                #  ClusterVariable('Cal1NumXtals/log10(CalEnergyRaw)', 0, 150,
                #                  label = 'NumXtals')
                  ]
"""
VARIABLE_LIST  = [ClusterVariable('CalTransRms', 0, 100),
                  ClusterVariable('CalLRmsAsym', 0, 0.25),
                  ClusterVariable('log10(%s/CalTransRms)'%CAL_LONG_EXP,
                                  0, 2.5, label = 'MomentRatio'),
                  ClusterVariable('CalNumXtals/log10(CalEnergyRaw)', 0, 150,
                                  label = 'NumXtals')
                  ]


FILE_PATH_DICT = {'gam': '/data/Pass8/mc/allGamma-GR-v15r39-Lyon/allGamma-GR-v15r39-Lyon_merit.root',
                  'had': '/data/Pass8/mc/allPro-GR-v15r39p1-Lyon/allPro-GR-v15r39p1-Lyon_merit.root',
                  'mip': '/data/Pass8/mc/allPro-GR-v15r39p1-Lyon/allPro-GR-v15r39p1-Lyon_merit.root',
                  #'ghost': '/data/Pass8/mc/allGamma-GR-v18r4p2-PT/skimPhilippePTEvents_merit.root'
                  }
"""
mainFilePath = '/data/work/recon/v18r8p2/work/firstlook'


FILE_PATH_DICT = {'gam': '%s/all_gamma-v18r8p2-1*.root'%mainFilePath,
                  'had': '%s/CrProtonMix-v18r8p2-1*.root'%mainFilePath,
                  'mip': '%s/high_e_surface_muons-v18r8p2-1*.root'%mainFilePath}



TOPOLOGY_DICT = {'gam': 0,
                 'had': 1,
                 'mip': 2
             #    'ghost':3
                 }

PRE_CUT_DICT = {'gam': None,
                'had': 'abs(CalMIPRatio - 1) > 0.75',
                'mip': 'abs(CalMIPRatio - 1) < 0.75'
              #  'ghost': None
                }
COLORS_DICT = {'gam': ROOT.kRed,
               'had': ROOT.kBlue,
               'mip': ROOT.kBlack
              # 'ghost': ROOT.kGray + 2
               }


def hname(label, topology):
    return 'fPdf_%s_%s' % (label, topology)


def hVarname(label,topology,Ebin):
    return 'fPdf_%s_%s_varBin_%s' % (label, topology,Ebin)

def getColor(topology):
    try:
        return COLORS_DICT[topology]
    except KeyError:
        return ROOT.kBlack

def getCut(topology):
    try:
        cut = PRE_CUT_DICT[topology]
    except KeyError:
        return PRE_CUT
    if cut is None or cut.strip() == '':
        return PRE_CUT
    return '(%s) && (%s)' % (PRE_CUT, cut)

def getTopologyIndex(topology):
    try:
        return TOPOLOGY_DICT[topology]
    except KeyError:
        return -1

def normalizeHist(hist1d, miny = 0, maxy = 1):
    numBinsX = hist1d.GetNbinsX()
    sum = 0.0
    for i in xrange(1, numBinsX + 1):
        sum += hist1d.GetBinContent(i)
    for i in xrange(1, numBinsX + 1):
        try:
            value = hist1d.GetBinContent(i)/sum
            hist1d.SetBinContent(i, value)
        except ZeroDivisionError:
            pass
    hist1d.SetMinimum(miny)
    hist1d.SetMaximum(maxy)

def normalizeSlices(hist2d, minz = 1e-3, maxz = 1):
    numBinsX = hist2d.GetNbinsX()
    numBinsY = hist2d.GetNbinsY()
    for i in xrange(1, numBinsX + 1):
        sum = 0.0
        for j in xrange(0, numBinsY + 2):
            sum += hist2d.GetBinContent(i, j)
        for j in xrange(0, numBinsY + 2):
            try:
                value = hist2d.GetBinContent(i, j)/sum
                hist2d.SetBinContent(i, j, value)
            except ZeroDivisionError:
                pass
    hist2d.SetMinimum(minz)
    hist2d.SetMaximum(maxz)



OBJECT_POOL = {}

def toPool(rootObject):
    OBJECT_POOL[rootObject.GetName()] = rootObject


class ClusterClassifier:

    def __init__(self, varBins=True):
        print 'Opening files...'
        self.RootTreeDict = {}
        for (topology, filePath) in FILE_PATH_DICT.items():
            self.RootTreeDict[topology] = ROOT.TChain('MeritTuple')
            self.RootTreeDict[topology].Add(filePath)
        print 'Creating histograms for pdfs...'
        self.PdfHistDict = {}
        self.PdfHistSliceDict = {}
        self.PdfVarBinsDict = {}
        for topology in FILE_PATH_DICT.keys():
            self.PdfHistDict[topology] = {}
            self.PdfHistSliceDict[topology] = {}
            self.PdfVarBinsDict[topology] = {}
            for var in VARIABLE_LIST:
                if varBins:
                    print 'Processing %s for %s with varBins' % (var, topology)
                    self.__createHistSlice(var, topology)
                else:
                    print 'Processing %s for %s with fixedBins'% (var, topology)
                    self.__createPdfHist(var, topology)
        print 'Done.'


    def __createHistSlice(self, var, topology):
        # Create the equal sized bin histogram...
        hName = hname(var.Label, topology)
        hTitle = '%s P.D.F. (%s)' % (var.Expression, topology)
        h = ROOT.TH1F(hName, hTitle,
                      INI_NUM_BINS, var.MinValue, var.MaxValue)
        h.SetTitle(var.Expression)
        self.PdfHistDict[topology][var.Label] = h
        cut = getCut(topology)
        
        self.PdfHistSliceDict[topology][var.Label] = {}
        self.PdfVarBinsDict[topology][var.Label] = {}
        #Project info to hist for each energy range and add to dict.
        #Call variableBinning method and re-project to histo.

        for i,(Emin,Emax) in enumerate(ENERGY_BINS):
            Ecut = "&&(log10(CalEnergyRaw)>=%s&&log10(CalEnergyRaw)<=%s)"\
            %(Emin,Emax)
            fullCut = cut+Ecut
            self.RootTreeDict[topology].Project(hName, var.Expression, fullCut)

            h.SetTitle('Projection_equalBins')
            h.SetLineColor(getColor(topology))
            h.GetXaxis().SetLabelSize(0.06)
            h.GetYaxis().SetLabelSize(0.06)
            h.GetXaxis().SetTitleSize(0.06)
            h.GetXaxis().SetTitleOffset(0.80)
        #    h.Draw()
        #    ROOT.gPad.Update()
        #    raw_input()
            self.PdfHistSliceDict[topology][var.Label][i] = h
            self.getVariableBinning(var,topology,i)
            Varbinning = self.PdfVarBinsDict[topology][var.Label][i]
            print Varbinning
            numBins = len(Varbinning) - 1
            NewhName  = "%s_varBin_%s"%(hName,i)
            NewhTitle = "%s_varBin_%s"%(hTitle,i)
            PdfSlice = ROOT.TH1F(NewhName, NewhTitle,numBins,Varbinning)
            self.RootTreeDict[topology].Project(NewhName,var.Expression, fullCut)
            PdfSlice.Draw()
            ROOT.gPad.Update()
           # raw_input()
            self.PdfVarBinsDict[topology][var.Label][i] = PdfSlice



    def __createPdfHist(self, var, topology):
        # Create the two-dimensional histogram...
        hName = hname(var.Label, topology)
        hTitle = '%s P.D.F. (%s)' % (var.Expression, topology)
        h = ROOT.TH2F(hName, hTitle, NUM_E_BINS, LOG_E_MIN, LOG_E_MAX,
                      var.NumBins, var.MinValue, var.MaxValue)
        h.SetXTitle('log10(CalEnergyRaw)')
        h.SetYTitle(var.Expression)
        self.PdfHistDict[topology][var.Label] = h
        expr = '%s:log10(CalEnergyRaw)' % var.Expression
        cut = getCut(topology)
        self.RootTreeDict[topology].Project(hName, expr, cut)
        # ... then normalize the vertical slices.
        normalizeSlices(h)
        # ... finally create a TH1 for each slice.
        self.PdfHistSliceDict[topology][var.Label] = {}
        for i in range(NUM_E_BINS):
            hSlice = h.ProjectionY('%s_slice%d' % (h.GetName(), i), i+1, i+1)
            hSlice.SetTitle('P.D.F.')
            hSlice.SetLineColor(getColor(topology))
            hSlice.GetXaxis().SetLabelSize(0.06)
            hSlice.GetYaxis().SetLabelSize(0.06)
            hSlice.GetXaxis().SetTitleSize(0.06)
            hSlice.GetXaxis().SetTitleOffset(0.80)
            self.PdfHistSliceDict[topology][var.Label][i] = hSlice



    def getNumBins(self,var,topology,i):
        #Decide the number of bins to use based on the statistics in the 
        #histo. 
        histo =  self.PdfHistSliceDict[topology][var.Label][i]
        TotalEntries = histo.GetEntries()

        if TotalEntries <= var.NumBins and TotalEntries>10.0:
            Numbins = 0.25*TotalEntries
            print "****Too few entries!"
            print "Tot:%s, going to use %s bins instead of %s"%\
                (TotalEntries,Numbins,var.NumBins)
        elif TotalEntries<=10.0:
            print "** Total entries less than 10! Going to use 1 bin!"
            Numbins = 1.0
        else:
            Numbins = var.NumBins
        
        return Numbins


    def getVariableBinning(self,var,topology,i):
         binList = []
         counter = 0
         tot = 0
         histo        =  self.PdfHistSliceDict[topology][var.Label][i]
         TotalEntries = histo.GetEntries()
         Numbins      = self.getNumBins(var,topology,i)
         #If there are less than 10 entries in histo, use one single bin over
         #entire range, otherwise calculate variable bins.
         if Numbins>1:
             try:
                 NeededEntries = TotalEntries/Numbins
             except ZeroDivisionError:
                 NeededEntries = TotalEntries
             PADDING       = 0.10*NeededEntries
             binList.append(var.MinValue)
             overFlow = histo.GetBinContent(INI_NUM_BINS + 2)
             underFlow = histo.GetBinContent(0)
             totExtra = underFlow+overFlow

             for bin in range(0,INI_NUM_BINS + 2):
                 entriesPerBin = histo.GetBinContent(bin)
                 binLowEdge    = histo.GetBinLowEdge(bin)
                 binWidth      = histo.GetBinWidth(bin)
                 binHighEdge   = binLowEdge + binWidth
             
                 counter+=entriesPerBin

                 if counter>=NeededEntries - PADDING and counter!=0:
                     tot += counter  
                     print "***tot = ",tot
                     if tot!=(TotalEntries):
                        print "%s.) %s %s %d %.4f"%\
                             (bin,counter,NeededEntries,tot,binHighEdge)
                        binList.append(binHighEdge)
                        counter = 0
                     else:
                         binList.append(var.MaxValue)
                         print "Ran out of statistics (%s/%s)!"%\
                             (tot,TotalEntries) 
                         print "Widenning the bin!"
                             
                         counter = 0

             if binList[-1]<var.MaxValue:
                 print "Previous value",binList[-1]
                 binList.append(var.MaxValue)
                 print "Adding max to bin list!"

         else:
             binList = [var.MinValue,var.MaxValue]
             print "NumBins is equal to 1, taking %s as bins!"%binList
         
         NewBins = array('f',binList)
         self.PdfVarBinsDict[topology][var.Label][i] = NewBins


    def drawAllPdfHists(self):
        for var in VARIABLE_LIST:
            self.drawPdfHists(var)

    def drawPdfHists(self, var):
        cName  = '%s_2d' % var.Label
        cTitle = '%s (2d)' % var.Expression
        c = ROOT.TCanvas(cName, cTitle, 1000, 800)
        toPool(c)
        c.Divide(2, 2)
        for (i, topology) in enumerate(FILE_PATH_DICT.keys()):
            c.cd(i + 1)
            ROOT.gPad.SetRightMargin(0.15)
            self.getPdfHist(topology, var).Draw('colz,text')
            ROOT.gPad.SetLogz(True)
        c.cd()
        c.Update()
        cName = '%s_slices' % var.Label
        cTitle = '%s (slices)' % var.Expression
        c = ROOT.TCanvas(cName, cTitle, 1000, 600)
        toPool(c)
        c.Divide(4, 3)
        for i in range(NUM_E_BINS):
            legend = ROOT.TLegend(0.65, 0.67, 0.90, 0.85)
            legend.SetName('%s_legend_slice%d' % (var.Label, i))
            legend.SetFillStyle(0)
            legend.SetLineStyle(0)
            legend.SetLineWidth(0)
            legend.SetBorderSize(0)
            legend.SetTextSize(0.08)
            toPool(legend)
            c.cd(i + 1)
            ymax = 0
            for (j, topology) in enumerate(FILE_PATH_DICT.keys()):
                hSlice = self.getPdfSliceHist(topology, var, i)
                y = hSlice.GetMaximum()
                if y > ymax:
                    ymax = y
            for (j, topology) in enumerate(FILE_PATH_DICT.keys()):
                hSlice = self.getPdfSliceHist(topology, var, i)
                hSlice.SetMaximum(1.2*ymax)
                hSlice.Draw('same'*(j!=0))
                legend.AddEntry(hSlice, topology)
            legend.Draw()
            logemin = LOG_E_MIN  +\
                      i*float(LOG_E_MAX - LOG_E_MIN)/float(NUM_E_BINS)
            logemax = LOG_E_MIN  +\
                      (i + 1)*float(LOG_E_MAX - LOG_E_MIN)/float(NUM_E_BINS)
            emin = (10**logemin)
            emax = (10**logemax)
            label = ROOT.TLatex(0.15, 0.8, '%.d--%d MeV' %\
                                (emin, emax))
            label.SetName('%s_label_slice%d' % (var.Label, i))
            label.SetTextSize(0.06)
            label.SetNDC()
            toPool(label)
            label.Draw()
        c.cd()
        c.Update()

    def getPdfHist(self, topology, var):
        return self.PdfHistDict[topology][var.Label]


    def getVarBinHistSlice(self,topology,var,i):
        return self.PdfVarBinsDict[topology][var.Label][i]


    def getPdfSliceHist(self, topology, var, i):
        return self.PdfHistSliceDict[topology][var.Label][i]


    def getSliceInfo(self,histo):
        infoList = []
        numBins = histo.GetNbinsX()
        numEntries = histo.GetEntries()
        print "Number of bins:",numBins, numEntries
        sumProb = 0
        for i in xrange(1, numBins):
            binWidth    = histo.GetBinWidth(i)
            binVal      = histo.GetBinContent(i)
            binLowEdge  = histo.GetBinLowEdge(i)
            binHighEdge = binLowEdge + binWidth
            try:
                prob = binVal/(float(numEntries)*binWidth)
                sumProb += prob
            except ZeroDivisionError:
                prob =  0.0
                print "Zero Division Error!"
            info = tuple(['%.5f'%binLowEdge,'%.5f'%binHighEdge,'%.5f'%prob])
            infoList.append(info)
            histo.SetBinContent(i, prob)
            print "%s.)Bin Width:%.2f\tBinContent:%d\t Prob:%.2f"\
                %(i,binWidth,binVal,(prob))
       
        return infoList


    def writeOutputFile(self, filePath,varBins=True):
        print 'Writing output file %s...' % filePath
        outputFile = ROOT.TFile(filePath, 'RECREATE')
        for topology in FILE_PATH_DICT.keys():
            for var in VARIABLE_LIST:
                if varBins:
                    for i,(Emin,Emax) in enumerate(ENERGY_BINS):
                        self.getVarBinHistSlice(topology,var,i).Write()
                else:
                    self.getPdfHist(topology, var).Write()
        outputFile.Close()
        print 'Done.'


    def getGRversion(self):
        GRversion = mainFilePath.split('/')[4]
        GRversion = GRversion.split('/')[0]
        return GRversion


 
    def writeXmlFile(self,filepath,varBins=True):
        print 'Writing output file %s...' % filepath
        writer = pXmlWriter('%s'%filepath)
        writer.writeComment('GR-%s used for training'%self.getGRversion())
        writer.writeComment('Generated by ClusterClassifier on %s'%\
                                time.asctime())
        writer.writeComment('Precut used in training:')
        writer.indent()
        for topology in FILE_PATH_DICT.keys():
             cut = getCut(topology)
             writer.writeComment('%s : %s'% (topology,cut))
        writer.openTag('VariableBinsInfo')
        writer.newLine()
        writer.indent()
        writer.writeComment('Energy intervals for the histograms log10(MeV).')
        writer.openTag('EnergyBins')
        writer.indent()
        for i,(Emin,Emax) in enumerate(ENERGY_BINS):
             writer.writeTag('Energy',{'bin':"%s"%i,'Emin':"%s"%Emin,'Emax':"%s"%Emax})
        writer.backup()
        writer.closeTag('EnergyBins')
        writer.newLine
        writer.writeComment('Histogram info for each topology considered (gam, had) and variables in equally populated bins. xmin, xmax are the bin low edge and hig edge, and prob is the probability in that bin.')
        for topology in FILE_PATH_DICT.keys():
            writer.openTag('Topology',{'name':"%s"%topology,})
            writer.newLine()
            writer.indent()
            for var in VARIABLE_LIST:
                writer.openTag('Variable',{'name':"%s"%var.Label,})
                writer.newLine()
                writer.indent()
                if varBins:
                    for i,(Emin,Emax) in enumerate(ENERGY_BINS):
                        writer.openTag('Energy',{'bin':"%s"%i,})
                        histo = self.getVarBinHistSlice(topology,var,i)
                        infoList = self.getSliceInfo(histo)
                        writer.indent()
                        for (xmin,xmax,prob) in infoList:
                            writer.writeTag('BinValues',{'xmin':"%s"%xmin,'xmax':"%s"%xmax,'pdv':"%s"%prob})
                        writer.backup()
                        writer.closeTag('Energy')
                    writer.backup()
                writer.closeTag('Variable')
                writer.newLine()
            writer.backup()
        
            writer.closeTag('Topology')
            writer.newLine()
        writer.backup()
        writer.closeTag('VariableBinsInfo')
        writer.closeFile()

if __name__ == '__main__':
    c = ClusterClassifier()
    c.writeOutputFile('cluclassTestingCode.root')
  #  c.drawAllPdfHists()
#    c.writeOutputFile('cluclassVarBins_dEdx.root')
    #c.writeXmlFile('xml_TestCode.xml')
