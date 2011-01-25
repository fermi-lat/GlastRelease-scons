import ROOT
import time
from array         import array
from pXmlWriter    import *
from ClusterConfig import *

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
        GRversion = MAIN_FILE_PATH.split('/')[4]
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
    c.writeXmlFile('xml_TestCode.xml')
