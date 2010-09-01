import ROOT

ROOT.gStyle.SetCanvasColor(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPaintTextFormat('.3f')
ROOT.gStyle.SetOptStat(0)


LOG_E_MIN = 1
LOG_E_MAX = 6
NUM_E_BINS = 10
PRE_CUT   = 'CalNumXtals >= 4'


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


CAL_LONG_EXP = 'CalLongRms*log(CalLATRLn - CalCntRLn)*((CalLATRLn - CalCntRLn) != 0)'

VARIABLE_LIST  = [ClusterVariable('CalTransRms', 0, 100),
                  ClusterVariable('CalLRmsAsym', 0, 0.25),
                  ClusterVariable('log10(%s/CalTransRms)' % CAL_LONG_EXP,
                                  0, 2.5, label = 'MomentRatio'),
                  #ClusterVariable('log10(CalLongRms/CalTransRms)',
                  #                0, 2.5, label = 'MomentRatio2'),
                  ClusterVariable('CalNumXtals/log10(CalEnergyRaw)', 0, 150,
                                  label = 'NumXtals')
                  ]
FILE_PATH_DICT = {'gam': '/data/mc/allGamma-GR-v15r39-Lyon/allGamma-GR-v15r39-Lyon_merit.root',
                  'had': '/data/mc/allPro-GR-v15r39p1-Lyon/allPro-GR-v15r39p1-Lyon_merit.root',
                  'mip': '/data/mc/allMuon-GR-v15r35/allMuon-GR-v15r35_merit.root'
                  }
TOPOLOGY_DICT = {'gam': 0,
                 'had': 1,
                 'mip': 2
                 }
PRE_CUT_DICT = {'gam': None,
                'had': 'abs(CalMIPRatio - 1) > 0.75',
                'mip': None
                }
COLORS_DICT = {'gam': ROOT.kRed,
               'had': ROOT.kBlue,
               'mip': ROOT.kBlack
               }


def hname(label, topology):
    return 'fPdf_%s_%s' % (label, topology)

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

    def __init__(self):
        print 'Opening files...'
        self.RootTreeDict = {}
        for (topology, filePath) in FILE_PATH_DICT.items():
            self.RootTreeDict[topology] = ROOT.TChain('MeritTuple')
            self.RootTreeDict[topology].Add(filePath)
        print 'Creating histograms for pdfs...'
        self.PdfHistDict = {}
        self.PdfHistSliceDict = {}
        for topology in FILE_PATH_DICT.keys():
            self.PdfHistDict[topology] = {}
            self.PdfHistSliceDict[topology] = {}
            for var in VARIABLE_LIST:
                print 'Processing %s for %s' % (var, topology)
                self.__createPdfHist(var, topology)
        print 'Done.'

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

    def getPdfSliceHist(self, topology, var, i):
        return self.PdfHistSliceDict[topology][var.Label][i]

    def writeOutputFile(self, filePath):
        print 'Writing output file %s...' % filePath
        outputFile = ROOT.TFile(filePath, 'RECREATE')
        for topology in FILE_PATH_DICT.keys():
            for var in VARIABLE_LIST:
                self.getPdfHist(topology, var).Write()
        outputFile.Close()
        print 'Done.'





if __name__ == '__main__':
    c = ClusterClassifier()
    c.drawAllPdfHists()
    c.writeOutputFile('cluclass.root')
