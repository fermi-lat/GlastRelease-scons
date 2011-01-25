import array

from ClusterConfig   import *
from math            import log10



class ClusterPredictor:

    def __init__(self, pdfFilePath, dataFilePath = None, cut = '1', varBins=True):
        print 'Retrieving histograms for pdfs...'
        self.PdfFile = ROOT.TFile(pdfFilePath)
        self.DataFilePath = dataFilePath
        self.Cut = cut
        self.PdfHistDict = {}
        for topology in FILE_PATH_DICT.keys():
            self.PdfHistDict[topology] = {}
            print self.PdfHistDict[topology]
            for var in VARIABLE_LIST:
                if varBins:
                    self.PdfHistDict[topology][var.Label] = {}
                    for i,(Emin,Emax) in enumerate(ENERGY_BINS):
                        print 'Processing %s for %s using variable binning.'\
                            % (var.Label, topology)
                        hName = hVarname(var.Label, topology,i)
                        print "Histo name is",hName
                        print self.PdfFile.Get(hName)
                        self.PdfHistDict[topology][var.Label][i] = \
                            self.PdfFile.Get(hName)


                else:
                    print 'Processing %s for %s' % (var, topology)
                    hName = hname(var.Label, topology)
                    self.PdfHistDict[topology][var.Label] = \
                        self.PdfFile.Get(hName)
      
        print 'Done.'
        if self.DataFilePath is not None:
            self.__setupDataTree()

    def __setupDataTree(self):
        print 'Creating the TTreeFormula objects...'
        self.RootTree = ROOT.TChain('MeritTuple')
        self.RootTree.Add(self.DataFilePath)
        print "Added %s files to tree" %self.DataFilePath
        self.numEvents = self.RootTree.GetEntries()
        self.TreeFormulaDict = {}
        self.addTreeFormula('cut', self.Cut)
        self.addTreeFormula('CalEnergyRaw', 'CalEnergyRaw')
        for var in VARIABLE_LIST:
            self.addTreeFormula(var.Label, var.Expression)
        print 'Done.'


    def getNumEntries(self):
        return self.numEvents

    def addTreeFormula(self, name, expr):
        formula = ROOT.TTreeFormula(name, expr, self.RootTree)
        self.TreeFormulaDict[name] = formula

    def getPdfHist(self, topology, var):
        return self.PdfHistDict[topology][var.Label]

    def getHistSlice(self,topology,var,Energybin):
        return self.PdfHistDict[topology][var.Label][Energybin]

    def getPdfValue(self, topology, var):
        h = self.getPdfHist(topology, var)
        try:
            logE = log10(self.getVar('CalEnergyRaw'))
        except ValueError:
            return 0.0
        binx = 1 + int(NUM_E_BINS*(float(logE - LOG_E_MIN)/\
                                       (LOG_E_MAX - LOG_E_MIN)))
        value = self.getVar(var.Label)
        if value < var.MinValue:
            biny = 0
        elif value >= var.MaxValue:
            biny = var.NumBins + 1
        else:
            biny = 1 + int(var.NumBins*(float(value - var.MinValue)/\
                                        (var.MaxValue - var.MinValue))) 
        pdfValue = h.GetBinContent(binx, biny)
      #  print 'PDF value for %s (logE = %.3f, %s = %.3f) = %.3f' %\
      #        (topology, logE, var.Label, value, pdfValue)
        return pdfValue

    def getPdfVarBinValue(self,topology,var):
        #For each event you need the topology, variable and energy to get
        #the prob. The energy is needed in order to select the right
        #slice from the root file.
        try:
            logE = log10(self.getVar('CalEnergyRaw'))
        except ValueError:
            return 0.0

        Energybin =  int(NUM_E_BINS*(float(logE - LOG_E_MIN)/\
                                            (LOG_E_MAX - LOG_E_MIN)))


        h = self.getHistSlice(topology, var,Energybin)
        numEntries = h.GetEntries()
       
        value = self.getVar(var.Label)
        bin = self.getBin(h,value,var)
        binVal = h.GetBinContent(bin)
        binWidth = h.GetBinWidth(bin)
        try:
            pdfValue = binVal/(float(numEntries)*binWidth)
        except ZeroDivisionError:
            pdfValue =  0.0
           # print "Zero Division Error!"
           # print 'PDF value for %s (logE = %.3f, %s = %.3f) = %.3f' %\
           #   (topology, logE, var.Label, value, pdfValue)
       
        return pdfValue
        

    def getBin(self,histogram,value,var):
        numBins = histogram.GetNbinsX()
        if value < var.MinValue:
            bin = 0
        elif value >= var.MaxValue:
            bin = numBins + 1
        else:
            for i in range(1, numBins + 1):
                lowEdge = histogram.GetBinLowEdge(i)
                highEdge = lowEdge +  histogram.GetBinWidth(i)
                if value>=lowEdge and value<=highEdge:
                    bin = i
        return bin

    def getEntry(self, entry):
        self.RootTree.GetEntry(entry)

    def getTreeFormula(self, name):
        return self.TreeFormulaDict[name]

    def getTreeFormulaValue(self, name):
        return self.getTreeFormula(name).EvalInstance()
        
    def cutPassed(self):
        self.getTreeFormula('cut').UpdateFormulaLeaves()
        return  self.getTreeFormulaValue('cut')


    def getVar(self, varLabel):
        self.getTreeFormula(varLabel).UpdateFormulaLeaves()
        return self.getTreeFormulaValue(varLabel)

    def classifyEvent(self, entry = None,varBins = True):
        if entry is not None:
            self.getEntry(entry)
        self.ProbDict = {}
        if not self.cutPassed():
            return
        for topology in FILE_PATH_DICT.keys():
            self.ProbDict[topology] = 1.0
        for topology in FILE_PATH_DICT.keys():
            for var in VARIABLE_LIST:
                if varBins:
                    pdfValue = self.getPdfVarBinValue(topology, var)
                    
                else:
                    pdfValue = self.getPdfValue(topology, var)
                
                self.ProbDict[topology] *= pdfValue
        probSum = sum(self.ProbDict.values())
        for topology in FILE_PATH_DICT.keys():
            try:
                self.ProbDict[topology] /= probSum
            except ZeroDivisionError:
                pass
        maxProb = 0
        classTopology = None
        for topology in FILE_PATH_DICT.keys():
            if self.ProbDict[topology] > maxProb:
                maxProb = self.ProbDict[topology]
                classTopology = topology
        return getTopologyIndex(classTopology)



if __name__ == '__main__':
    dataFilePath = FILE_PATH_DICT['mip']
    cut = PRE_CUT
    p = ClusterPredictor('cluclassTestingCode.root', dataFilePath, cut)

    hDict = {}
    for topology in FILE_PATH_DICT.keys():
        hName  = 'hProb_%s' % topology
        hTitle = '%s probability' % topology
        h = ROOT.TH2F(hName, hTitle, NUM_E_BINS, LOG_E_MIN, LOG_E_MAX,
                      55, 0, 1.01)
        h.SetXTitle('log10(CalEnergyRaw)')
        h.GetXaxis().SetTitleOffset(1.55)
        h.SetYTitle(hTitle)
        h.GetYaxis().SetTitleOffset(1.55)
        h.SetZTitle('Normalized # events')
        h.GetZaxis().SetTitleOffset(1.25)
        hDict[topology] = h

    hClass = ROOT.TH1F('hClass', 'hClass', 3, 0, 3)
    hClass.SetXTitle('Predicted topology')
    hClass.SetYTitle('Normalized # events')
    for (key, value) in TOPOLOGY_DICT.items():
        hClass.GetXaxis().SetBinLabel(value + 1, key)

    #numEvents = p.getNumEntries()
    numEvents = 80000
    histoList = []    
    for i in xrange(numEvents):
        classTopologyIndex = p.classifyEvent(i)
        if p.cutPassed():
            hClass.Fill(classTopologyIndex)
            logE = log10(p.getVar('CalEnergyRaw'))
            for topology in FILE_PATH_DICT.keys():
                hDict[topology].Fill(logE, p.ProbDict[topology])
    for topology in FILE_PATH_DICT.keys():
        normalizeSlices(hDict[topology])
    cName = 'cClass'
    cTitle = 'Classification of test dataset' 
    c = ROOT.TCanvas(cName, cTitle, 1000, 800)
    c.Divide(2, 2)
   
    for (i, topology) in enumerate(FILE_PATH_DICT.keys()):
        c.cd(i + 1)
        hDict[topology].Draw('lego')
        histoList.append(hDict[topology])
    c.cd(4)
    saveToFile(histoList,"TestHistoFile.root")
    normalizeHist(hClass)
    hClass.Draw()
    c.cd()
    c.Update()
 #   c.SaveAs('TestVarBinProbs_3Classes.pdf')
