import array

from ClusterClassifier import *
from math              import log10



class ClusterPredictor:

    def __init__(self, pdfFilePath, dataFilePath = None, cut = '1'):
        print 'Retrieving histograms for pdfs...'
        self.PdfFile = ROOT.TFile(pdfFilePath)
        self.PdfHistDict = {}
        for topology in FILE_PATH_DICT.keys():
            self.PdfHistDict[topology] = {}
            for var in VARIABLE_LIST:
                print 'Processing %s for %s' % (var, topology)
                hName = hname(var.Label, topology)
                self.PdfHistDict[topology][var.Label] = self.PdfFile.Get(hName)
        print 'Done.'
        if dataFilePath is not None:
            self.__setupDataTree()

    def __setupDataTree(self):
        print 'Creating the TTreeFormula objects...'
        self.RootTree = ROOT.TChain('MeritTuple')
        self.RootTree.Add(dataFilePath)
        self.TreeFormulaDict = {}
        self.addTreeFormula('cut', cut)
        self.addTreeFormula('CalEnergyRaw', 'CalEnergyRaw')
        for var in VARIABLE_LIST:
            self.addTreeFormula(var.Label, var.Expression)
        print 'Done.'

    def addTreeFormula(self, name, expr):
        formula = ROOT.TTreeFormula(name, expr, self.RootTree)
        self.TreeFormulaDict[name] = formula

    def getPdfHist(self, topology, var):
        return self.PdfHistDict[topology][var.Label]

    def getPdfValue(self, topology, var):
        h = self.getPdfHist(topology, var)
        try:
            logE = log10(self.getVar('CalEnergyRaw'))
        except ValueError:
            return 0.0
        binx = 1 + int(NUM_E_BINS*(float(logE - LOG_E_MIN)/\
                                       (LOG_E_MAX - LOG_E_MIN)))
        biny = 1 + int(var.NumBins*(float(self.getVar(var.Label) -\
                                          var.MinValue)/\
                                    (var.MaxValue - var.MinValue))) 
        return h.GetBinContent(binx, biny)

    def getEntry(self, entry):
        self.RootTree.GetEntry(entry)

    def getTreeFormula(self, name):
        return self.TreeFormulaDict[name]

    def getTreeFormulaValue(self, name):
        return self.getTreeFormula(name).EvalInstance()
        
    def cutPassed(self):
        return self.getTreeFormulaValue('cut')

    def getVar(self, varLabel):
        return self.getTreeFormulaValue(varLabel)

    def classifyEvent(self, entry):
        self.getEntry(entry)
        self.ProbDict = {}
        if not self.cutPassed():
            return
        for topology in FILE_PATH_DICT.keys():
            self.ProbDict[topology] = 1.0
        for topology in FILE_PATH_DICT.keys():
            for var in VARIABLE_LIST:
                self.ProbDict[topology] *= self.getPdfValue(topology, var)
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
    numEvents = 20000
    cut = PRE_CUT
    p = ClusterPredictor('cluclass.root', dataFilePath, cut)
    hDict = {}
    for topology in FILE_PATH_DICT.keys():
        hName  = 'hProb_%s' % topology
        hTitle = '%s probability' % topology
        h = ROOT.TH2F(hName, hTitle, NUM_E_BINS, LOG_E_MIN, LOG_E_MAX,
                      50, 0, 1.01)
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
    c.cd(4)
    normalizeHist(hClass)
    hClass.Draw()
    c.cd()
    c.Update()
