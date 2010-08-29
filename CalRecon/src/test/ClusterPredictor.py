import array

from ClusterClassifier import *
from math              import log10



class ClusterPredictor:

    def __init__(self, pdfFilePath, dataFilePath, cut = '1'):
        self.PdfFile = ROOT.TFile(pdfFilePath)
        self.RootTree = ROOT.TChain('MeritTuple')
        self.RootTree.Add(dataFilePath)
        print 'Retrieving histograms for pdfs...'
        self.PdfHistDict = {}
        for topology in FILE_PATH_DICT.keys():
            self.PdfHistDict[topology] = {}
            for var in VARIABLE_LIST:
                print 'Processing %s for %s' % (var, topology)
                hName = hname(var.Label, topology)
                self.PdfHistDict[topology][var.Label] = self.PdfFile.Get(hName)
        print 'Done.'
        print 'Creating the TTreeFormula objects...'
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




if __name__ == '__main__':
    dataFilePath = FILE_PATH_DICT['had']
    cut = PRE_CUT
    p = ClusterPredictor('cluclass.root', dataFilePath, cut)
    h = ROOT.TH2F('pgamma', 'gamma prob.', NUM_E_BINS, LOG_E_MIN, LOG_E_MAX,
                  50, 0, 1.01)
    for i in xrange(10000):
        p.classifyEvent(i)
        if p.cutPassed():
            h.Fill(log10(p.getVar('CalEnergyRaw')), p.ProbDict['gamma'])
    h.Draw('lego')
    ROOT.gPad.Update()
