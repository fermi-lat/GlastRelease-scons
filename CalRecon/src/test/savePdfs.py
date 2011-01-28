import ROOT
ROOT.gROOT.SetStyle('Plain')

from ClusterConfig     import *


class ClusterPdfs:
    
    def __init__(self,pdfFilePath):
        self.PdfFile = ROOT.TFile(pdfFilePath)        
        self.PdfHistDict = {}
        self.PdfESliceHistDict = {}
        
        for topology in CLASS_FILE_PATH_DICT.keys():
            self.PdfHistDict[topology]       = {}
            self.PdfESliceHistDict[topology] = {}
             
            for var in VARIABLE_LIST:
                self.PdfHistDict[topology][var.Label]       = {}
                self.PdfESliceHistDict[topology][var.Label] = {}
                
                for i,(Emin,Emax) in enumerate(ENERGY_BINS):
                    print 'Processing %s for %s using variable binning.'\
                          % (var.Label, topology)
                    hName = hVarname(var.Label, topology,i)
                    print "Histo name is",hName
                    print self.PdfFile.Get(hName)
                    self.PdfHistDict[topology][var.Label][i] = \
                                                      self.PdfFile.Get(hName)





    def getHistSlice(self,topology,var,Energybin):
        return self.PdfHistDict[topology][var.Label][Energybin]


    def getpdfHist(self,topology,var,i):
        return self.PdfESliceHistDict[topology][var.Label][i]



    def getPdf(self,topology,var,Energybin):
   
   
        h = self.getHistSlice(topology, var,Energybin)
        numEntries = h.GetEntries()
        numBins = h.GetNbinsX()
        
        hTitle = '%s P.D.F.Distribution (%s)' % (var.Expression, topology)
        hName = hname("%s_PDF_%s"%(var.Label,Energybin), topology)
        pdfh = ROOT.TH1F()
        h.Copy(pdfh)
        pdfh.SetName(hName)
       

        for bin in range(1,numBins +1):
        
            binVal = h.GetBinContent(bin)
            binWidth = h.GetBinWidth(bin)

            try:
                pdfValue = binVal/(float(numEntries)*binWidth)
            except ZeroDivisionError:
                pdfValue =  0.0

            pdfh.SetBinContent(bin,pdfValue)
            
       
        self.PdfESliceHistDict[topology][var.Label][Energybin] = pdfh
       



    def writeOutputFile(self, filePath):
        print 'Writing output file %s...' % filePath
        outputFile = ROOT.TFile(filePath, 'RECREATE')
        for topology in CLASS_FILE_PATH_DICT.keys():
            for var in VARIABLE_LIST:
                for i,(Emin,Emax) in enumerate(ENERGY_BINS):
                    self.getpdfHist(topology,var,i).Write()

        outputFile.Close()
        print 'Done.'



        
if __name__ == '__main__':
    c = ClusterPdfs('cluclass_MomNumXtalsdEdx_largeStat.root')
    for energyBin, (emin,emax) in enumerate(ENERGY_BINS):
        for topology in CLASS_FILE_PATH_DICT.keys():
            for var in VARIABLE_LIST:
                c.getPdf(topology,var,energyBin)
    c.writeOutputFile("cluclass_pdfHists.root")
