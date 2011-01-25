import ROOT
from ClusterConfig import *
ROOT.gROOT.SetStyle('Plain')


def getYProjection(Histo2D,name,Emin,Emax):
    nYbins = Histo2D.GetNbinsY()
    nXbins = Histo2D.GetNbinsX()
  
    ymin   = Histo2D.GetYaxis().GetBinLowEdge(1)
    ymax   = Histo2D.GetYaxis().GetBinUpEdge(nYbins)

    Histo1D = ROOT.TH1F(name,name,nYbins,ymin,ymax)
 
    for binx in range(1,nXbins + 1):
        if Histo2D.GetBinCenter(binx) > Emin and Histo2D.GetBinCenter(binx)<Emax:
          
            for biny in range(1,nYbins + 1):
                content = Histo2D.GetBinContent(binx,biny)
                Histo1D.AddBinContent(biny,content)
             
    return Histo1D


def drawProbSlices(topology,ClassifierFileName):
    
    cName = 'cClass'
    cTitle = 'Probability of topology' 
    c = ROOT.TCanvas(cName, cTitle, 1000, 800)
    c.Divide(3,3)
   
  #  ClassifierFileName = 'TestNormHistoFileTrainSample_%s.root'%topology
    histfile = ROOT.TFile(ClassifierFileName)
    pool = []
    pool.append(histfile)
    for i,(Emin,Emax) in enumerate(ENERGY_BINS):
        c.cd(i+1)
        hist2d = histfile.Get("hNormProb_%s"%topology)
        hist1d = getYProjection(hist2d,"%s probability logE(%s-%s)"%(topology,Emin,Emax),Emin,Emax)
        hist1d.GetXaxis().SetTitle("%s probability"%topology)
        hist1d.GetYaxis().SetTitle("Normalized # of events")
        hist1d.GetYaxis().SetRangeUser(0,1)
        hist1d.Draw()
        pool.append(hist1d)
        ROOT.gPad.SetGridy(True)
        ROOT.gPad.SetGridx(True)
        c.Update()
    return c, pool

def draw2dHists(topology,NormClassFileName,ClassFileName):
    cName = '2dHists'
    cTitle = 'Normalized vs standard 2d hist'
    c2d = ROOT.TCanvas(cName, cTitle, 850, 400)
    pool2 = []
    c2d.Divide(2,1)
    
    NormClassFile = ROOT.TFile(NormClassFileName)
    ClassFile     = ROOT.TFile(ClassFileName)

    c2d.cd(1)
    Norm2dHist = NormClassFile.Get("hNormProb_%s"%topology)
    Norm2dHist.Draw("lego")
    pool2.append(Norm2dHist)
    pool2.append(NormClassFile)
    pool2.append(ClassFile)
    c2d.cd(2)
    Stand2dHist = ClassFile.Get("hProb_%s"%topology)
    Stand2dHist.Draw("lego")
    pool2.append(Stand2dHist)
    c2d.Update()

    return c2d,pool2
    
if __name__== '__main__':
    topology = "gam"
    NormClassFileName = 'TestNormHistoFileTrainSample_%s.root'%topology
    ClassFileName     = 'TestHistoFileTrainSample_%s.root'%topology
    (c,pool) = drawProbSlices(topology,NormClassFileName)
    (c2d,pool2) = draw2dHists(topology,NormClassFileName,ClassFileName)
    raw_input()
