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


def drawProbSlices(topology,ClassifierFileName,norm=True):
    
    cName = 'cClass'
    cTitle = 'Probability of topology' 
    c = ROOT.TCanvas(cName, cTitle, 1000, 800)
    c.Divide(3,3)
   
    histfile = ROOT.TFile(ClassifierFileName)
    pool = []
    pool.append(histfile)
    for i,(Emin,Emax) in enumerate(ENERGY_BINS):
        c.cd(i+1)
        if norm:
            hist2d = histfile.Get("hNormProb_%s"%topology)
        else:
            hist2d = histfile.Get("hProb_%s"%topology)
        hist1d = getYProjection(hist2d,"%s probability logE(%s-%s)"%(topology,Emin,Emax),Emin,Emax)
        hist1d.GetXaxis().SetTitle("%s probability"%topology)
        if norm:
            hist1d.GetYaxis().SetTitle("Normalized # of events")
            hist1d.GetYaxis().SetRangeUser(0,1)
        else:
            hist1d.GetYaxis().SetTitle("# of events")
        
        hist1d.Draw()
        hist1d.SetLineWidth(2)
        hist1d.SetFillColor(COLORS_DICT[topology])
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
   # Norm2dHist.SetLineColor(ROOT.kRed)
    Norm2dHist.SetTitle("Normalized %s probability"%topology)
    pool2.append(Norm2dHist)
    pool2.append(NormClassFile)
    pool2.append(ClassFile)
    c2d.cd(2)
    Stand2dHist = ClassFile.Get("hProb_%s"%topology)
    Stand2dHist.Draw("lego")
    #Stand2dHist.SetLineColor(ROOT.kRed)
    pool2.append(Stand2dHist)
    c2d.Update()

    return c2d,pool2


def drawPanel(topologyList,FileName):
    cName = 'PanelHists'
    cTitle = 'Classification'
    c3d = ROOT.TCanvas(cName, cTitle, 1000, 800)
    pool3 = []
    c3d.Divide(2,3)
    rootFile = ROOT.TFile(FileName)
    for i,topology in enumerate(topologyList):
        c3d.cd(i+1)
        hist = rootFile.Get("hNormProb_%s"%topology)
        hist.Draw("lego")
        pool3.append(hist)
    c3d.cd(5)
    hclass = rootFile.Get("hClass")
    hclass.Draw()
    pool3.append(hclass)
    pool3.append(rootFile)
    c3d.Update()
    return c3d,pool3



if __name__== '__main__':
    #BaseFileName = 'cluclassTestBinning_3Classes'
    #BaseFileName = 'cluclassTestBinning_AllVars_Mst_ClassData'
    BaseFileName = 'cluclassTestBinning_AllVars_Mst'
    #BaseFileName = 'cluclassTestBinning'
   # BaseFileName = 'cluclassTestBinning_MstVar'

    topology = "gam"
    topologyList = ['mip','gam','had','ghost']
    
    NormClassFileName = '%s_NormHisto_%s.root'%(BaseFileName,topology)
    ClassFileName     = '%s_Histo_%s.root'%(BaseFileName,topology)
    
   # NormClassFileName = 'TestNormHistoFile3Classes_largeStat_%s.root'%topology
   # ClassFileName     = 'TestHistoFile3Classes_largeStat_%s.root'%topology

   
    (c,pool)    = drawProbSlices(topology,NormClassFileName)
    (c2d,pool2) = draw2dHists(topology,NormClassFileName,ClassFileName)
    (c3d,pool3) = drawPanel(topologyList,NormClassFileName)
    raw_input()


   
