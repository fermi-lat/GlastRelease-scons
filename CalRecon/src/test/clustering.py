
from ReconReader import *

ROOT.gStyle.SetOptStat(111111)

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-c', '--skim-cut', type = str, dest = 'c',
                  default = '1',
                  help = 'a cut to filter the events')
parser.add_option('-s', '--save-canvas', type = str, dest = 's',
                  default = None,
                  help = 'a path to save all canvas in .png format')
(opts, args) = parser.parse_args()
if len(args) == 0:
    sys.exit('Please provide a recon input root file.')
elif len(args) > 2:
    sys.exit('Too many arguments.')
reconFilePath = args[0]
try:
    meritFilePath = args[1]
except IndexError:
    meritFilePath = None

SaveAllCanvas = False
if opts.s != None:
    SaveAllCanvas = True
    savedCanvasPath = opts.s
                                              
ANALYSIS_BIN_LIST = ['McEnergy < 100',
                     'McEnergy >= 100 && McEnergy<500',
                     'McEnergy >= 500 && McEnergy<1000',
                     'McEnergy >= 1000 && McEnergy<5000',
                     'McEnergy >= 5000 && McEnergy<20000',
                     'McEnergy >= 20000'
                     ]
ANALYSIS_BIN_LIST = ['McEnergy > 100'
                     ]

reader = ReconReader(reconFilePath, meritFilePath, None, opts.c)
numEvents = 10000

# BAD MAIN Clster Solution TTreeFourmula
BAD_ANGLE_VALUE = -0.1
UberSolutionIsBetter = ROOT.TTreeFormula('UberSolutionIsBetter',
                                          '(acos(-(Tkr1ZDir*CalUberZDir + Tkr1YDir*CalUberYDir + Tkr1XDir*CalUberXDir))-CalTrackAngle) < %f' % BAD_ANGLE_VALUE ,
                                     reader.MeritChain)


# Create TTreeFormulas for the analysis bins.
treeFormulaList = []
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    treeFormula =  ROOT.TTreeFormula('analysisBin%d' % i, cut,
                                     reader.MeritChain)
    treeFormulaList.append(treeFormula)

# Create the histograms.
hNumCluList = []
hNumIsolatedCluList = []
hDistIsolatedCluList = []
hDistSecondCluList = []
hFirstCluFracEneList = []
hSecondCluFracEneList = []
hFirstAndSecondCluFracEneList = []
hSecondCluEneList = []
hFirstCluFracXtalsList = []
hDist_vs_EnergyList = []
hCluAngle_vs_UberAngle = []
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    hTitle = cut
    hName = 'NumClu_%d' % i
    h = ROOT.TH1F(hName, hTitle, 20, 0, 20)
    h.GetXaxis().SetTitle('Number of clusters')
    hNumCluList.append(h)
    hName = 'NumIsolatedClu_%d' % i
    h = ROOT.TH1F(hName, hTitle, 20, 0, 20)
    h.GetXaxis().SetTitle('Number of isolated clusters')
    hNumIsolatedCluList.append(h)
    hName = 'DistIsolatedClu_%d' % i
    h = ROOT.TH1F(hName, hTitle, 100, 0, 500)
    h.GetXaxis().SetTitle('Distance of isolated clusters from the main axis')
    hDistIsolatedCluList.append(h)
    hName = 'DistSecondClu_%d' % i
    h = ROOT.TH1F(hName, hTitle, 100, 0, 500)
    h.GetXaxis().SetTitle('Distance of the 2nd cluster from the main axis')
    hDistSecondCluList.append(h)

    hName = 'FirstCluFracEne_%d' % i
    h = ROOT.TH1F(hName, hTitle, 100, 0.0, 1.0)
    h.GetXaxis().SetTitle('Fraction of energy in the first cluster')
    hFirstCluFracEneList.append(h)

    hName = 'SecondCluFracEne_%d' % i
    h = ROOT.TH1F(hName, hTitle, 100, 0.0, 1.0)
    h.GetXaxis().SetTitle('Fraction of energy in the second cluster')
    hSecondCluFracEneList.append(h)
    
    hName = 'FirstAndSecondCluFracEne_%d' % i
    h = ROOT.TH1F(hName, hTitle, 100, 0.0, 1.0)
    h.GetXaxis().SetTitle('Fraction of energy in the 1st+2nd (non isolated) cluster')
    hFirstAndSecondCluFracEneList.append(h)
    hName = 'SecondCluEne_%d' % i
    h = ROOT.TH1F(hName, hTitle, 100, 0.0, 100)
    h.GetXaxis().SetTitle('Energy in the second cluster')
    hSecondCluEneList.append(h)
    hName = 'FirstCluFracXtals_%d' % i
    h = ROOT.TH1F(hName, hTitle, 100, 0.0, 1.0)
    h.GetXaxis().SetTitle('Fraction of xtals in the first cluster')
    hFirstCluFracXtalsList.append(h)
    hName = 'Dist_vs_Energy_%d' % i
    h = ROOT.TH2F(hName, hTitle, 100, 0.0, 500, 100, 0, 30)
    hDist_vs_EnergyList.append(h)
    hName = 'CalAngle_vs_UberAngle_%d' % i
    h = ROOT.TH2F(hName, hTitle, 50, -3, 2, 50, -3, 2)
    h.GetXaxis().SetTitle('log10(First cluster - McDir angle)')
    h.GetYaxis().SetTitle('log10(Uber cluster - McDir angle)')
    hCluAngle_vs_UberAngle.append(h)




# Start the event loop.
for event in xrange(numEvents):
    if reader.getEntry(event):
        
        numClusters = reader.getNumClusters()
        numIsolatedClusters = 0
        if numClusters == 0:
            firstCluFracEne = -1
            secondCluFracEne = -1
            firstAndSecondCluFracEne = -1
            secondCluEne = -1
            firstCluFracXtals = -1
            distSecondClu = -1
        else:
            clusterList = reader.getCalClusterList()
            for cluster in clusterList:
                if cluster.getTotNumXtals() == 1:
                    numIsolatedClusters += 1
            uberCluster = clusterList[0]
            if numClusters == 1:
                firstCluFracEne = 2
                secondCluFracEne = -1
                firstAndSecondCluFracEne = -1
                secondCluEne = -1
                firstCluFracXtals = 2
                distSecondClu = -1
            else:
                cluster1 = clusterList[1]
                cluster2 = clusterList[2]
                firstCluFracEne = cluster1.getEnergy()/uberCluster.getEnergy()
                secondCluFracEne = cluster2.getEnergy()/uberCluster.getEnergy()
                if cluster2.getTotNumXtals() == 1:
                    firstAndSecondCluFracEne = 2
                else:
                    firstAndSecondCluFracEne = (cluster1.getEnergy() + cluster2.getEnergy())/uberCluster.getEnergy()
                secondCluEne = cluster2.getEnergy()
                firstCluFracXtals = float(cluster1.getTotNumXtals())/\
                                    uberCluster.getTotNumXtals()
                distSecondClu = cluster2.distToAxis(cluster1)

                RefXDir = reader.getMeritVariable("McXDir")
                RefYDir = reader.getMeritVariable("McYDir")
                RefZDir = reader.getMeritVariable("McZDir")
                CalXDir = reader.getMeritVariable("CalXDir")
                CalYDir = reader.getMeritVariable("CalYDir")
                CalZDir = reader.getMeritVariable("CalZDir")
                CalUberXDir = reader.getMeritVariable("CalUberXDir")
                CalUberYDir = reader.getMeritVariable("CalUberYDir")
                CalUberZDir = reader.getMeritVariable("CalUberZDir")
                CluAngle = log10(acos(-(RefZDir*CalZDir+\
                                        RefYDir*CalYDir+\
                                        RefXDir*CalXDir)))
                UberAngle = log10(acos(-(RefZDir*CalUberZDir+\
                                         RefYDir*CalUberYDir+\
                                         RefXDir*CalUberXDir)))
                
        for (i, treeFormula) in enumerate(treeFormulaList):
            if treeFormula.EvalInstance():
                hNumCluList[i].Fill(numClusters - int(numClusters > 1))
                hNumIsolatedCluList[i].Fill(numIsolatedClusters)
                hFirstCluFracEneList[i].Fill(firstCluFracEne)
                hSecondCluFracEneList[i].Fill(secondCluFracEne)

                hFirstAndSecondCluFracEneList[i].Fill(firstAndSecondCluFracEne)
                hSecondCluEneList[i].Fill(secondCluEne)
                hFirstCluFracXtalsList[i].Fill(firstCluFracXtals)
                hDistSecondCluList[i].Fill(distSecondClu)
                if numClusters > 0:
                    for cluster in clusterList:
                        if cluster.getTotNumXtals() == 1:
                            dist = cluster.distToAxis(cluster1)
                            energy = cluster.getEnergy()
                            
                            hDistIsolatedCluList[i].Fill(dist)
                            hDist_vs_EnergyList[i].Fill(dist, energy)
                            # if UberSolutionIsBetter.EvalInstance():
                            hCluAngle_vs_UberAngle[i].Fill(CluAngle, UberAngle)
                            
# And eventually draw stuff.
cFirstCluFracEne = ROOT.TCanvas('FirstCluFracEne', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cFirstCluFracEne.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cFirstCluFracEne.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hFirstCluFracEneList[i].Draw()
cFirstCluFracEne.cd()
cFirstCluFracEne.Update()
if SaveAllCanvas:
    cFirstCluFracEne.Print(os.path.join(savedCanvasPath,
                                              cFirstCluFracEne.GetName()+".png"))

cSecondCluFracEne = ROOT.TCanvas('SecondCluFracEne', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cSecondCluFracEne.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cSecondCluFracEne.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hSecondCluFracEneList[i].Draw()
cSecondCluFracEne.cd()
cSecondCluFracEne.Update()
if SaveAllCanvas:
    cSecondCluFracEne.Print(os.path.join(savedCanvasPath,
                                         cSecondCluFracEne.GetName()+".png"))

cFirstAndSecondCluFracEne = ROOT.TCanvas('FirstAndSecondCluFracEne', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cFirstAndSecondCluFracEne.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cFirstAndSecondCluFracEne.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hFirstAndSecondCluFracEneList[i].Draw()
cFirstAndSecondCluFracEne.cd()
cFirstAndSecondCluFracEne.Update()
if SaveAllCanvas:
    cFirstAndSecondCluFracEne.Print(os.path.join(savedCanvasPath,
                                                 cFirstAndSecondCluFracEne.GetName()+".png"))

cSecondCluEne = ROOT.TCanvas('SecondCluEne', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cSecondCluEne.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cSecondCluEne.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hSecondCluEneList[i].Draw()
cSecondCluEne.cd()
cSecondCluEne.Update()
if SaveAllCanvas:
    cSecondCluEne.Print(os.path.join(savedCanvasPath,
                                     cSecondCluEne.GetName()+".png"))
    
cFirstCluFracXtals = ROOT.TCanvas('FirstCluFracXtals', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cFirstCluFracXtals.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cFirstCluFracXtals.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hFirstCluFracXtalsList[i].Draw()
cFirstCluFracXtals.cd()
cFirstCluFracXtals.Update()
if SaveAllCanvas:
    cFirstCluFracXtals.Print(os.path.join(savedCanvasPath,
                                          cFirstCluFracXtals.GetName()+".png"))
    
cNumClu = ROOT.TCanvas('NumClu', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cNumClu.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cNumClu.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hNumCluList[i].Draw()
    hNumIsolatedCluList[i].SetLineColor(ROOT.kRed)
    hNumIsolatedCluList[i].Draw('sames')
cNumClu.cd()
cNumClu.Update()
if SaveAllCanvas:
    cNumClu.Print(os.path.join(savedCanvasPath,
                               cNumClu.GetName()+".png"))
    
cDistIsolatedClu = ROOT.TCanvas('DistIsolatedClu', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cDistIsolatedClu.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cDistIsolatedClu.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hDistIsolatedCluList[i].Draw()
cDistIsolatedClu.cd()
cDistIsolatedClu.Update()
if SaveAllCanvas:
    cDistIsolatedClu.Print(os.path.join(savedCanvasPath,
                                        cDistIsolatedClu.GetName()+".png"))
    
cDistSecondClu = ROOT.TCanvas('DistSecondClu', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cDistSecondClu.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cDistSecondClu.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hDistSecondCluList[i].Draw()
cDistSecondClu.cd()
cDistSecondClu.Update()
if SaveAllCanvas:
    cDistSecondClu.Print(os.path.join(savedCanvasPath,
                                      cDistSecondClu.GetName()+".png"))
    
cDist_vs_Energy = ROOT.TCanvas('Dist_vs_Energy', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cDist_vs_Energy.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cDist_vs_Energy.cd(i + 1)
    hDist_vs_EnergyList[i].Draw()
cDist_vs_Energy.cd()
cDist_vs_Energy.Update()
if SaveAllCanvas:
    cDist_vs_Energy.Print(os.path.join(savedCanvasPath,
                                       cDist_vs_Energy.GetName()+".png"))
    
#f = ROOT.TF1("f", "x", -5, 5)
cCluAngle_vs_UberAngle = ROOT.TCanvas('CluAngle_vs_UberAngle', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cCluAngle_vs_UberAngle.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cCluAngle_vs_UberAngle.cd(i + 1)
    hCluAngle_vs_UberAngle[i].Draw()
    #f.Draw("same")
cCluAngle_vs_UberAngle.cd()
cCluAngle_vs_UberAngle.Update()
if SaveAllCanvas:
    cCluAngle_vs_UberAngle.Print(os.path.join(savedCanvasPath,
                                              cCluAngle_vs_UberAngle.GetName()+".png"))
