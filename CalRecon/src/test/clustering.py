
from ReconReader import *

ROOT.gStyle.SetOptStat(111111)

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-c', '--skim-cut', type = str, dest = 'c',
                  default = '1',
                  help = 'a cut to filter the events')
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


ANALYSIS_BIN_LIST = ['McEnergy < 100',
                     'McEnergy >= 100 && McEnergy<500',
                     'McEnergy >= 500 && McEnergy<1000',
                     'McEnergy >= 1000 && McEnergy<5000',
                     'McEnergy >= 5000 && McEnergy<20000',
                     'McEnergy >= 20000'
                     ]


reader = ReconReader(reconFilePath, meritFilePath, None, opts.c)
numEvents = 1000

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
hSecondCluEneList = []
hFirstCluFracXtalsList = []
hDist_vs_EnergyList = []
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

# Start the event loop.
for event in xrange(numEvents):
    if reader.getEntry(event):
        numClusters = reader.getNumClusters()
        numIsolatedClusters = 0
        if numClusters == 0:
            firstCluFracEne = -1
            secondCluFracEne = -1
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
                secondCluEne = -1
                firstCluFracXtals = 2
                distSecondClu = -1
            else:
                cluster1 = clusterList[1]
                cluster2 = clusterList[2]
                firstCluFracEne = cluster1.getEnergy()/uberCluster.getEnergy()
                secondCluFracEne = cluster2.getEnergy()/uberCluster.getEnergy()
                secondCluEne = cluster2.getEnergy()
                firstCluFracXtals = float(cluster1.getTotNumXtals())/\
                                    uberCluster.getTotNumXtals()
                distSecondClu = cluster2.distToAxis(cluster1)
        for (i, treeFormula) in enumerate(treeFormulaList):
            if treeFormula.EvalInstance():
                hNumCluList[i].Fill(numClusters - int(numClusters > 1))
                hNumIsolatedCluList[i].Fill(numIsolatedClusters)
                hFirstCluFracEneList[i].Fill(firstCluFracEne)
                hSecondCluFracEneList[i].Fill(secondCluFracEne)
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

cSecondCluFracEne = ROOT.TCanvas('SecondCluFracEne', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cSecondCluFracEne.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cSecondCluFracEne.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hSecondCluFracEneList[i].Draw()
cSecondCluFracEne.cd()
cSecondCluFracEne.Update()

cSecondCluEne = ROOT.TCanvas('SecondCluEne', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cSecondCluEne.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cSecondCluEne.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hSecondCluEneList[i].Draw()
cSecondCluEne.cd()
cSecondCluEne.Update()

cFirstCluFracXtals = ROOT.TCanvas('FirstCluFracXtals', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cFirstCluFracXtals.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cFirstCluFracXtals.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hFirstCluFracXtalsList[i].Draw()
cFirstCluFracXtals.cd()
cFirstCluFracXtals.Update()

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

cDistIsolatedClu = ROOT.TCanvas('DistIsolatedClu', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cDistIsolatedClu.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cDistIsolatedClu.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hDistIsolatedCluList[i].Draw()
cDistIsolatedClu.cd()
cDistIsolatedClu.Update()

cDistSecondClu = ROOT.TCanvas('DistSecondClu', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cDistSecondClu.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cDistSecondClu.cd(i + 1)
    ROOT.gPad.SetLogy(True)
    hDistSecondCluList[i].Draw()
cDistSecondClu.cd()
cDistSecondClu.Update()

cDist_vs_Energy = ROOT.TCanvas('Dist_vs_Energy', '', 1100, 600)
ROOT.gPad.SetTitle(ROOT.gPad.GetName())
cDist_vs_Energy.Divide(3, 2)
for (i, cut) in enumerate(ANALYSIS_BIN_LIST):
    cDist_vs_Energy.cd(i + 1)
    hDist_vs_EnergyList[i].Draw()
cDist_vs_Energy.cd()
cDist_vs_Energy.Update()
