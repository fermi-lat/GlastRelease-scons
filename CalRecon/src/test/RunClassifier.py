import ROOT
ROOT.gROOT.SetStyle('Plain')

#Main configuration is specified in the ClusterConfig file.

from ClusterConfig     import *
from ClusterClassifier import *
from ClusterPredictor  import *


#Test script to run the classifier and the predictor for a given topology.
#Write the output to a root file with the 2d histograms of the classified
#topologies as well as the xml file.

#BaseName = 'cluclassTestBinning_AllVars_Mst_ClassData'
BaseName = 'cluclass_testBatchMode'

ClassifierFileName = '%s.root'%BaseName

#Specify the name of the output root and xml files here.

#classifier = ClusterClassifier(True,False)
#classifier.writeOutputFile(ClassifierFileName)
#classifier.writeXmlFile('xml_MomNumXtalsdEdx_largeStat.xml')



#Specify the topology you want to classify here.
ClassType = "mip"
dataFilePathList = CLASS_FILE_PATH_DICT[ClassType]
cut = PRE_CUT

predictor = ClusterPredictor(ClassifierFileName,dataFilePathList,cut)

hDict = {}
hNormDict = {}

def setUp2dHist(norm=True):
    if norm:
        hName  = 'hNormProb_%s' % topology
    else:
        hName  = 'hProb_%s' % topology
    hTitle = '%s probability' % topology
    h = ROOT.TH2F(hName, hTitle, NUM_E_BINS, LOG_E_MIN, LOG_E_MAX,
                  55, 0, 1.01)
    h.SetXTitle('log10(CalEnergyRaw)')
    h.GetXaxis().SetTitleOffset(1.55)
    h.SetYTitle(hTitle)
    h.GetYaxis().SetTitleOffset(1.55)
    if norm:
        h.SetZTitle('Normalized # events')
    else:
        h.SetZTitle('# events')
    h.GetZaxis().SetTitleOffset(1.25)
    return h



for topology in CLASS_FILE_PATH_DICT.keys():
    hNorm = setUp2dHist()

    hNormDict[topology] = hNorm
    h     = setUp2dHist(False)
    hDict[topology] = h

hClass = ROOT.TH1F('hClass', 'hClass', 4, 0, 4)
hClass.SetXTitle('Predicted topology')
hClass.SetYTitle('Normalized # events')
for (key, value) in TOPOLOGY_DICT.items():
    hClass.GetXaxis().SetBinLabel(value + 1, key)


#Pick whether you want to classify using the full training sample or just a
#small number of events.

numEvents = predictor.getNumEntries()
#numEvents = 8000
histoList = []    
NormhistoList = []    
for i in xrange(numEvents):
    classTopologyIndex = predictor.classifyEvent(i)
    if predictor.cutPassed():
        hClass.Fill(classTopologyIndex)
        logE = log10(predictor.getVar('CalEnergyRaw'))
        for topology in CLASS_FILE_PATH_DICT.keys():
            hDict[topology].Fill(logE, predictor.ProbDict[topology])
            hNormDict[topology].Fill(logE, predictor.ProbDict[topology])
            
            
for topology in CLASS_FILE_PATH_DICT.keys():
    normalizeSlices(hNormDict[topology])

    histoList.append(hDict[topology])
    NormhistoList.append(hNormDict[topology])


cName = 'cClass'
cTitle = 'Classification of test dataset' 
c = ROOT.TCanvas(cName, cTitle, 1000, 800)
c.Divide(2, 3)
   
for (i, topology) in enumerate(CLASS_FILE_PATH_DICT.keys()):
    c.cd(i + 1)
    hNormDict[topology].Draw('lego')
   
c.cd(5)
normalizeHist(hClass)
hClass.Draw()
c.cd()
c.Update()
NormhistoList.append(hClass)
#Save the canvas as pdf.
c.SaveAs('%s_%s.pdf'%(BaseName,ClassType))
#Save both the normalized 2d histos and well as those which are not normalized.
saveToFile(histoList,"%s_Histo_%s.root"%(BaseName,ClassType))
saveToFile(NormhistoList,"%s_NormHisto_%s.root"%(BaseName,ClassType))


#NewVar: log10(Cal1dEdxAve/Cal1FullLength), and log10(Cal1FitChiSquare)
#NewVar1:Cal1MomNumCoreXtals and not log10(Cal1FitChiSquare)
#NewVar2: log10(Cal1FitChiSquare) and not Cal1MomNumCoreXtals --Bad
