//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Fri Oct 17 16:23:02 2003 by ROOT version3.04/02)
//   from TTree 1/Glast tuple
//   found on file: fulltup.root
//////////////////////////////////////////////////////////


#ifndef tree_class_h
#define tree_class_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class tree_class {
   public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
//Declaration of leaves types
   Double_t        Run;
   Double_t        Event_ID;
   Double_t        MC_src_Id;
   Double_t        elapsed_time;
   Double_t        FilterStatus_HI;
   Double_t        FilterStatus_LO;
   Double_t        McId;
   Double_t        McCharge;
   Double_t        McEnergy;
   Double_t        McLogEnergy;
   Double_t        McX0;
   Double_t        McY0;
   Double_t        McZ0;
   Double_t        McXDir;
   Double_t        McYDir;
   Double_t        McZDir;
   Double_t        McXErr;
   Double_t        McYErr;
   Double_t        McZErr;
   Double_t        McXDirErr;
   Double_t        McYDirErr;
   Double_t        McZDirErr;
   Double_t        McDirErr;
   Double_t        McTkr1DirErr;
   Double_t        McTkr2DirErr;
   Double_t        GltWord;
   Double_t        GltTower;
   Double_t        GltXTower;
   Double_t        GltYTower;
   Double_t        GltLayer;
   Double_t        GltTotal;
   Double_t        GltNumTowers;
   Double_t        GltType;
   Double_t        GltMoment;
   Double_t        GltZDir;
   Double_t        TkrNumTracks;
   Double_t        TkrSumKalEne;
   Double_t        TkrSumConEne;
   Double_t        TkrEnergy;
   Double_t        TkrEnergySum;
   Double_t        TkrEnergyCorr;
   Double_t        TkrEdgeCorr;
   Double_t        TkrHDCount;
   Double_t        TkrTotalHits;
   Double_t        TkrThinHits;
   Double_t        TkrThickHits;
   Double_t        TkrBlankHits;
   Double_t        TkrRadLength;
   Double_t        TkrTwrEdge;
   Double_t        TkrTrackLength;
   Double_t        Tkr1Chisq;
   Double_t        Tkr1FirstChisq;
   Double_t        Tkr1Hits;
   Double_t        Tkr1FirstHits;
   Double_t        Tkr1FirstLayer;
   Double_t        Tkr1DifHits;
   Double_t        Tkr1Gaps;
   Double_t        Tkr1FirstGaps;
   Double_t        Tkr1Qual;
   Double_t        Tkr1Type;
   Double_t        Tkr1TwrEdge;
   Double_t        Tkr1PrjTwrEdge;
   Double_t        Tkr1DieEdge;
   Double_t        Tkr1TwrGap;
   Double_t        Tkr1KalEne;
   Double_t        Tkr1ConEne;
   Double_t        Tkr1KalThetaMS;
   Double_t        Tkr1XDir;
   Double_t        Tkr1YDir;
   Double_t        Tkr1ZDir;
   Double_t        Tkr1Phi;
   Double_t        Tkr1Theta;
   Double_t        Tkr1X0;
   Double_t        Tkr1Y0;
   Double_t        Tkr1Z0;
   Double_t        Tkr1ThetaErr;
   Double_t        Tkr1PhiErr;
   Double_t        Tkr1ErrAsym;
   Double_t        Tkr1CovDet;
   Double_t        Tkr1SXX;
   Double_t        Tkr1SXY;
   Double_t        Tkr1SYY;
   Double_t        Tkr1ToTFirst;
   Double_t        Tkr1ToTAve;
   Double_t        Tkr1ToTTrAve;
   Double_t        Tkr1ToTAsym;
   Double_t        Tkr1ChisqAsym;
   Double_t        Tkr1SSDVeto;
   Double_t        Tkr2Chisq;
   Double_t        Tkr2FirstChisq;
   Double_t        Tkr2Hits;
   Double_t        Tkr2FirstHits;
   Double_t        Tkr2FirstLayer;
   Double_t        Tkr2DifHits;
   Double_t        Tkr2Gaps;
   Double_t        Tkr2FirstGaps;
   Double_t        Tkr2Qual;
   Double_t        Tkr2Type;
   Double_t        Tkr2TwrEdge;
   Double_t        Tkr2PrjTwrEdge;
   Double_t        Tkr2DieEdge;
   Double_t        Tkr2KalEne;
   Double_t        Tkr2ConEne;
   Double_t        Tkr2KalThetaMS;
   Double_t        Tkr2XDir;
   Double_t        Tkr2YDir;
   Double_t        Tkr2ZDir;
   Double_t        Tkr2Phi;
   Double_t        Tkr2X0;
   Double_t        Tkr2Y0;
   Double_t        Tkr2Z0;
   Double_t        VtxXDir;
   Double_t        VtxYDir;
   Double_t        VtxZDir;
   Double_t        VtxPhi;
   Double_t        VtxTheta;
   Double_t        VtxX0;
   Double_t        VtxY0;
   Double_t        VtxZ0;
   Double_t        VtxAngle;
   Double_t        VtxDOCA;
   Double_t        VtxHeadSep;
   Double_t        VtxS1;
   Double_t        VtxS2;
   Double_t        VtxDocaWgt;
   Double_t        VtxHSWgt;
   Double_t        VtxS1Wgt;
   Double_t        VtxT12Wgt;
   Double_t        VtxT2QWgt;
   Double_t        VtxTotalWgt;
   Double_t        CalEnergySum;
   Double_t        CalEnergyCorr;
   Double_t        CalEneSumCorr;
   Double_t        Cal_Energy_LLCorr;
   Double_t        CalLeakCorr2;
   Double_t        CalEdgeSumCorr;
   Double_t        CalTotSumCorr;
   Double_t        CalCsIRLn;
   Double_t        CalTotRLn;
   Double_t        CalCntRLn;
   Double_t        CalDeadTotRat;
   Double_t        CalDeadCntRat;
   Double_t        CalTPred;
   Double_t        CalDeltaT;
   Double_t        CalTwrEdge;
   Double_t        CalLATEdge;
   Double_t        CalTENrm;
   Double_t        CalTrackSep;
   Double_t        CalTrackDoca;
   Double_t        CalTwrGap;
   Double_t        CalELayer0;
   Double_t        CalELayer1;
   Double_t        CalELayer2;
   Double_t        CalELayer3;
   Double_t        CalELayer4;
   Double_t        CalELayer5;
   Double_t        CalELayer6;
   Double_t        CalELayer7;
   Double_t        CalLyr0Ratio;
   Double_t        CalLyr7Ratio;
   Double_t        CalBkHalfRatio;
   Double_t        CalXtalsTrunc;
   Double_t        CalXtalRatio;
   Double_t        CalLongRms;
   Double_t        CalLRmsRatio;
   Double_t        CalTransRms;
   Double_t        CalMIPDiff;
   Double_t        CalMIPRatio;
   Double_t        CalXEcntr;
   Double_t        CalYEcntr;
   Double_t        CalZEcntr;
   Double_t        CalXDir;
   Double_t        CalYDir;
   Double_t        CalZDir;
   Double_t        CalX0;
   Double_t        CalY0;
   Double_t        CalZ0;
   Double_t        AcdTotalEnergy;
   Double_t        AcdTileCount;
   Double_t        AcdDoca;
   Double_t        AcdActiveDist;
   Double_t        AcdGammaDoca;
   Double_t        AcdActDistTop;
   Double_t        AcdActDistSideRow0;
   Double_t        AcdActDistSideRow1;
   Double_t        AcdActDistSideRow2;
   Double_t        AcdNoSideRow0;
   Double_t        AcdNoSideRow1;
   Double_t        AcdNoSideRow2;
   Double_t        EvtEnergySumOpt;
   Double_t        EvtEnergyTracker;
   Double_t        EvtEnergyRaw;
   Double_t        EvtMcEnergySigma;
   Double_t        EvtCalEdgeAngle;
   Double_t        EvtTkrEdgeAngle;
   Double_t        EvtLogESum;
   Double_t        EvtTkr1EFrac;
   Double_t        EvtVtxKin;
   Double_t        EvtVtxEAngle;
   Double_t        EvtTkrComptonRatio;
   Double_t        EvtTkrEComptonRatio;
   Double_t        EvtPSFModel;
   Double_t        EvtTkr1EChisq;
   Double_t        EvtTkr1EFirstChisq;
   Double_t        EvtTkr1EQual;
   Double_t        EvtTkr1PSFMdRat;
   Double_t        EvtTkr1ECovDet;
   Double_t        EvtTkr2EChisq;
   Double_t        EvtTkr2EFirstChisq;
   Double_t        EvtTkr2EQual;
   Double_t        EvtCalETLRatio;
   Double_t        EvtCalEXtalRatio;
   Double_t        EvtCalEXtalTrunc;
   Double_t        EvtCalETrackDoca;
   Double_t        EvtCalETrackSep;
   Double_t        EvtVtxEEAngle;
   Double_t        EvtVtxEDoca;
   Double_t        EvtVtxEHeadSep;
   Double_t        IMgoodCalProb;
   Double_t        IMvertexProb;
   Double_t        IMcoreProb;
   Double_t        IMpsfErrPred;
   Double_t        IMgammaProb;

//List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event_ID;   //!
   TBranch        *b_MC_src_Id;   //!
   TBranch        *b_elapsed_time;   //!
   TBranch        *b_FilterStatus_HI;   //!
   TBranch        *b_FilterStatus_LO;   //!
   TBranch        *b_McId;   //!
   TBranch        *b_McCharge;   //!
   TBranch        *b_McEnergy;   //!
   TBranch        *b_McLogEnergy;   //!
   TBranch        *b_McX0;   //!
   TBranch        *b_McY0;   //!
   TBranch        *b_McZ0;   //!
   TBranch        *b_McXDir;   //!
   TBranch        *b_McYDir;   //!
   TBranch        *b_McZDir;   //!
   TBranch        *b_McXErr;   //!
   TBranch        *b_McYErr;   //!
   TBranch        *b_McZErr;   //!
   TBranch        *b_McXDirErr;   //!
   TBranch        *b_McYDirErr;   //!
   TBranch        *b_McZDirErr;   //!
   TBranch        *b_McDirErr;   //!
   TBranch        *b_McTkr1DirErr;   //!
   TBranch        *b_McTkr2DirErr;   //!
   TBranch        *b_GltWord;   //!
   TBranch        *b_GltTower;   //!
   TBranch        *b_GltXTower;   //!
   TBranch        *b_GltYTower;   //!
   TBranch        *b_GltLayer;   //!
   TBranch        *b_GltTotal;   //!
   TBranch        *b_GltNumTowers;   //!
   TBranch        *b_GltType;   //!
   TBranch        *b_GltMoment;   //!
   TBranch        *b_GltZDir;   //!
   TBranch        *b_TkrNumTracks;   //!
   TBranch        *b_TkrSumKalEne;   //!
   TBranch        *b_TkrSumConEne;   //!
   TBranch        *b_TkrEnergy;   //!
   TBranch        *b_TkrEnergySum;   //!
   TBranch        *b_TkrEnergyCorr;   //!
   TBranch        *b_TkrEdgeCorr;   //!
   TBranch        *b_TkrHDCount;   //!
   TBranch        *b_TkrTotalHits;   //!
   TBranch        *b_TkrThinHits;   //!
   TBranch        *b_TkrThickHits;   //!
   TBranch        *b_TkrBlankHits;   //!
   TBranch        *b_TkrRadLength;   //!
   TBranch        *b_TkrTwrEdge;   //!
   TBranch        *b_TkrTrackLength;   //!
   TBranch        *b_Tkr1Chisq;   //!
   TBranch        *b_Tkr1FirstChisq;   //!
   TBranch        *b_Tkr1Hits;   //!
   TBranch        *b_Tkr1FirstHits;   //!
   TBranch        *b_Tkr1FirstLayer;   //!
   TBranch        *b_Tkr1DifHits;   //!
   TBranch        *b_Tkr1Gaps;   //!
   TBranch        *b_Tkr1FirstGaps;   //!
   TBranch        *b_Tkr1Qual;   //!
   TBranch        *b_Tkr1Type;   //!
   TBranch        *b_Tkr1TwrEdge;   //!
   TBranch        *b_Tkr1PrjTwrEdge;   //!
   TBranch        *b_Tkr1DieEdge;   //!
   TBranch        *b_Tkr1TwrGap;   //!
   TBranch        *b_Tkr1KalEne;   //!
   TBranch        *b_Tkr1ConEne;   //!
   TBranch        *b_Tkr1KalThetaMS;   //!
   TBranch        *b_Tkr1XDir;   //!
   TBranch        *b_Tkr1YDir;   //!
   TBranch        *b_Tkr1ZDir;   //!
   TBranch        *b_Tkr1Phi;   //!
   TBranch        *b_Tkr1Theta;   //!
   TBranch        *b_Tkr1X0;   //!
   TBranch        *b_Tkr1Y0;   //!
   TBranch        *b_Tkr1Z0;   //!
   TBranch        *b_Tkr1ThetaErr;   //!
   TBranch        *b_Tkr1PhiErr;   //!
   TBranch        *b_Tkr1ErrAsym;   //!
   TBranch        *b_Tkr1CovDet;   //!
   TBranch        *b_Tkr1SXX;   //!
   TBranch        *b_Tkr1SXY;   //!
   TBranch        *b_Tkr1SYY;   //!
   TBranch        *b_Tkr1ToTFirst;   //!
   TBranch        *b_Tkr1ToTAve;   //!
   TBranch        *b_Tkr1ToTTrAve;   //!
   TBranch        *b_Tkr1ToTAsym;   //!
   TBranch        *b_Tkr1ChisqAsym;   //!
   TBranch        *b_Tkr1SSDVeto;   //!
   TBranch        *b_Tkr2Chisq;   //!
   TBranch        *b_Tkr2FirstChisq;   //!
   TBranch        *b_Tkr2Hits;   //!
   TBranch        *b_Tkr2FirstHits;   //!
   TBranch        *b_Tkr2FirstLayer;   //!
   TBranch        *b_Tkr2DifHits;   //!
   TBranch        *b_Tkr2Gaps;   //!
   TBranch        *b_Tkr2FirstGaps;   //!
   TBranch        *b_Tkr2Qual;   //!
   TBranch        *b_Tkr2Type;   //!
   TBranch        *b_Tkr2TwrEdge;   //!
   TBranch        *b_Tkr2PrjTwrEdge;   //!
   TBranch        *b_Tkr2DieEdge;   //!
   TBranch        *b_Tkr2KalEne;   //!
   TBranch        *b_Tkr2ConEne;   //!
   TBranch        *b_Tkr2KalThetaMS;   //!
   TBranch        *b_Tkr2XDir;   //!
   TBranch        *b_Tkr2YDir;   //!
   TBranch        *b_Tkr2ZDir;   //!
   TBranch        *b_Tkr2Phi;   //!
   TBranch        *b_Tkr2X0;   //!
   TBranch        *b_Tkr2Y0;   //!
   TBranch        *b_Tkr2Z0;   //!
   TBranch        *b_VtxXDir;   //!
   TBranch        *b_VtxYDir;   //!
   TBranch        *b_VtxZDir;   //!
   TBranch        *b_VtxPhi;   //!
   TBranch        *b_VtxTheta;   //!
   TBranch        *b_VtxX0;   //!
   TBranch        *b_VtxY0;   //!
   TBranch        *b_VtxZ0;   //!
   TBranch        *b_VtxAngle;   //!
   TBranch        *b_VtxDOCA;   //!
   TBranch        *b_VtxHeadSep;   //!
   TBranch        *b_VtxS1;   //!
   TBranch        *b_VtxS2;   //!
   TBranch        *b_VtxDocaWgt;   //!
   TBranch        *b_VtxHSWgt;   //!
   TBranch        *b_VtxS1Wgt;   //!
   TBranch        *b_VtxT12Wgt;   //!
   TBranch        *b_VtxT2QWgt;   //!
   TBranch        *b_VtxTotalWgt;   //!
   TBranch        *b_CalEnergySum;   //!
   TBranch        *b_CalEnergyCorr;   //!
   TBranch        *b_CalEneSumCorr;   //!
   TBranch        *b_Cal_Energy_LLCorr;   //!
   TBranch        *b_CalLeakCorr2;   //!
   TBranch        *b_CalEdgeSumCorr;   //!
   TBranch        *b_CalTotSumCorr;   //!
   TBranch        *b_CalCsIRLn;   //!
   TBranch        *b_CalTotRLn;   //!
   TBranch        *b_CalCntRLn;   //!
   TBranch        *b_CalDeadTotRat;   //!
   TBranch        *b_CalDeadCntRat;   //!
   TBranch        *b_CalTPred;   //!
   TBranch        *b_CalDeltaT;   //!
   TBranch        *b_CalTwrEdge;   //!
   TBranch        *b_CalLATEdge;   //!
   TBranch        *b_CalTENrm;   //!
   TBranch        *b_CalTrackSep;   //!
   TBranch        *b_CalTrackDoca;   //!
   TBranch        *b_CalTwrGap;   //!
   TBranch        *b_CalELayer0;   //!
   TBranch        *b_CalELayer1;   //!
   TBranch        *b_CalELayer2;   //!
   TBranch        *b_CalELayer3;   //!
   TBranch        *b_CalELayer4;   //!
   TBranch        *b_CalELayer5;   //!
   TBranch        *b_CalELayer6;   //!
   TBranch        *b_CalELayer7;   //!
   TBranch        *b_CalLyr0Ratio;   //!
   TBranch        *b_CalLyr7Ratio;   //!
   TBranch        *b_CalBkHalfRatio;   //!
   TBranch        *b_CalXtalsTrunc;   //!
   TBranch        *b_CalXtalRatio;   //!
   TBranch        *b_CalLongRms;   //!
   TBranch        *b_CalLRmsRatio;   //!
   TBranch        *b_CalTransRms;   //!
   TBranch        *b_CalMIPDiff;   //!
   TBranch        *b_CalMIPRatio;   //!
   TBranch        *b_CalXEcntr;   //!
   TBranch        *b_CalYEcntr;   //!
   TBranch        *b_CalZEcntr;   //!
   TBranch        *b_CalXDir;   //!
   TBranch        *b_CalYDir;   //!
   TBranch        *b_CalZDir;   //!
   TBranch        *b_CalX0;   //!
   TBranch        *b_CalY0;   //!
   TBranch        *b_CalZ0;   //!
   TBranch        *b_AcdTotalEnergy;   //!
   TBranch        *b_AcdTileCount;   //!
   TBranch        *b_AcdDoca;   //!
   TBranch        *b_AcdActiveDist;   //!
   TBranch        *b_AcdGammaDoca;   //!
   TBranch        *b_AcdActDistTop;   //!
   TBranch        *b_AcdActDistSideRow0;   //!
   TBranch        *b_AcdActDistSideRow1;   //!
   TBranch        *b_AcdActDistSideRow2;   //!
   TBranch        *b_AcdNoSideRow0;   //!
   TBranch        *b_AcdNoSideRow1;   //!
   TBranch        *b_AcdNoSideRow2;   //!
   TBranch        *b_EvtEnergySumOpt;   //!
   TBranch        *b_EvtEnergyTracker;   //!
   TBranch        *b_EvtEnergyRaw;   //!
   TBranch        *b_EvtMcEnergySigma;   //!
   TBranch        *b_EvtCalEdgeAngle;   //!
   TBranch        *b_EvtTkrEdgeAngle;   //!
   TBranch        *b_EvtLogESum;   //!
   TBranch        *b_EvtTkr1EFrac;   //!
   TBranch        *b_EvtVtxKin;   //!
   TBranch        *b_EvtVtxEAngle;   //!
   TBranch        *b_EvtTkrComptonRatio;   //!
   TBranch        *b_EvtTkrEComptonRatio;   //!
   TBranch        *b_EvtPSFModel;   //!
   TBranch        *b_EvtTkr1EChisq;   //!
   TBranch        *b_EvtTkr1EFirstChisq;   //!
   TBranch        *b_EvtTkr1EQual;   //!
   TBranch        *b_EvtTkr1PSFMdRat;   //!
   TBranch        *b_EvtTkr1ECovDet;   //!
   TBranch        *b_EvtTkr2EChisq;   //!
   TBranch        *b_EvtTkr2EFirstChisq;   //!
   TBranch        *b_EvtTkr2EQual;   //!
   TBranch        *b_EvtCalETLRatio;   //!
   TBranch        *b_EvtCalEXtalRatio;   //!
   TBranch        *b_EvtCalEXtalTrunc;   //!
   TBranch        *b_EvtCalETrackDoca;   //!
   TBranch        *b_EvtCalETrackSep;   //!
   TBranch        *b_EvtVtxEEAngle;   //!
   TBranch        *b_EvtVtxEDoca;   //!
   TBranch        *b_EvtVtxEHeadSep;   //!
   TBranch        *b_IMgoodCalProb;   //!
   TBranch        *b_IMvertexProb;   //!
   TBranch        *b_IMcoreProb;   //!
   TBranch        *b_IMpsfErrPred;   //!
   TBranch        *b_IMgammaProb;   //!

   tree_class(TTree *tree=0);
   ~tree_class();
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop();
   Bool_t Notify();
   void   Show(Int_t entry = -1);
};

#endif

#ifdef tree_class_cxx
tree_class::tree_class(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("fulltup.root");
      if (!f) {
         f = new TFile("fulltup.root");
      }
      tree = (TTree*)gDirectory->Get("1");

   }
   Init(tree);
}

tree_class::~tree_class()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree_class::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t tree_class::LoadTree(Int_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tree_class::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run",&Run);
   fChain->SetBranchAddress("Event_ID",&Event_ID);
   fChain->SetBranchAddress("MC_src_Id",&MC_src_Id);
   fChain->SetBranchAddress("elapsed_time",&elapsed_time);
   fChain->SetBranchAddress("FilterStatus_HI",&FilterStatus_HI);
   fChain->SetBranchAddress("FilterStatus_LO",&FilterStatus_LO);
   fChain->SetBranchAddress("McId",&McId);
   fChain->SetBranchAddress("McCharge",&McCharge);
   fChain->SetBranchAddress("McEnergy",&McEnergy);
   fChain->SetBranchAddress("McLogEnergy",&McLogEnergy);
   fChain->SetBranchAddress("McX0",&McX0);
   fChain->SetBranchAddress("McY0",&McY0);
   fChain->SetBranchAddress("McZ0",&McZ0);
   fChain->SetBranchAddress("McXDir",&McXDir);
   fChain->SetBranchAddress("McYDir",&McYDir);
   fChain->SetBranchAddress("McZDir",&McZDir);
   fChain->SetBranchAddress("McXErr",&McXErr);
   fChain->SetBranchAddress("McYErr",&McYErr);
   fChain->SetBranchAddress("McZErr",&McZErr);
   fChain->SetBranchAddress("McXDirErr",&McXDirErr);
   fChain->SetBranchAddress("McYDirErr",&McYDirErr);
   fChain->SetBranchAddress("McZDirErr",&McZDirErr);
   fChain->SetBranchAddress("McDirErr",&McDirErr);
   fChain->SetBranchAddress("McTkr1DirErr",&McTkr1DirErr);
   fChain->SetBranchAddress("McTkr2DirErr",&McTkr2DirErr);
   fChain->SetBranchAddress("GltWord",&GltWord);
   fChain->SetBranchAddress("GltTower",&GltTower);
   fChain->SetBranchAddress("GltXTower",&GltXTower);
   fChain->SetBranchAddress("GltYTower",&GltYTower);
   fChain->SetBranchAddress("GltLayer",&GltLayer);
   fChain->SetBranchAddress("GltTotal",&GltTotal);
   fChain->SetBranchAddress("GltNumTowers",&GltNumTowers);
   fChain->SetBranchAddress("GltType",&GltType);
   fChain->SetBranchAddress("GltMoment",&GltMoment);
   fChain->SetBranchAddress("GltZDir",&GltZDir);
   fChain->SetBranchAddress("TkrNumTracks",&TkrNumTracks);
   fChain->SetBranchAddress("TkrSumKalEne",&TkrSumKalEne);
   fChain->SetBranchAddress("TkrSumConEne",&TkrSumConEne);
   fChain->SetBranchAddress("TkrEnergy",&TkrEnergy);
   fChain->SetBranchAddress("TkrEnergySum",&TkrEnergySum);
   fChain->SetBranchAddress("TkrEnergyCorr",&TkrEnergyCorr);
   fChain->SetBranchAddress("TkrEdgeCorr",&TkrEdgeCorr);
   fChain->SetBranchAddress("TkrHDCount",&TkrHDCount);
   fChain->SetBranchAddress("TkrTotalHits",&TkrTotalHits);
   fChain->SetBranchAddress("TkrThinHits",&TkrThinHits);
   fChain->SetBranchAddress("TkrThickHits",&TkrThickHits);
   fChain->SetBranchAddress("TkrBlankHits",&TkrBlankHits);
   fChain->SetBranchAddress("TkrRadLength",&TkrRadLength);
   fChain->SetBranchAddress("TkrTwrEdge",&TkrTwrEdge);
   fChain->SetBranchAddress("TkrTrackLength",&TkrTrackLength);
   fChain->SetBranchAddress("Tkr1Chisq",&Tkr1Chisq);
   fChain->SetBranchAddress("Tkr1FirstChisq",&Tkr1FirstChisq);
   fChain->SetBranchAddress("Tkr1Hits",&Tkr1Hits);
   fChain->SetBranchAddress("Tkr1FirstHits",&Tkr1FirstHits);
   fChain->SetBranchAddress("Tkr1FirstLayer",&Tkr1FirstLayer);
   fChain->SetBranchAddress("Tkr1DifHits",&Tkr1DifHits);
   fChain->SetBranchAddress("Tkr1Gaps",&Tkr1Gaps);
   fChain->SetBranchAddress("Tkr1FirstGaps",&Tkr1FirstGaps);
   fChain->SetBranchAddress("Tkr1Qual",&Tkr1Qual);
   fChain->SetBranchAddress("Tkr1Type",&Tkr1Type);
   fChain->SetBranchAddress("Tkr1TwrEdge",&Tkr1TwrEdge);
   fChain->SetBranchAddress("Tkr1PrjTwrEdge",&Tkr1PrjTwrEdge);
   fChain->SetBranchAddress("Tkr1DieEdge",&Tkr1DieEdge);
   fChain->SetBranchAddress("Tkr1TwrGap",&Tkr1TwrGap);
   fChain->SetBranchAddress("Tkr1KalEne",&Tkr1KalEne);
   fChain->SetBranchAddress("Tkr1ConEne",&Tkr1ConEne);
   fChain->SetBranchAddress("Tkr1KalThetaMS",&Tkr1KalThetaMS);
   fChain->SetBranchAddress("Tkr1XDir",&Tkr1XDir);
   fChain->SetBranchAddress("Tkr1YDir",&Tkr1YDir);
   fChain->SetBranchAddress("Tkr1ZDir",&Tkr1ZDir);
   fChain->SetBranchAddress("Tkr1Phi",&Tkr1Phi);
   fChain->SetBranchAddress("Tkr1Theta",&Tkr1Theta);
   fChain->SetBranchAddress("Tkr1X0",&Tkr1X0);
   fChain->SetBranchAddress("Tkr1Y0",&Tkr1Y0);
   fChain->SetBranchAddress("Tkr1Z0",&Tkr1Z0);
   fChain->SetBranchAddress("Tkr1ThetaErr",&Tkr1ThetaErr);
   fChain->SetBranchAddress("Tkr1PhiErr",&Tkr1PhiErr);
   fChain->SetBranchAddress("Tkr1ErrAsym",&Tkr1ErrAsym);
   fChain->SetBranchAddress("Tkr1CovDet",&Tkr1CovDet);
   fChain->SetBranchAddress("Tkr1SXX",&Tkr1SXX);
   fChain->SetBranchAddress("Tkr1SXY",&Tkr1SXY);
   fChain->SetBranchAddress("Tkr1SYY",&Tkr1SYY);
   fChain->SetBranchAddress("Tkr1ToTFirst",&Tkr1ToTFirst);
   fChain->SetBranchAddress("Tkr1ToTAve",&Tkr1ToTAve);
   fChain->SetBranchAddress("Tkr1ToTTrAve",&Tkr1ToTTrAve);
   fChain->SetBranchAddress("Tkr1ToTAsym",&Tkr1ToTAsym);
   fChain->SetBranchAddress("Tkr1ChisqAsym",&Tkr1ChisqAsym);
   fChain->SetBranchAddress("Tkr1SSDVeto",&Tkr1SSDVeto);
   fChain->SetBranchAddress("Tkr2Chisq",&Tkr2Chisq);
   fChain->SetBranchAddress("Tkr2FirstChisq",&Tkr2FirstChisq);
   fChain->SetBranchAddress("Tkr2Hits",&Tkr2Hits);
   fChain->SetBranchAddress("Tkr2FirstHits",&Tkr2FirstHits);
   fChain->SetBranchAddress("Tkr2FirstLayer",&Tkr2FirstLayer);
   fChain->SetBranchAddress("Tkr2DifHits",&Tkr2DifHits);
   fChain->SetBranchAddress("Tkr2Gaps",&Tkr2Gaps);
   fChain->SetBranchAddress("Tkr2FirstGaps",&Tkr2FirstGaps);
   fChain->SetBranchAddress("Tkr2Qual",&Tkr2Qual);
   fChain->SetBranchAddress("Tkr2Type",&Tkr2Type);
   fChain->SetBranchAddress("Tkr2TwrEdge",&Tkr2TwrEdge);
   fChain->SetBranchAddress("Tkr2PrjTwrEdge",&Tkr2PrjTwrEdge);
   fChain->SetBranchAddress("Tkr2DieEdge",&Tkr2DieEdge);
   fChain->SetBranchAddress("Tkr2KalEne",&Tkr2KalEne);
   fChain->SetBranchAddress("Tkr2ConEne",&Tkr2ConEne);
   fChain->SetBranchAddress("Tkr2KalThetaMS",&Tkr2KalThetaMS);
   fChain->SetBranchAddress("Tkr2XDir",&Tkr2XDir);
   fChain->SetBranchAddress("Tkr2YDir",&Tkr2YDir);
   fChain->SetBranchAddress("Tkr2ZDir",&Tkr2ZDir);
   fChain->SetBranchAddress("Tkr2Phi",&Tkr2Phi);
   fChain->SetBranchAddress("Tkr2X0",&Tkr2X0);
   fChain->SetBranchAddress("Tkr2Y0",&Tkr2Y0);
   fChain->SetBranchAddress("Tkr2Z0",&Tkr2Z0);
   fChain->SetBranchAddress("VtxXDir",&VtxXDir);
   fChain->SetBranchAddress("VtxYDir",&VtxYDir);
   fChain->SetBranchAddress("VtxZDir",&VtxZDir);
   fChain->SetBranchAddress("VtxPhi",&VtxPhi);
   fChain->SetBranchAddress("VtxTheta",&VtxTheta);
   fChain->SetBranchAddress("VtxX0",&VtxX0);
   fChain->SetBranchAddress("VtxY0",&VtxY0);
   fChain->SetBranchAddress("VtxZ0",&VtxZ0);
   fChain->SetBranchAddress("VtxAngle",&VtxAngle);
   fChain->SetBranchAddress("VtxDOCA",&VtxDOCA);
   fChain->SetBranchAddress("VtxHeadSep",&VtxHeadSep);
   fChain->SetBranchAddress("VtxS1",&VtxS1);
   fChain->SetBranchAddress("VtxS2",&VtxS2);
   fChain->SetBranchAddress("VtxDocaWgt",&VtxDocaWgt);
   fChain->SetBranchAddress("VtxHSWgt",&VtxHSWgt);
   fChain->SetBranchAddress("VtxS1Wgt",&VtxS1Wgt);
   fChain->SetBranchAddress("VtxT12Wgt",&VtxT12Wgt);
   fChain->SetBranchAddress("VtxT2QWgt",&VtxT2QWgt);
   fChain->SetBranchAddress("VtxTotalWgt",&VtxTotalWgt);
   fChain->SetBranchAddress("CalEnergySum",&CalEnergySum);
   fChain->SetBranchAddress("CalEnergyCorr",&CalEnergyCorr);
   fChain->SetBranchAddress("CalEneSumCorr",&CalEneSumCorr);
   fChain->SetBranchAddress("Cal_Energy_LLCorr",&Cal_Energy_LLCorr);
   fChain->SetBranchAddress("CalLeakCorr2",&CalLeakCorr2);
   fChain->SetBranchAddress("CalEdgeSumCorr",&CalEdgeSumCorr);
   fChain->SetBranchAddress("CalTotSumCorr",&CalTotSumCorr);
   fChain->SetBranchAddress("CalCsIRLn",&CalCsIRLn);
   fChain->SetBranchAddress("CalTotRLn",&CalTotRLn);
   fChain->SetBranchAddress("CalCntRLn",&CalCntRLn);
   fChain->SetBranchAddress("CalDeadTotRat",&CalDeadTotRat);
   fChain->SetBranchAddress("CalDeadCntRat",&CalDeadCntRat);
   fChain->SetBranchAddress("CalTPred",&CalTPred);
   fChain->SetBranchAddress("CalDeltaT",&CalDeltaT);
   fChain->SetBranchAddress("CalTwrEdge",&CalTwrEdge);
   fChain->SetBranchAddress("CalLATEdge",&CalLATEdge);
   fChain->SetBranchAddress("CalTENrm",&CalTENrm);
   fChain->SetBranchAddress("CalTrackSep",&CalTrackSep);
   fChain->SetBranchAddress("CalTrackDoca",&CalTrackDoca);
   fChain->SetBranchAddress("CalTwrGap",&CalTwrGap);
   fChain->SetBranchAddress("CalELayer0",&CalELayer0);
   fChain->SetBranchAddress("CalELayer1",&CalELayer1);
   fChain->SetBranchAddress("CalELayer2",&CalELayer2);
   fChain->SetBranchAddress("CalELayer3",&CalELayer3);
   fChain->SetBranchAddress("CalELayer4",&CalELayer4);
   fChain->SetBranchAddress("CalELayer5",&CalELayer5);
   fChain->SetBranchAddress("CalELayer6",&CalELayer6);
   fChain->SetBranchAddress("CalELayer7",&CalELayer7);
   fChain->SetBranchAddress("CalLyr0Ratio",&CalLyr0Ratio);
   fChain->SetBranchAddress("CalLyr7Ratio",&CalLyr7Ratio);
   fChain->SetBranchAddress("CalBkHalfRatio",&CalBkHalfRatio);
   fChain->SetBranchAddress("CalXtalsTrunc",&CalXtalsTrunc);
   fChain->SetBranchAddress("CalXtalRatio",&CalXtalRatio);
   fChain->SetBranchAddress("CalLongRms",&CalLongRms);
   fChain->SetBranchAddress("CalLRmsRatio",&CalLRmsRatio);
   fChain->SetBranchAddress("CalTransRms",&CalTransRms);
   fChain->SetBranchAddress("CalMIPDiff",&CalMIPDiff);
   fChain->SetBranchAddress("CalMIPRatio",&CalMIPRatio);
   fChain->SetBranchAddress("CalXEcntr",&CalXEcntr);
   fChain->SetBranchAddress("CalYEcntr",&CalYEcntr);
   fChain->SetBranchAddress("CalZEcntr",&CalZEcntr);
   fChain->SetBranchAddress("CalXDir",&CalXDir);
   fChain->SetBranchAddress("CalYDir",&CalYDir);
   fChain->SetBranchAddress("CalZDir",&CalZDir);
   fChain->SetBranchAddress("CalX0",&CalX0);
   fChain->SetBranchAddress("CalY0",&CalY0);
   fChain->SetBranchAddress("CalZ0",&CalZ0);
   fChain->SetBranchAddress("AcdTotalEnergy",&AcdTotalEnergy);
   fChain->SetBranchAddress("AcdTileCount",&AcdTileCount);
   fChain->SetBranchAddress("AcdDoca",&AcdDoca);
   fChain->SetBranchAddress("AcdActiveDist",&AcdActiveDist);
   fChain->SetBranchAddress("AcdGammaDoca",&AcdGammaDoca);
   fChain->SetBranchAddress("AcdActDistTop",&AcdActDistTop);
   fChain->SetBranchAddress("AcdActDistSideRow0",&AcdActDistSideRow0);
   fChain->SetBranchAddress("AcdActDistSideRow1",&AcdActDistSideRow1);
   fChain->SetBranchAddress("AcdActDistSideRow2",&AcdActDistSideRow2);
   fChain->SetBranchAddress("AcdNoSideRow0",&AcdNoSideRow0);
   fChain->SetBranchAddress("AcdNoSideRow1",&AcdNoSideRow1);
   fChain->SetBranchAddress("AcdNoSideRow2",&AcdNoSideRow2);
   fChain->SetBranchAddress("EvtEnergySumOpt",&EvtEnergySumOpt);
   fChain->SetBranchAddress("EvtEnergyTracker",&EvtEnergyTracker);
   fChain->SetBranchAddress("EvtEnergyRaw",&EvtEnergyRaw);
   fChain->SetBranchAddress("EvtMcEnergySigma",&EvtMcEnergySigma);
   fChain->SetBranchAddress("EvtCalEdgeAngle",&EvtCalEdgeAngle);
   fChain->SetBranchAddress("EvtTkrEdgeAngle",&EvtTkrEdgeAngle);
   fChain->SetBranchAddress("EvtLogESum",&EvtLogESum);
   fChain->SetBranchAddress("EvtTkr1EFrac",&EvtTkr1EFrac);
   fChain->SetBranchAddress("EvtVtxKin",&EvtVtxKin);
   fChain->SetBranchAddress("EvtVtxEAngle",&EvtVtxEAngle);
   fChain->SetBranchAddress("EvtTkrComptonRatio",&EvtTkrComptonRatio);
   fChain->SetBranchAddress("EvtTkrEComptonRatio",&EvtTkrEComptonRatio);
   fChain->SetBranchAddress("EvtPSFModel",&EvtPSFModel);
   fChain->SetBranchAddress("EvtTkr1EChisq",&EvtTkr1EChisq);
   fChain->SetBranchAddress("EvtTkr1EFirstChisq",&EvtTkr1EFirstChisq);
   fChain->SetBranchAddress("EvtTkr1EQual",&EvtTkr1EQual);
   fChain->SetBranchAddress("EvtTkr1PSFMdRat",&EvtTkr1PSFMdRat);
   fChain->SetBranchAddress("EvtTkr1ECovDet",&EvtTkr1ECovDet);
   fChain->SetBranchAddress("EvtTkr2EChisq",&EvtTkr2EChisq);
   fChain->SetBranchAddress("EvtTkr2EFirstChisq",&EvtTkr2EFirstChisq);
   fChain->SetBranchAddress("EvtTkr2EQual",&EvtTkr2EQual);
   fChain->SetBranchAddress("EvtCalETLRatio",&EvtCalETLRatio);
   fChain->SetBranchAddress("EvtCalEXtalRatio",&EvtCalEXtalRatio);
   fChain->SetBranchAddress("EvtCalEXtalTrunc",&EvtCalEXtalTrunc);
   fChain->SetBranchAddress("EvtCalETrackDoca",&EvtCalETrackDoca);
   fChain->SetBranchAddress("EvtCalETrackSep",&EvtCalETrackSep);
   fChain->SetBranchAddress("EvtVtxEEAngle",&EvtVtxEEAngle);
   fChain->SetBranchAddress("EvtVtxEDoca",&EvtVtxEDoca);
   fChain->SetBranchAddress("EvtVtxEHeadSep",&EvtVtxEHeadSep);
   fChain->SetBranchAddress("IMgoodCalProb",&IMgoodCalProb);
   fChain->SetBranchAddress("IMvertexProb",&IMvertexProb);
   fChain->SetBranchAddress("IMcoreProb",&IMcoreProb);
   fChain->SetBranchAddress("IMpsfErrPred",&IMpsfErrPred);
   fChain->SetBranchAddress("IMgammaProb",&IMgammaProb);
   Notify();
}

Bool_t tree_class::Notify()
{
   // Called when loading a new file.
   // Get branch pointers.
   b_Run = fChain->GetBranch("Run");
   b_Event_ID = fChain->GetBranch("Event_ID");
   b_MC_src_Id = fChain->GetBranch("MC_src_Id");
   b_elapsed_time = fChain->GetBranch("elapsed_time");
   b_FilterStatus_HI = fChain->GetBranch("FilterStatus_HI");
   b_FilterStatus_LO = fChain->GetBranch("FilterStatus_LO");
   b_McId = fChain->GetBranch("McId");
   b_McCharge = fChain->GetBranch("McCharge");
   b_McEnergy = fChain->GetBranch("McEnergy");
   b_McLogEnergy = fChain->GetBranch("McLogEnergy");
   b_McX0 = fChain->GetBranch("McX0");
   b_McY0 = fChain->GetBranch("McY0");
   b_McZ0 = fChain->GetBranch("McZ0");
   b_McXDir = fChain->GetBranch("McXDir");
   b_McYDir = fChain->GetBranch("McYDir");
   b_McZDir = fChain->GetBranch("McZDir");
   b_McXErr = fChain->GetBranch("McXErr");
   b_McYErr = fChain->GetBranch("McYErr");
   b_McZErr = fChain->GetBranch("McZErr");
   b_McXDirErr = fChain->GetBranch("McXDirErr");
   b_McYDirErr = fChain->GetBranch("McYDirErr");
   b_McZDirErr = fChain->GetBranch("McZDirErr");
   b_McDirErr = fChain->GetBranch("McDirErr");
   b_McTkr1DirErr = fChain->GetBranch("McTkr1DirErr");
   b_McTkr2DirErr = fChain->GetBranch("McTkr2DirErr");
   b_GltWord = fChain->GetBranch("GltWord");
   b_GltTower = fChain->GetBranch("GltTower");
   b_GltXTower = fChain->GetBranch("GltXTower");
   b_GltYTower = fChain->GetBranch("GltYTower");
   b_GltLayer = fChain->GetBranch("GltLayer");
   b_GltTotal = fChain->GetBranch("GltTotal");
   b_GltNumTowers = fChain->GetBranch("GltNumTowers");
   b_GltType = fChain->GetBranch("GltType");
   b_GltMoment = fChain->GetBranch("GltMoment");
   b_GltZDir = fChain->GetBranch("GltZDir");
   b_TkrNumTracks = fChain->GetBranch("TkrNumTracks");
   b_TkrSumKalEne = fChain->GetBranch("TkrSumKalEne");
   b_TkrSumConEne = fChain->GetBranch("TkrSumConEne");
   b_TkrEnergy = fChain->GetBranch("TkrEnergy");
   b_TkrEnergySum = fChain->GetBranch("TkrEnergySum");
   b_TkrEnergyCorr = fChain->GetBranch("TkrEnergyCorr");
   b_TkrEdgeCorr = fChain->GetBranch("TkrEdgeCorr");
   b_TkrHDCount = fChain->GetBranch("TkrHDCount");
   b_TkrTotalHits = fChain->GetBranch("TkrTotalHits");
   b_TkrThinHits = fChain->GetBranch("TkrThinHits");
   b_TkrThickHits = fChain->GetBranch("TkrThickHits");
   b_TkrBlankHits = fChain->GetBranch("TkrBlankHits");
   b_TkrRadLength = fChain->GetBranch("TkrRadLength");
   b_TkrTwrEdge = fChain->GetBranch("TkrTwrEdge");
   b_TkrTrackLength = fChain->GetBranch("TkrTrackLength");
   b_Tkr1Chisq = fChain->GetBranch("Tkr1Chisq");
   b_Tkr1FirstChisq = fChain->GetBranch("Tkr1FirstChisq");
   b_Tkr1Hits = fChain->GetBranch("Tkr1Hits");
   b_Tkr1FirstHits = fChain->GetBranch("Tkr1FirstHits");
   b_Tkr1FirstLayer = fChain->GetBranch("Tkr1FirstLayer");
   b_Tkr1DifHits = fChain->GetBranch("Tkr1DifHits");
   b_Tkr1Gaps = fChain->GetBranch("Tkr1Gaps");
   b_Tkr1FirstGaps = fChain->GetBranch("Tkr1FirstGaps");
   b_Tkr1Qual = fChain->GetBranch("Tkr1Qual");
   b_Tkr1Type = fChain->GetBranch("Tkr1Type");
   b_Tkr1TwrEdge = fChain->GetBranch("Tkr1TwrEdge");
   b_Tkr1PrjTwrEdge = fChain->GetBranch("Tkr1PrjTwrEdge");
   b_Tkr1DieEdge = fChain->GetBranch("Tkr1DieEdge");
   b_Tkr1TwrGap = fChain->GetBranch("Tkr1TwrGap");
   b_Tkr1KalEne = fChain->GetBranch("Tkr1KalEne");
   b_Tkr1ConEne = fChain->GetBranch("Tkr1ConEne");
   b_Tkr1KalThetaMS = fChain->GetBranch("Tkr1KalThetaMS");
   b_Tkr1XDir = fChain->GetBranch("Tkr1XDir");
   b_Tkr1YDir = fChain->GetBranch("Tkr1YDir");
   b_Tkr1ZDir = fChain->GetBranch("Tkr1ZDir");
   b_Tkr1Phi = fChain->GetBranch("Tkr1Phi");
   b_Tkr1Theta = fChain->GetBranch("Tkr1Theta");
   b_Tkr1X0 = fChain->GetBranch("Tkr1X0");
   b_Tkr1Y0 = fChain->GetBranch("Tkr1Y0");
   b_Tkr1Z0 = fChain->GetBranch("Tkr1Z0");
   b_Tkr1ThetaErr = fChain->GetBranch("Tkr1ThetaErr");
   b_Tkr1PhiErr = fChain->GetBranch("Tkr1PhiErr");
   b_Tkr1ErrAsym = fChain->GetBranch("Tkr1ErrAsym");
   b_Tkr1CovDet = fChain->GetBranch("Tkr1CovDet");
   b_Tkr1SXX = fChain->GetBranch("Tkr1SXX");
   b_Tkr1SXY = fChain->GetBranch("Tkr1SXY");
   b_Tkr1SYY = fChain->GetBranch("Tkr1SYY");
   b_Tkr1ToTFirst = fChain->GetBranch("Tkr1ToTFirst");
   b_Tkr1ToTAve = fChain->GetBranch("Tkr1ToTAve");
   b_Tkr1ToTTrAve = fChain->GetBranch("Tkr1ToTTrAve");
   b_Tkr1ToTAsym = fChain->GetBranch("Tkr1ToTAsym");
   b_Tkr1ChisqAsym = fChain->GetBranch("Tkr1ChisqAsym");
   b_Tkr1SSDVeto = fChain->GetBranch("Tkr1SSDVeto");
   b_Tkr2Chisq = fChain->GetBranch("Tkr2Chisq");
   b_Tkr2FirstChisq = fChain->GetBranch("Tkr2FirstChisq");
   b_Tkr2Hits = fChain->GetBranch("Tkr2Hits");
   b_Tkr2FirstHits = fChain->GetBranch("Tkr2FirstHits");
   b_Tkr2FirstLayer = fChain->GetBranch("Tkr2FirstLayer");
   b_Tkr2DifHits = fChain->GetBranch("Tkr2DifHits");
   b_Tkr2Gaps = fChain->GetBranch("Tkr2Gaps");
   b_Tkr2FirstGaps = fChain->GetBranch("Tkr2FirstGaps");
   b_Tkr2Qual = fChain->GetBranch("Tkr2Qual");
   b_Tkr2Type = fChain->GetBranch("Tkr2Type");
   b_Tkr2TwrEdge = fChain->GetBranch("Tkr2TwrEdge");
   b_Tkr2PrjTwrEdge = fChain->GetBranch("Tkr2PrjTwrEdge");
   b_Tkr2DieEdge = fChain->GetBranch("Tkr2DieEdge");
   b_Tkr2KalEne = fChain->GetBranch("Tkr2KalEne");
   b_Tkr2ConEne = fChain->GetBranch("Tkr2ConEne");
   b_Tkr2KalThetaMS = fChain->GetBranch("Tkr2KalThetaMS");
   b_Tkr2XDir = fChain->GetBranch("Tkr2XDir");
   b_Tkr2YDir = fChain->GetBranch("Tkr2YDir");
   b_Tkr2ZDir = fChain->GetBranch("Tkr2ZDir");
   b_Tkr2Phi = fChain->GetBranch("Tkr2Phi");
   b_Tkr2X0 = fChain->GetBranch("Tkr2X0");
   b_Tkr2Y0 = fChain->GetBranch("Tkr2Y0");
   b_Tkr2Z0 = fChain->GetBranch("Tkr2Z0");
   b_VtxXDir = fChain->GetBranch("VtxXDir");
   b_VtxYDir = fChain->GetBranch("VtxYDir");
   b_VtxZDir = fChain->GetBranch("VtxZDir");
   b_VtxPhi = fChain->GetBranch("VtxPhi");
   b_VtxTheta = fChain->GetBranch("VtxTheta");
   b_VtxX0 = fChain->GetBranch("VtxX0");
   b_VtxY0 = fChain->GetBranch("VtxY0");
   b_VtxZ0 = fChain->GetBranch("VtxZ0");
   b_VtxAngle = fChain->GetBranch("VtxAngle");
   b_VtxDOCA = fChain->GetBranch("VtxDOCA");
   b_VtxHeadSep = fChain->GetBranch("VtxHeadSep");
   b_VtxS1 = fChain->GetBranch("VtxS1");
   b_VtxS2 = fChain->GetBranch("VtxS2");
   b_VtxDocaWgt = fChain->GetBranch("VtxDocaWgt");
   b_VtxHSWgt = fChain->GetBranch("VtxHSWgt");
   b_VtxS1Wgt = fChain->GetBranch("VtxS1Wgt");
   b_VtxT12Wgt = fChain->GetBranch("VtxT12Wgt");
   b_VtxT2QWgt = fChain->GetBranch("VtxT2QWgt");
   b_VtxTotalWgt = fChain->GetBranch("VtxTotalWgt");
   b_CalEnergySum = fChain->GetBranch("CalEnergySum");
   b_CalEnergyCorr = fChain->GetBranch("CalEnergyCorr");
   b_CalEneSumCorr = fChain->GetBranch("CalEneSumCorr");
   b_Cal_Energy_LLCorr = fChain->GetBranch("Cal_Energy_LLCorr");
   b_CalLeakCorr2 = fChain->GetBranch("CalLeakCorr2");
   b_CalEdgeSumCorr = fChain->GetBranch("CalEdgeSumCorr");
   b_CalTotSumCorr = fChain->GetBranch("CalTotSumCorr");
   b_CalCsIRLn = fChain->GetBranch("CalCsIRLn");
   b_CalTotRLn = fChain->GetBranch("CalTotRLn");
   b_CalCntRLn = fChain->GetBranch("CalCntRLn");
   b_CalDeadTotRat = fChain->GetBranch("CalDeadTotRat");
   b_CalDeadCntRat = fChain->GetBranch("CalDeadCntRat");
   b_CalTPred = fChain->GetBranch("CalTPred");
   b_CalDeltaT = fChain->GetBranch("CalDeltaT");
   b_CalTwrEdge = fChain->GetBranch("CalTwrEdge");
   b_CalLATEdge = fChain->GetBranch("CalLATEdge");
   b_CalTENrm = fChain->GetBranch("CalTENrm");
   b_CalTrackSep = fChain->GetBranch("CalTrackSep");
   b_CalTrackDoca = fChain->GetBranch("CalTrackDoca");
   b_CalTwrGap = fChain->GetBranch("CalTwrGap");
   b_CalELayer0 = fChain->GetBranch("CalELayer0");
   b_CalELayer1 = fChain->GetBranch("CalELayer1");
   b_CalELayer2 = fChain->GetBranch("CalELayer2");
   b_CalELayer3 = fChain->GetBranch("CalELayer3");
   b_CalELayer4 = fChain->GetBranch("CalELayer4");
   b_CalELayer5 = fChain->GetBranch("CalELayer5");
   b_CalELayer6 = fChain->GetBranch("CalELayer6");
   b_CalELayer7 = fChain->GetBranch("CalELayer7");
   b_CalLyr0Ratio = fChain->GetBranch("CalLyr0Ratio");
   b_CalLyr7Ratio = fChain->GetBranch("CalLyr7Ratio");
   b_CalBkHalfRatio = fChain->GetBranch("CalBkHalfRatio");
   b_CalXtalsTrunc = fChain->GetBranch("CalXtalsTrunc");
   b_CalXtalRatio = fChain->GetBranch("CalXtalRatio");
   b_CalLongRms = fChain->GetBranch("CalLongRms");
   b_CalLRmsRatio = fChain->GetBranch("CalLRmsRatio");
   b_CalTransRms = fChain->GetBranch("CalTransRms");
   b_CalMIPDiff = fChain->GetBranch("CalMIPDiff");
   b_CalMIPRatio = fChain->GetBranch("CalMIPRatio");
   b_CalXEcntr = fChain->GetBranch("CalXEcntr");
   b_CalYEcntr = fChain->GetBranch("CalYEcntr");
   b_CalZEcntr = fChain->GetBranch("CalZEcntr");
   b_CalXDir = fChain->GetBranch("CalXDir");
   b_CalYDir = fChain->GetBranch("CalYDir");
   b_CalZDir = fChain->GetBranch("CalZDir");
   b_CalX0 = fChain->GetBranch("CalX0");
   b_CalY0 = fChain->GetBranch("CalY0");
   b_CalZ0 = fChain->GetBranch("CalZ0");
   b_AcdTotalEnergy = fChain->GetBranch("AcdTotalEnergy");
   b_AcdTileCount = fChain->GetBranch("AcdTileCount");
   b_AcdDoca = fChain->GetBranch("AcdDoca");
   b_AcdActiveDist = fChain->GetBranch("AcdActiveDist");
   b_AcdGammaDoca = fChain->GetBranch("AcdGammaDoca");
   b_AcdActDistTop = fChain->GetBranch("AcdActDistTop");
   b_AcdActDistSideRow0 = fChain->GetBranch("AcdActDistSideRow0");
   b_AcdActDistSideRow1 = fChain->GetBranch("AcdActDistSideRow1");
   b_AcdActDistSideRow2 = fChain->GetBranch("AcdActDistSideRow2");
   b_AcdNoSideRow0 = fChain->GetBranch("AcdNoSideRow0");
   b_AcdNoSideRow1 = fChain->GetBranch("AcdNoSideRow1");
   b_AcdNoSideRow2 = fChain->GetBranch("AcdNoSideRow2");
   b_EvtEnergySumOpt = fChain->GetBranch("EvtEnergySumOpt");
   b_EvtEnergyTracker = fChain->GetBranch("EvtEnergyTracker");
   b_EvtEnergyRaw = fChain->GetBranch("EvtEnergyRaw");
   b_EvtMcEnergySigma = fChain->GetBranch("EvtMcEnergySigma");
   b_EvtCalEdgeAngle = fChain->GetBranch("EvtCalEdgeAngle");
   b_EvtTkrEdgeAngle = fChain->GetBranch("EvtTkrEdgeAngle");
   b_EvtLogESum = fChain->GetBranch("EvtLogESum");
   b_EvtTkr1EFrac = fChain->GetBranch("EvtTkr1EFrac");
   b_EvtVtxKin = fChain->GetBranch("EvtVtxKin");
   b_EvtVtxEAngle = fChain->GetBranch("EvtVtxEAngle");
   b_EvtTkrComptonRatio = fChain->GetBranch("EvtTkrComptonRatio");
   b_EvtTkrEComptonRatio = fChain->GetBranch("EvtTkrEComptonRatio");
   b_EvtPSFModel = fChain->GetBranch("EvtPSFModel");
   b_EvtTkr1EChisq = fChain->GetBranch("EvtTkr1EChisq");
   b_EvtTkr1EFirstChisq = fChain->GetBranch("EvtTkr1EFirstChisq");
   b_EvtTkr1EQual = fChain->GetBranch("EvtTkr1EQual");
   b_EvtTkr1PSFMdRat = fChain->GetBranch("EvtTkr1PSFMdRat");
   b_EvtTkr1ECovDet = fChain->GetBranch("EvtTkr1ECovDet");
   b_EvtTkr2EChisq = fChain->GetBranch("EvtTkr2EChisq");
   b_EvtTkr2EFirstChisq = fChain->GetBranch("EvtTkr2EFirstChisq");
   b_EvtTkr2EQual = fChain->GetBranch("EvtTkr2EQual");
   b_EvtCalETLRatio = fChain->GetBranch("EvtCalETLRatio");
   b_EvtCalEXtalRatio = fChain->GetBranch("EvtCalEXtalRatio");
   b_EvtCalEXtalTrunc = fChain->GetBranch("EvtCalEXtalTrunc");
   b_EvtCalETrackDoca = fChain->GetBranch("EvtCalETrackDoca");
   b_EvtCalETrackSep = fChain->GetBranch("EvtCalETrackSep");
   b_EvtVtxEEAngle = fChain->GetBranch("EvtVtxEEAngle");
   b_EvtVtxEDoca = fChain->GetBranch("EvtVtxEDoca");
   b_EvtVtxEHeadSep = fChain->GetBranch("EvtVtxEHeadSep");
   b_IMgoodCalProb = fChain->GetBranch("IMgoodCalProb");
   b_IMvertexProb = fChain->GetBranch("IMvertexProb");
   b_IMcoreProb = fChain->GetBranch("IMcoreProb");
   b_IMpsfErrPred = fChain->GetBranch("IMpsfErrPred");
   b_IMgammaProb = fChain->GetBranch("IMgammaProb");
   return kTRUE;
}

void tree_class::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree_class::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree_class_cxx


