#ifndef RootTreeAnalysis_cxx
#define RootTreeAnalysis_cxx 1

#include "MuonAnalysis.h"
#include <string>

//TObjArray tkrLayerHistArr;
UInt_t digiEventId, reconEventId, mcEventId;
UInt_t digiRunNum, reconRunNum, mcRunNum;
  
Short_t pdl_mpv[768];    
Float_t pdl_width[768];    
Float_t logPoverN_slope[96];
Float_t logPoverN_norm[192];
int range=0;
Float_t cde_h= 21.35;
Float_t cde_w= 27.35;
Float_t cde_l= 326.0;
Float_t cde_r= cde_w/cde_h;
Float_t x_or= -724.76;
Float_t y_or= -728.22;
Float_t z_or= 57;
TCanvas *canvas= new TCanvas("usefull_canvas", "");
TH2F *th2= new TH2F("usefull_th2", "", 8, 0, 7, 12, 0, 11 );
TObjArray *hist_array;


// PDL= 0, CORR_PDL = 1, CORR_PDL_NOMUON=2, LOG_R_OVER_L = 3, GAIN = 4
// MUONS_VS_GAIN=5

//histdefine
  void RootTreeAnalysis::McHistDefine() {
    // Purpose and Method:  Monte Carlo histogram definitions
  }

  void RootTreeAnalysis::DrawHist( TH2F* hist, char* Xname, char* Yname ){
  }

  void RootTreeAnalysis::HistDefine() {
    // Purpose and Method:  Setup Histograms
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    switch (job_id){    
      case 0: //Pdl
        char name[]="Pdl00000";
        histFile = new TFile(m_histFileName,"RECREATE");
        histFile->cd(); 
        hist_array = new TObjArray(768);
        for (int layer=0;layer <8;layer++)  
          for(int col=0;col<12;col++)
            for(int fc=0;fc<2;fc++)
              for(int rg=0;rg<4;rg++){
                int id =rg*192+layer*24+col*2+fc;
                sprintf( name, "Pdl%1d%1d%02d%1d",rg, layer,col, fc);
                hist_array->AddAt(new TH1F(name, "Pedestal in counts per ADC", 
                      1001, 0, 1000), id);
                ((TH1F*) hist_array->At(id))->SetXTitle("ADC");
                ((TH1F*) hist_array->At(id))->SetYTitle("Counts");
              }
        break;

      case 1: //Correlated Pedestals
        char name[]="Pdl_L81_X000X";
        histFile = new TFile(m_histFileName,"UPDATE");
        Float_t mpv[4], width[4];
        TNtuple *tree= (TNtuple *) histFile->Get("Pdl_tree");
        if( !tree ){
          printf( "ERROR: no Pdl_tree\n Run jobid=0 first\n");
          return;
        };
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("MPV_LEX8", 1);
        tree->SetBranchStatus("MPV_LEX1", 1);
        tree->SetBranchStatus("MPV_HEX8", 1);
        tree->SetBranchStatus("MPV_HEX1", 1);
        tree->SetBranchStatus("Width_LEX8", 1);
        tree->SetBranchStatus("Width_LEX1", 1);
        tree->SetBranchStatus("Width_HEX8", 1);
        tree->SetBranchStatus("Width_HEX1", 1);

        tree->SetBranchAddress("MPV_LEX8", mpv);
        tree->SetBranchAddress("MPV_LEX1", mpv+1);
        tree->SetBranchAddress("MPV_HEX8", mpv+2);
        tree->SetBranchAddress("MPV_HEX1", mpv+3);
        tree->SetBranchAddress("Width_LEX8", width);
        tree->SetBranchAddress("Width_LEX1", width+1);
        tree->SetBranchAddress("Width_HEX8", width+2);
        tree->SetBranchAddress("Width_HEX1", width+3);
        int axes[4];      
        
        histFile->cd(); 
        hist_array = new TObjArray(384);
        for (int layer=0;layer <8;layer++)           
          for(int col=0;col<12;col++)
            for(int fc=0;fc<2;fc++){
              sprintf( name, "Pdl_L81_X%1d%02d%01d", layer, col, fc);
              int id= layer*24+col*2+fc;
              tree->GetEntry(id);
              axes[0]= mpv[0]-5*width[0];
              axes[1]= mpv[0]+10*width[0];
              axes[2]= mpv[1]-5*width[1];
              axes[3]= mpv[1]+10*width[1];
              hist_array->AddAt( new TH2F(name, "Pedestal LEX8 versus LEX1", 
                    axes[1]-axes[0], axes[0], axes[1], 
                    axes[3]-axes[2], axes[2], axes[3] ), id );
              ((TH2F*)hist_array->At(id))->SetXTitle("LEX8 (ADC)");
              ((TH2F*)hist_array->At(id))->SetYTitle("LEX1 (ADC)");
              sprintf( name, "Pdl_H81_X%1d%02d%01d", layer, col, fc);
              axes[0]= mpv[2]-5*width[2];
              axes[1]= mpv[2]+10*width[2];
              axes[2]= mpv[3]-5*width[3];
              axes[3]= mpv[3]+10*width[3];
              id+=192;
              hist_array->AddAt( new TH2F(name, "Pedestal HEX8 versus HEX1", 
                    axes[1]-axes[0], axes[0], axes[1], 
                    axes[3]-axes[2], axes[2], axes[3] ), id );
              ((TH2F*)hist_array->At(id))->SetXTitle("HEX8 (ADC)");
              ((TH2F*)hist_array->At(id))->SetYTitle("HEX1 (ADC)");
          }
        break;

      case 2: // Correlated Pedestals, 
              // cut away from muon path(~random trigger)
        char name[]="CleanPdl_L81_X000X";
        histFile = new TFile(m_histFileName,"UPDATE");
        Float_t mpv[4], width[4];
        TNtuple *tree= (TNtuple *) histFile->Get("Pdl_tree");
        if( !tree ){
          printf( "ERROR: no Pdl_tree\n Run jobid=0 first\n");
          return;
        };
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("MPV_LEX8", 1);
        tree->SetBranchStatus("MPV_LEX1", 1);
        tree->SetBranchStatus("MPV_HEX8", 1);
        tree->SetBranchStatus("MPV_HEX1", 1);
        tree->SetBranchStatus("Width_LEX8", 1);
        tree->SetBranchStatus("Width_LEX1", 1);
        tree->SetBranchStatus("Width_HEX8", 1);
        tree->SetBranchStatus("Width_HEX1", 1);

        tree->SetBranchAddress("MPV_LEX8", mpv);
        tree->SetBranchAddress("MPV_LEX1", mpv+1);
        tree->SetBranchAddress("MPV_HEX8", mpv+2);
        tree->SetBranchAddress("MPV_HEX1", mpv+3);
        tree->SetBranchAddress("Width_LEX8", width);
        tree->SetBranchAddress("Width_LEX1", width+1);
        tree->SetBranchAddress("Width_HEX8", width+2);
        tree->SetBranchAddress("Width_HEX1", width+3);
        int axes[4];      
        
        histFile->cd(); 
        hist_array = new TObjArray(384);
        for (int layer=0;layer <8;layer++)           
          for(int col=0;col<12;col++)
            for(int fc=0;fc<2;fc++){
              sprintf( name, "CleanPdl_L81_X%1d%02d%01d", layer, col, fc);
              int id= layer*24+col*2+fc;
              tree->GetEntry(id);
              axes[0]= mpv[0]-5*width[0];
              axes[1]= mpv[0]+5*width[0];
              axes[2]= mpv[1]-5*width[1];
              axes[3]= mpv[1]+5*width[1];
              hist_array->AddAt( new TH2F(name, 
                    "Pedestal LEX8 versus LEX1, cut on muon path", 
                    axes[1]-axes[0], axes[0], axes[1], 
                    axes[3]-axes[2], axes[2], axes[3] ), id );
              ((TH2F*)hist_array->At(id))->SetXTitle("LEX8 (ADC)");
              ((TH2F*)hist_array->At(id))->SetYTitle("LEX1 (ADC)");
              sprintf( name, "CleanPdl_H81_X%1d%02d%01d", layer, col, fc);
              axes[0]= mpv[2]-5*width[2];
              axes[1]= mpv[2]+10*width[2];
              axes[2]= mpv[3]-5*width[3];
              axes[3]= mpv[3]+10*width[3];
              id+=192;
              hist_array->AddAt( new TH2F(name, 
                    "Pedestal LEX8 versus LEX1, cut on muon path", 
                    axes[1]-axes[0], axes[0], axes[1], 
                    axes[3]-axes[2], axes[2], axes[3] ), id );
              ((TH2F*)hist_array->At(id))->SetXTitle("LEX8 (ADC)");
              ((TH2F*)hist_array->At(id))->SetYTitle("LEX1 (ADC)");
          }
        break;

      case 3: //logPoverN
        char name[]="logPoverN0000X";
        histFile = new TFile(m_histFileName,"UPDATE");
        histFile->cd(); 
        hist_array= new TObjArray(96);
        //for(int rg=0;range<4;range++)  
        for(int col=0;col<12;col++)
          for (int layer=0;layer <8;layer++){
            int id= layer*12+col;
            sprintf( name, "logPoverN%1d%1d%02dX", range, layer, col);
            hist_array->AddAt(new TProfile(name, 
                  "Log(#frac{POS}{NEG}) versus position in Cristal", 
                  20, -cde_l*0.6, cde_l*0.6 ), id);
            ((TProfile*)hist_array->At(id))->SetXTitle("Position");
            ((TProfile*)hist_array->At(id))->SetYTitle("Log(#frac{POS}{NEG})");
          }  
      break;
      
      case 4://Gain
        char name[]="Gain_00000";
        char name2[]="non_lin00000";
        histFile = new TFile(m_histFileName,"UPDATE");
        histFile->cd(); 
        hist_array = new TObjArray(384);
        //for(int range=0;range<4;range++) 
          for (int layer=0;layer <8;layer++)           
            for(int col=0;col<12;col++)
                for(int fc=0;fc<2;fc++){
                  int id= layer*24+col*2+fc;
                  sprintf( name, "Gain_%1d%1d%02d%1d", range, layer, col, fc);
                  sprintf( name2, "non_lin%1d%1d%02d%1d", range, layer, col, fc);
                  hist_array->AddAt( new TH1F(name, "Muon Signal", 4096, 0, 4095),
                      id );
                  ((TH1F*)hist_array->At(id))->SetXTitle("Muon (ADC)");
                  ((TH1F*)hist_array->At(id))->SetYTitle("Counts");
                
                }
        break;

      case 5://muons vs position
        char name[]="muons_vs_pos00000";
        histFile = new TFile(m_histFileName,"UPDATE");
        histFile->cd(); 
        hist_array = new TObjArray(192);
        //for(int range=0;range<4;range++) 
          for (int layer=0;layer <8;layer++)           
            for(int col=0;col<12;col++)
                for(int fc=0;fc<2;fc++){
                  int id= layer*24+col*2+fc;
                  sprintf( name, "muons_vs_pos%1d%1d%02d%1d", 
                      range, layer, col, fc);
                  hist_array->AddAt( new TH2F(name, "Muon Signal versus Position",
                        800, 0, 799, 40, -cde_l*0.5, cde_l*0.5 ), id);
                  ((TH2F*)hist_array->At(id))->SetXTitle("Muon (ADC)");
                  ((TH2F*)hist_array->At(id))->SetYTitle("Position");
                }
        break;
      default:
        histFile= new TFile( m_histFileName, "UPDATE");
    }
    gROOT->cd();
  }


//get data
  void RootTreeAnalysis::Pdl( void ) {
    CalDigi *caldigi;
    short int lr, col, fc;
    Short_t rg; 
    Int_t cde_nb;
    
    for( cde_nb=0; caldigi=(CalDigi*) evt->getCalDigiCol()->At(cde_nb); 
                                                              cde_nb++ ){
      CalXtalId packedid= caldigi->getPackedId();
      lr=packedid.getLayer();
      col=packedid.getColumn();
      for( fc=0; fc<2; fc++){
        for( rg=0; rg<4; rg++){
          Short_t value= caldigi->getAdcSelectedRange( rg, fc );
          ((TH1F*)hist_array->At(rg*192+lr*24+col*2+fc))->Fill( value );
        }
      }
    }
  }
    
  void RootTreeAnalysis::Pdl_Corr_NoMuon( void ){
    Short_t *pdl;
    CalDigi *caldigi;
    short int lr, col, fc, rg, columns[9]={0,0,0,0,0,0,0,0,0};
    Int_t cde_nb;
    Float_t *pdl_sigma;
    Int_t XtalAdc[4][8][12][2];
    for(rg=0; rg<4; rg++)
      for(lr=0; lr<8; lr++)
        for(col=0; col<8; col++)
          for(fc=0; fc<2; fc++)
            XtalAdc[rg][lr][col][fc]=0;
    
   
    //cut: +- 1 on the sides of +5 sigma pulses
    for( cde_nb=0; caldigi=(CalDigi*) evt->getCalDigiCol()->At(cde_nb); 
                                                              cde_nb++ ){
      CalXtalId packedid= caldigi->getPackedId();    
      lr= packedid.getLayer();
      col= packedid.getColumn();
      pdl= pdl_mpv + lr*24 + col*2;
      pdl_sigma= pdl_width + lr*24 + col*2;
      
      for( fc=0; fc<2; fc++ , pdl_sigma++, pdl++ ){
        if( XtalAdc[0][lr][col][fc]==-1) continue;

        XtalAdc[0][lr][col][fc]= caldigi->getAdcSelectedRange( 0, fc);
        if  ( XtalAdc[0][lr][col][0] - *(pdl) >  (*pdl_sigma)  *5){          
          //neigbours cut
          for( rg=(col==0)?0:col-1; rg<=col+1 && rg<12; rg++)
            XtalAdc[0][lr][rg][0]= -1;
          //8 hits cut
          columns[lr]++;
          columns[8]+=(columns[lr]==1); 
          fc=2;
        }
        else 
          for ( rg=1; rg<4; rg++, pdl_sigma+= 190, pdl+= 190 )
            XtalAdc[rg][lr][col][fc]= caldigi->getAdcSelectedRange( rg ,fc);
      }
    }
    //8 hits cut
    if (columns[8]<8) return;


    //fill hists
    for(lr=0; lr<8; lr++)
      for(col=0; col<12; col++){
        for(fc=0; fc<2 && XtalAdc[0][lr][col][fc]!=-1; fc++){
          ((TH2F*)hist_array->At(lr*24+col*2+fc))->Fill( 
                                                  XtalAdc[0][lr][col][fc],
                                                  XtalAdc[1][lr][col][fc] );
          ((TH2F*)hist_array->At(192+lr*24+col*2+fc))->Fill( 
                                                  XtalAdc[2][lr][col][fc],
                                                  XtalAdc[3][lr][col][fc] );
        }
      }
  }

  void RootTreeAnalysis::Pdl_Correlated( void ){
      CalDigi *caldigi;
    Float_t enem, enep;
    short int histid, fc, cur_col;
    Int_t cde_nb;

    for( cde_nb=0; caldigi=(CalDigi*) evt->getCalDigiCol()->At( cde_nb ); 
         cde_nb++){
      CalXtalId packedid= caldigi->getPackedId();
      histid= packedid.getLayer()*24+packedid.getColumn()*2;
      for( short int fc=0; fc<2; fc++ )
        for( short int rg=0; rg<4; rg+=2 )
          ((TH2F*) hist_array->At(histid+fc+rg*96))->Fill( 
                                       caldigi->getAdcSelectedRange(rg, fc),
                                       caldigi->getAdcSelectedRange(rg+1,fc) );
    }
  }

  Bool_t RootTreeAnalysis::Make_Cut( Double_t *position,  Int_t *XtalAdc, 
                                                          Short_t *col ) {
    //CalDigi: make cut
		Float_t *pdl_sigma;
    Short_t *pdl;
    CalDigi *caldigi;
    short int lr, fc, cur_col;
    Int_t cur_XtalAdc[2];
    Short_t nhits=8;
    Int_t cde_nb;

    //cut for 1 hit for every layer and no more
    col[0]=col[1]=col[2]=col[3]=-1;
    col[4]=col[5]=col[6]=col[7]=-1;
    for( cde_nb=0; caldigi=(CalDigi*) evt->getCalDigiCol()->At( cde_nb ); 
                                                                    cde_nb++){
      //get Xtal charateristics
      CalXtalId packedid= caldigi->getPackedId();
      lr= packedid.getLayer();
      cur_col= packedid.getColumn();
      pdl= pdl_mpv + range*192 + lr*24 + cur_col*2;
      pdl_sigma= pdl_width + range*192 + lr*24 + cur_col*2;
       //get cur_XtalAdc
      cur_XtalAdc[0]= caldigi->getAdcSelectedRange( range, 0) - pdl[0];
      cur_XtalAdc[1]= caldigi->getAdcSelectedRange( range, 1) - pdl[1];
      //fill XtalAdc and col
      if  (    (cur_XtalAdc[0] < (pdl_sigma[0]*5))
            || (cur_XtalAdc[1] < (pdl_sigma[1]*5)) ) continue;
      if ( col[lr]!=-1 ) return kFALSE;// 1 per layer cut 
      XtalAdc[lr*2]=    cur_XtalAdc[0];
      XtalAdc[lr*2+1]=  cur_XtalAdc[1];
      col[lr]= cur_col;
      nhits--;
    }
    if ( nhits!=0 ) return kFALSE;//8 in all cut
    
    
    //Recon: if recon exists get direction given by tracker
    if (rec){
			TkrRecon *tkrRec = rec->getTkrRecon();
      if (!tkrRec)  return;
      short int ii;  
      Float_t added[2];
      TVector3 tvect3;
      const TObjArray *vertexCol = tkrRec->getVertexCol();
      TkrVertex *vert = (TkrVertex*)vertexCol->At(0);
      

      if(!vertexCol->GetEntries()) return kFALSE;

			tvect3= vert->getDirection();
			added[0]= -tvect3.X()/tvect3.Z();
			added[1]=  -tvect3.Y()/tvect3.Z();
			position[16]= -tvect3.CosTheta();
			if( position[16]==0 ) return kFALSE;

			tvect3 = vert->getPosition();
			position[0]= tvect3.X() + ( z_or+tvect3.Z())* added[0];
			position[1]= tvect3.Y() + ( z_or+tvect3.Z())* added[1];
			added[0]*=  cde_h; added[1]*=  cde_h;

			for( lr=1; lr<8; lr++){
				fc= lr%2;
				position[ lr*2]=  position[ lr*2-1] + added[(fc==0)];
				position[ lr*2+1]=  position[ lr*2-2] + added[fc];
			}
      return kTRUE;
    }

    // otherwise get direction, position from cal
		Double_t trajectory[2][2];
    TGraphErrors *tpos;
    TF1 *tline= new TF1( "tline", "[0]*x+[1]", 0, 8 );
    canvas->cd();
    th2->Draw();

    //get direction
    for( fc=0, position[16]=0; fc<2; fc++){ 
    //  fill TGraphErrors
      tpos= new TGraphErrors(); 
      for( lr=0; lr<8; lr+=2 ){
        tpos->SetPoint( lr/2, lr, col[lr+fc]);
        tpos->SetPointError( lr/2, 2, 2*cde_r);
      }
    //  get trajectory parameters
      tpos->Draw("*");
      tpos->Fit( tline, "WNQ" );
      trajectory[1-fc][0]= tline->GetParameter(0);
      trajectory[1-fc][1]= tline->GetParameter(1);
      position[16]+= (trajectory[1-fc][0])**2;
      tpos->~TGraphErrors();
    }
    position[16]= 1/sqrt(1 + cde_r**2 * position[16] );// position[16]=cos(Theta)

    //get position
    position[0]= (trajectory[0][1]-5.5)*cde_w/12;
    position[1]= (trajectory[1][1]-5.5)*cde_w/12;
    if( position[0]<-cde_l/2 || position[0]>cde_l/2 ) trajectory[1][0]=-6;
    if( position[1]<-cde_l/2 || position[1]>cde_l/2) trajectory[0][0]=-6;
    for( lr=1; lr<8; lr++ ){
      position[lr*2]= position[lr*2-1] + trajectory[lr%2][0]*cde_w/12;
      position[lr*2+1]= position[(lr-1)*2] +  trajectory[lr%2==0][0]*cde_w/12;
      if( position[lr*2+1]<-cde_l/2 || position[lr*2]<-cde_l/2 ) break;
      if( position[lr*2+1]>cde_l/2 || position[lr*2]>cde_l/2 ) break;
    }

    tline->~TF1();
    canvas->Clear("D");
    caldigi=0;

    return (lr==8); // lr==8 if and only if positions have been correctly rebuilt
  }

  void RootTreeAnalysis::Make_Hist( Double_t *position, Int_t *XtalAdc,
                                                          Short_t *col  ){
    short int lr, fc;
    int histid;
    Float_t Xtalval[2], value_raw;
    Float_t slcoef, value;

    //fill hists
    histFile->cd();
      switch( job_id ){
        case 3://fill logPoverN
          for( lr=0; lr<8; lr++){
            value= XtalAdc[2*lr];
            value/= XtalAdc[2*lr+1];
            value= log(value);
            histid= lr*12 + col[lr];
            ((TH2F*) hist_array->At(histid))->Fill( position[2*lr], value );
          }
          break;; 
          
        case 4://fill muon hist
          for( lr=0; lr<8; lr++ ){
            histid= lr*12 + col[lr];
            slcoef = logPoverN_slope[histid]*(position[2*lr]-5.5)/2;
            slcoef *= position[16]; //position[16]= cos(theta_direction)
            for( fc=0, histid *= 2; fc<2; fc++, histid++ ){
              value= XtalAdc[lr*2+fc];
              value+=  (fc==0)?   -value*slcoef
                              :   value*slcoef;
              ((TH1F*) hist_array->At(histid))->Fill( value );
            }
          }
        
          break;
        
        case 5://fill muons vs position
          if ( (position[16]<0.98) && (position[16]>-0.98) ) return;
          for( lr=0; lr<8; lr++ ){
            histid= lr*12 + col[lr];
            for( fc=0, histid *= 2; fc<2; fc++, histid++){
              value= XtalAdc[lr*2+fc]*position[16];
              ((TH2F*) hist_array->At( histid ))->Fill( value,
									position[2*lr] );
            }
          }
          break;
      }   
    gROOT->cd();
  }


//treat data
  void RootTreeAnalysis::FitPdl( void ){
		float value, mean, rms;
    Float_t results_values[12];
    short int ii, col, lr, fc, rg;
    int length, iii;
    char name[]="Pdl00000";
    TH1F* th;
    TF1* gaus;
    histFile = new TFile(m_histFileName, "UPDATE");
    histFile->cd();
      TNtuple *results= new TNtuple("Pdl_tree", "Uncorrelated pedestals",
          "cde_id:MPV_LEX8:MPV_LEX1:MPV_HEX8:MPV_HEX1:Width_LEX8:Width_LEX1:Width_HEX8:Width_HEX1");
    gROOT->cd();
    
    for( lr=0; lr<8; lr++ )
      for( col=0; col<12; col++ )
        for( fc=0; fc<2; fc++ ){
          results_values[0]= lr*1000 + col*10 + fc;
          for( rg=0; rg<4; rg++ ){
            //get th1
            sprintf( name, "Pdl%1d%1d%02d%1d", rg, lr, col, fc);
            th= (TH1F*) (TH1F*)histFile->Get( name );
            if (!th){
              printf("ERROR NO TH\n cde_id:\n");
              printf("layer: %d, column: %d, face: %d, range: %d\n", 
                  lr, col, fc, rg);
              return;
            }
            
            
            //fit th1
            for ( ii=0; ii<3; ii++){
              mean= th->GetMean();
              rms= th->GetRMS();
              th->SetAxisRange(mean-2.5*rms, mean+2.5*rms);
            }
            th->Fit("gaus", "", "", th->GetMean()-2.5*th->GetRMS(),  
                th->GetMean()+ 2.5*th->GetRMS() );
            
            gaus= (TF1*)(th->GetFunction("gaus"));
            results_values[rg+1]= gaus->GetParameter(1);
            results_values[rg+5]=gaus->GetParameter(2);
          }        
          results->Fill(results_values);
        }

    histFile->Write( "", TObject::kOverwrite );
    histFile->Close();
    gROOT->cd();
  }

  void RootTreeAnalysis::FitCorrPdl( char* option ){
      int id;
    TF1 *fits= new TF1( "line_fit", "[0]+[1]*x", 0, 800);
    TProfile *axes;
    Double_t fit_ranges[4];
    char name[]="CleanPdl_H81_X0000";
    Float_t results[17];
    TH2F* pdl;
    TF2 *gauss= new TF2("gauss", 
        "[0]*exp( -0.5/([1]*[1])*((x-[4])*cos([3])+(y-[5])*sin([3]))**2 ) * exp( -0.5/([2]*[2])*((y-[5])*cos([3])-(x-[4])*sin([3]))**2)" );
    
    gauss->SetParName( 0, "area" );
    gauss->SetParName( 1, "sigma_X" );
    gauss->SetParName( 2, "sigma_Y" );
    gauss->SetParName( 3, "angle" );
    gauss->SetParName( 4, "X_MPV" );
    gauss->SetParName( 5, "Y_MPV" );
    TH1D *proj;
    histFile= new TFile( m_histFileName, "UPDATE");
    histFile->cd();
			TNtuple *tree= new TNtuple("Pdl_tree", "Correlated Pedestal Characteristics", 
      "cde_id:MPV_LEX8:MPV_LEX1:MPV_HEX8:MPV_HEX1:Width_LEX8:Width_LEX1:Width_HEX8:Width_HEX1:Std_Dev_LEX:Std_Dev_HEX:Orth_Dev_LEX:Orth_Dev_HEX:angle_LEX:angle_HEX:Gain_LEX:Gain_HEX");
    gROOT->cd();
    gStyle->SetOptFit(1111);

    for( short int lr=0; lr<8; lr++ )
      for( short int col=0; col<12; col++ )
        for( short int fc=0; fc<2; fc++ ){
          results[0]= lr*1000+col*10+fc;
          for( short int rg=0; rg<4; rg++ ){ 
            sprintf( name, "%sPdl_%c81_X%01d%02d%01d", 
                option, (Char_t ) (76-4*rg), lr, col, fc );
            pdl= (TH2F*) histFile->Get(name);
            if ( !pdl ) continue;

            axes= pdl->ProfileX( "prof", -1, -1);
            for( int jj=0; jj<3; jj++ ){
              fit_ranges[0]=   pdl->GetMean(2)-pdl->GetRMS(2)*3;
              fit_ranges[1]=   pdl->GetMean(2)+pdl->GetRMS(2)*3;
              fit_ranges[2]=   pdl->GetMean(1)-pdl->GetRMS(1)*3;
              fit_ranges[3]=   pdl->GetMean(1)+pdl->GetRMS(1)*3;
              pdl->GetXaxis()->SetRangeUser( fit_ranges[2], fit_ranges[3] );
              pdl->GetYaxis()->SetRangeUser( fit_ranges[0], fit_ranges[1] );
            }

            fits->SetRange( fit_ranges[2], fit_ranges[3] );
            axes->Fit( fits, "QR" );

            //gausss
            gauss->SetParameter( "area", 1.5*pdl->Integral(
                  fit_ranges[2], fit_ranges[3], fit_ranges[0], fit_ranges[1]) );
            gauss->SetParameter( "sigma_X", pdl->GetRMS(1) );
            gauss->SetParameter( "sigma_Y", pdl->GetRMS(2) );
            gauss->SetParameter( "angle", atan(fits->GetParameter(1)) );
            gauss->SetParameter( "X_MPV", pdl->GetMean(1) );
            gauss->SetParameter( "Y_MPV", pdl->GetMean(2) );
            gauss->SetRange( fit_ranges[2], fit_ranges[3], 
                fit_ranges[0], fit_ranges[1] );

            gauss->SetParLimits( 0, 1, pdl->GetEntries() );
            gauss->SetParLimits( 1, gauss->GetParameter("sigma_X")/20, 
                3*(fit_ranges[3]-fit_ranges[2]) );
            gauss->SetParLimits( 2, gauss->GetParameter("sigma_Y")/20, 
                3*(fit_ranges[1]-fit_ranges[0]) );
            gauss->SetParLimits( 3, -acos(-1), acos(-1) );
            gauss->SetParLimits( 4, fit_ranges[2], fit_ranges[3] );
            gauss->SetParLimits( 5, fit_ranges[0], fit_ranges[1] );

            pdl->Fit( gauss, "" );
          
            //save results
            results[1+2*rg]= gauss->GetParameter("X_MPV");
            results[2+2*rg]= gauss->GetParameter("Y_MPV");
            results[9+rg]= gauss->GetParameter("sigma_X");
            results[11+rg]= gauss->GetParameter("sigma_Y");
            results[13+rg]= gauss->GetParameter("angle");

            proj= pdl->ProjectionX("proj", -1, -1);
            proj->SetAxisRange(fit_ranges[2], fit_ranges[3]);
            proj->Fit("gaus", "Q");
            results[5+2*rg]= ((TF1*) proj->GetFunction("gaus"))->GetParameter(2);
            proj->~TH1D();

            proj= pdl->ProjectionY("proj", -1, -1);
            proj->SetAxisRange(fit_ranges[0], fit_ranges[1]);
            proj->Fit("gaus", "Q");
            results[6+2*rg]= ((TF1*) proj->GetFunction("gaus"))->GetParameter(2);
            proj->~TH1D();

            fits->SetRange( fit_ranges[3], 800 );
            axes->Fit( fits, "QR" );
            results[15+rg]= fits->GetParameter(1);
              
            if( results[9+rg]<results[11+rg]){
              Float_t val= results[11+rg]; 
              results[11+rg]= results[9+rg]; 
              results[9+rg]= val;
            }
            for( ; results[13+rg]<0; results[13+rg]+=acos(-1)/2 );
            for( ; results[13+rg]>acos(-1)/2; results[13+rg]-=acos(-1)/2 );

            if( results[9+rg]==0 ) results[11+rg]= results[9+rg]= -1000;
            axes->~TProfile();
          }

          tree->Fill(results);
        }
    histFile->Write("", TObject::kOverwrite );
    histFile->Close();
    return;
  }

  void RootTreeAnalysis::FitLogPoverN( void ){
		int histid;
    short int lr, col;
    TProfile* tp;
    TF1 *func;
    Double_t logratio, slope;
    Float_t results_values[4];
    char name[200];
    histFile= new TFile( m_histFileName, "UPDATE" );
    sprintf( name, "LogPoverN_tree" );
    histFile->cd();
      TNtuple *results= new TNtuple( name, "fit values of Log(Left/Rigth)", 
          "cde_id:norm_POS:norm_NEG:slope");
		gROOT->cd();
    
      
    for( lr=0; lr<8; lr++)
      for( col=0;col<12;col++){
        //get tprofile
        histid = lr*12+col;      
        sprintf( name, "logPoverN%1d%1d%02dX", range, lr, col);
        tp=  (TProfile*) histFile->Get( name );
        if (!tp){
          printf("ERROR NO TProfile\n cde_id:\n");
          printf("layer: %d, column: %d, face: %d, range: %d\n", 
                lr, col, fc, rg);
          return;
        }
        

        //fit tprofile
        tp->SetAxisRange(-cde_l*0.3,cde_l*0.3);
        tp->Fit("pol1","Q");
        func= (TF1*) tp->GetFunction("pol1");
        logratio= func->GetParameter( 0 );
        slope= func->GetParameter( 1 );
      
        //write results
        results_values[0]= lr*1000 + col*10;
        results_values[1]= exp(-logratio/2);
        results_values[2]= exp(logratio/2);
        results_values[3]= slope;
        results->Fill( results_values );
        tp->SetAxisRange(-7, 7);
     }

    histFile->Write( "", TObject::kOverwrite );
    gROOT->cd();
    histFile->Close();
  }

  void RootTreeAnalysis::FitGain( int rebinning  ){
    int histid;
    short int lr, col, fc, ii, ll;
    TH1 *th1;
    TF1 *landau, *convoluted;
    Double_t constant, sigma, mpv, range_val[2];
    Float_t results_values[3];
    char name[]="Gain_00000";
    char newname[]="Gain_00000_fitted";
    canvas->cd();
    histFile= new TFile( m_histFileName, "UPDATE" );
    histFile->cd();
			hist_array= new TObjArray(193);
			hist_array->AddAt( new TNtuple("Gain_tree", 
			"pedestal substracted, attenuation and path correct muon signal in ADC units",
			"cde_id:MPV_LEX8:Width_LEX8"), 0);
		gROOT->cd();
    
    for( lr=0; lr<8; lr++)
      for( col=0;col<12;col++)
        for( fc=0;fc<2;fc++){
          //get th1
          histid = lr*12+col;      
          sprintf( name, "Gain_%1d%1d%02d%1d", range, lr, col, fc);
          sprintf( newname, "Gain_%1d%1d%02d%1d_fitted", range, lr, col, fc);
          hist_array->AddAt((TH1F*) ((TH1F*) histFile->Get( name ))->Clone( newname ), 
              lr*24+col*2+fc+1);
          th1= (TH1F*) hist_array->At( lr*24+col*2+fc+1);
          if (!th1){
            printf("ERROR NO TH1\n cde_id:\n");
            printf("layer: %d, column: %d, face: %d, range: %d\n", 
                lr, col, fc, rg);
            return;
          }
          

          //fit th1
          th1->Rebin( rebinning );
          for( ll=0; ll<5; ll++ )
            th1->SetAxisRange(   th1->GetMean() - 2* th1->GetRMS(),
                                            th1->GetMean() + 3* th1->GetRMS() );
          range_val[0]=   th1->GetMean() - 2* th1->GetRMS();
          range_val[1]=   th1->GetMean() + 3* th1->GetRMS();

          printf("Landau fit: \n");;
          //  fit landau
          //      set tf1 names
          sprintf( name, "muons_ldau%1d%1d%02d%1d", range, lr, col, fc);
          landau= new TF1( name, "landau",  th1->GetMean()-2*th1->GetRMS(), 
																						th1->GetMean()+3*th1->GetRMS());
          landau->SetParNames("Constant", "MPV", "Sigma");
          //      set tf1  start values
          th1->SetAxisRange( th1->GetMean()-th1->GetRMS()/2, th1->GetMean());
          mpv=  th1->GetMean();
          constant=  th1->GetEntries();
          th1->SetAxisRange( 0, 4095 );
          sigma=  th1->GetRMS()/5;
          landau->SetParameter(0, constant );
          landau->SetParameter( 1, mpv );
          landau->SetParameter( 2, sigma );
          landau->SetRange( range_val[0], range_val[1] );
          //      set tf1 parameter limits
          landau->SetParLimits( 2, sigma/10, range_val[1]-range_val[0] );
          landau->SetParLimits( 1, range_val[0], range_val[1] );
          landau->SetParLimits( 0, constant/100, constant );
          //    fit
          th1->Fit( landau, "R" );
          
          landau= (TF1*) th1->GetFunction( name );
          mpv= landau->GetParameter(1); 
          sigma= landau->GetParameter(2);
          constant= 50* landau->GetParameter(0);
          

          printf("Gauss-convoluted Landau fit: \n");;
          //  fit convoluted
          sprintf( name, "muons_LGau%1d%1d%02d%1d", range, lr, col, fc);
          //      set tf1 names
          convoluted= new TF1( name, langaufun, 
																						th1->GetMean()-2*th1->GetRMS(), 
																						th1->GetMean()+3*th1->GetRMS(), 4);
          convoluted->SetParNames("Width","MPV","Area","GSigma");
          convoluted->SetLineColor(2);
          //      set tf1  start values
          convoluted->SetParameter( 0, sigma);
          convoluted->SetParameter( 1, mpv);
          convoluted->SetParameter( 2, constant );
          convoluted->SetParameter( 3, 70 );
          convoluted->SetRange( range_val[0], range_val[1] );
          //      set tf1 parameter limits
          convoluted->SetParLimits( 0, 1, 50 );
          convoluted->SetParLimits( 1, range_val[0], range_val[1] );
          convoluted->SetParLimits( 2, constant/5, constant *5 );
          convoluted->SetParLimits( 3, 20, 150 );
          
          //    fit
          th1->Fit( convoluted,"R+");
          convoluted= (TF1*) th1->GetFunction( name );
          
          

          results_values[0]= lr*1000 + col*10 + fc;
          results_values[1]= (Float_t) get_MPV( convoluted, 
                                                th1->GetMean()-th1->GetRMS()/2,
                                                th1->GetMean()+th1->GetRMS()/2,
                                      0.000001 );
          results_values[2]= th1->GetRMS();
          ((TNtuple*) hist_array->At(0))->Fill( results_values );
        }

    histFile->Write( "", TObject::kOverwrite );
    gROOT->cd();
    histFile->Close();
  }  
  
  void RootTreeAnalysis::FitTapering( void ){
    short int lr, col, fc;
    char name[]="muons_vs_pos00000";
    Double_t results[3];
    TH1D *proj;
    TH2F *muons;
    TProfile *gr;
    TF1 *line= new TF1("line", "[0]+[1]/40*x", 0, 40 );
    histFile= new TFile( m_histFileName, "UPDATE");
    try_file->cd();
      TProfile *pr[192];
      TNtuple *tree= new TNtuple("Taper_tree", "", "cde_id:slope:norm");
      for( lr=0; lr<8; lr++ )
        for( col=0; col<12; col++ )
          for( fc=0; fc<2; fc++ ){
          sprintf( name, "taper0%01d%02d%01d", lr, col, fc);
          pr[xx*192+lr*24+col*2+fc]= new TProfile( name, name, 15, 0, 1 );
        }
    gROOT->cd();

    for( lr=0; lr<8; lr++ )
      for( col=0; col<12; col++ )
        for( fc=0; fc<2; fc++ ){
          results[0]= lr*1000+ col*10+fc;
          sprintf( name, "taper0%01d%02d%01d", lr, col, fc);
          muons= (TH2F*) rec_file->Get( name );
          gr= pr[lr*24+col*2+fc];
          for( xx=0; xx<40; xx++ ){
            proj= muons->ProjectionX( "proj", xx, xx );
            if ( proj->GetEntries()<5 )
              gr->Fill(xx/40.,  proj->GetMaximumBin());
            else{
              proj->Fit("landau", "Q");
              gr->Fill(xx/40., ((TF1*) proj->GetFunction("landau"))->GetParameter(1));
            }
            proj->~TH1D();
          }
          gr->Fit("expo", "", "", .2, .8 );
          results[1]= ((TF1*) gr->GetFunction("expo"))->GetParameter(0);
          results[2]= ((TF1*) gr->GetFunction("expo"))->GetParameter(1);

          tree->Fill(results);
        }
    histFile->Write("", TObject::kOverwrite );
    histFile->Close();
  }        

  Double_t RootTreeAnalysis::get_MPV( TF1* function, 
                     Double_t min, Double_t max,  Double_t precision ){
      Double_t steps, step_interval, max_range, min_range, value, old_value;
    

    //init variables
    step_interval= (max_range -steps)/10;
    max_range= max;
    min_range= min;
    if (min_range<0) min_range=0;
    old_value= function->Eval( min_range );
    

    //get peak
    for( step_interval= (max-min)/3; step_interval > precision/2; ){
      old_value= function->Eval( min_range );
      for(  steps= min_range+step_interval; steps<= max_range+1; 
                                          steps+= step_interval ){
        value= function->Eval( steps );
        if ( old_value>=value ) {
          max_range= steps;
          if (min_range < steps-2*step_interval) min_range= steps-2*step_interval;
          step_interval= (max_range-min_range)/3;
          break;
        }
        old_value= value;
      }
    }

    return (old_value>value)? steps-step_interval: steps;
  }

  Double_t langaufun(Double_t* x, Double_t *par){
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation), 
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

    // variables
    // Numeric constants
    Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    Double_t mpshift  = -0.22278298;       // Landau maximum location
    
    // Control constants
    Double_t np = 100.0;      // number of convolution steps
    Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0]; 

    // Range of convolution integral
    xlow =*x - sc * par[3];
    xupp = *x + sc * par[3];

    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx, mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(*x,xx,par[3]);

      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(*x,xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
  }


//retrieve data
  Bool_t RootTreeAnalysis::getPdl( Short_t *Pdlmpv, Float_t *Pdlsigma ){
      TNtuple *tree= (TNtuple*) histFile->Get("Pdl_tree");
    if (!tree){
      printf("ERROR: no Pdl_tree in file:\n\t%s\n", m_histFileName );
      return kFALSE;
    }
    Float_t mpv[4], width[4], cde_id;
    Short_t lr, col, rg, entry;
    Int_t fc; 
    
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("MPV_LEX8", 1);
    tree->SetBranchStatus("MPV_LEX1", 1);
    tree->SetBranchStatus("MPV_HEX8", 1);
    tree->SetBranchStatus("MPV_HEX1", 1);
    tree->SetBranchStatus("Width_LEX8", 1);
    tree->SetBranchStatus("Width_LEX1", 1);
    tree->SetBranchStatus("Width_HEX8", 1);
    tree->SetBranchStatus("Width_HEX1", 1);
    tree->SetBranchStatus("cde_id", 1);

    tree->SetBranchAddress("MPV_LEX8", mpv);
    tree->SetBranchAddress("MPV_LEX1", mpv+1);
    tree->SetBranchAddress("MPV_HEX8", mpv+2);
    tree->SetBranchAddress("MPV_HEX1", mpv+3);
    tree->SetBranchAddress("Width_LEX8", width);
    tree->SetBranchAddress("Width_LEX1", width+1);
    tree->SetBranchAddress("Width_HEX8", width+2);
    tree->SetBranchAddress("Width_HEX1", width+3);
    tree->SetBranchAddress("cde_id", &cde_id);
    
    for ( entry=0; tree->GetEntry( entry ); entry++){
      fc= cde_id; 
      lr= fc/1000; fc-= lr*1000;
      col= fc/10; fc-= col*10;
      for(rg = 0; rg<4; rg++ ){ 
        *(Pdlmpv + rg*192 + lr*24 + col*2 + fc)= ( Short_t )mpv[rg];
        *(Pdlsigma + rg*192 + lr*24 + col*2 + fc )= width[rg];
      }
    }
		return kTRUE;
  }

  Bool_t RootTreeAnalysis::getLogPoverN( Float_t *logPoverN_s ){
		Float_t norm_POS, norm_NEG, slope, cde_id;
    Short_t lr, rg, entry;
    Int_t col;
    TNtuple *tree= (TNtuple*) histFile->Get("Taper_tree");
    Bool_t taper=(tree==0);
		if (taper){
      printf("no Taper_tree in file:\n\t%s\n", m_histFileName );
			tree= (TNtuple*) histFile->Get("LogPoverN_tree");
			if (!tree){
				printf("ERROR: no LogPoverN_tree in file:\n\t%s\n", m_histFileName );
				return kFALSE;
			}
		}
    tree->SetBranchStatus( "*", 0);
    tree->SetBranchStatus( "slope", 1);
    tree->SetBranchAddress( "slope", &slope );
    tree->SetBranchStatus("cde_id", 1);
    tree->SetBranchAddress("cde_id", &cde_id);
    
    for( entry=0; tree->GetEntry(entry); entry++){
      col= cde_id;
      rg= col/10000; col-= rg*10000;
      lr= col/1000; col-= lr*1000;
      col/= 10;
      logPoverN_s[rg*192 + lr*24 + col*2+0]= slope;
			if( !taper ) tree->GetEntry(entry);
			logPoverN_s[rg*192 + lr*24 + col*2+1]= slope;
    }
		return kTRUE;
  }

  
//write to xml
  void RootTreeAnalysis::PdlToXml( char* filename ){
      std::ofstream xml( filename );
    histFile= new TFile( m_histFileName );
    Float_t mpv[4], width[4];
    TNtuple *tree;
      tree= (TNtuple*) histFile->Get("Pdl_tree");
    if (!(tree)){
      printf("ERROR:NO PEDESTALS TREE\n");
      return;
    }
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("MPV_LEX8", 1);
    tree->SetBranchStatus("MPV_LEX1", 1);
    tree->SetBranchStatus("MPV_HEX8", 1);
    tree->SetBranchStatus("MPV_HEX1", 1);
    tree->SetBranchStatus("Width_LEX8", 1);
    tree->SetBranchStatus("Width_LEX1", 1);
    tree->SetBranchStatus("Width_HEX8", 1);
    tree->SetBranchStatus("Width_HEX1", 1);

    tree->SetBranchAddress("MPV_LEX8", mpv);
    tree->SetBranchAddress("MPV_LEX1", mpv+1);
    tree->SetBranchAddress("MPV_HEX8", mpv+2);
    tree->SetBranchAddress("MPV_HEX1", mpv+3);
    tree->SetBranchAddress("Width_LEX8", width);
    tree->SetBranchAddress("Width_LEX1", width+1);
    tree->SetBranchAddress("Width_HEX8", width+2);
    tree->SetBranchAddress("Width_HEX1", width+3);
    short int rg, layer, column, face;
    
    //header
    xml<<"<?xml version=\"1.0\" ?>"<<endl<<endl;
    xml<<"<!-- $Header: my Pdl-->"<<endl<<endl;
    xml<<"<!-- Approximately real  CAL ped XML file for EM, according to calCalib.dtd -->"<<endl<<endl;

    xml<<"<!DOCTYPE calCalib SYSTEM \"../calCalib_v3.dtd\" [] >"<<endl<<endl;

    //calCalib
    xml<<"<calCalib>"<<endl;
    xml<<"\t<generic\tinstrument=\"EM\" timestamp=\"2003-11-2-12:56\""<<endl;
    xml<<"\t\t\t\tcalibType=\"CAL_Ped\" fmtVersion=\"v2r0\">"<<endl<<endl;
    xml<<"\t\t<inputSample\tstartTime=\"2003-2-21-05:49:12\" stopTime=\"2003-2-24-07:07:02\""<<endl; 
    xml<<"\t\t\t\t\ttriggers=\"random\" mode=\"normal\" source=\"stuff\" >"<<endl<<endl;
    
    xml<<"\t\tTimes are start and stop time of calibration run."<<endl;
    xml<<"\t\tOther attributes are just made up for code testing."<<endl;    
    
    xml<<"\t\t</inputSample>"<<endl;
    xml<<"\t</generic>"<<endl<<endl<<endl;
    
    xml<<"<!-- EM instrument: 8 layers, 12 columns -->"<<endl<<endl;
    xml<<"\t<dimension nRow=\"1\" nCol=\"1\" nLayer=\"8\" nXtal=\"12\" nFace=\"2\" />"<<endl<<endl;
    
    xml<<"\t<tower iRow=\"0\" iCol=\"0\">"<<endl;
    
    for( layer=0; layer<8; layer++ ){
      xml<<"\t\t<layer iLayer=\""<<layer<<"\">"<<endl;
      
      for( column=0; column<12; column++ ){
        xml<<"\t\t\t<xtal iXtal=\""<<column<<"\">"<<endl;
        
        for( face=1; face>-1; face-- ){
          xml<<"\t\t\t\t<face end=\""<<((face==0)?"POS":"NEG")<<"\">"<<endl;

          tree->GetEntry( layer*24 + column*2 + face );
          for( rg=0; rg<4; rg++ ){
            xml<<"\t\t\t\t\t<calPed avg=\""<<mpv[rg]
              <<"\" sig=\""<<width<<"\" range=\""
              <<((rg<2)?"L":"H")<<"EX"<<((rg%2==0)?"8":"1")<<"\" />"<<endl;
          }
    
          xml<<"\t\t\t\t</face>"<<endl;  
        }

        xml<<"\t\t\t</xtal>"<<endl;
      }

      xml<<"\t\t</layer>"<<endl;
    }
    
    //end file
    xml<<"\t</tower>"<<endl;
    xml<<"</calCalib>"<<endl;
    tree->SetBranchStatus("*", 1);
		histFile->Close();
  }

  void RootTreeAnalysis::CorrelatedPdlToXml( char* filename ){
    std::ofstream xml( filename );
    histFile= new TFile( m_histFileName );
    short int diode, layer, column, face;
    Float_t mpv[4], sigma[2], angle[2], ratio[2];
    TNtuple *tree;
      tree= (TNtuple*) histFile->Get("Pdl_tree");
    if (!(tree)){
      printf("ERROR:NO PEDESTALS TREE\n");
      return;
    }
    tree->SetBranchAddress("*", 0);
    tree->SetBranchAddress("MPV_LEX8", 1);
    tree->SetBranchAddress("MPV_LEX1", 1);
    tree->SetBranchAddress("MPV_HEX8", 1);
    tree->SetBranchAddress("MPV_HEX1", 1);
    tree->SetBranchAddress("angle_LEX", 1);
    tree->SetBranchAddress("angle_HEX", 1);
    tree->SetBranchAddress("Std_Dev_LEX", 1);
    tree->SetBranchAddress("Std_Dev_HEX", 1);
    tree->SetBranchAddress("Orth_Dev_LEX", 1);
    tree->SetBranchAddress("Orth_Dev_HEX", 1);

    tree->SetBranchAddress("MPV_LEX8", mpv);
    tree->SetBranchAddress("MPV_LEX1", mpv+1);
    tree->SetBranchAddress("MPV_HEX8", mpv+2);
    tree->SetBranchAddress("MPV_HEX1", mpv+3);
    tree->SetBranchAddress("angle_LEX", angle);
    tree->SetBranchAddress("angle_HEX", angle+1);
    tree->SetBranchAddress("Std_Dev_LEX", sigma);
    tree->SetBranchAddress("Std_Dev_HEX", sigma+1);
    tree->SetBranchAddress("Orth_Dev_LEX", ratio);
    tree->SetBranchAddress("Orth_Dev_HEX", ratio+1);
    
    
    xml<<"<?xml version=\"1.0\" ?>"<<endl<<endl
      <<"<!-- $Header: my Pdl-->"<<endl<<endl
      <<"<!-- Approximately real  CAL ped XML file for EM,"
      <<"according to calCalib.dtd -->"<<endl<<endl
      <<"<!DOCTYPE calCalib SYSTEM \"../calCalib_v3.dtd\" [] >"<<endl<<endl;

    xml<<"<calCalib>"<<endl
      <<"\t<generic\tinstrument=\"EM\" timestamp=\"2003-11-2-12:56\""<<endl
      <<"\t\t\t\tcalibType=\"CAL_Ped\" fmtVersion=\"v2r0\">"<<endl<<endl
      <<"\t\t<inputSample\tstartTime=\"2003-2-21-05:49:12\" "
      <<"stopTime=\"2003-2-24-07:07:02\""<<endl
      <<"\t\t\t\t\ttriggers=\"random\" mode=\"normal\" source=\"stuff\" >"<<endl
      <<endl;
    
    xml<<"\t\tTimes are start and stop time of calibration run."<<endl
      <<"\t\tOther attributes are just made up for code testing."<<endl
      <<"\t\t</inputSample>"<<endl
      <<"\t</generic>"<<endl<<endl<<endl;
    
    xml<<"<!-- EM instrument: 8 layers, 12 columns -->"<<endl<<endl
      <<"\t<dimension nRow=\"1\" nCol=\"1\" nLayer=\"8\" "
      <<"nXtal=\"12\" nFace=\"2\" />"<<endl<<endl
      <<"\t<tower iRow=\"0\" iCol=\"0\">"<<endl;
    
    for( layer=0; layer<8; layer++ ){
      xml<<"\t\t<layer iLayer=\""<<layer<<"\">"<<endl;
      
      for( column=0; column<12; column++ ){
        xml<<"\t\t\t<xtal iXtal=\""<<column<<"\">"<<endl;
        
        for( face=1; face>-1; face-- ){
          xml<<"\t\t\t\t<face end=\""<<((face==0)?"POS":"NEG")<<"\">"<<endl;

          tree->GetEntry( layer*24 + column*2 + face );
          for( diode=0; diode<2; diode++ ){
            xml<<"\t\t\t\t\t<calPed avg=\""<<mpv[2*diode]
              <<"\" sig=\""<<sigma[diode]
              <<"\" cos=\""<<cos(angle[diode])
              <<"\" range=\""<<((diode==1)?"H":"L")<<"EX8\" />"<<endl
              <<"\t\t\t\t\t<calPed avg=\""<<mpv[2*diode+1]
              <<"\" sig=\""<<ratio[diode]
              <<"\" cos=\""<<-1000
              <<"\" range=\""<<((diode)?"H":"L")<<"EX1\" />"<<endl;
          }
    
          xml<<"\t\t\t\t</face>"<<endl;  
        }

        xml<<"\t\t\t</xtal>"<<endl;
      }

      xml<<"\t\t</layer>"<<endl;
    }
    
    xml<<"\t</tower>"<<endl;
    xml<<"</calCalib>"<<endl;
    
    xml.close();
    tree->SetBranchAddress("*", 0);
    histFile->Close();
  }

  void RootTreeAnalysis::GainToXml( char* filename ){
    Float_t Nb_MeV= 12.5;
      std::ofstream xml( filename );
    histFile= new TFile( m_histFileName );
    Float_t peak, gain, rms, mgain[2], norm;
    TNtuple *tree_muons= (TNtuple*) histFile->Get("Gain_tree");
    if (!tree_muons){
      printf("ERROR:No Gain_tree TREE\n");
      return;
    }
    tree_muons->SetBranchStatus("*", 0);
    tree_muons->SetBranchStatus("MPV_LEX8", 1);
    tree_muons->SetBranchAddress("MPV_LEX8", &peak);
    tree_muons->SetBranchStatus("Width_LEX8", 1);
    tree_muons->SetBranchAddress("Width_LEX8", &rms);
    TNtuple *tree_pdl= (TNtuple*) histFile->Get("Pdl_tree");
    if (!tree_pdl){
      printf("ERROR:NO Pdl_tree TREE\n");
      return;
    }
    tree_pdl->SetBranchStatus("*", 0);
    tree_pdl->SetBranchStatus("Gain_LEX", 1);
    tree_pdl->SetBranchAddress("Gain_LEX", mgain);
    tree_pdl->SetBranchStatus("Gain_HEX", 1);
    tree_pdl->SetBranchAddress("Gain_HEX", mgain+1);
    TNtuple *tree_att= (TNtuple*) histFile->Get("Taper_tree");
    Bool_t taper= tree_att!=0;
    if ( taper ){
      tree_att= (TNtuple*) histFile->Get("LogPoverN_tree");
      if (!tree_att){
        printf("ERROR:NO LogPoverN_tree or Taper_tree TREE\n");
        return;
      }
      tree_att->SetBranchStatus("*", 0);
      tree_att->SetBranchStatus("norm_POS", 1);
      tree_att->SetBranchAddress("norm_POS", &norm );
      printf("ERROR:NO LogPoverN_tree TREE\n");
    }
    else { 
      tree_att->SetBranchStatus("*", 0);
      tree_att->SetBranchStatus("norm", 1);
      tree_att->SetBranchAddress("norm", &norm );
    }
    short int rg, layer, column, face;
    
    //header
    xml<<"<?xml version=\"1.0\" ?>"<<endl<<endl;
    xml<<"<!-- $Header: my Pdl-->"<<endl<<endl;
    xml<<"<!-- Approximately real  CAL gain XML file for EM, according to calCalibv3.dtd -->"<<endl<<endl;

    xml<<"<!DOCTYPE calCalib SYSTEM \"../calCalib_v3.dtd\" [] >"<<endl<<endl;

    //calCalib
    xml<<"<calCalib>"<<endl;
    xml<<"\t<generic\tinstrument=\"EM\" timestamp=\"2003-11-2-12:56\""<<endl;
    xml<<"\t\t\t\tcalibType=\"CAL_ElecGain\" fmtVersion=\"v2r0\">"<<endl<<endl;
    xml<<"\t\t<inputSample\tstartTime=\"2003-2-21-05:49:12\" stopTime=\"2003-2-24-07:07:02\""<<endl; 
    xml<<"\t\t\t\t\ttriggers=\"random\" mode=\"normal\" source=\"stuff\" >"<<endl<<endl;
    
    xml<<"\t\tTimes are start and stop time of calibration run."<<endl;
    xml<<"\t\tOther attributes are just made up for code testing."<<endl;    
    
    xml<<"\t\t</inputSample>"<<endl;
    xml<<"\t</generic>"<<endl<<endl<<endl;
    
    xml<<"<!-- EM instrument: 8 layers, 12 columns -->"<<endl<<endl;
    xml<<"\t<dimension nRow=\"1\" nCol=\"1\" nLayer=\"8\" nXtal=\"12\" nFace=\"2\" />"<<endl<<endl;
    
    xml<<"\t<tower iRow=\"0\" iCol=\"0\">"<<endl;
    
    for( layer=0; layer<8; layer++ ){
      xml<<"\t\t<layer iLayer=\""<<layer<<"\">"<<endl;
      
      for( column=0; column<12; column++ ){
        xml<<"\t\t\t<xtal iXtal=\""<<column<<"\">"<<endl;
        
        for( face=1; face>-1; face-- ){
          xml<<"\t\t\t\t<face end=\""<<((face==0)?"POS":"NEG")<<"\">"<<endl;
          tree_pdl->GetEntry( layer*24 + column*2 + face );
          tree_muons->GetEntry( layer*24 + column*2 + face );
          if( !taper ) tree_att->GetEntry( layer*24 + column*2 + face );
          for( rg=0; rg<4; rg++ ){
            if( taper ) 
              tree_att->GetEntry( 192*rg+layer*24 + column*2 + face );
            switch( rg ){
              case 0:
                gain= peak/Nb_MeV; 
                if( taper || gain ) gain *= norm;
                else gain/= norm;
                rms/= peak; break;
              case 2:  
                gain /= 8;
                rms /= 8; break;
              default:  
                gain *= mgain[rg/2];
                rms *= mgain[rg/2]; break;
            }

            xml<<"\t\t\t\t\t<calGain avg=\""<<gain
              <<"\" sig=\""<<rms<<"\" range=\""
              <<((rg<2)?"L":"H")<<"EX"<<((rg%2==0)?"8":"1")<<"\" />"<<endl;
          }
    
          xml<<"\t\t\t\t</face>"<<endl;  
        }

        xml<<"\t\t\t</xtal>"<<endl;
      }

      xml<<"\t\t</layer>"<<endl;
    }
    
    //end file
    xml<<"\t</tower>"<<endl;
    xml<<"</calCalib>"<<endl;
    xml.close();

    tree_pdl->SetBranchStatus("*", 1);
    tree_muons->SetBranchStatus("*", 1);
		tree_att->SetBranchStatus("*", 1);
    histFile.Close();
  }

  void RootTreeAnalysis::MuonSlopesToXml( char* filename ){
      std::ofstream xml( filename );
    histFile= new TFile( m_histFileName );
    Float_t slope;
    TNtuple *tree= (TNtuple*) histFile->Get("logPoverN");
    if (!tree){
      printf("ERROR:NO CAL CORR TREE\n");
      return;
    }
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("slope", 1);
    tree->SetBranchAddress("slope", &slope);
    short int rg, layer, column, face;
    
    //header
    xml<<"<?xml version=\"1.0\" ?>"<<endl<<endl;
    xml<<"<!-- $Header: my Pdl-->"<<endl<<endl;
    xml<<"<!-- Approximately real  CAL muon slope XML file for EM, according to calCalibv3.dtd -->"<<endl<<endl;

    xml<<"<!DOCTYPE calCalib SYSTEM \"../calCalib_v3.dtd\" [] >"<<endl<<endl;

    //calCalib
    xml<<"<calCalib>"<<endl
      <<"\t<generic\tinstrument=\"EM\" timestamp=\"2003-11-2-12:56\""<<endl
      <<"\t\t\t\tcalibType=\"CAL_MuSlope\" fmtVersion=\"v2r0\">"<<endl<<endl
      <<"\t\t<inputSample\tstartTime=\"2003-2-21-05:49:12\""
      <<" stopTime=\"2003-2-24-07:07:02\""<<endl
      <<"\t\t\t\t\ttriggers=\"random\" mode=\"normal\" source=\"stuff\" >"<<endl
      <<endl;
    
    xml<<"\t\tTimes are start and stop time of calibration run."<<endl
      <<"\t\tOther attributes are just made up for code testing."<<endl;    
    
    xml<<"\t\t</inputSample>"<<endl
      <<"\t</generic>"<<endl<<endl<<endl;
    
    xml<<"<!-- EM instrument: 8 layers, 12 columns -->"<<endl<<endl
      <<"\t<dimension nRow=\"1\" nCol=\"1\" "
      <<"nLayer=\"8\" nXtal=\"12\" nFace=\"1\" />"<<endl<<endl;
    
    xml<<"\t<tower iRow=\"0\" iCol=\"0\">"<<endl;
    
    for( layer=0; layer<8; layer++ ){
      xml<<"\t\t<layer iLayer=\""<<layer<<"\">"<<endl;
      
      for( column=0; column<12; column++ ){
        xml<<"\t\t\t<xtal iXtal=\""<<column<<"\">"<<endl;
        xml<<"\t\t\t\t<face end=\"NA\">"<<endl;
        tree->GetEntry( layer*24 + column*2 + face );
        for( rg=0; rg<4; rg++ ){
          xml<<"\t\t\t\t\t<muSlope slope=\""<<slope<<"\" range=\"";
          xml<<((rg<2)?"L":"H")<<"EX"<<((rg%2==0)?"8":"1")<<"\" />"<<endl;
        }
        xml<<"\t\t\t\t</face>"<<endl;  
        xml<<"\t\t\t</xtal>"<<endl;
      }
      xml<<"\t\t</layer>"<<endl;
    }
    
    //end file
    xml<<"\t</tower>"<<endl;
    xml<<"</calCalib>"<<endl;
    xml.close();
		tree->SetBranchStatus("*", 0);
		histFile->Close();
  }


void RootTreeAnalysis::Go(Int_t numEvents){
  //init    
  // Purpose and Method:  Event Loop
  //   All analysis goes here
  //  To read only selected branches - saves processing time
  //  Comment out any branches you are not interested in.

  //branch status 
  if (digiTree) {
      digiTree->SetBranchStatus("*",0);  // disable all branches
      // activate desired brances
      digiTree->SetBranchStatus("m_cal*",1);  
      digiTree->SetBranchStatus("m_tkr*",1);  
      digiTree->SetBranchStatus("m_eventId", 1); 
      digiTree->SetBranchStatus("m_runId", 1);
  }
  
  // leaving all recon branches activated for now
  if (reconTree) {
      reconTree->SetBranchStatus("*",0);  // disable all branches
      // activate desired branches
      reconTree->SetBranchStatus("m_cal", 1);  
      reconTree->SetBranchStatus("m_tkr", 1);
  }
  
    
  Int_t nentries = GetEntries();
  std::cout << "\nNum Events in File is: " << nentries << std::endl;
  Int_t curI=0;
  Int_t nMax = TMath::Min(numEvents+m_StartEvent,nentries);
  if (m_StartEvent == nentries) {
      std::cout << " all events in file read" << std::endl;
      return;
  }
  if (nentries <= 0) return;
  
  // Keep track of how many bytes we have read in from the data files
  Int_t nbytes = 0, nb = 0;

  HistDefine();
  if( (job_id>1) && !getPdl( pdl_mpv, pdl_width ) ) return;
	if( (job_id==4) && !getLogPoverN( logPoverN_slope ) ) return;
  

  // BEGINNING OF EVENT LOOP
  Double_t position[19];
  Int_t XtalAdc[16];
  Short_t col[8];
  for (Int_t ievent= m_StartEvent; ievent<nMax; ievent++, curI=ievent) {
    //init 
    if (evt) evt->Clear();
    if (rec) rec->Clear();
    
    digiEventId = 0; reconEventId = 0; mcEventId = 0;
    digiRunNum = 0; reconRunNum = 0; mcRunNum = 0;
    
    nb = GetEvent(ievent);
    nbytes += nb;

    if ( evt ) {
      switch(job_id){
        case 0:
          Pdl(); break;
        case 1:
          Pdl_Correlated(); break;
        case 2:
          Pdl_Corr_NoMuon(); break;
        default:
          if ( (Make_Cut( position, XtalAdc, col ))==0 ) continue;
          Make_Hist( position, XtalAdc, col );
      }
      if (ievent%1000==0) {
        printf("[%ld]", ievent); fflush(stdout);
      }
      if (ievent>=50000 && ievent%50000==0) 
        histFile->Write("", TObject::kOverwrite);
    }
  }  // end analysis code in event loop

  //end
  m_StartEvent = curI;
  histFile->Write("", TObject::kOverwrite);
  histFile->Close();
}
#endif
