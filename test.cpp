// compile with : g++ test.cpp EGammaSF.cc `root-config --cflags --glibs --libs --evelibs`
           
#include "TTree.h"
#include "TROOT.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2D.h"
#include "TF1.h"
#include "EGammaSF.h"
#include "TRandom.h"
//#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd5/include/TTreeReaderValueBase.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <exception>


void DrawTest(ScaleFactorHelper* sfhelper,float eta)
{
       TGraph * g = new TGraph(1);
       TGraph* gup = new TGraph(1);
       TGraph* gdown = new TGraph(1);
       for (int j=0; j<800;j++){
           float pt = 10.+j*1.5;
           float sf = sfhelper->GetSFSmooth(pt,eta);
           float unc = sfhelper->GetUncertaintySmooth(pt,eta);
           g->SetPoint(j,pt,sf);
           std::cout << j << " pt " << pt << " sf " << sf << std::endl;
           gup->SetPoint(j,pt,sf+unc);
           gdown->SetPoint(j,pt,sf-unc);
       }   
 TCanvas* c = new TCanvas("test","test",400,400);
 c->SetLogx();
 g->SetLineColor(kBlack);
 g->SetLineWidth(2);
 gup->SetLineColor(kGreen);
 gup->SetLineWidth(2);
//gup->Draw("LA");
 g->Draw("LA");
 gdown->SetLineColor(kGreen);
 gdown->SetLineColor(2);
 //gdown->Draw("Lsame");
 c->SaveAs("test.pdf");
 
 
}


int main(int argc, char** argv)
{
     try{
        //initialize scale factor helper : 
//          ScaleFactorHelper* bla0  = new ScaleFactorHelper(EGammaInput::electronCutBasedVetoID,1);
//          ScaleFactorHelper* bla1 = new ScaleFactorHelper(EGammaInput::electronMedium,1);
//          ScaleFactorHelper* bla2 = new ScaleFactorHelper(EGammaInput::electronLoose,1);
         ScaleFactorHelper* bla3 = new ScaleFactorHelper(EGammaInput::electronTight,0);
//          ScaleFactorHelper* bla4 = new ScaleFactorHelper(EGammaInput::electronMVA80,1);
//          ScaleFactorHelper* bla5 = new ScaleFactorHelper(EGammaInput::electronMVA90,1);
//          
//           
//          ScaleFactorHelper* bla0p = new ScaleFactorHelper(EGammaInput::photonMedium,1);
//          ScaleFactorHelper* bla1p = new ScaleFactorHelper(EGammaInput::photonLoose,1);
//          ScaleFactorHelper* bla2p = new ScaleFactorHelper(EGammaInput::photonTight,1);
//         //ScaleFactorHelper* bla = new ScaleFactorHelper(EGammaInput::electronVetoPixelSeedEndcap,1);
        
        
        
        
//         TRandom *r = new TRandom(42);
//         for(int i=0;i<80;i+=1)
//         {
//         float pt = TMath::Abs(r->Gaus(130,80));
//         float unc = bla->GetUncertainty(pt,2.1);
//         float sf = bla->GetSF(pt,2.1);
//         float sf_smooth = bla->GetSFSmooth(pt, 2.1);
//         float un_rel = bla->GetUncertaintySmooth(pt, 2.1);
//         std::cout <<"pt : " << pt << "  sf : " << sf << " +- " << unc << " smoothed : "<< sf_smooth << " unc smoothed : " << un_rel << std::endl;
//         }
//         
//         float eff = bla->GetEfficiency(600,-1.2,0);
//         float effdata = bla->GetEfficiency(600,-1.2,1);
//         std::cout << bla->GetSF(600,-1.2) << std::endl;
//         std::cout << "GetEfficiency mc " << eff << " efficiency data : " << effdata << std::endl;
//         std::cout << " effdata/effmc  " << effdata/eff << std::endl;
        DrawTest(bla3,-0.4); 
         
      }
      catch(const std::exception& ex) {
        std::cout << "EXCEPTION!!!" << std::endl;
      }

    
    
   return 0; 
    
}
