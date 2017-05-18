#include "EGammaSF.h"


#include "TH1F.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <iomanip>



ScaleFactorHelper::ScaleFactorHelper(EGammaInput what, bool debugging )
 {
   debug_flag_ = debugging;
   input_ = what;
   ConfigParser* parser = new ConfigParser("config.txt",what);
   std::string filename = parser->NameFile();
        PrintDebug("opening file "+filename);
   TFile file(filename.c_str(),"READ");
    if (file.IsZombie()) throw file_not_found();
        std::cout << &file << std::endl;
   egm2d_ = *dynamic_cast<TH2F*>(file.Get(parser->NameSF().c_str()));          
   efficiency_mc_ = *dynamic_cast<TH2F*>(file.Get(parser->NameEffMC().c_str()));
   efficiency_data_ = *dynamic_cast<TH2F*>(file.Get(parser->NameEffData().c_str()));

   TCanvas* c = new TCanvas();
   egm2d_.Draw("COLZ");
   c->SaveAs("histo.pdf");
  
   for(int i = 1; i< egm2d_.GetXaxis() -> GetNbins() +1;i++)
   {
        TGraphErrors g = GetTGraph(i);
        DrawSF(g, egm2d_.GetXaxis() -> GetBinCenter(i));
   }
 }
 
 
ScaleFactorHelper::~ScaleFactorHelper(void)
 {
     
 }
 
float ScaleFactorHelper::GetSF(float pT, float superClusterEta)
{
  float sf;
  if( pT!=input_pt_ or superClusterEta != input_eta_)
  {
    SetEtaBin(superClusterEta);
    SetPtBin(pT);
  }
  sf = egm2d_.GetBinContent(etaBin_,ptBin_);
  
  std::cout << "pt : " <<pT << " pt bin : " << ptBin_ << " eta: " << superClusterEta << " eta bin : " << etaBin_ << std::endl;
  
  return sf;  
}

float ScaleFactorHelper::GetUncertainty(float pT, float superClusterEta)
{
 float sf_unc;
  if( pT!=input_pt_ or superClusterEta != input_eta_)
  {
    SetEtaBin(superClusterEta);
    SetPtBin(pT);
  }
  sf_unc = egm2d_.GetBinError(etaBin_,ptBin_);
  // add additional uncertainty to electron reconstruction scale-factor:
  if(input_ == EGammaInput::electronRecoSF)
  {
      //std::cout << " add additional uncertainty for electron REco " << std::endl;
        if(pT > 20. or pT > 80.)
        {
          //std::cout << "pt is smaller 20 or larger 80 GeV :  " <<  0.01*egm2d_.GetBinContent(etaBin_,ptBin_) << std::endl; 
          sf_unc = TMath::Sqrt(pow(sf_unc,2) + pow(0.01*egm2d_.GetBinContent(etaBin_,ptBin_),2));   
        }
  }
  
  return sf_unc;
}


void ScaleFactorHelper::SetEtaBin(float superClusterEta)
{
  input_eta_ = superClusterEta;  
  etaBin_ = egm2d_.GetXaxis() -> FindBin(superClusterEta);
  maxBineta_ = egm2d_.GetXaxis() -> GetNbins();
  if (etaBin_==0 or etaBin_==maxBineta_+1) 
  {std::string err = "tried to evaluate scale factor for |eta| >"; err.append(std::to_string(TMath::Abs(egm2d_.GetXaxis()->GetBinLowEdge(1)))); throw my_range_error(err);}
}

void ScaleFactorHelper::SetPtBin(float pT)
{
  input_pt_ = pT;  
  ptBin_ = egm2d_.GetYaxis() -> FindBin(pT);
  maxBinpt_ = egm2d_.GetYaxis() -> GetNbins();
  if(ptBin_ ==0) {std::string err = "tried to evaluate scale factor for pt < "; err.append(std::to_string(egm2d_.GetYaxis()->GetBinLowEdge(ptBin_+1)));
      std::cout << err << std::endl;
      ptBin_=1;
      /*throw my_range_error(err);*/}
  if(ptBin_ >= maxBinpt_ and maxBinpt_ > 1 ) ptBin_ = maxBinpt_-1;  // last pt bin is only used as control but the scale-factor here should not be used because of limited statistics
}


TGraphErrors ScaleFactorHelper::GetTGraph(int etaBin)
{
    int max = egm2d_.GetYaxis()-> GetNbins();
    float eta = egm2d_.GetXaxis()-> GetBinCenter(etaBin);
    TGraphErrors *g = new TGraphErrors(max);
    
    std::cout << max << std::endl;
    for(int i=1; i< max+1;i++)
    {
        std::cout << i << std::endl;
        float pt = egm2d_.GetYaxis() -> GetBinCenter(i);
        float sf = GetSF(pt,eta);
        float pt_unc = egm2d_.GetYaxis() -> GetBinWidth(i) /2.;
        std::cout << " pt " << pt << " pt unc " << pt_unc << std::endl;
        float sf_unc = GetUncertainty(pt,eta);
        g->SetPoint(i-1,pt,sf);
        g->SetPointError(i-1,pt_unc,sf_unc);
        std::cout << " bla " << i << std::endl;
    }
    return *g;
}


void ScaleFactorHelper::DrawSF(TGraphErrors g, float eta)
{
  TCanvas* cg = new TCanvas();
  g.GetXaxis()->SetTitle("pT (GeV)");
  g.GetYaxis()->SetTitle("scale factor");
  g.SetLineColor(kBlue);
  g.SetLineWidth(2);
  g.SetMarkerColor(kBlack);
  g.SetMarkerStyle(8);
  g.SetMaximum(1.1);
  g.SetMinimum(0.8);
  g.Draw("ALP");
  TPaveText* addInfo = new TPaveText(0.9,0.02,0.64,0.3,"NDC");
  addInfo->SetFillColor(0);
  addInfo->SetLineColor(0);
  addInfo->SetFillStyle(0);
  addInfo->SetBorderSize(0);
  addInfo->SetTextFont(42);
  addInfo->SetTextSize(0.040);
  addInfo->SetTextAlign(12);
  addInfo->AddText(Form("#eta =  %.1f ", eta));
  addInfo->Draw("same");
  cg->SaveAs(Form("graph_%.1f.pdf",eta));  
    
    
}

    
