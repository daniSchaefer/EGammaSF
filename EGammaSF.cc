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
   fit_flag_ = parser->FitFlag();
        PrintDebug(Form("initialising histogram names using : %s, %s, %s", parser->NameSF().c_str() ,parser->NameEffMC().c_str(), parser->NameEffData().c_str()));
        PrintDebug(Form("fitting flag set to : %i",fit_flag_));
   //set fit range:
        int nBins = egm2d_.GetYaxis()->GetNbins();
        rangelow_= egm2d_.GetYaxis()->GetBinLowEdge(1);
        rangeup_ = egm2d_.GetYaxis()->GetBinLowEdge(nBins)+ egm2d_.GetYaxis()->GetBinWidth(nBins);
        PrintDebug(Form("set range to %.1f - %.1f",rangelow_,rangeup_));
   TCanvas* c = new TCanvas();
   egm2d_.Draw("COLZ");
   c->SaveAs("histo.pdf");
        PrintDebug("===============================================================");
        PrintDebug("starting fit for scale factors ");
        PrintDebug("===============================================================");
   for(int i = 1; i< egm2d_.GetXaxis() -> GetNbins() +1;i++)
   {
        TGraphErrors g = GetTGraph(i);
        TF1 fit = FitScaleFactor(g,i);
        smooth_sf_.push_back(fit);
        PrintDebug(Form("fitted scale factor for eta bin %i using function %i",i,fit_flag_));
        DrawSF(g,fit,egm2d_.GetXaxis() -> GetBinCenter(i));
   }
   for(int i=1;i<egm2d_.GetXaxis()->GetNbins()+1;i++)
   {
      TGraphErrors gunc = GetTGraph(i,1);
      TF1 fitunc = GetUncertaintyFunction(i);
      fitunc.SetParameter(2,gunc.GetMinimum());
      TFitResult* fitres = (gunc.Fit(&fitunc,"RW")).Get();
      //fitres->Print();
      DrawSF(gunc,egm2d_.GetXaxis() -> GetBinCenter(i));
       
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
  
  //std::cout << "pt : " <<pT << " pt bin : " << ptBin_ << " eta: " << superClusterEta << " eta bin : " << etaBin_ << std::endl;
  
  return sf;  
}

float ScaleFactorHelper::GetSFSmooth(float pT, float superClusterEta)
{
   SetEtaBin(superClusterEta);
   if(pT > rangeup_) return smooth_sf_.at(etaBin_-1).Eval(rangeup_);
   if(pT < rangelow_) return smooth_sf_.at(etaBin_-1).Eval(rangelow_);
   float sf = smooth_sf_.at(etaBin_-1).Eval(pT);
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


TGraphErrors ScaleFactorHelper::GetTGraph(int etaBin,bool forUncertainty)
{
    int max = egm2d_.GetYaxis()-> GetNbins();
    float eta = egm2d_.GetXaxis()-> GetBinCenter(etaBin);
    TGraphErrors *g = new TGraphErrors(max);
    
    for(int i=1; i< max+1;i++)
    {
        float pt = egm2d_.GetYaxis() -> GetBinCenter(i);
        float sf = GetSF(pt,eta);
        float pt_unc = egm2d_.GetYaxis() -> GetBinWidth(i) /2.;
        //std::cout << " pt " << pt << " pt unc " << pt_unc << std::endl;
        float sf_unc = GetUncertainty(pt,eta);
        if (forUncertainty)
        {
           g->SetPoint(i-1,pt,sf_unc/sf);
           g->SetPointError(i-1,pt_unc,0.);
        }
        else
        {
        g->SetPoint(i-1,pt,sf);
        g->SetPointError(i-1,pt_unc,sf_unc);
        }
    }
    return *g;
}

TF1 ScaleFactorHelper::FitScaleFactor(TGraphErrors g, int etaBin)
{
  TF1 fit = GetFitFunction(etaBin);
  TFitResult* res = (g.Fit(&fit,"MSR")).Get();
  res->Print();  
  return fit;  
}


void ScaleFactorHelper::DrawSF(TGraphErrors g, TF1 f, float eta)
{
    TCanvas* cg = new TCanvas();
    g.GetXaxis()->SetTitle("pT (GeV)");
    g.GetYaxis()->SetTitle("scale factor");
    g.SetLineColor(kBlue);
    g.SetLineWidth(2);
    g.SetMarkerColor(kBlack);
    g.SetMarkerStyle(8);
    f.SetLineColor(kRed);
    f.SetLineWidth(2);
    //g.SetMaximum(1.1);
    //g.SetMinimum(0.8);
    g.Draw("ALP");
    f.Draw("same");
    TPaveText* addInfo = new TPaveText(0.9,0.54,0.64,0.4,"NDC");
    addInfo->SetFillColor(0);
    addInfo->SetLineColor(0);
    addInfo->SetFillStyle(0);
    addInfo->SetBorderSize(0);
    addInfo->SetTextFont(42);
    addInfo->SetTextSize(0.040);
    addInfo->SetTextAlign(12);
    float Chi2 = g.Chisquare(&f,"R");
    int ndof   = maxBinpt_ -2;
    std::cout << " chi2 " <<  Chi2 << " ndof " << ndof << std::endl;
    addInfo->AddText(Form("#eta =  %.1f ", eta));
    addInfo->AddText(f.GetExpFormula());
    addInfo->AddText(Form("#chi^{2}/ndof = %.3f ",Chi2/ndof));
    addInfo->Draw("same");
    cg->SaveAs(Form("graph_%.1f_fitfunction%i.pdf",eta,fit_flag_));  
      
    
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
    g.Draw("ALP");
    TPaveText* addInfo = new TPaveText(0.9,0.54,0.64,0.4,"NDC");
    addInfo->SetFillColor(0);
    addInfo->SetLineColor(0);
    addInfo->SetFillStyle(0);
    addInfo->SetBorderSize(0);
    addInfo->SetTextFont(42);
    addInfo->SetTextSize(0.040);
    addInfo->SetTextAlign(12);
    addInfo->AddText(Form("#eta =  %.1f ", eta));
    addInfo->Draw("same");
    cg->SaveAs(Form("graph_%.1f_uncertainty.pdf",eta));  
}



float ScaleFactorHelper::GetEfficiency(float pT, float superClusterEta,bool isData)
{
    if( pT!=input_pt_ or superClusterEta != input_eta_)
    {
        SetEtaBin(superClusterEta);
        SetPtBin(pT);
    } 
    if(isData) return efficiency_data_.GetBinContent(etaBin_,ptBin_); 
    return efficiency_mc_.GetBinContent(etaBin_,ptBin_); 
    
}

 TF1 ScaleFactorHelper::GetFitFunction(int etaBin)
 {
     TF1* f;
     if(fit_flag_==0)
     {
         f = new TF1("f","[0]+ atan([1]*x)",rangelow_,rangeup_);
         return *f;
     }
     if(fit_flag_==1)
     {
         f = new TF1("f","[0]+ [1]*(1/x)",rangelow_,rangeup_);
         f->SetParLimits(1,-10,0.);
         return *f;
     }
     if(fit_flag_==2)
     {
         f = new TF1("f","[0]+ [1]*(1/(x*x))",rangelow_,rangeup_);
         f->SetParameter(0,egm2d_.GetBinContent(maxBinpt_,etaBin));
         f->SetParLimits(1,-600,0.);
         return *f;
     }
     if(fit_flag_==3)
     {
         f = new TF1("f","[0] + [1]*(1/x^[2])",rangelow_,rangeup_);
         //f->SetParameter(0,egm2d_.GetBinContent(maxBinpt_,etaBin));
         f->SetParLimits(1,-600,0.);
         return *f;
     }
     
     if(fit_flag_==11)
     {
         f = new TF1("f","[0]",rangelow_,rangeup_);
         //f->SetParameter(0,egm2d_.GetBinContent(maxBinpt_,etaBin));
         //f->SetParLimits(1,-600,0.);
         return *f;
     }
     
     
     throw my_range_error("wrong fitting flag used in scale factor fit");
 }

 
 
TF1 ScaleFactorHelper::GetUncertaintyFunction(int etaBin)
 {
     TF1* f;
     if(uncertainty_flag_==0)
     {
         f = new TF1("f","[2]+ [1]*(x-[0])^(2)",rangelow_,150.);
         f->SetParLimits(0,20.,150.);
         return *f;
     }     
     
     throw my_range_error("wrong fitting flag used in uncertainty smoothing");
 }
 
    
