#include "EGammaSF.h"


#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TGraphSmooth.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TRint.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <iomanip>

// convert enum class to string -> used for naming 

std::string GetString(EGammaInput sf){
       if (sf == EGammaInput::electronRecoSF) return "electronRecoSF";
       if (sf == EGammaInput::electronLoose) return "electronLoose";
       if (sf == EGammaInput::electronTight) return "electronTight";
       if (sf == EGammaInput::electronMedium) return "electronMedium";
       if (sf == EGammaInput::electronMVA80) return "electronMVA80";
       if (sf == EGammaInput::electronMVA90) return "electronMVA90";
       
       if (sf == EGammaInput::photonLoose) return "photonLoose";
       if (sf == EGammaInput::photonTight) return "photonTight";
       if (sf == EGammaInput::photonMedium) return "photonMedium";
       if (sf == EGammaInput::photonMVA90) return "photonMVA90";
       if (sf == EGammaInput::electronVetoPixelSeedEndcap) return "electronVetoPixelSeedEndcap";
       if (sf == EGammaInput::electronVetoPixelSeedBarrel) return "electronVetoPixelSeedBarrel";
       if (sf == EGammaInput::electronVetoPixelSeed) return "electronVetoPixelSeed";
       if (sf == EGammaInput::electronVetoCSEV) return "electronVetoCSEV";
       if (sf == EGammaInput::electronCutBasedVetoID ) return "electronCutBasedVetoID";
       
       
       throw my_range_error("invalid input parameters");
}


ScaleFactorHelper::ScaleFactorHelper(EGammaInput what, bool debugging )
 {
   logfile_.open("out.log",ios::trunc);
   logfile_ << "Begin scale factor fitting for "<< GetString(what) <<std::endl;
   logfile_.close();
   debug_flag_ = debugging;
   input_ = what;
   ConfigParser* parser = new ConfigParser("config.txt",what);
   std::string filename = parser->NameFile();
        PrintDebug("opening file "+filename);
   std::string openAs = "READ";
   TFile file(filename.c_str(),openAs.c_str());
        if (file.IsZombie()) throw file_not_found();
   std::string outfilename = "OUT"+filename;
   file_ = outfilename;
   egm2d_ = *dynamic_cast<TH2F*>(file.Get(parser->NameSF().c_str()));
   
   if (parser->NameEffData().find("NONE")==std::string::npos)
   {
        efficiency_mc_ = *dynamic_cast<TH2F*>(file.Get(parser->NameEffMC().c_str()));
        efficiency_data_ = *dynamic_cast<TH2F*>(file.Get(parser->NameEffData().c_str()));
   }
   efficiency_data_.SetName(parser->NameEffData().c_str());
   efficiency_mc_.SetName(parser->NameEffMC().c_str());
   file.Close();
   
   // set maximum of bins in input histograms  
   egm2d_.SetName(parser->NameSF().c_str());
   maxBineta_ = egm2d_.GetXaxis()->GetNbins();
   maxBinpt_  = egm2d_.GetYaxis()->GetNbins();
   
   
   //set fit range:
   int nBins = egm2d_.GetYaxis()->GetNbins();
   rangelow_= egm2d_.GetYaxis()->GetBinLowEdge(1);
   rangeup_ = egm2d_.GetYaxis()->GetBinLowEdge(nBins)+ egm2d_.GetYaxis()->GetBinWidth(nBins);
   
   
   // set flags for which fit function should be used as default and globally -> only important for the debugging mode where functions are chosen for the user
   fit_flag_ = parser->FitFlag();
   uncertainty_flag_ = parser->UncFlag();
   local_flag_ = parser->LocalFitFlag();
   local_unc_flag_ = parser->LocalUncFlag();
   unc_range_local_ = parser->GetFitRangeUnc();
   for (int e=1;e<= maxBineta_;e++)
   {
    if(local_flag_.find(e)==  local_flag_.end()) local_flag_.insert(std::pair<int, int>(e,fit_flag_));
    if(local_unc_flag_.find(e)==  local_unc_flag_.end()) local_unc_flag_.insert(std::pair<int, int>(e,uncertainty_flag_));
    if(unc_range_local_.find(e) == unc_range_local_.end()) unc_range_local_.insert(std::pair<int,int>(e,rangeup_));
   }
   

        PrintDebug(Form("initialising histogram names using : %s, %s, %s", parser->NameSF().c_str() ,parser->NameEffMC().c_str(), parser->NameEffData().c_str()));
        PrintDebug(Form("fitting flag set to : %i",fit_flag_));
        PrintDebug(Form("uncertainty flag set to %i",uncertainty_flag_));
        PrintDebug(Form("set range to %.1f - %.1f",rangelow_,rangeup_));
        PrintDebug(Form("Draw input sf histogram %s ",("histo_"+GetString(input_)+".pdf").c_str()));
        // draw input sf -histo
        if(debug_flag_){
            TCanvas* c = new TCanvas();
            egm2d_.Draw("COLZ");
            c->SaveAs(("histo_"+GetString(input_)+".pdf").c_str());
        }
   
   // if in debug mode open output file to write the fit functions into
   PrintDebug(Form("create output file %s ",outfilename.c_str())); 
   if(debug_flag_) openAs = "RECREATE";
   TFile out(outfilename.c_str(),openAs.c_str());
   // do fits /smoothing in these functions or load the predifined fit functions for user usage:
   if (input_ != EGammaInput::electronRecoSF){
        std::vector<int> fitflags = {4,11,6};
        std::vector<double> Chi2 = TryFits(fitflags,1, maxBineta_);
        for(int c=0;c<Chi2.size();c++)
        {
            if (Chi2.at(c) >3.0) {
            PrintDebug(Form( "ATTENTION: fit failed! chi2  = %.2f  for eta bin %i",Chi2.at(c),egm2d_.GetXaxis()->GetBinCenter(c+1)));
            fitflags = {6,0,1,3,9,10,4};
            fit_flag_= -99;
            std::cout << "use alternative fits !"<<std::endl;
            TryFits(fitflags,c+1,c+1);
            }
        }
        //SetFitFlagManually(9,1);
        InitializeSF();
        InitializeUnc();
        // draw fits for cross checking
        DrawAll();
        DrawAllUncertainty();
   }
   out.Close();   
 }
 
 
 
 
// destructor 
ScaleFactorHelper::~ScaleFactorHelper(void) { }



 // set local fit flags corresponding to the etaBin in use
 void ScaleFactorHelper::SetFitFlag(int etaBin) { fit_flag_         = local_flag_.at(etaBin);
                                                  uncertainty_flag_ = local_unc_flag_.at(etaBin);
                                                  unc_range_ =  unc_range_local_.at(etaBin);
}


 
TGraphErrors ScaleFactorHelper::GetGraphNumSmoothed(TH1F h)
{
    h.Smooth();
    TGraphErrors *g = new TGraphErrors(0);
    for ( int i =1; i< h.GetNbinsX()+1; i++)
    {
     g->SetPoint(i-1,h.GetXaxis()->GetBinCenter(i),h.GetBinContent(i));  
     g->SetPointError(i-1,h.GetBinWidth(i)/2.,h.GetBinError(i));
    }
    return *g;
}
 
 std::vector<double> ScaleFactorHelper::TryFits(std::vector<int> fitflags,int binmin,int binmax)
 {
   if(fit_flag_==-99 and debug_flag_)
   {
      std::vector<double> result={}; 
      for(int eb = binmin; eb < binmax+1;eb++)
      //for(int eb = 1; eb < 1+1;eb++)
      {    
    // try which function fits best
          PrintDebug(Form( "=========================================== eta bin %i ================================ ",eb));
          TCanvas* tc = new TCanvas("tc","tc",400,400);
          std::vector<double> chi2={};
          std::vector<int> ndof={};
          //std::vector<int> fitflags = {0,1,2,3,4,11,6};
           
          std::vector<TGraphErrors> graphs ={};
          std::vector<TF1> fits={};
          TGraphErrors gS = GetTGraph(eb);
          graphs.push_back(gS);
          for(int f=0; f< fitflags.size(); f++)
          {
            
            fit_flag_ = fitflags.at(f);
            TF1 fit = FitScaleFactor(gS,eb);
            chi2.push_back(fit.GetChisquare());
            ndof.push_back(ndof_[fitflags.at(f)]);
            std::cout << chi2[eb] << " "<<ndof[eb] << std::endl;
            fits.push_back(fit);
        
          }
          TLegend* leg = new TLegend(0.4304361,0.263552,0.808124,0.5624967);
           gS.SetLineColor(kBlue);
           gS.SetLineWidth(2);
           gS.Draw("alp");
           
           
          for(int i=0;i<fits.size();i++)
          {
           fits.at(i).SetLineColor(kRed+i);
           if (i==5) fits.at(i).SetLineColor(kMagenta);
           fits.at(i).SetLineWidth(2);
           fits.at(i).Draw("same");
           leg->AddEntry(&fits.at(i),fits.at(i).GetName(),"l");
          }
          leg->Draw("same");
          tc->SaveAs(("fitTestsSF_bin"+std::to_string(eb)+".pdf").c_str());
        int index=0;
        double Lmin =10000;
        for ( int i=0;i<chi2.size();i++)
        {
          double L = chi2.at(i) + ndof.at(i);
          PrintDebug(Form( "Chi2 %.2f fit function %i",chi2.at(i), fitflags.at(i)));
          std::cout << " chi2 " << chi2.at(i) << std::endl;
          if (L < Lmin and !(TMath::Abs(chi2.at(i))== std::numeric_limits<double>::infinity()) )
          {
           Lmin = L;
           index = i;   
          }
        }
        PrintDebug(Form("Lmin %.3f  index %i   ",Lmin,index));
        if ( int(chi2.at(index)*10) < 10 and int(chi2.at(index)*10) !=0)
        {
          for(unsigned int c = 0; c < chi2.size();c++)
          {
            if(chi2.at(c) < 1 and chi2.at(c) > chi2.at(index))
                index = c;
          }
        }
    // set local fit flag to the best function!
        PrintDebug(Form( " taking fit function : %i index %i",fitflags.at(index),index));
        local_flag_[eb] = fitflags.at(index);
        std::cout << local_flag_[eb] <<std::endl;
        result.push_back(chi2.at(index));
      } 
      return result;
   }
   
 }
 
 // inialize smooth scale factors :
 void ScaleFactorHelper::InitializeSF()
 {
        PrintDebug("===============================================================");
        PrintDebug("starting fit for scale factors ");
        PrintDebug("==============================================================="); 
        if(debug_flag_ ) {egm2d_.Write(); efficiency_data_.Write(); efficiency_mc_.Write();
            for(int i = 1; i< maxBineta_ +1;i++)
            {    
                SetFitFlag(i);
                TF1 fit;
                // don't smooth graph before fitting:
                TH1F h;
                TGraphErrors g = GetTGraph(i);
                //smooth graph before fitting:
//                  TH1F h = GetHisto(i);
//                  h.Smooth();
//                  TGraphErrors g = GetGraphNumSmoothed(h);
                
                ///////////////////////////////////
                
                graph_.push_back(g);
                if(fit_flag_!=100)
                {
                    fit = FitScaleFactor(g,i);
                    if(i==4 )
                        {
                        TCanvas* tc = new TCanvas("tc","tc",400,400);
                        g.Draw("alp");
                        fit.Draw("same");
                        tc->SaveAs("debug2.pdf");
                       std::cout << "debug2 "<<fit_flag_ << std::endl;
                       std::cout << rangelow_ << " "<< rangeup_ <<std::endl;
                       std::cout << fit.GetName() << std::endl;
                            
                        }
                    
                    fit.SetName(Form("smooth_sf_etabin%i",i));
                    PrintDebug(Form("fitted scale factor for eta bin %i using function %i",i,fit_flag_));
                    fit.Write();
                    smooth_sf_.push_back(fit);
                    num_smooth_sf_.push_back(h);
                     
                    
                }
                else
                {
                    h = GetHisto(i);
                    h.Smooth();
                    h.SetName(Form("smooth_sf_etabin%i",i));
                    h.Write();
                    num_smooth_sf_.push_back(h);
                    smooth_sf_.push_back(fit);
                }
            }
        }
        else{
            for(int i = 1; i< maxBineta_ +1;i++)
            {    
                SetFitFlag(i);
                graph_.push_back(GetTGraph(i));
                if(fit_flag_!=100){ 
                    TH1F h;
                    smooth_sf_.push_back(SetSFFunction(i));
                    num_smooth_sf_.push_back(h);
                }
                else {
                    TF1 fit;
                    num_smooth_sf_.push_back(SetSFHisto(i));
                    smooth_sf_.push_back(fit);
            }
        }
   }
   
 }
 
 TH1F ScaleFactorHelper::GetHisto(int etaBin)
 {
   TH1F* testhisto = new TH1F("test","test",maxBinpt_,egm2d_.GetYaxis()->GetXbins()->GetArray());
        for(int j=1;j<maxBinpt_+1;j++)
        {
         float pt = egm2d_.GetYaxis() -> GetBinCenter(j);
         float eta = egm2d_.GetXaxis() -> GetBinCenter(etaBin);
         std::cout << "test " << std::endl;
         std::cout << pt << std::endl;
         testhisto->Fill(pt, GetSF(pt,eta));
         testhisto->SetBinError(j,GetUncertainty(pt,eta));
         std::cout << "test2 " << std::endl;
        }   
     
   return *testhisto;  
 }
 
 
 void ScaleFactorHelper::InitializeUnc()
 {
        PrintDebug("===============================================================");
        PrintDebug("starting fit for scale factor uncertainties ");
        PrintDebug("==============================================================="); 
   for(int i=1;i<egm2d_.GetXaxis()->GetNbins()+1;i++)
   {
      SetFitFlag(i);
      TF1 fitunc;
      TH1F hunc;
      if(debug_flag_){
      TGraphErrors gunc = GetTGraph(i,1);
      graph_unc_.push_back(gunc);
      if( uncertainty_flag_ >= 1000)
      {    
         if( (input_ ==EGammaInput::electronLoose or input_ ==EGammaInput::electronTight or input_ ==EGammaInput::electronMedium or input_ ==EGammaInput::electronMVA80 or input_ == EGammaInput::electronMVA90 or input_ ==EGammaInput::electronCutBasedVetoID))
        {
            
        float eta = egm2d_.GetXaxis()-> GetBinCenter(i);
        float ymax = FixMaxUncPoint(i,45.);
        float xmax = egm2d_.GetYaxis() -> GetBinCenter(ptBin_);
        float ymin = GetUncertainty(45.,eta);
        float parb = (ymin);
        SetPtBin(45.);
        float parc = egm2d_.GetYaxis() -> GetBinCenter(ptBin_);
        float para = (ymax-parb)/pow(xmax-parc,2);
        
 
        TF1 f = GetUncertaintyFunction(i);
        
    
        f.SetParameter(2,parb);
        f.SetParameter(1,para);
        f.SetParameter(0,parc);
        

        f.SetName(Form("smooth_unc_etabin%i",i));
        smooth_unc_.push_back(f);
        f.Write();
            
        }
        if( (input_ == EGammaInput::photonLoose or input_ == EGammaInput::photonTight or input_ == EGammaInput::photonMedium or input_==EGammaInput::photonMVA90))
        {
            
        float eta = egm2d_.GetXaxis()-> GetBinCenter(i);
        SetPtBin(25.);
        float xmin = egm2d_.GetYaxis() -> GetBinCenter(ptBin_);
        float ymax = FixMaxUncPoint(i,xmin);
        float xmax = egm2d_.GetYaxis() -> GetBinCenter(ptBin_);
        float ymin = GetUncertainty(xmin,eta);
        float para = (ymax-ymin)/(xmax-xmin);
        float parb = ymin;
        
        std::cout << "xmin " << xmin << " ymin " << ymin << " xmax "<< xmax << " ymax " << ymax << std::endl;
        
        TF1 f = GetUncertaintyFunction(i);
        
        f.SetParameter(2,xmin);
        f.SetParameter(1,para);
        f.SetParameter(0,parb);
        

        f.SetName(Form("smooth_unc_etabin%i",i));
        smooth_unc_.push_back(f);
        f.Write();
            
        }
        
      }
      else
      {
          
        fitunc = GetUncertaintyFunction(i);
            fitunc.SetName(Form("smooth_unc_etabin%i",i));
            if (uncertainty_flag_ ==0) fitunc.SetParLimits(2,minsfunc_,10);
            if (uncertainty_flag_ ==1) fitunc.SetParLimits(0,maxsfunc_,1.2*maxsfunc_);
            TFitResult* fitres = (gunc.Fit(&fitunc,"REX0")).Get();
            fitunc.Write();
            smooth_unc_.push_back(fitunc);
     
      }
      }
      else {
          
           SetFitFlag(i);
           graph_.push_back(GetTGraph(i,1));
           fitunc  = SetUncertaintyFunction(i);
           smooth_unc_.push_back(fitunc);
          }
    }
 }
 
 float ScaleFactorHelper::FixMaxUncPoint(int etaBin, float pt_minUnc)
 {
     vector<float> pt;
     vector<float> unc;
     int index = 1;
     float max =-100000;
     float minUnc=0;
     float eta = egm2d_.GetXaxis()-> GetBinCenter(etaBin);
     
     
     SetPtBin(25.);
     float xmin = egm2d_.GetYaxis() -> GetBinCenter(ptBin_);
     float ymin = GetUncertainty(xmin,eta);
     
     
     for(int i=1;i<maxBinpt_+1;i++)
     {
       pt.push_back(egm2d_.GetYaxis() -> GetBinCenter(i));   
       unc.push_back(GetUncertainty(pt.at(i-1),eta));  
       if (TMath::Abs(pt.at(i-1)-pt_minUnc)<4. ){ minUnc= unc.at(i-1);}
           
     }

     
     
     for(int j=0;j<unc.size()-1;j++)
     {
        if (TMath::Abs(pt.at(j)-pt_minUnc)<4. ) continue; 
        float tmp = unc.at(j)/TMath::Abs(pt.at(j)-pt_minUnc); 
        if(input_ == EGammaInput::photonLoose or input_ == EGammaInput::photonMedium or input_ == EGammaInput::photonTight or input_ == EGammaInput::photonMVA90)
        {    
            tmp = ((unc.at(j)-ymin)/TMath::Abs((pt.at(j)-pt_minUnc)));
        }
        if(tmp>max){index=j; max = tmp;}
     }
     SetPtBin(pt.at(index));
     return unc.at(index);
 }
 
 
 
 
 
 
 float ScaleFactorHelper::GetSF(float pT, float superClusterEta)
{
  float sf;
  SetEtaBin(superClusterEta);
  SetPtBin(pT);
  sf = egm2d_.GetBinContent(etaBin_,ptBin_);
  return sf;  
}

float ScaleFactorHelper::GetSFSmooth(float pT, float superClusterEta)
{
   SetEtaBin(superClusterEta);
   float upperRange = egm2d_.GetYaxis()->GetBinLowEdge(maxBinpt_-1)+ egm2d_.GetYaxis()->GetBinWidth(maxBinpt_-1);
   SetPtBin(pT);
   SetFitFlag(etaBin_);
   float sf;
   if(input_ == EGammaInput::electronRecoSF or input_ == EGammaInput::electronVetoPixelSeedEndcap or input_ == EGammaInput::electronVetoPixelSeedBarrel or input_ == EGammaInput::electronVetoPixelSeed or input_ == EGammaInput::electronVetoCSEV) return GetSF(pT,superClusterEta);
   if( fit_flag_ !=100)
   {
        if(pT > upperRange) { return smooth_sf_.at(etaBin_-1).Eval(upperRange);}
        if(pT < rangelow_) return smooth_sf_.at(etaBin_-1).Eval(rangelow_);
        sf = smooth_sf_.at(etaBin_-1).Eval(pT);
        return sf; 
   }
   if (fit_flag_ ==100)
   {
        if(pT > upperRange) return num_smooth_sf_.at(etaBin_-1).GetBinContent(maxBinpt_);
        if(pT < rangelow_) return num_smooth_sf_.at(etaBin_-1).GetBinContent(1);
        SetPtBin(pT);
        sf = num_smooth_sf_.at(etaBin_-1).GetBinContent(ptBin_);
        return sf; 
   }
   throw my_range_error("smoothed scale factor functions not correctly loaded");   
}

TF1 ScaleFactorHelper::SetSFFunction(int etaBin)
{
    TFile file(file_.c_str(),"READ");
    if (file.IsZombie()) throw file_not_found();
    TF1 func = *dynamic_cast<TF1*>(file.Get(Form("smooth_sf_etabin%i",etaBin)));
    if( func.IsZombie()) throw my_range_error("smoothed sf fit function could not be opened");
    file.Close();
    return func;
}

TF1 ScaleFactorHelper::SetUncertaintyFunction(int etaBin)
{
   TFile file(file_.c_str(),"READ");
    if (file.IsZombie()) throw file_not_found();
    TF1 func = *dynamic_cast<TF1*>(file.Get(Form("smooth_unc_etabin%i",etaBin)));
    if( func.IsZombie()) throw my_range_error("smoothed sf fit function could not be opened");
    file.Close();
    return func; 
}

TH1F ScaleFactorHelper::SetSFHisto(int etaBin)
{
    TFile file(file_.c_str(),"READ");
    if (file.IsZombie()) throw file_not_found();
    TH1F func = *dynamic_cast<TH1F*>(file.Get(Form("smooth_sf_etabin%i",etaBin)));
    if( func.IsZombie()) throw my_range_error("smoothed sf histogram could not be opened");
    file.Close();
    return func;
}


float ScaleFactorHelper::GetUncertaintySmooth(float pT, float superClusterEta)
{
  SetEtaBin(superClusterEta);
  SetPtBin(pT);
  SetFitFlag(etaBin_);
  if(uncertainty_flag_==100){ std::cout << "warning : scale factor must be numerically smoothed in order to use this uncertainty option" << std::endl; return num_smooth_sf_.at(etaBin_).GetBinError(ptBin_);} 
 
  if(input_ == EGammaInput::electronRecoSF or input_ == EGammaInput::electronVetoPixelSeedEndcap or input_ == EGammaInput::electronVetoPixelSeedBarrel or input_ == EGammaInput::electronVetoPixelSeed or input_ == EGammaInput::electronVetoCSEV) return GetUncertainty(pT,superClusterEta);
  float unc = smooth_unc_.at(etaBin_-1).Eval(pT);
  float minimal_uncertainty =minUnc_[etaBin_]; 
  if(unc < minimal_uncertainty){ unc = minimal_uncertainty;}
  if((input_ ==EGammaInput::electronLoose or input_ ==EGammaInput::electronTight or input_ ==EGammaInput::electronMedium or input_ ==EGammaInput::electronMVA80 or input_ == EGammaInput::electronMVA90 or input_ == EGammaInput::electronCutBasedVetoID)){
   if(pT < 45. and unc > maxUnc_sm45_[etaBin_]){unc = maxUnc_sm45_[etaBin_];}   
   if(pT > 45. and unc > maxUnc_l45_[etaBin_]){unc = maxUnc_l45_[etaBin_];}    
  }
  if(unc > maxUnc_[etaBin_]) {unc = maxUnc_[etaBin_];} // maximal uncertainty given by the maximal uncertainty measured
  return unc;  
}


float ScaleFactorHelper::GetUncertainty(float pT, float superClusterEta)
{
 float sf_unc;
 SetEtaBin(superClusterEta);
 SetPtBin(pT);
 sf_unc = egm2d_.GetBinError(etaBin_,ptBin_);
  // add additional uncertainty to electron reconstruction scale-factor:
  if(input_ == EGammaInput::electronRecoSF)
  {
      //std::cout << " add additional uncertainty for electron REco " << std::endl;
        if(pT < 20. or pT > 80.)
        {
          //std::cout << "pt is smaller 20 or larger 80 GeV :  " <<  0.01*egm2d_.GetBinContent(etaBin_,ptBin_) << std::endl; 
          sf_unc = TMath::Sqrt(pow(sf_unc,2) + pow(0.01*egm2d_.GetBinContent(etaBin_,ptBin_),2));   
        }
  }
  
  return sf_unc;
}


void ScaleFactorHelper::SetEtaBin(float superClusterEta)
{
  if (superClusterEta != input_eta_){  
  input_eta_ = superClusterEta;  
  etaBin_ = egm2d_.GetXaxis() -> FindBin(superClusterEta);
  maxBineta_ = maxBineta_;
  if (etaBin_==0 or etaBin_==maxBineta_+1) 
  {std::string err = "tried to evaluate scale factor for |eta| >"; err.append(std::to_string(TMath::Abs(egm2d_.GetXaxis()->GetBinLowEdge(1)))); throw my_range_error(err);}
  }
}

void ScaleFactorHelper::SetPtBin(float pT)
{
  if (pT < 0 ) throw my_range_error("pT < 0 not allowed");  
  if (pT != input_pt_){  
  input_pt_ = pT;  
  ptBin_ = egm2d_.GetYaxis() -> FindBin(pT);
  maxBinpt_ = egm2d_.GetYaxis() -> GetNbins();
  if(ptBin_ ==0) {std::string err = "tried to evaluate scale factor for pt < "; err.append(std::to_string(egm2d_.GetYaxis()->GetBinLowEdge(ptBin_+1)));
      ptBin_=1;
      /*throw my_range_error(err);*/}
  if(ptBin_ >= maxBinpt_ and maxBinpt_ > 1 ) ptBin_ = maxBinpt_-1;  // last pt bin is only used as control but the scale-factor here should not be used because of limited statistics
  if(maxBinpt_ ==1) ptBin_ = 1; // this histogram has only one bin, the setting from above does not apply here!
  }
}


TGraphErrors ScaleFactorHelper::GetTGraph(int etaBin,bool forUncertainty)
{
    int max = egm2d_.GetYaxis()-> GetNbins();
    float eta = egm2d_.GetXaxis()-> GetBinCenter(etaBin);
    TGraphErrors *g = new TGraphErrors(max);
    float findmin=10000;
    float findmax=-100000;
    float findmaxL45=-100000;
    float findmaxS45=-100000;
    float ptmin =0;
    float ptmax =0;
    for(int i=1; i< max+1;i++)
    {
        float pt = egm2d_.GetYaxis() -> GetBinCenter(i);
        float sf = GetSF(pt,eta);
        float pt_unc = egm2d_.GetYaxis() -> GetBinWidth(i) /2.;
        float sf_unc = GetUncertainty(pt,eta);
        if (sf_unc < findmin) {findmin =sf_unc; ptmin = pt;}
        if (sf_unc > findmax) {findmax =sf_unc; ptmax = pt;}
        
        if (sf_unc > findmaxS45 and pt <45.) {findmaxS45 =sf_unc; }
        if (sf_unc > findmaxL45 and pt >45.) {findmaxL45 =sf_unc; }
        
        
        
        if (forUncertainty)
        {
           //float sf_unc_unc = sf_unc/TMath::Sqrt(2*(maxBinpt_*maxBineta_-2)); 
           //float sf_unc_unc = TMath::Sqrt( sf_unc*sf_unc*TMath::Sqrt( 2/(TMath::Abs(sf_unc-1))));
            g->SetPoint(i-1,pt,sf_unc);
           //g->SetPointError(i-1,pt_unc,sf_unc_unc);
        }
        else
        {
            g->SetPoint(i-1,pt,sf);
            g->SetPointError(i-1,pt_unc,sf_unc);
        }
    }
    minsfunc_=findmin;
    maxsfunc_=findmax;
    maxUnc_.insert ( std::pair<int,float>( etaBin_, findmax ));
    minUnc_.insert ( std::pair<int,float>( etaBin_, findmin ));
    maxUnc_sm45_.insert ( std::pair<int,float>( etaBin_, findmaxS45 ));
    maxUnc_l45_.insert ( std::pair<int,float>( etaBin_, findmaxL45 ));
    
    ptmaxUnc_.insert ( std::pair<int,float>( etaBin_, ptmax ));
    ptminUnc_.insert ( std::pair<int,float>( etaBin_, ptmin ));
    
    
    return *g;
}

TF1 ScaleFactorHelper::FitScaleFactor(TGraphErrors g, int etaBin)
{
  TF1 fit = GetFitFunction(etaBin);
  TFitResult* res = (g.Fit(&fit,"SR")).Get();
  res->Print();  
  return fit;  
}


void ScaleFactorHelper::DrawAll()
{
  if (debug_flag_){  
  TCanvas* drawAll = new TCanvas("drawAll","drawAll",600,1000);
  drawAll->Divide(2,maxBineta_/2.);
  std::vector<TGraph*> gs={};
  std::vector<TGraph*> gs_UncUp={};
  std::vector<TGraph*> gs_UncDown={};
   for(int i = 1; i< maxBineta_ +1;i++)
   {   float max =-10;
       float min =1000;
       SetFitFlag(i);
       TPad* pad =(TPad*) drawAll->cd(i);
       
       graph_.at(i-1).GetXaxis()->SetTitle("pT (GeV)");
       graph_.at(i-1).GetYaxis()->SetTitle("scale factor");
       graph_.at(i-1).SetLineColor(kBlue);
       graph_.at(i-1).SetLineWidth(2);
       graph_.at(i-1).SetMarkerColor(kBlack);
       graph_.at(i-1).SetMarkerStyle(8);
       
       float eta = egm2d_.GetXaxis() -> GetBinCenter(i);
       TGraph * g = new TGraph(1);
       TGraph* gup = new TGraph(1);
       TGraph* gdown = new TGraph(1);
       for (int j=0; j<800;j++){
           float pt = rangelow_+j*1.5;
           float sf = GetSFSmooth(pt,eta);
           float unc = GetUncertaintySmooth(pt,eta);
           g->SetPoint(j,pt,sf);
           gup->SetPoint(j,pt,sf+unc);
           gdown->SetPoint(j,pt,sf-unc);
           if (sf+unc > max) max = sf+unc;
           if (sf-unc < min) min = sf-unc;


       }
       gs.push_back(g);
       gs_UncUp.push_back(gup);
       gs_UncDown.push_back(gdown);
    pad->SetLogx();
    if (graph_.at(i-1).Eval(500.) > max) max = graph_.at(i-1).Eval(500.);
    if (graph_.at(i-1).Eval(25.) < min) min = graph_.at(i-1).Eval(25.);
    
    graph_.at(i-1).SetMaximum(1.05*max);
     graph_.at(i-1).SetMinimum(0.95*min);
    graph_.at(i-1).Draw("ALP");
    gs[i-1]->SetLineColor(kMagenta);
    gs[i-1]->SetLineWidth(2);
    gs[i-1]->Draw("Lsame");
    gs_UncUp[i-1]->SetLineColor(kGreen);
    gs_UncDown[i-1]->SetLineColor(kGreen);
    gs_UncUp[i-1]->SetLineWidth(2);
    gs_UncDown[i-1]->SetLineWidth(2);
    gs_UncUp[i-1]->Draw("Lsame");
    gs_UncDown[i-1]->Draw("Lsame");
    
    TPaveText* addInfo = new TPaveText(0.9,0.54,0.64,0.4,"NDC");
    addInfo->SetFillColor(0);
    addInfo->SetLineColor(0);
    addInfo->SetFillStyle(0);
    addInfo->SetBorderSize(0);
    addInfo->SetTextFont(42);
    addInfo->SetTextSize(0.040);
    addInfo->SetTextAlign(12);
    addInfo->AddText(Form("#eta =  %.1f ", eta));
    if(fit_flag_!=100) addInfo->AddText(smooth_sf_.at(i-1).GetExpFormula());
    addInfo->Draw("same");

        
   }
   std::string modifier = GetString(input_);
   drawAll->Modified();
   drawAll->SaveAs(Form("All_%s.pdf",modifier.c_str()));
   std::cout << " draw all scale factors All_"<< modifier <<std::endl;
 
  }
}

void ScaleFactorHelper::DrawAllUncertainty()
{
  if (debug_flag_){  
  TCanvas* drawAll = new TCanvas("drawAll","drawAll",600,1000);
  drawAll->Divide(2,maxBineta_/2.); 
  std::vector<TGraph*> gs={};
   for(int i = 1; i< maxBineta_ +1;i++)
   {  
       SetFitFlag(i);
       drawAll->cd(i);
       graph_unc_.at(i-1).GetXaxis()->SetTitle("pT (GeV)");
       graph_unc_.at(i-1).GetYaxis()->SetTitle("(unc sf)/sf");
       graph_unc_.at(i-1).SetLineColor(kBlue);
       graph_unc_.at(i-1).SetLineWidth(2);
       graph_unc_.at(i-1).SetMarkerColor(kBlack);
       graph_unc_.at(i-1).SetMarkerStyle(8);
       
        float max = -100000;
        int bin = 0;
        bool which =0;
        float eta = egm2d_.GetXaxis() -> GetBinCenter(i);
        TGraph * g = new TGraph(1);
        for (int j=0; j<800;j++){
           float pt = rangelow_+j*1.5;
           float u = GetUncertaintySmooth(pt,eta);
           g->SetPoint(j,pt,u);
           if(u > max) {max = u;}

       }  
    double* y = graph_unc_.at(i-1).GetY();
    for( int j= 0;j< sizeof(y)/sizeof(y[0]);j++)
    {
     if(y[j]+graph_unc_.at(i-1).GetErrorY(j) > max) {max = y[j]; bin = j; which =1;}   
    }
    if(which)
    {
     max = max+ 1.1* graph_unc_.at(i-1).GetErrorY(bin);  
    }
    else {max = max*1.05;}
    gs.push_back(g);
    gs[i-1]->SetLineColor(kMagenta);
    gs[i-1]->SetLineWidth(2);
    graph_unc_.at(i-1).SetMaximum(max);
    graph_unc_.at(i-1).Draw("ALP");
    gs[i-1]->Draw("Lsame");
    TPaveText* addInfo = new TPaveText(0.9,0.54,0.64,0.4,"NDC");
    addInfo->SetFillColor(0);
    addInfo->SetLineColor(0);
    addInfo->SetFillStyle(0);
    addInfo->SetBorderSize(0);
    addInfo->SetTextFont(42);
    addInfo->SetTextSize(0.040);
    addInfo->SetTextAlign(12);
    addInfo->AddText(Form("#eta =  %.1f ", eta));
    addInfo->AddText(smooth_unc_.at(i-1).GetExpFormula());
    addInfo->Draw("same");

        
    }
   std::string modifier = GetString(input_);
   drawAll->Modified();
   drawAll->SaveAs(Form("AllUnc_%s.pdf",modifier.c_str()));
 
  }
}



float ScaleFactorHelper::GetEfficiency(float pT, float superClusterEta,bool isData)
{
    if( pT!=input_pt_ or superClusterEta != input_eta_)
    {
        SetEtaBin(superClusterEta);
        SetPtBin(pT);
    }
    std::string name(efficiency_data_.GetName());
    if(name.compare("NONE")==0) throw my_range_error("efficiency of this quantity is undefined");
    if(isData) return efficiency_data_.GetBinContent(etaBin_,ptBin_); 
    return efficiency_mc_.GetBinContent(etaBin_,ptBin_); 
    
}

 TF1 ScaleFactorHelper::GetFitFunction(int etaBin)
 {
     TF1* f;
     if(fit_flag_==0)
     {
         f = new TF1("atan","[0]+ atan([1]*x)",rangelow_,rangeup_);
         ndof_.insert ( std::pair<int,int>( fit_flag_, 2 ));
         return *f;
     }
     if(fit_flag_==1)
     {
         f = new TF1("overX","[0]+ [1]*(1/x)",rangelow_,rangeup_);
         f->SetParLimits(1,-10000.,0.);
         ndof_.insert ( std::pair<int,int>( fit_flag_, 2 ));
         return *f;
     }
     if(fit_flag_==2)
     {
         f = new TF1("overX2","[0]+ [1]*(1/(x*x))",rangelow_,rangeup_);
         ndof_.insert ( std::pair<int,int>( fit_flag_, 2 ));
         f->SetParLimits(1,-6000,0.);
         return *f;
     }
     if(fit_flag_==3)
     {
         f = new TF1("overXa","[0] + [1]*(1/x^[2])",rangelow_,rangeup_);
         ndof_.insert ( std::pair<int,int>( fit_flag_, 3 ));
         return *f;
     }
     
     
     if(fit_flag_==4)
     {
       f = new TF1("erf","TMath::Erf(([1]*x-[0]))",rangelow_,rangeup_);
       ndof_.insert ( std::pair<int,int>( fit_flag_, 2 ));
       return *f;  
     }
     
     if(fit_flag_==5)
     {
       f = new TF1("overexp","[2]/(1+exp(-[1]*(x-[0])))",rangelow_,rangeup_);
       ndof_.insert ( std::pair<int,int>( fit_flag_, 3 ));
       return *f;  
     }
     if(fit_flag_==6)
     {
       f = new TF1("pol2","[0]+[1]*x ",rangelow_,rangeup_);
       ndof_.insert ( std::pair<int,int>( fit_flag_, 2 ));
       return *f;  
     }
     if(fit_flag_==7)
     {
       f = new TF1("pol3","[0]+[1]*x +[2]*x*x ",rangelow_,rangeup_);
       ndof_.insert ( std::pair<int,int>( fit_flag_, 3 ));
       return *f;  
     }
     if(fit_flag_==8)
     {
       f = new TF1("pol4","[0]+[1]*x +[2]*x*x + [3]*x*x*x",rangelow_,rangeup_);
       ndof_.insert ( std::pair<int,int>( fit_flag_, 4 ));
       return *f;  
     }
     if(fit_flag_==9)
     {
       f = new TF1("1+1/sqrt(x)","[0]+[1]/TMath::Sqrt(x)",rangelow_,rangeup_);
       ndof_.insert ( std::pair<int,int>( fit_flag_, 2 ));
       return *f;  
     }
     
     if(fit_flag_==10)
     {
       f = new TF1("test","([1]+ [0]*TMath::Sqrt(x)) ",rangelow_,rangeup_);
       ndof_.insert ( std::pair<int,int>( fit_flag_, 2 ));
       return *f;  
     }
     
     
     if(fit_flag_==11)
     {
         f = new TF1("line","[0]",rangelow_,rangeup_);
         ndof_.insert ( std::pair<int,int>( fit_flag_, 1 ));
         return *f;
     }
     
     std::cout << fit_flag_ <<std::endl;
     throw my_range_error("wrong fitting flag used in scale factor fit");
 }

 
 
TF1 ScaleFactorHelper::GetUncertaintyFunction(int etaBin)
 {
     TF1* f;
     if(uncertainty_flag_==0 or uncertainty_flag_==1000)
     {
         f = new TF1("func","[2] + [1]*(x-[0])^(2)",rangelow_,unc_range_);
         f->SetParLimits(0,30,55);
         f->SetParLimits(1,0,10000);
         return *f;
     }
     if(uncertainty_flag_==1)
     {
         f = new TF1("func","[0]+ ([1]/x^2)",rangelow_,unc_range_);
         return *f;
     }  
     
     if(uncertainty_flag_==6 or uncertainty_flag_==1006)
     {
       f = new TF1("line2","[0]+[1]*(x -[2] )",rangelow_,unc_range_);
       return *f;  
     }
     
     
     throw my_range_error("wrong fitting flag used in uncertainty smoothing");
 }
 
    
