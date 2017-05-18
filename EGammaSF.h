#ifndef EGAMMASF_H
#define EGAMMASF_H

#include "TH2F.h"
#include "TGraphErrors.h"
#include<iostream>
#include<fstream>
#include<exception>
#include<string>
using namespace std;

// class of allowed input values:

enum class EGammaInput : int {electronRecoSF = 1 , electronMVA90 =2, electronMVA80 =3 , electronTight =4, electronMedium =5, electronLoose = 6, photonMVA90 = 7, photonLoose=8, photonTight=9, photonMedium = 10};



// some operator overloads for debugging purposes

template<typename T>
std::ostream& operator<<(typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e)
{
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
}



// exceptions:

class file_not_found: public std::exception
{
  virtual const char* what() const throw()
  {
    return "File could not be openened. Please check if it exists.";
  }
};

class my_range_error: public std::exception
{
  public:
  my_range_error(std::string message) {
    std::cout << message << std::endl;
  }
};



// small parser class to include the config file config.txt
// this file only has to be changed if the filenames of the input root-files or the contained histogram names need to be changed
 class ConfigParser {
     
 public:
     ConfigParser(std::string config_name, EGammaInput sf){
      std::ifstream  config(config_name);
      if(!config.good()) throw file_not_found();
      ReadFile(config, sf);
         
     };
     ~ConfigParser(void);
     std::string NameEffMC(){ return effmc_;};
     std::string NameEffData(){return effdata_;};
     std::string NameSF(){ return sf_;};
     std::string NameFile(){ return filename_;};
     
 private:
     void ReadFile(std::ifstream &config, EGammaInput sf) {
      std::string cat = GetString(sf);
      for(;;) {
        std::string line;
        std::getline(config, line);
        if(!config) break;
        // read from stream
        //std::cout << line << std::endl;
        if(line.find(cat)!=std::string::npos)
        {
            for(;;){
                    std::string line2;
                    std::getline(config,line2);
            if( line2.find("name")!=std::string::npos)
            {
                 std::size_t n1 = line2.find_first_of("=");
                 std::size_t n2 = line2.find_first_of("\n");
                 filename_ = line2.substr(n1+1,n2);
            }
             if( line2.find("EffMC")!=std::string::npos)
            {
                 std::size_t n1 = line2.find_first_of("=");
                 std::size_t n2 = line2.find_first_of("\n");
                 effmc_ = line2.substr(n1+1,n2);
            }
             if( line2.find("EffData")!=std::string::npos)
            {
                 std::size_t n1 = line2.find_first_of("=");
                 std::size_t n2 = line2.find_first_of("\n");
                 effdata_ = line2.substr(n1+1,n2);
            }
             if( line2.find("SF")!=std::string::npos)
            {
                 std::size_t n1 = line2.find_first_of("=");
                 std::size_t n2 = line2.find_first_of("\n");
                 sf_ = line2.substr(n1+1,n2);
            }
            
                if( line2.find("]")!=std::string::npos) break;
            }
            
        }    
        }   
      //std::cout << filename_ << " " << effmc_ << " "<< effdata_ << " "<< sf_ << " "<< cat << std::endl;   
     };
     
     std::string effmc_ = "EGamma_EffMC2D";
     std::string effdata_ = "EGamma_EffData2D";
     std::string sf_ = "EGamma_SF2D";
     std::string filename_ ="";
     
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
       throw my_range_error("invalid input parameters");
     };
 
 };














class ScaleFactorHelper {

public:
    
//============================================ 
// Constuctor -> should be used once to initialize variables
// should initialize:
//      open right root file
//      extract TH2D containing scale factors
//      set fit_flag_ for which function should be used for the fit 
//      set uncertainty_flag for which function should be fixed for uncertainties    
//      fit pT for each eat bin
//      fit uncertainty for each eta bin
//      -> initialize functions (TF1) for each eta bin   (name should contain number of eta bin)
// should contain:
//      mode for debugging -> save canvas with fits etc.
//      mode for analysers    

   ScaleFactorHelper( EGammaInput what, bool debugging = 0 );
//============================================



   float GetEfficiency(float pT, float superClusterEta,bool isData);     // => return efficiency in MC or data 
// //============================================
// // return right bin content of TH2D histo   
    float GetSF(float pT, float superClusterEta);
// //============================================   
   
                                                                            // => functions have to handle pT/eta values in underflow or overflow
//============================================                               // => they should initialize pT,eta,sf, unc private variables and only recalculate when pT/eta given is different than in
 // get error of right bin of TH2D histo                                            iteration before! (saves time)
 // add additional uncertainties                                              // => should throw and error if the class wasn't properly initialized in the constructor
 // return uncertainty
       
    float GetUncertainty(float pT, float superClusterEta);
// //============================================   
//       
//    
//    float GetSFSmooth(float pT, float superClusterEta);                    // => get smoothed out SF and uncertainty using the functions fitted during construction
//    float GetUncertaintySmooth(float pT, float superClusterEta);           // => should do something sensible with under/overflow
//                                                                              // => throw error for uninitialized class
//    
// //============================================
// // close all files    
    ~ScaleFactorHelper( void );
// //============================================
//    
 private:

//===============================================
// sets private variable pTBin_ / etaBin_ i.e. the bin number of given pT/eta value in TH2D
// this function has to handle the underflow/overflow properly
void SetEtaBin(float superClusterEta);
// last bin in pT is not used as SF (low statistics) but just as controll -> take pT bin before instead!
void SetPtBin( float pT);
//===============================================    


//===============================================
// get all pT values (and uncertainties) corresponding to the bin in eta etaBin and put them into a TGraphErrors
// this function has to handle additional uncertainties applied depending on pT (for example 1% extra uncertainty for pT < 20 GeV) see:  https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
 TGraphErrors GetTGraph(int etaBin);
//===============================================    
 
// //===============================================
// // fit a predefined function to input graph and return a clone of the fit function with fixed values
// TF1 FitScaleFactor(TGraphErrors g_eta);
// //===============================================
// 
// 
// //===============================================
// // fix uncertainty function from input graph => start as in  https://indico.cern.ch/event/482669/contributions/1991571/attachments/1253608/1849618/20160406CharlotScaleFactors.pdf slide 6
// // add additional functions later /fit?
// // this assumes symetric error bars!
// TF1 SmoothUncertainty(TGraphErrors g_eta);
// //===============================================
// 
// 
// //===============================================
// // return a fit function according to fit_flag_ i.e. 0 simple factor, 1 polynomial degree 1 or other functions that seem resonable and should be tried out
// TF1 GetFitFunction();
// //===============================================
// 
// 
// //===============================================
// // return a fit function according to uncertainty_flag_ (functions that seem resonable and should be tried out )
// // => SmoothUncertainty needs to react differently to the funcitons (most likely)
// TF1 GetUncertaintyFunction();
// //===============================================
 
 
 //===============================================
 // if debug_flag_ set draw a canvas with TGraph of SF and uncertainties (and functions) belonging to etaBin_
 void DrawSF(TGraphErrors g, float eta);
 //===============================================
 
 
void PrintDebug(std::string stuff)
{
  if(debug_flag_)
  {
    std::cout << stuff << std::endl;   
  }
};
// 
// 
// 
// //               private variables
// //===============================================
// //===============================================
        bool debug_flag_ =0; // if set to 1 print/draw additional information
//        int fit_flag_ =0;    // use different functions for the fit of SF (smoothing)
//        int uncertainty_flag_ =0; // use different assumptions to estimate the scale factor uncertainties ( smoothing )
//        

        EGammaInput input_;
        int etaBin_ = -99;                 // safe bin number of etaBin currently used
        int ptBin_  = -99;                  // safe bin number of pT bin currently used
        int maxBinpt_;
        int maxBineta_;
//        std::string particle_;
//        std::string identification_;     // safe initialization values
//        std::string working_point_;
//        
        float input_pt_;
        float input_eta_;              // safe eta/pt values SF are currently evaluated for
//        
        TH2F egm2d_;                    // safe 2d histo containing scale factors
        TH2F efficiency_mc_;            // safe 2d histo containing efficiency in mc
        TH2F efficiency_data_;          // safe 2d histo containing efficiency in data
//        
//        std::vector<TGraphErrors> graphs_; // safe graphs containing SF+Unc per eta etaBin
//        std::vector<TF1> smooth_sf_; // safe fit functins for smooth scale factors for each eta bin 
//        std::vector<TF1> smooth_unc_; // safe uncertainty functions for smooth uncertianties and for each eta bin
// 
// //===============================================
// //===============================================

};




#endif
