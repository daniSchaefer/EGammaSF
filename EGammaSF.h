#ifndef EGAMMASF_H
#define EGAMMASF_H

#include <iostream>
#include <exception>

// problem input values => c++ standard is to typesafe i.e.:

enum class EGammaInput : int {electronRecoSF = 1 , electronMVA90 =2, electronMVA80 =3 , electronTight =4, electronMedium =5, electronLoose = 6, photonMVA90 = 7, photonLoose=8, photonTight=9, photonMedium = 10};























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

   ScaleFactorHelper(std::string particle, std::string ID, std::string workingPoint, bool debugging = 0 );
//============================================


   ScaleFactorHelper( void); // if function is used like this just print the help and throw a warning that initialization is missing!
   void Initialize(std::string particle, std::string ID, std::string workingPoint ); // does same initialization as constructor with arguments -> i.e. constructor just uses Initialize function
   
   
//============================================   
// print the values that are used for initialization  (for debugging) 
// also print which functions will be used for fitting !   
   void Print()
   {
     if( !is_initialized_)
     {
        std::cout << "Warning: class is not initialized properly!" << std::endl; 
     }
     else
     {
       std::cout << "Calculating scale factors for "<< particle_ << " using " << identification_ << " working point " << working_point_ << std::endl;  
       if(debug_flag_)
       {
        std::cout << "debugging mode on" << std::end;   
       }
     }
   };
//============================================   
   
   
//============================================
// print quick overview on how to use the function, how to initialize it etc.   
   void PrintHelp();
//============================================
  
   
   double GetEfficiency(double pT, double superClusterEta, std::string mode);     // => return efficiency in MC or data 
//============================================
// return right bin content of TH2D histo   
   double GetSF(double pT, double superClusterEta);
//============================================   
  
                                                                           // => functions have to handle pT/eta values in underflow or overflow
//============================================                               // => they should initialize pT,eta,sf, unc private variables and only recalculate when pT/eta given is different than in
// get error of right bin of TH2D histo                                            iteration before! (saves time)
// add additional uncertainties                                              // => should throw and error if the class wasn't properly initialized in the constructor
// return uncertainty
      
   double GetUncertainty(double pT, double superClusterEta);
//============================================   
      
   
   double GetSFSmooth(double pT, double superClusterEta);                    // => get smoothed out SF and uncertainty using the functions fitted during construction
   double GetUncertaintySmooth(double pT, double superClusterEta);           // => should do something sensible with under/overflow
                                                                             // => throw error for uninitialized class
   
//============================================
// close all files    
   ~ScaleFactorHelper( void );
//============================================
   
private:
    
//===============================================
// should handle cases for particle
// same for ID
// handle error for misspelled/wrong particle 
// => set the initialization variables i.e. particle_ etc.
//should handle each wrong argument seperately:
//      for neither electron nor photon message should contain the two possibilities (a kind of how to use message)
//      wrong WP -> should deliver list of all possible working points    
void HandleInput(std::string particle, std::string ID, std::string workingPoint );
//===============================================    
    
//============================================
// pick right root-file to extract SF from
    // initialized variables:
    //      electron/photon
    //      MVA/cut based
    //      Working point
    //has to handle:
    //      missing files 
void OpenFile( );                                  // root TFile: if constuctor fails IsZombie() will return true

                                                  // try{ }
   //============================================


//===============================================
// sets private variable pTBin_ / etaBin_ i.e. the bin number of given pT/eta value in TH2D
// this function has to handle the underflow/overflow properly
void SetEtaBin(double superClusterEta);
// last bin in pT is not used as SF (low statistics) but just as controll -> take pT bin before instead!
void SetPtBin( double pT);
//===============================================    


//===============================================
// get all pT values (and uncertainties) corresponding to the bin in eta etaBin and put them into a TGraphErrors
// this function has to handle additional uncertainties applied depending on pT (for example 1% extra uncertainty for pT < 20 GeV) see:  https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
TGraphErrors GetTGraph(int etaBin);
//===============================================    

//===============================================
// fit a predefined function to input graph and return a clone of the fit function with fixed values
TF1 FitScaleFactor(TGraphErrors g_eta);
//===============================================


//===============================================
// fix uncertainty function from input graph => start as in  https://indico.cern.ch/event/482669/contributions/1991571/attachments/1253608/1849618/20160406CharlotScaleFactors.pdf slide 6
// add additional functions later /fit?
// this assumes symetric error bars!
TF1 SmoothUncertainty(TGraphErrors g_eta);
//===============================================


//===============================================
// return a fit function according to fit_flag_ i.e. 0 simple factor, 1 polynomial degree 1 or other functions that seem resonable and should be tried out
TF1 GetFitFunction();
//===============================================


//===============================================
// return a fit function according to uncertainty_flag_ (functions that seem resonable and should be tried out )
// => SmoothUncertainty needs to react differently to the funcitons (most likely)
TF1 GetUncertaintyFunction();
//===============================================


//===============================================
// if debug_flag_ set draw a canvas with TGraph of SF and uncertainties (and functions) belonging to etaBin_
void DrawSF();
//===============================================


void PrintDebug(std::string stuff)
{
  if(debug_flag_)
  {
    std::cout << stuff << std::endl;   
  }
};



//               private variables
//===============================================
//===============================================
       bool debug_flag_ =0; // if set to 1 print/draw additional information
       int fit_flag_ =0;    // use different functions for the fit of SF (smoothing)
       int uncertainty_flag_ =0; // use different assumptions to estimate the scale factor uncertainties ( smoothing )
       bool is_initialized_ =0;  // if initialization worked this is set to 1 => if 0 all Get* function throw error
       
       
       int etaBin_ = -99;                 // safe bin number of etaBin currently used
       int pTBin_  = -99;                  // safe bin number of pT bin currently used
       std::string particle_;
       std::string identification_;     // safe initialization values
       std::string working_point_;
       
       double input_pt_;
       double input_eta_;              // safe eta/pt values SF are currently evaluated for
       
       TH2D* egm2d_;                    // safe 2d histo containing scale factors
       TH2D* efficiency_mc_;            // safe 2d histo containing efficiency in mc
       TH2D* efficiency_data_;          // safe 2d histo containing efficiency in data
       
       std::vector<TGraphErrors> graphs_; // safe graphs containing SF+Unc per eta etaBin
       std::vector<TF1> smooth_sf_; // safe fit functins for smooth scale factors for each eta bin 
       std::vector<TF1> smooth_unc_; // safe uncertainty functions for smooth uncertianties and for each eta bin

//===============================================
//===============================================

};




#endif
