# EGammaSF


provides simple class to help analysist extract the scale-factors provided by EGamma POG.
use makefile and test.cpp to get minimal working example


## Setup for users:    
    For the general setup you just need the two files EGammaSF.cc and EGammaSF.h as well as the root file containing the scale factors provided by the EGamma POG, for example for the electron tight ID scale factors you need the file "egammaEffi_EGM2D_electron_cutbased_tightID.root".  For an example on how to include the class into a short c++ program see test.cpp.
    The SF can then be simply calculated by creating an instance of the sf helper in your code:
        ScaleFactorHelper* sf = new ScaleFactorHelper(EGammaInput::electronTight);
        sf->GetSF(pt,eta); -> returns the binned sf for a given pt and eta
        sf->GetSFSmooth(pt, eta); -> returns a smoothed sf for a given pt and eta
        sf->GetUncertaintySmooth(pt, eta); -> returns the binned uncertainty for given pt, eta
        sf->GetUncertainty(pt, eta);       -> return smoothed uncertainty for a given pt, eta
        
    NOTE: The iniatialization of the sf class needs to be done only once per scale factor!
    NOTE: Check the naming convention of the input ROOT-files in config.txt -> name=... if you want to change the name this line in the config file has to be adjusted
    
    
    
    
## Class properties:
     needs config file in which the filename for the corresponding root file, and the names of the histograms are specified 
     construct with ScaleFactorHelper( EGammaInput, bool )
        - first argument determines which scale factors are used. possibilities are:
                - electronCutBasedVetoID
                - electronRecoSF
                - electronLoose
                - electronTight
                - electronMedium
                - electronMVA80
                - electronMVA90
                - photonLoose
                - photonTight
                - photonMedium 
                - photonMVA90 
         - second argument set to one for more output for debugging, default 0       

     after initialization scale factors/ uncertainties / efficiencies can be calculated with the method GetSF(pt,eta), GetUncertainty(pt,eta), GetEfficiency(pt,eta,isData)
     
 

 
## For Egamma Developers:
    only in the debugging stage, fits are made for the smoothing of sf and the smoothed uncertainties are calculated and saved to the same file that contains the sf histogram
    if the class is used in user mode, those functions are read from a ROOT-file that has to be provided by Egamma 
 
     
### structure of config-file:
    all arguments for one input have to started with "[<name of EGammaInput>" and stopped with "]"
    general syntax is "<name of variable to set>=<variable>" without putting any spaces in between
    minimal example:
        [electronLoose
            name=nameOfRootfile.root
        ]
        the name of the root-file that contains the scale factors has to be stated here! every other variable has a default
        
### setting the fit-flag:
    when using the FitFunction=<number> option in the config file the default fit function for the sf is changed
    functions that are available are the following:
    0 -> [0]+ atan([1]*pt)
    1 -> [0] + [1]/x
    2 -> [0] + [1]/x^2
    3 -> [0] + [1]/x^[2]
    4 -> Erf([1]x-[0])
    5 -> [1]/exp(-[2]*(x-[0]))
    6 -> [0] + [1] x
    7 -> [0] + [1] x + [2]*x^2
    8 -> [0] + [1] x + [2]*x^2 +[3]*x^3
    9 -> [0] +[1]/sqrt(x)
    10-> [0] +[1]*sqrt(x)
    11-> [0]
    
    added option to use numerical smoothing instead of fitting a function using:
    FitFunction=100
    
    different fit functions/ smoothing can be used in different eta regions:
    
          using : LocalFitFunction<1<3=1
          this will use a 1/x fit for the 1,2nd and third bin in eta
          all bins not covered by this are fitted with the function specified in FitFunction or the default 0 if no fit function is specified
          an arbitrary amount of local fit functions can be used
          
    alternatively the option findBestFit can be used in order to loop over all possible fit functions, and set the function with the smallest chi2 value (for similiar chi2 the function with fewer parameters is favored) as fit function
    For FitFunction=findBestFit function 11,6 and 4 are fit and the option with relative chi2 closest to 1 kept, except if the relative chi2 value is larger than 3 in that case all functions are fit and the best fit kept
    
    
### evaluation of sf fit range:
    The fitted sf uses the fit function for evaluation only up to the end of the second to last pT bin of the sf histogram. This ensures that high pt scale factors are applied correctly. To change this default range to something else "SetRangeUser<1<1=120" can be specified in the config (binwise as for the local fit flag. The argument is the value of electron pT at which the evaluation of sf using the fit function stops. All sf for pT larger than this values, are set to the sf value at this cutoff pT.
