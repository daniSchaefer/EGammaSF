# EGammaSF


provides simple class to help analysist extract the scale-factors provided by EGamma POG.


## While in debugging stage:
    use makefile and test.cpp to get minimal working example
    
    
## Class properties:
     needs config file in which the filename for the corresponding root file, and the names of the histograms are specified 
     construct with ScaleFactorHelper( EGammaInput, bool )
        - first argument determines which scale factors are used. possibilities are:
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
     
     
     
## structure of config-file:
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
    
