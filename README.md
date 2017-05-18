# EGammaSF


provides simple class to help analysist extract the scale-factors provided by EGamma POG.


While in debugging stage:
    use makefile and test.cpp to get minimal working example
    
    
Class properties:
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
     
