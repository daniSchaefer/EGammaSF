// compile with : g++ test.cpp EGammaSF.cc `root-config --cflags --glibs --libs --evelibs`
           
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2D.h"
#include "TF1.h"
#include "EGammaSF.h"
//#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd5/include/TTreeReaderValueBase.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>
#include <exception>



int main(int argc, char** argv)
{
     try{
        // initialize scale factor helper : 
        ScaleFactorHelper* bla = new ScaleFactorHelper(EGammaInput::electronTight,1);
         
        for(int i=0;i<10;i+=1)
        {
        float unc = bla->GetUncertainty(30+i*10,2.1);
        float sf = bla->GetSF(30.0+i*10,2.1);
        std::cout << "sf : " << sf << " +- " << unc << std::endl;
        }
      }
      catch(const std::exception& ex) {
        std::cout << "EXCEPTION!!!" << std::endl;
      }

    
    
   return 0; 
    
}