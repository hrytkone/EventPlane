#ifndef EVENTPLANE_H
#define EVENTPLANE_H

#include "Const.h"

#include <iostream>

#include "TComplex.h"
#include "TFile.h"
#include "TTree.h"

#include <DataFormatsFV0/BCData.h>
#include <DataFormatsFV0/ChannelData.h>
#include <DataFormatsFV0/MCLabel.h>
#include <DataFormatsFT0/Digit.h>
#include <DataFormatsFT0/ChannelData.h>
#include <DataFormatsFT0/MCLabel.h>
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

class Eventplane {

public:
    Eventplane();
    Eventplane(double bmin, double bmax);
    virtual ~Eventplane() {;}

    int OpenFiles(TString nameKineFile, TString nameFV0DigitFile, TString nameFT0DigitFile);
    void CloseFiles();

    void SetB(double bmin, double bmax) { fBmin = bmin; fBmax = bmax; }

    std::vector<TComplex> GetQvecA(TString det);
    std::vector<std::vector<TComplex>> GetQvecBC();
    
    void CalculateQvecFV0(std::vector<TComplex> &QvecContainer);   
    void CalculateQvecFT0(std::vector<TComplex> &QvecContainer, TString side);   

    void SumQvec(TComplex &Qvec, double nch, int chno, TString det);
    double GetFV0Phi(int chno);
    double GetFT0APhi(int chno);
    double GetFT0CPhi(int chno);
    double GetEventPlane(TComplex Qvec);

private:
    double fBmin;
    double fBmax;
    
    TFile *fInKine;
    TFile *fInFV0Digit;
    TFile *fInFT0Digit;

    TTree *fKineTree; 
    TTree *fFV0DigitTree; 
    TTree *fFT0DigitTree; 

    o2::dataformats::MCEventHeader *mcheader;
    TBranch *hdrbr;
};

#endif
