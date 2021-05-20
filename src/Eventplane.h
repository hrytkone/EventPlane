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
#include <SimulationDataFormat/MCTruthContainer.h>
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

class Eventplane {

public:
    Eventplane();
    virtual ~Eventplane() {;}

    std::vector<TComplex> GetQvecA(TString det, TString sInFile);
    std::vector<std::vector<TComplex>> GetQvecBC(TString sInFile);
    
    void CalculateQvecFV0(TTree *digitTree, std::vector<TComplex> &QvecContainer);   
    void CalculateQvecFT0(TTree *digitTree, std::vector<TComplex> &QvecContainer, TString side);   

    void SumQvec(TComplex &Qvec, double nch, int chno, TString det);
    double GetFV0Phi(int chno);
    double GetFT0APhi(int chno);
    double GetFT0CPhi(int chno);
    double GetEventPlane(TComplex Qvec);
};

#endif
