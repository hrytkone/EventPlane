#ifndef HISTOS_H
#define HISTOS_H

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

class Histos {

public:

    Histos(TString oufile, TFile *fout);
    virtual ~Histos() {;}

    // Q-vectors 
    TH2D *hQvecAfv0;
    TH2D *hQvecAft0c;
    TH2D *hQvecAft0a;
    TH2D *hQvecB;
    TH2D *hQvecC;

    // Resolution components
    TH1D *hRabFV0;
    TH1D *hRacFV0;
    TH1D *hRbcFV0;
    
    TH1D *hRabFT0A;
    TH1D *hRacFT0A;
    TH1D *hRbcFT0A;
    
    TH1D *hRabFT0C;
    TH1D *hRacFT0C;
    TH1D *hRbcFT0C;

    // Event plane distributions
    TH1D *hEPAfv0;
    TH1D *hEPAft0a;
    TH1D *hEPAft0c;
    TH1D *hEPB;
    TH1D *hEPC;
};

#endif
