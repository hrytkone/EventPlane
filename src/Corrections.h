#ifndef CORRECTIONS_H
#define CORRECTIONS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

const int nCorrections = 8;

class Corrections {

public:
    Corrections(TString fileName);
    virtual ~Corrections() {;}
    
    double *CalculateCorrections(TH2D *hQ);
    void GetRecenteringCorrection(TH2D *hQ, double &xmean, double &ymean);
    void GetWidthCorrection(TH2D *hQ, double &xdev, double &ydev);
    void GetTwistAndRescaleCorrection(TH2D *hQ, double &aplus, double &aminus, double &lambdaplus, double &lambdaminus);

    void SaveCorrections(TH2D *hQvec); // Save corrections to the text file
    void LoadCorrections(TString filename); // Load corrections from the text file   
    void DoCorrections(double &qx, double &qy);
    void Print();

private:
    double corrections[nCorrections];
    TString saveFileName;
    std::ofstream outputFile;
    std::ifstream inputFile;

};

#endif
