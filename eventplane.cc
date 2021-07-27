#include <iostream>
#include <vector>
#include <stdlib.h>

#include "src/Histos.h"
#include "src/DataManager.h"
#include "src/Eventplane.h"
#include "src/Corrections.h"

#include "TH2D.h"

#include "TString.h"

int main(int argc, char **argv) {

    TString fOutName = argc > 1 ? argv[1] : "eventplane.root";
    TString sKineFile = argc > 2 ? argv[2] : "";
    TString sFV0File = argc > 3 ? argv[3] : "";
    TString sFT0File = argc > 4 ? argv[4] : "";
    double bmin = argc > 5 ? std::stof(argv[5]) : 0.;
    double bmax = argc > 6 ? std::stof(argv[6]) : 20.;
    bool bDoCorrections = argc > 7 ? std::stoi(argv[7]) : 0;

    TFile *fOut;

    Histos *histos = new Histos(fOutName, fOut); // This opents the output file
    DataManager *dm = new DataManager();
    Eventplane *ep = new Eventplane(bmin, bmax);
    Corrections *corrFV0 = new Corrections("corrections_fv0.txt");
    Corrections *corrFT0A = new Corrections("corrections_ft0a.txt");
    Corrections *corrFT0C = new Corrections("corrections_ft0c.txt");
    
    if (bDoCorrections) {
        corrFV0->LoadCorrections("corrections_fv0.txt");
        corrFV0->Print();
        corrFT0A->LoadCorrections("corrections_ft0a.txt");
        corrFT0A->Print();
        corrFT0C->LoadCorrections("corrections_ft0c.txt");
        corrFT0C->Print();
    }
   
    std::cout << "Calculating event plane to events between b = [ " << bmin << " " << bmax << " ]" << std::endl;
    
    int filesOpen = ep->OpenFiles(sKineFile, sFV0File, sFT0File);
    if (!filesOpen) continue;

    std::vector<TComplex> QvecAfv0 = ep->GetQvecA("FV0");
    std::vector<TComplex> QvecAft0a = ep->GetQvecA("FT0A");
    std::vector<TComplex> QvecAft0c = ep->GetQvecA("FT0C");
    std::vector<std::vector<TComplex>> QvecBC = ep->GetQvecBC();

    std::cout << "Q-vec sizes (FV0 FT0A FT0C TPC)" << std::endl;
    std::cout << QvecAfv0.size() << "  " << QvecAft0a.size() << "  " << QvecAft0c.size() << "  " << QvecBC.size()<< "  " << std::endl;

    int nev = QvecAfv0.size();
    for (int iev = 0; iev < nev; iev++) {
        double qx = QvecAfv0[iev].Re(); double qy = QvecAfv0[iev].Im();
        
        if (bDoCorrections) corrFV0->DoCorrections(qx, qy);
        
        double epA = ep->GetEventPlane(TComplex(qx, qy));
        double epB = ep->GetEventPlane(QvecBC[iev][0]);
        double epC = ep->GetEventPlane(QvecBC[iev][1]);

        histos->hQvecAfv0->Fill(qx, qy);
        histos->hEPAfv0->Fill(epA);
        histos->hEPB->Fill(epB);
        histos->hEPC->Fill(epC);
        histos->hRabFV0->Fill(TMath::Cos(2*(epA - epB)));
        histos->hRacFV0->Fill(TMath::Cos(2*(epA - epC)));
        histos->hRbcFV0->Fill(TMath::Cos(2*(epB - epC)));
    }
	
    int nevft0a = QvecAft0a.size();
    for (int iev = 0; iev < nevft0a; iev++) {
        double qx = QvecAft0a[iev].Re(); double qy = QvecAft0a[iev].Im();
        
        if (bDoCorrections) corrFT0A->DoCorrections(qx, qy);
        
        double epA = ep->GetEventPlane(TComplex(qx, qy));
        double epB = ep->GetEventPlane(QvecBC[iev][0]);
        double epC = ep->GetEventPlane(QvecBC[iev][1]);
        
        histos->hQvecAft0a->Fill(qx, qy);
        histos->hEPAft0a->Fill(epA);
        histos->hRabFT0A->Fill(TMath::Cos(2*(epA - epB)));
        histos->hRacFT0A->Fill(TMath::Cos(2*(epA - epC)));
        histos->hRbcFT0A->Fill(TMath::Cos(2*(epB - epC)));
    }
        
    int nevft0c = QvecAft0c.size();
    for (int iev = 0; iev < nevft0c; iev++) {
        double qx = QvecAft0c[iev].Re(); double qy = QvecAft0c[iev].Im();
        
        if (bDoCorrections) corrFT0C->DoCorrections(qx, qy);
        
        double epA = ep->GetEventPlane(TComplex(qx, qy));
        double epB = ep->GetEventPlane(QvecBC[iev][0]);
        double epC = ep->GetEventPlane(QvecBC[iev][1]); 
        
        histos->hQvecAft0c->Fill(qx, qy);
        histos->hEPAft0c->Fill(epA);
        histos->hRabFT0C->Fill(TMath::Cos(2*(epA - epB)));
        histos->hRacFT0C->Fill(TMath::Cos(2*(epA - epC)));
        histos->hRbcFT0C->Fill(TMath::Cos(2*(epB - epC)));
    }

    ep->CloseFiles();

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    return 0;
}
