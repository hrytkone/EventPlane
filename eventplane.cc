#include <iostream>
#include <vector>
#include <stdlib.h>

#include "src/Histos.h"
#include "src/DataManager.h"
#include "src/Eventplane.h"
#include "src/Corrections.h"

#include "TStopwatch.h"
#include "TH2D.h"

int main(int argc, char **argv) {

    TString fOutName = argc > 1 ? argv[1] : "eventplane.root";
    TString sDirName = argc > 2 ? argv[2] : "";
    bool bDoCorrections = argc > 3 ? std::stoi(argv[3]) : 0;
    bool bCalculateCorrections = argc > 4 ? std::stoi(argv[4]) : 0; // if this is true only calculate corrections and save then to a file

    TStopwatch timer;
    timer.Start();

    TFile *fOut = TFile::Open(fOutName, "RECREATE");
    
    Histos *histos = new Histos();
    DataManager *dm = new DataManager();
    Corrections *corrFV0 = new Corrections("corrections_fv0.txt");
    Corrections *corrFT0A = new Corrections("corrections_ft0a.txt");
    Corrections *corrFT0C = new Corrections("corrections_ft0c.txt");
   
    TString cmdfv0(Form("ls %s*/fv0*.root", sDirName.Data()));
    std::vector<TString> fv0Files = dm->GetFileNames(cmdfv0);

    TString cmdft0(Form("ls %s*/ft0*.root", sDirName.Data()));
    std::vector<TString> ft0Files = dm->GetFileNames(cmdft0);

    TString cmdmc(Form("ls %s*/output_nockov.root", sDirName.Data()));
    std::vector<TString> mcFiles = dm->GetFileNames(cmdmc);

    if (bCalculateCorrections) {
        TH2D *hQvecfv0 = new TH2D("hQvecfv0", "hQvecfv0", 251, -0.5, 0.5, 251, -0.5, 0.5);
        TH2D *hQvecft0a = new TH2D("hQvecft0a", "hQvecft0a", 251, -0.5, 0.5, 251, -0.5, 0.5);
        TH2D *hQvecft0c = new TH2D("hQvecft0c", "hQvecft0c", 251, -0.5, 0.5, 251, -0.5, 0.5);

        int nfiles = fv0Files.size();
        for (int ifile = 0; ifile < nfiles; ifile++) {
            Eventplane ep;
            std::vector<TComplex> QvecAfv0 = ep.GetQvecA("FV0", fv0Files[ifile]);
            for (int i = 0; i < (int)QvecAfv0.size(); i++) hQvecfv0->Fill(QvecAfv0[i].Re(), QvecAfv0[i].Im());            
            
            std::vector<TComplex> QvecAft0a = ep.GetQvecA("FT0A", ft0Files[ifile]);
            for (int i = 0; i < (int)QvecAft0a.size(); i++) hQvecft0a->Fill(QvecAft0a[i].Re(), QvecAft0a[i].Im());            
            
            std::vector<TComplex> QvecAft0c = ep.GetQvecA("FT0C", ft0Files[ifile]);
            for (int i = 0; i < (int)QvecAft0c.size(); i++) hQvecft0c->Fill(QvecAft0c[i].Re(), QvecAft0c[i].Im());            
        }

        corrFV0->SaveCorrections(hQvecfv0);
        corrFT0A->SaveCorrections(hQvecft0a);
        corrFT0C->SaveCorrections(hQvecft0c);

        return 0;
    }

    if (bDoCorrections) {
        corrFV0->LoadCorrections("corrections_fv0.txt");
        corrFV0->Print();
        corrFT0A->LoadCorrections("corrections_ft0a.txt");
        corrFT0A->Print();
        corrFT0C->LoadCorrections("corrections_ft0c.txt");
        corrFT0C->Print();
    }
   
    // loop over files
    int nfiles = fv0Files.size();
    //for (int ifile = 0; ifile < nfiles-3; ifile++) {
    for (int ifile = 0; ifile < 1; ifile++) {
        Eventplane ep;
        std::vector<TComplex> QvecAfv0 = ep.GetQvecA("FV0", fv0Files[ifile]);
        std::vector<TComplex> QvecAft0a = ep.GetQvecA("FT0A", ft0Files[ifile]);
        std::vector<TComplex> QvecAft0c = ep.GetQvecA("FT0C", ft0Files[ifile]);
        std::vector<std::vector<TComplex>> QvecBC = ep.GetQvecBC(mcFiles[ifile]);

        int nev = QvecAfv0.size();
        for (int iev = 0; iev < nev; iev++) {
            double qx = QvecAfv0[iev].Re(); double qy = QvecAfv0[iev].Im();
            
            if (bDoCorrections) corrFV0->DoCorrections(qx, qy);
            
            double epA = ep.GetEventPlane(TComplex(qx, qy));
            double epB = ep.GetEventPlane(QvecBC[iev][0]);
            double epC = ep.GetEventPlane(QvecBC[iev][1]);

            histos->hQvecAfv0->Fill(qx, qy);
            histos->hEPAfv0->Fill(epA);
            histos->hRabFV0->Fill(TMath::Cos(2*(epA - epB)));
            histos->hRacFV0->Fill(TMath::Cos(2*(epA - epC)));
            histos->hRbcFV0->Fill(TMath::Cos(2*(epB - epC)));
        }
	
	    int nevft0a = QvecAft0a.size();
        for (int iev = 0; iev < nevft0a; iev++) {
            double qx = QvecAft0a[iev].Re(); double qy = QvecAft0a[iev].Im();
            
            if (bDoCorrections) corrFT0A->DoCorrections(qx, qy);
            
            double epA = ep.GetEventPlane(TComplex(qx, qy));
            double epB = ep.GetEventPlane(QvecBC[iev][0]);
            double epC = ep.GetEventPlane(QvecBC[iev][1]);
            
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
            
            double epA = ep.GetEventPlane(TComplex(qx, qy));
            double epB = ep.GetEventPlane(QvecBC[iev][0]);
            double epC = ep.GetEventPlane(QvecBC[iev][1]); 
            
            histos->hQvecAft0c->Fill(qx, qy);
            histos->hEPAft0c->Fill(epA);
            histos->hRabFT0C->Fill(TMath::Cos(2*(epA - epB)));
            histos->hRacFT0C->Fill(TMath::Cos(2*(epA - epC)));
            histos->hRbcFT0C->Fill(TMath::Cos(2*(epB - epC)));
        }
    }
    
    fOut->Write("", TObject::kOverwrite);
    fOut->Close();
    timer.Print();

    return 0;
}
