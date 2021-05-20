#include "Histos.h"
#include "TMath.h"

Histos::Histos() {

    hQvecAfv0 = new TH2D("hQvecAfv0", "hQvecAfv0", 401, -1.0, 1.0, 401, -1.0, 1.0);
    hQvecAft0c = new TH2D("hQvecAft0c", "hQvecAft0c", 401, -1.0, 1.0, 401, -1.0, 1.0);
    hQvecAft0a = new TH2D("hQvecAft0a", "hQvecAft0a", 401, -1.0, 1.0, 401, -1.0, 1.0);
    hQvecB = new TH2D("hQvecB", "hQvecB", 401, -1.0, 1.0, 401, -1.0, 1.0);
    hQvecC = new TH2D("hQvecC", "hQvecC", 401, -1.0, 1.0, 401, -1.0, 1.0);
    
    hRabFV0 = new TH1D("hRabFV0", "hRabFV0", 401, -1.0, 1.0);
    hRacFV0 = new TH1D("hRacFV0", "hRacFV0", 401, -1.0, 1.0);
    hRbcFV0 = new TH1D("hRbcFV0", "hRbcFV0", 401, -1.0, 1.0);
    
    hRabFT0A = new TH1D("hRabFT0A", "hRabFT0A", 401, -1.0, 1.0);
    hRacFT0A = new TH1D("hRacFT0A", "hRacFT0A", 401, -1.0, 1.0);
    hRbcFT0A = new TH1D("hRbcFT0A", "hRbcFT0A", 401, -1.0, 1.0);
    
    hRabFT0C = new TH1D("hRabFT0C", "hRabFT0C", 401, -1.0, 1.0);
    hRacFT0C = new TH1D("hRacFT0C", "hRacFT0C", 401, -1.0, 1.0);
    hRbcFT0C = new TH1D("hRbcFT0C", "hRbcFT0C", 401, -1.0, 1.0);

    hEPAfv0 = new TH1D("hEPAfv0", "hEPAfv0", 201, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hEPAft0a = new TH1D("hEPAft0a", "hEPAft0a", 201, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hEPAft0c = new TH1D("hEPAft0c", "hEPAft0c", 201, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hEPB = new TH1D("hEPB", "hEPB", 201, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hEPC = new TH1D("hEPC", "hEPC", 201, -TMath::Pi()/2.0, TMath::Pi()/2.0);
}
