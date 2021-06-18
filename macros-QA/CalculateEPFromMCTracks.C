#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TComplex.h>
#include <TDatabasePDG.h>

#include <TStopwatch.h>
#include <memory>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>
#include <unordered_map>

#include "DataFormatsFV0/Hit.h"
#include "DataFormatsFV0/BCData.h"
#include "DataFormatsFV0/ChannelData.h"
#include "DataFormatsFV0/MCLabel.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsFT0/ChannelData.h"
#include "DataFormatsFT0/MCLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "FV0Simulation/DigitizationConstant.h"
//#include "FV0Simulation/MCLabel.h"
#include "FairLogger.h"

std::vector<TComplex> GetQvecAFromHits(TString sHitsFile);
std::vector<std::vector<TComplex>> GetQvecABC(TString sKineFile);

double GetFV0Phi(double phi, double eta, double startAngle = 0.);
void SumQvec(TComplex &Qvec, double phi, double eta);

double GetEventPlane(TComplex Qvec);

bool isSumZero(double i) { return i==0.0 ? 1 : 0; }
bool isQvecEmpty(TComplex i) { return (i.Re()==0.0 && i.Im()==0.0) ? 1 : 0; }

//void CalculateEPFromMCTracks(TString sOutName = "", TString sDirName = "/scratch/project_2003583/simO2_outputs")
void CalculateEPFromMCTracks(TString sOutName = "", TString sDirName = "/scratch/project_2003583/2021-06-09_string-melting")
{
    TStopwatch timer;
    timer.Start();

    TFile *fOut = TFile::Open(sOutName, "RECREATE");

    TH1D *hRab = new TH1D("hRabFV0", "hRabFV0", 401, -1.0, 1.0);
    TH1D *hRac = new TH1D("hRacFV0", "hRacFV0", 401, -1.0, 1.0);
    TH1D *hRbc = new TH1D("hRbcFV0", "hRbcFV0", 401, -1.0, 1.0);
    
    TH2D *hEPAB = new TH2D("hEPABfv0", "hEPABfv0", 101, -TMath::Pi()/2.0, TMath::Pi()/2.0, 101, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    TH2D *hEPAC = new TH2D("hEPACfv0", "hEPACfv0", 101, -TMath::Pi()/2.0, TMath::Pi()/2.0, 101, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    TH2D *hEPBC = new TH2D("hEPBCfv0", "hEPBCfv0", 101, -TMath::Pi()/2.0, TMath::Pi()/2.0, 101, -TMath::Pi()/2.0, TMath::Pi()/2.0);

    TH1D *hEPA = new TH1D("hEPAfv0", "hEPAfv0", 201, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    TH1D *hEPB = new TH1D("hEPBfv0", "hEPBfv0", 100, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    TH1D *hEPC = new TH1D("hEPCfv0", "hEPCfv0", 100, -TMath::Pi()/2.0, TMath::Pi()/2.0);

    //TString cmd(Form("ls %s/run_5p5TeV_midcent_job*/o2sim_Kine_nockov.root", sDirName.Data()));
    TString cmd(Form("ls %s/run_*/o2sim_Kine.root", sDirName.Data()));

    std::vector<TString> mctrackFiles;

    TString tok;
    Ssiz_t from = 0;
    TString output = gSystem->GetFromPipe(cmd);
    while (output.Tokenize(tok, from, "\n")) {
        mctrackFiles.push_back(tok);
    }
    
    int nFiles = mctrackFiles.size();

    std::cout << "\nGo through files" << std::endl;
    for (int iFile = 0; iFile < nFiles; iFile++) {
    //for (int iFile = 0; iFile < 40; iFile++) {

        std::cout << "\r\tProcessing job " << iFile+1 << " / " << nFiles << std::flush;
        std::vector<std::vector<TComplex>> QvecABC = GetQvecABC(mctrackFiles[iFile]);

        for (int i=0; i<QvecABC.size(); i++) {

            double xa, ya, xb, yb, xc, yc;
            xa = QvecABC[i][0].Re(); ya = QvecABC[i][0].Im();
            xb = QvecABC[i][1].Re(); yb = QvecABC[i][1].Im();
            xc = QvecABC[i][2].Re(); yc = QvecABC[i][2].Im();

            if (xa==0 && ya==0) continue;

            double epA = GetEventPlane(TComplex(xa, ya));
            double epB = GetEventPlane(TComplex(xb, yb));
            double epC = GetEventPlane(TComplex(xc, yc));

            hRab->Fill(TMath::Cos(2*(epA - epB)));
            hRac->Fill(TMath::Cos(2*(epA - epC)));
            hRbc->Fill(TMath::Cos(2*(epB - epC)));

            hEPA->Fill(epA);
            hEPB->Fill(epB);
            hEPC->Fill(epC);

            hEPAB->Fill(epA, epB);
            hEPAC->Fill(epA, epC);
            hEPBC->Fill(epB, epC);
        }
 
        QvecABC.clear();
    }

    fOut->cd();
    fOut->Write("",TObject::kOverwrite);
    fOut->Close();

    std::cout << "\n";
    timer.Print();
}

std::vector<std::vector<TComplex>> GetQvecABC(TString sKineFile)
{
    TFile *fIn = TFile::Open(sKineFile);
    TTree* kineTree = (TTree*)fIn->Get("o2sim");

    std::vector<o2::MCTrack>* mctrack = nullptr;
    auto mcbr = kineTree->GetBranch("MCTrack");
    mcbr->SetAddress(&mctrack);

    std::vector<std::vector<TComplex>> QvecCont;
    TComplex QvecA(0), QvecB(0), QvecC(0);

    UInt_t nEntries = kineTree->GetEntries();
    for (UInt_t ient = 0; ient < nEntries; ient++) {

        mcbr->GetEntry(ient);

        int nTracksA = 0, nTracksB = 0, nTracksC = 0;
        for (auto &t : *mctrack) {

            if (t.GetPt() < 0.2) continue;

            Int_t pid = t.GetPdgCode();
            if (TMath::Abs(pid)>1000000000) continue;
            Double_t charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge();
            if (charge==0.0) continue;
                        
            double len = TMath::Sqrt(t.Vx()*t.Vx() + t.Vy()*t.Vy() + t.Vz()*t.Vz());
            if (len > 0.5) continue;

            double eta = t.GetEta();
            
            //if (t.leftTrace(o2::detectors::DetID::FV0)) {
            if ( eta > 2.2 && eta < 5.0 ) {
                SumQvec(QvecA, TMath::ATan2(t.Py(), t.Px()), eta);
                nTracksA++;
            }

            if ( eta > -0.8 && eta < -0.1 ) {
                SumQvec(QvecB, TMath::ATan2(t.Py(), t.Px()), -100.0);
                nTracksB++;
            }

            if ( eta > 0.1 && eta < 0.8 ) {
                SumQvec(QvecC, TMath::ATan2(t.Py(), t.Px()), -100.0);
                nTracksC++;
            }
        }

        QvecA /= (double)nTracksA;
        QvecB /= (double)nTracksB;
        QvecC /= (double)nTracksC;

        QvecCont.push_back(std::vector<TComplex>());
        QvecCont[ient].push_back(QvecA);
        QvecCont[ient].push_back(QvecB);
        QvecCont[ient].push_back(QvecC);

        QvecA = TComplex(0, 0);
        QvecB = TComplex(0, 0);
        QvecC = TComplex(0, 0);
    }

    fIn->Close();

    return QvecCont;
}

void SumQvec(TComplex &Qvec, double phi, double eta)
{
    if (eta != -100.) phi = GetFV0Phi(phi, eta, -TMath::Pi());
    Qvec += TComplex(TMath::Cos(2.0 * phi), TMath::Sin(2.0 * phi));
}

double GetFV0Phi(double phi, double eta, double startAngle)
{   
    double lower, upper;
    int nSec = 8;
    if (eta<2.8) nSec = 16;
    for (int i=0; i<nSec; i++) {
        lower = startAngle+2*TMath::Pi()*i/(float)nSec;
        upper = startAngle+2*TMath::Pi()*(i+1)/(float)nSec;
        if (lower < phi && phi < upper) {
            return lower + (upper-lower)/2.0;
        }
    }
    return 0;
}


double GetEventPlane(TComplex Qvec)
{
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/2.0;
}
