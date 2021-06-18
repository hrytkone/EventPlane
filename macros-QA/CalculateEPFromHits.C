//#if !defined(__CLING__) || defined(__ROOTCLING__)

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
std::vector<std::vector<TComplex>> GetQvecBC(TString sKineFile);

double GetFV0PhiHits(int pmtNumber, double phi0);
void SumQvec(TComplex &Qvec, float nch, int chNumber, double phi0);
void SumQvec(TComplex &Qvec, double phi);

double GetEventPlane(TComplex Qvec);

bool isSumZero(double i) { return i==0.0 ? 1 : 0; }
bool isQvecEmpty(TComplex i) { return (i.Re()==0.0 && i.Im()==0.0) ? 1 : 0; }

void CalculateEPFromHits(TString sOutName = "", TString sDirName = "/scratch/project_2003583/simO2_outputs")
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

    TString cmd1(Form("ls %s/run_5p5TeV_midcent_job*/o2sim_Kine_nockov.root", sDirName.Data()));
    TString cmd2(Form("ls %s/run_5p5TeV_midcent_job*/o2sim_HitsFV0.root", sDirName.Data()));

    std::vector<TString> hitFiles, mctrackFiles;

    TString tok;
    Ssiz_t from = 0;
    TString output = gSystem->GetFromPipe(cmd1);
    while (output.Tokenize(tok, from, "\n")) {
        mctrackFiles.push_back(tok);
    }
    
    from = 0;
    output = gSystem->GetFromPipe(cmd2);
    while (output.Tokenize(tok, from, "\n")) {
        hitFiles.push_back(tok);
    }

    int nFiles = mctrackFiles.size();

    std::cout << "\nGo through files" << std::endl;
    int nEvents = 0;
    for (int iFile = 40; iFile < 80; iFile++) {

        std::cout << "\r\tProcessing job " << iFile+1 << " / " << nFiles << std::flush;

        std::vector<TComplex> QvecA = GetQvecAFromHits(hitFiles[iFile]);
        std::vector<std::vector<TComplex>> QvecBC = GetQvecBC(mctrackFiles[iFile]);

        for (int i=0; i<QvecA.size(); i++) {

            double xa, ya, xb, yb, xc, yc;
            xa = QvecA[i].Re(); ya = QvecA[i].Im();
            xb = QvecBC[i][0].Re(); yb = QvecBC[i][0].Im();
            xc = QvecBC[i][1].Re(); yc = QvecBC[i][1].Im();

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
 
        QvecA.clear();
        QvecBC.clear();
    }

    fOut->cd();
    fOut->Write("",TObject::kOverwrite);
    fOut->Close();

    std::cout << "\n";
    timer.Print();
}

std::vector<TComplex> GetQvecAFromHits(TString sHitsFile)
{
    TFile *fIn = TFile::Open(sHitsFile);
    TTree* hitsTree = (TTree*)fIn->Get("o2sim");

    std::vector<o2::fv0::Hit> fv0hits, *fv0hitsPtr = &fv0hits;
    hitsTree->SetBranchAddress("FV0Hit", &fv0hitsPtr);

    std::vector<TComplex> QvecCont;
    TComplex Qvec(0);

    UInt_t nEntries = hitsTree->GetEntries();
    for (UInt_t ient = 0; ient < nEntries; ient++) {
      hitsTree->GetEntry(ient);

      double elossSum = 0.0;

      int nhits = fv0hits.size();
      for (int ihit = 0; ihit < nhits; ihit++) {

          o2::fv0::Hit* hit = &(fv0hits.at(ihit));

          int idet = hit->GetDetectorID();
          double eloss = hit->GetEnergyLoss();
          double phi0 = TMath::ATan2(hit->GetStartY(), hit->GetStartX());

          SumQvec(Qvec, eloss, idet, phi0);
          elossSum += eloss;
      }

      Qvec /= elossSum;
      QvecCont.push_back(Qvec);
      Qvec = TComplex(0, 0);

    }

    fIn->Close();

    return QvecCont;
}

std::vector<std::vector<TComplex>> GetQvecBC(TString sKineFile)
{
    TFile *fIn = TFile::Open(sKineFile);
    TTree* kineTree = (TTree*)fIn->Get("o2sim");

    std::vector<o2::MCTrack>* mctrack = nullptr;
    auto mcbr = kineTree->GetBranch("MCTrack");
    mcbr->SetAddress(&mctrack);

    std::vector<std::vector<TComplex>> QvecCont;
    TComplex QvecB(0), QvecC(0);

    UInt_t nEntries = kineTree->GetEntries();
    for (UInt_t ient = 0; ient < nEntries; ient++) {

        mcbr->GetEntry(ient);

        int nTracksB = 0, nTracksC = 0, counter = 0;
        for (auto &t : *mctrack) {

            if (t.GetPt() < 0.2) continue;

            Int_t pid = t.GetPdgCode();
            if (TMath::Abs(pid)>1000000000) continue;
            Double_t charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge();
            if (charge==0.0) continue;
            double len = TMath::Sqrt(t.Vx()*t.Vx() + t.Vy()*t.Vy() + t.Vz()*t.Vz());
            if (len > 0.5) continue;
            //if (!isHadron(pid)) continue;

            // why phi shifted?
            //std::cout << "phi = " << t.GetPhi() << "\tp : [ " << t.Px() << " " << t.Py() << " ]" << std::endl;

            double eta = t.GetEta();

            if ( eta>-0.8 && eta<-0.1 ) {
                SumQvec(QvecB, TMath::ATan2(t.Py(), t.Px()));
                nTracksB++;
            }

            if ( eta>0.1 && eta<0.8 ) {
                SumQvec(QvecC, TMath::ATan2(t.Py(), t.Px()));
                nTracksC++;
            }
        }

        QvecB /= (double)nTracksB;
        QvecC /= (double)nTracksC;

        QvecCont.push_back(std::vector<TComplex>());
        QvecCont[ient].push_back(QvecB);
        QvecCont[ient].push_back(QvecC);

        QvecB = TComplex(0, 0);
        QvecC = TComplex(0, 0);
    }

    fIn->Close();

    return QvecCont;
}

double GetFV0PhiHits(int pmtNumber, double phi0)
{
    const double pi = TMath::Pi();
    int pmtMapSmallCh[8] = {5,4,3,2,6,7,0,1};

    double phi = pi/8.0 + (pmtMapSmallCh[pmtNumber%8])*pi/4.0 - pi;
    if (pmtNumber>31) { // 5th ring has 16 segments
        if (phi0>phi) {
            return phi + pi/16.0;
        } else {
            return phi - pi/16.0;
        }
    } else {
        return phi;
    }
    return 0;
}

void SumQvec(TComplex &Qvec, float nch, int chNumber, double phi0)
{
    double phi = GetFV0PhiHits(chNumber, phi0);
    Qvec += TComplex(nch*TMath::Cos(2.0 * phi), nch*TMath::Sin(2.0 * phi));
}

void SumQvec(TComplex &Qvec, double phi)
{
    Qvec += TComplex(TMath::Cos(2.0 * phi), TMath::Sin(2.0 * phi));
}

double GetEventPlane(TComplex Qvec)
{
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/2.0;
}

//#endif
