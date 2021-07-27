#include "Eventplane.h"

Eventplane::Eventplane()
{
    fBmin = 0.;
    fBmax = 20.;
}

Eventplane::Eventplane(double bmin, double bmax)
{
    fBmin = bmin;
    fBmax = bmax;
}

int Eventplane::OpenFiles(TString nameKineFile, TString nameFV0DigitFile, TString nameFT0DigitFile)
{
    fInKine = TFile::Open(nameKineFile, "READ");
    fInFV0Digit = TFile::Open(nameFV0DigitFile, "READ");
    fInFT0Digit = TFile::Open(nameFT0DigitFile, "READ");

    if (!fInKine || !fInFV0Digit || !fInFT0Digit) {
        std::cout << "Could not open the files, skipping them.." << std::endl;
        return 0;
    }

    fKineTree = (TTree*)fInKine->Get("o2sim");
    fFV0DigitTree = (TTree*)fInFV0Digit->Get("o2sim");
    fFT0DigitTree = (TTree*)fInFT0Digit->Get("o2sim");

    return 1;
}

void Eventplane::CloseFiles()
{
    fInKine->Close();
    fInFV0Digit->Close();
    fInFT0Digit->Close();
}

std::vector<TComplex> Eventplane::GetQvecA(TString det)
{
    std::vector<TComplex> QvecContainer;
    std::cout << "\nQvec A : " << std::endl;
    if (det == "FV0") { 
        CalculateQvecFV0(QvecContainer);
    } else if (det == "FT0A" || det == "FT0C") {
        CalculateQvecFT0(QvecContainer, det);
    } else {
        std::cout << "Eventplane::GetQvecA : No detector specified!" << std::endl;
        return QvecContainer;
    }
  
    return QvecContainer;
}

std::vector<std::vector<TComplex>> Eventplane::GetQvecBC()
{
    std::vector<std::vector<TComplex>> QvecCont;
    TComplex QvecB(0), QvecC(0);

    std::vector<o2::MCTrack>* mctrack = nullptr;
    auto mcbr = fKineTree->GetBranch("MCTrack");
    mcbr->SetAddress(&mctrack);

    o2::dataformats::MCEventHeader* mcheader = nullptr;
    auto hdrbr = fKineTree->GetBranch("MCEventHeader.");
    hdrbr->SetAddress(&mcheader);

    int iqvec = 0;
    UInt_t nEntries = fKineTree->GetEntries();
    std::cout << "\nQvec B and C : " << std::endl;
    for (UInt_t ient = 0; ient < nEntries; ient++) {

        fKineTree->GetEntry(ient);
        
        double b = mcheader->GetB();
        std::cout << "\tEvent " << ient << "  b : " << b << std::endl;
        if (b < fBmin || b > fBmax) continue;
        int nTracksB = 0, nTracksC = 0;
        for (auto &t : *mctrack) {

            if (t.GetPt() < 0.2) continue;

            Int_t pid = t.GetPdgCode();
            if (TMath::Abs(pid)>1000000000) continue;
            if (TDatabasePDG::Instance()->GetParticle(pid)==NULL) continue;
            Double_t charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge();
            if (charge==0.0) continue;
            double len = TMath::Sqrt(t.Vx()*t.Vx() + t.Vy()*t.Vy() + t.Vz()*t.Vz());
            if (len > 0.5) continue;
            //if (!isHadron(pid)) continue;

            double eta = t.GetEta();
            double phi = TMath::ATan2(t.Py(), t.Px());

            if ( eta>-0.8 && eta<-0.1 ) {
                QvecB += TComplex(TMath::Cos(2.0 * phi), TMath::Sin(2.0 * phi));
                nTracksB++;
            }

            if ( eta>0.1 && eta<0.8 ) {
                QvecC += TComplex(TMath::Cos(2.0 * phi), TMath::Sin(2.0 * phi));
                nTracksC++;
            }
        }

        QvecB /= (double)nTracksB;
        QvecC /= (double)nTracksC;

        QvecCont.push_back(std::vector<TComplex>());
        QvecCont[iqvec].push_back(QvecB);
        QvecCont[iqvec++].push_back(QvecC);
        std::cout << "\t\t --> event kept" << std::endl;

        QvecB = TComplex(0, 0);
        QvecC = TComplex(0, 0);
    }
 
    fKineTree->ResetBranchAddresses();

    delete mctrack;
    delete mcheader;
   
    return QvecCont;
}

void Eventplane::CalculateQvecFV0(std::vector<TComplex> &QvecContainer)
{
    std::vector<o2::fv0::BCData> bcdata, *bcdataPtr = &bcdata;
    std::vector<o2::fv0::ChannelData> chdata, *chdataPtr = &chdata;
    o2::dataformats::MCTruthContainer<o2::fv0::MCLabel> labels, *labelsPtr = &labels;

    fFV0DigitTree->SetBranchAddress("FV0DigitBC", &bcdataPtr);
    fFV0DigitTree->SetBranchAddress("FV0DigitCh", &chdataPtr);
    fFV0DigitTree->SetBranchAddress("FV0DigitLabels", &labelsPtr);
 
    o2::dataformats::MCEventHeader* mcheader = nullptr;
    auto hdrbr = fKineTree->GetBranch("MCEventHeader.");
    hdrbr->SetAddress(&mcheader);
   
    TComplex Qvec(0);
    int nEntries = fFV0DigitTree->GetEntries();
    for (int ient = 0; ient < nEntries; ient++) {
        
        fFV0DigitTree->GetEntry(ient);

        double signalSum = 0;
        int prevLabelEvent = -1;
        int nbc = bcdata.size();
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            const auto lb = labels.getLabels(ibc);
            int labelEvent = 0;
            if (lb.size() > 0) {
                labelEvent = lb[0].getEventID();
            } else {
                //std::cout << "Problem with labels" << std::endl;
                continue;
            }
           
	    	if (prevLabelEvent >= labelEvent) continue;

            // Add empty vector if there is a gap in event numbers
            if (labelEvent - prevLabelEvent > 1) {
                int ievmiss = labelEvent - 1;
                fKineTree->GetEntry(ievmiss);
                double b = mcheader->GetB();
                std::cout << "\tEvent " << ievmiss << "  b : " << b << std::endl;
                if (fBmin < b && fBmax > b) {
                    std::cout << "Event " << ievmiss << " missing, add empty vector" << std::endl;
                    QvecContainer.push_back(TComplex(0, 0));
                }
            }

            prevLabelEvent = labelEvent;
            
            fKineTree->GetEntry(labelEvent);

            double b = mcheader->GetB();
            std::cout << "\tEvent " << labelEvent << "  b : " << b << std::endl;
            if (fBmin > b || fBmax < b) continue;

            double charge;
            int channel;
            const auto &bcd = bcdata[ibc];
            int chEnt = bcd.ref.getFirstEntry();
            for (int ich = 0; ich < bcd.ref.getEntries(); ich++) { // Go through FV0 channels
                const auto &chd = chdata[chEnt++];
                charge = chd.chargeAdc;
                channel = chd.pmtNumber;
                SumQvec(Qvec, charge, channel, "FV0");
                signalSum += charge;
            }
            
            if (signalSum != 0) {
                Qvec /= signalSum;
                QvecContainer.push_back(TComplex(Qvec));
            } else {
                QvecContainer.push_back(TComplex(0, 0));
            }
            std::cout << "\t\t --> event kept" << std::endl;

            Qvec = TComplex(0, 0);
            signalSum = 0;
            
        }
    }

    fKineTree->ResetBranchAddresses();
    fFV0DigitTree->ResetBranchAddresses();

    delete bcdataPtr;
    delete chdataPtr;
    delete labelsPtr;
    delete mcheader;
}

void Eventplane::CalculateQvecFT0(std::vector<TComplex> &QvecContainer, TString det)
{
    std::vector<o2::ft0::Digit> bcdata, *bcdataPtr = &bcdata;
    std::vector<o2::ft0::ChannelData> chdata, *chdataPtr = &chdata;
    o2::dataformats::MCTruthContainer<o2::ft0::MCLabel> labels, *labelsPtr = &labels;

    fFT0DigitTree->SetBranchAddress("FT0DIGITSBC", &bcdataPtr);
    fFT0DigitTree->SetBranchAddress("FT0DIGITSCH", &chdataPtr);
    fFT0DigitTree->SetBranchAddress("FT0DIGITSMCTR", &labelsPtr);
 
    o2::dataformats::MCEventHeader* mcheader = nullptr;
    auto hdrbr = fKineTree->GetBranch("MCEventHeader.");
    hdrbr->SetAddress(&mcheader);
   
    TComplex Qvec(0);
    int nEntries = fFT0DigitTree->GetEntries();
    for (int ient = 0; ient < nEntries; ient++) {
        
        fFT0DigitTree->GetEntry(ient);

        double signalSum = 0;
        int prevLabelEvent = -1;
        int nbc = bcdata.size();
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            const auto lb = labels.getLabels(ibc);
            int labelEvent = 0;
            if (lb.size() > 0) {
                labelEvent = lb[0].getEventID();
            } else {
                std::cout << "Problem with labels" << std::endl;
                continue;
            }
            
            if (labelEvent == prevLabelEvent) continue;

            // Add empty vector if there is a gap in event numbers
            if (labelEvent - prevLabelEvent > 1) {
                int ievmiss = labelEvent - 1;
                fKineTree->GetEntry(ievmiss);
                double b = mcheader->GetB();
                std::cout << "\tEvent " << ievmiss << "  b : " << b << std::endl;
                if (fBmin < b && fBmax > b) {
                    std::cout << "Event " << ievmiss << " missing, add empty vector" << std::endl;
                    QvecContainer.push_back(TComplex(0, 0));
                }
            }

            prevLabelEvent = labelEvent;
            
            fKineTree->GetEntry(labelEvent);
            
            double b = mcheader->GetB();
            std::cout << "\tEvent " << labelEvent << "  b : " << b << std::endl;
            if (fBmin > b || fBmax < b) continue;

            double charge;
            int channel;
            const auto &bcd = bcdata[ibc];
            int chEnt = bcd.ref.getFirstEntry();
            for (int ich = 0; ich < bcd.ref.getEntries(); ich++) { // Go through FT0 channels
                const auto &chd = chdata[chEnt++];
                charge = chd.QTCAmpl;
                channel = chd.ChId;
                if (det == "FT0A" && channel < 96) {
                    SumQvec(Qvec, charge, channel, det);
                    signalSum += charge;
                }
                if (det == "FT0C" && channel > 95) {
                    SumQvec(Qvec, charge, channel, det);
                    signalSum += charge;
                }
            }
 
            if (signalSum != 0) {
                Qvec /= signalSum;
                QvecContainer.push_back(TComplex(Qvec));
            } else {
                QvecContainer.push_back(TComplex(0, 0));
            }
            std::cout << "\t\t --> event kept" << std::endl;

            Qvec = TComplex(0, 0);
            signalSum = 0;
        }
    }

    fKineTree->ResetBranchAddresses();
    fFT0DigitTree->ResetBranchAddresses();

    delete bcdataPtr;
    delete chdataPtr;
    delete labelsPtr;
    delete mcheader;
}

void Eventplane::SumQvec(TComplex &Qvec, double nch, int chno, TString det)
{
    double phi = 0.0;
    if (det == "FV0") {
        phi = GetFV0Phi(chno);
    } else if (det == "FT0A") {
        phi = GetFT0APhi(chno);
    } else if (det == "FT0C") {
        phi = GetFT0CPhi(chno-96);
    } else {
        std::cout << "Eventplane::SumQvec : warning, phi not taken from any detector" << std::endl;
    }
    Qvec += TComplex(nch*TMath::Cos(2.0 * phi), nch*TMath::Sin(2.0 * phi));
}

double Eventplane::GetFV0Phi(int chno)
{
    const double pi = TMath::Pi();
    int pmtMapSmallCh[8] = {5,4,3,2,6,7,0,1};
    std::map<int, int> pmtMapBigCh = {
        {46, 0}, {38, 1}, {47, 2}, {39, 3}, {43, 4}, {35, 5}, {42, 6}, {34, 7},
        {41, 8}, {33, 9}, {40, 10}, {32, 11}, {44, 12}, {36, 13}, {45, 14}, {37, 15}
    };

    if (chno>31) { // 5th ring has 16 segments
        return pi/16.0 + (pmtMapBigCh[chno])*pi/8.0 - pi;
    } else {
        return pi/8.0 + (pmtMapSmallCh[chno%8])*pi/4.0 - pi;
    }
}

double Eventplane::GetFT0APhi(int chno)
{
    return TMath::ATan2(ft0ay[chno], ft0ax[chno]);
}

double Eventplane::GetFT0CPhi(int chno)
{
    return TMath::ATan2(ft0cy[chno], ft0cx[chno]);
}

double Eventplane::GetEventPlane(TComplex Qvec)
{
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/2.0;
}
