#include "Eventplane.h"

Eventplane::Eventplane()
{
}

std::vector<TComplex> Eventplane::GetQvecA(TString det, TString sInFile)
{
    TFile *fIn = TFile::Open(sInFile);
    TTree *digitTree = (TTree*)fIn->Get("o2sim");

    std::vector<TComplex> QvecContainer;
    if (det == "FV0") { 
        CalculateQvecFV0(digitTree, QvecContainer);
    } else if (det == "FT0A" || det == "FT0C") {
        CalculateQvecFT0(digitTree, QvecContainer, det);
    } else {
        std::cout << "Eventplane::GetQvecA : No detector specified!" << std::endl;
        return QvecContainer;
    }
  
    fIn->Close();

    return QvecContainer;
}

std::vector<std::vector<TComplex>> Eventplane::GetQvecBC(TString sInFile)
{
    TFile *fIn = TFile::Open(sInFile);
    TTree* kineTree = (TTree*)fIn->Get("o2sim");

    std::vector<o2::MCTrack>* mctrack = nullptr;
    auto mcbr = kineTree->GetBranch("MCTrack");
    mcbr->SetAddress(&mctrack);

    std::vector<std::vector<TComplex>> QvecCont;
    TComplex QvecB(0), QvecC(0);

    UInt_t nEntries = kineTree->GetEntries();
    for (UInt_t ient = 0; ient < nEntries; ient++) {

        mcbr->GetEntry(ient);

        int nTracksB = 0, nTracksC = 0;
        for (auto &t : *mctrack) {

            if (t.GetPt() < 0.2) continue;

            Int_t pid = t.GetPdgCode();
            if (TMath::Abs(pid)>1000000000) continue;
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
        QvecCont[ient].push_back(QvecB);
        QvecCont[ient].push_back(QvecC);

        QvecB = TComplex(0, 0);
        QvecC = TComplex(0, 0);
    }

    fIn->Close();

    return QvecCont;
}

void Eventplane::CalculateQvecFV0(TTree *digitTree, std::vector<TComplex> &QvecContainer)
{
    std::vector<o2::fv0::BCData> bcdata, *bcdataPtr = &bcdata;
    std::vector<o2::fv0::ChannelData> chdata, *chdataPtr = &chdata;
    o2::dataformats::MCTruthContainer<o2::fv0::MCLabel> labels, *labelsPtr = &labels;

    digitTree->SetBranchAddress("FV0DigitBC", &bcdataPtr);
    digitTree->SetBranchAddress("FV0DigitCh", &chdataPtr);
    digitTree->SetBranchAddress("FV0DigitLabels", &labelsPtr);

    TComplex Qvec(0);
    int nEntries = digitTree->GetEntries();
    for (int ient = 0; ient < nEntries; ient++) {
        
        digitTree->GetEntry(ient);

        double signalSum = 0;
        int prevLabelEvent = -1;
        int nbc = bcdata.size();
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            const auto lb = labels.getLabels(ibc+1);
            int labelEvent = 0;
            if (lb.size() > 0) {
                labelEvent = lb[0].getEventID();
            } else {
                //std::cout << "Problem with labels" << std::endl;
                continue;
            }
           
	    	if (prevLabelEvent >= labelEvent) continue;

            prevLabelEvent = labelEvent;

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

            Qvec = TComplex(0, 0);
            signalSum = 0;
            
        }
    }
}

void Eventplane::CalculateQvecFT0(TTree *digitTree, std::vector<TComplex> &QvecContainer, TString det)
{
    std::vector<o2::ft0::Digit> bcdata, *bcdataPtr = &bcdata;
    std::vector<o2::ft0::ChannelData> chdata, *chdataPtr = &chdata;
    o2::dataformats::MCTruthContainer<o2::ft0::MCLabel> labels, *labelsPtr = &labels;

    digitTree->SetBranchAddress("FT0DIGITSBC", &bcdataPtr);
    digitTree->SetBranchAddress("FT0DIGITSCH", &chdataPtr);
    digitTree->SetBranchAddress("FT0DIGITSMCTR", &labelsPtr);

    TComplex Qvec(0);
    int nEntries = digitTree->GetEntries();
    for (int ient = 0; ient < nEntries; ient++) {
        
        digitTree->GetEntry(ient);

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
            prevLabelEvent = labelEvent;
            
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

            Qvec = TComplex(0, 0);
            signalSum = 0;
        }
    }
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
