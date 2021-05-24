double GetFV0Phi(int chno);

void CompareToReactionPlane(TString sDirName)   
{
    TString digitfile(Form("%s/fv0digits.root", sDirName.Data()));
    TString mcfile(Form("%s/o2sim_Kine.root", sDirName.Data()));

    TFile *finkine = TFile::Open(mcfile);
    TTree *kineTree = (TTree*)finkine->Get("o2sim");
    
    kineTree->SetBranchStatus("*", 0);
    kineTree->SetBranchStatus("MCEventHeader*", 1);
    
    TFile *findigit = TFile::Open(digitfile);
    TTree *digitTree = (TTree*)findigit->Get("o2sim");
    
    /**std::vector<o2::MCTrack>* mctrack = nullptr;
    auto mcbr = kineTree->GetBranch("MCTrack");
    mcbr->SetAddress(&mctrack);**/

    o2::dataformats::MCEventHeader* mcheader = nullptr;
    auto hdrbr = kineTree->GetBranch("MCEventHeader.");
    hdrbr->SetAddress(&mcheader);

    std::vector<o2::fv0::BCData> bcdata, *bcdataPtr = &bcdata;
    std::vector<o2::fv0::ChannelData> chdata, *chdataPtr = &chdata;
    o2::dataformats::MCTruthContainer<o2::fv0::MCLabel> labels, *labelsPtr = &labels;

    digitTree->SetBranchAddress("FV0DigitBC", &bcdataPtr);
    digitTree->SetBranchAddress("FV0DigitCh", &chdataPtr);
    digitTree->SetBranchAddress("FV0DigitLabels", &labelsPtr);

    std::vector<double> rp, ep;

    int nent = kineTree->GetEntries();
    for (int ient = 0; ient < nent; ient++) {
        //std::cout << "mc event " << ient+1 << "/" << nent << std::endl;
        kineTree->GetEntry(ient);
        double reactionplane = mcheader->GetRotZ();
        if (reactionplane > TMath::Pi()) reactionplane -= 2*TMath::Pi();
        rp.push_back(reactionplane);
    }

    nent = digitTree->GetEntries();
    for (int ient = 0; ient < nent; ient++) {
        
        digitTree->GetEntry(ient);

        double signalsum = 0.;
        int prevLabelEvent = -1;
        int nbc = bcdata.size();
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            const auto lb = labels.getLabels(ibc+1);
            int labelEvent = 0;
            if (lb.size() > 0) {
                labelEvent = lb[0].getEventID();
            } else {
                continue;
            }

            if (prevLabelEvent >= labelEvent) continue;

            prevLabelEvent = labelEvent;

            double charge;
            int channel;
            const auto &bcd = bcdata[ibc];
            int chent = bcd.ref.getFirstEntry();
            double qx = 0., qy = 0.;
            for (int ich = 0; ich < bcd.ref.getEntries(); ich++) {
                const auto &chd = chdata[chent++];
                charge = chd.chargeAdc;
                channel = chd.pmtNumber;
                double phi = GetFV0Phi(channel);
                qx += charge*TMath::Cos(2.0 * phi); 
                qy += charge*TMath::Sin(2.0 * phi);               
            }
            ep.push_back(TMath::ATan2(qy, qx)/2.0);
        }
    }

    finkine->Close();
    findigit->Close();

    TFile *fout = TFile::Open("output.root", "UPDATE");

    TH1D *hDiff = new TH1D("hDiff", "hDiff", 100, 0., 2.*TMath::Pi());
    
    std::cout << "RP container size : " << rp.size() << std::endl;
    std::cout << "EP container size : " << ep.size() << std::endl;
    for (int i = 0; i < (int)rp.size(); i++) {
        //std::cout << "RP : " << rp[i] << "   EP : " << ep[i] << "   |RP-EP| : " << TMath::Abs(rp[i] - ep[i]) << std::endl;
        hDiff->Fill(TMath::Abs(rp[i] - ep[i]));
    }

    //TCanvas *c1 = new TCanvas("c1", "c1");
    //hDiff->Draw("HIST");

    fout->Write("", TObject::kOverwrite);
    fout->Close();
}

double GetFV0Phi(int chno)
{
    const double pi = TMath::Pi();
    int pmtMapSmallCh[8] = {5,4,3,2,6,7,0,1};
    std::map<int, int> pmtMapBigCh = {
        {45, 0}, {37, 1}, {44, 2}, {36, 3}, {43, 4}, {35, 5}, {42, 6}, {34, 7},
        {41, 8}, {33, 9}, {40, 10}, {32, 11}, {47, 12}, {39, 13}, {46, 14}, {38, 15}
    };
         
    if (chno>31) { // 5th ring has 16 segments
        return pi/16.0 + (pmtMapBigCh[chno])*pi/8.0 - pi;
    } else {
        return pi/8.0 + (pmtMapSmallCh[chno%8])*pi/4.0 - pi;
    }
}

