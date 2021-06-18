void CompareElossToCharge(TString sDirName)
{
    TString hitfile(Form("%s/o2sim_HitsFV0.root", sDirName.Data()));
    TString digitfile(Form("%s/fv0digits.root", sDirName.Data()));

    // Save eloss and charge for each event
    std::vector<std::vector<double>> vecEloss;
    std::vector<std::vector<double>> vecCharge;

    // Go through hits
    TFile *finhits = TFile::Open(hitfile);
    TTree *hitsTree = (TTree*)finhits->Get("o2sim");

    std::vector<o2::fv0::Hit> fv0hits, *fv0hitsPtr = &fv0hits;
    hitsTree->SetBranchAddress("FV0Hit", &fv0hitsPtr);

    int nent = hitsTree->GetEntries();
    for (int ient = 0; ient < nent; ient++) {
        
        std::vector<double> vecElossCh(40);
        
        hitsTree->GetEntry(ient);

        int nhits = fv0hits.size();
        for (int ihit = 0; ihit < nhits; ihit++) {
            o2::fv0::Hit* hit = &(fv0hits.at(ihit));

            int idet = hit->GetDetectorID();
            double eloss = hit->GetEnergyLoss();
            vecElossCh[idet] += eloss;
        }

        vecEloss.push_back(vecElossCh);
        vecElossCh.clear();
    }

    finhits->Close();

    // Go through digits
    TFile *findigits = TFile::Open(digitfile);
    TTree *digitTree = (TTree*)findigits->Get("o2sim");

    std::vector<o2::fv0::BCData> bcdata, *bcdataPtr = &bcdata;
    std::vector<o2::fv0::ChannelData> chdata, *chdataPtr = &chdata;
    o2::dataformats::MCTruthContainer<o2::fv0::MCLabel> labels, *labelsPtr = &labels;

    digitTree->SetBranchAddress("FV0DigitBC", &bcdataPtr);
    digitTree->SetBranchAddress("FV0DigitCh", &chdataPtr);
    digitTree->SetBranchAddress("FV0DigitLabels", &labelsPtr);
    
    nent = digitTree->GetEntries();
    for (int ient = 0; ient < nent; ient++) {
        
        std::vector<double> vecChargeCh(40);

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

            const auto &bcd = bcdata[ibc];
            int chent = bcd.ref.getFirstEntry();
            for (int ich = 0; ich < bcd.ref.getEntries(); ich++) {
                const auto &chd = chdata[chent++];
                double charge = chd.chargeAdc;
                int ch = chd.pmtNumber;
                if (ch < 40) {
                    vecChargeCh[ch] += charge;
                } else {
                    vecChargeCh[ch-8] += charge;
                }
            }
            vecCharge.push_back(vecChargeCh);
            vecChargeCh.clear();
        }
    }

    findigits->Close();

    TFile *fout;
    TH2D *hHitPerCharge;

    bool bOutputExists = gSystem->AccessPathName("hitpercharge.root");
    if (bOutputExists) {
        fout = TFile::Open("hitpercharge.root", "NEW");
        hHitPerCharge = new TH2D("hHitPerCharge", "hHitPerCharge", 40, 0.5, 40.5, 301, 0., 0.002);
    } else {
        fout = TFile::Open("hitpercharge.root", "UPDATE");
        hHitPerCharge = (TH2D*)fout->Get("hHitPerCharge");
    } 

    // Fill hits per charge histogram for each channel
    for (int i = 0; i < (int)vecCharge.size(); i++) {
        for (int j = 0; j < (int)vecCharge[i].size(); j++) {
            double div = vecEloss[i][j]/vecCharge[i][j];
            hHitPerCharge->Fill(j, div);
        }
    }

    //TCanvas *c1 = new TCanvas("c1", "c1");
    //hHitPerCharge->Draw("COLZ");

	//c1->SaveAs("hitpercharge.png");

    fout->Write("", TObject::kOverwrite);
    fout->Close();
}
