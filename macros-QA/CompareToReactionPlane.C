int GetBin(double val, double *bins, int nbins);
double GetFV0Phi(int chno);

void CompareToReactionPlane(TString sDirName="")   
{
    TString digitfile(Form("%sfv0digits.root", sDirName.Data()));
    TString mcfile(Form("%so2sim_Kine.root", sDirName.Data()));

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

    std::vector<double> rp, ep, amplsum;

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

            double sum = 0;
            const auto &bcd = bcdata[ibc];
            int chent = bcd.ref.getFirstEntry();
            double qx = 0., qy = 0.;
            for (int ich = 0; ich < bcd.ref.getEntries(); ich++) {
                const auto &chd = chdata[chent++];
                double charge = chd.chargeAdc;
                int channel = chd.pmtNumber;
                double phi = GetFV0Phi(channel);
				sum += charge;
                qx += charge*TMath::Cos(2.0 * phi); 
                qy += charge*TMath::Sin(2.0 * phi);               
            }
			qx /= sum;
			qy /= sum;
            ep.push_back(TMath::ATan2(qy, qx)/2.0);
            sum /= 1000.;
            amplsum.push_back(sum);
        }
    }

    finkine->Close();
    findigit->Close();

    // Save information to a root-file
    TFile *fout; 
	TH2D *hRpVsEp;
    TH2D *hDiffVsAmplSum;

    const int nPhiBin = 9;
    double pi = TMath::Pi();
    double binPhi[nPhiBin] = {-pi, -3.*pi/2., -pi/2., -pi/4., 0., pi/4., pi/2., 3.*pi/2., pi};
	TH1D *hDiff[nPhiBin-1];

	bool bFileExists = gSystem->AccessPathName("output.root");
	if (bFileExists) {
		fout = TFile::Open("output.root", "NEW");
    	hRpVsEp = new TH2D("hRpVsEp", "hRpVsEp", 301, -TMath::Pi(), TMath::Pi(), 301, -TMath::Pi(), TMath::Pi());
    	hDiffVsAmplSum = new TH2D("hDiffVsAmplSum", "hDiffVsAmplSum", 301, 0., 2.*TMath::Pi(), 301, 0., 300);
        for (int i = 0; i < nPhiBin-1; i++)
    	    hDiff[i] = new TH1D(Form("hDiff[%.2f;%.2f]", binPhi[i], binPhi[i+1]), Form("hDiff[%.2f;%.2f]", binPhi[i], binPhi[i+1]), 301, 0., 2.*TMath::Pi());
	} else {
		fout = TFile::Open("output.root", "UPDATE");
		hRpVsEp = (TH2D*)fout->Get("hRpVsEp");
		hDiffVsAmplSum = (TH2D*)fout->Get("hDiffVsAmplSum");
        for (int i = 0; i < nPhiBin-1; i++) 
		    hDiff[i] = (TH1D*)fout->Get(Form("hDiff[%.2f;%.2f]", binPhi[i], binPhi[i+1]));
	}
    
    for (int i = 0; i < (int)rp.size(); i++) {
        //std::cout << "RP : " << rp[i] << "   EP : " << ep[i] << "   |RP-EP| : " << TMath::Abs(rp[i] - ep[i]) << std::endl;
        double diff = TMath::Abs(rp[i] - ep[i]);
        int bin = GetBin(ep[i], binPhi, nPhiBin); 
        hDiff[bin]->Fill(diff);
        hRpVsEp->Fill(rp[i], ep[i]);
        hDiffVsAmplSum->Fill(diff, amplsum[i]);
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

int GetBin(double val, double *bins, int nbins)
{
    for (int i = 0; i < nbins-1; i++) {
        if (bins[i] <= val && bins[i+1] > val) return i;
    }
    return -1;
}
