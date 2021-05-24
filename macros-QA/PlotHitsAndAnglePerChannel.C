/**
 *	Macro to check if the phi angle in the code alings with the
 *	given channel.
 *
 *	Uses hits files for FV0 and FT0 from simulation directory
 *	given as argument.
 */

const int nfv0ch = 40;

double GetFV0Phi(int chno);

void PlotHitsAndAnglePerChannel(TString sDirName)
{

	TString fv0file(Form("%s/o2sim_HitsFV0.root", sDirName.Data()));

    TH2D *hHitsFV0[nfv0ch];
	TGraph *gAngle[nfv0ch];
    for (int ihist = 0; ihist < nfv0ch; ihist++) {
		if (ihist < 32) hHitsFV0[ihist] = new TH2D(Form("hHitsFV0%d", ihist), Form("Channel %d", ihist), 100, -30., 30, 100, -30., 30.);
		if (ihist >= 32) hHitsFV0[ihist] = new TH2D(Form("hHitsFV0%d", ihist), Form("Channel %d", ihist), 100, -100., 100., 100, -100., 100.);
		double phi;
		phi = GetFV0Phi(ihist);
		
		Double_t x[2], y[2];
		x[0] = TMath::Cos(phi); 
		x[1] = 100.*TMath::Cos(phi); 
		y[0] = TMath::Sin(phi); 
		y[1] = 100.*TMath::Sin(phi); 

		gAngle[ihist] = new TGraph(2, x, y);
		gAngle[ihist]->SetMarkerStyle(20);
		gAngle[ihist]->SetMarkerSize(0.5);
		gAngle[ihist]->SetMarkerColor(2);
		gAngle[ihist]->SetLineWidth(2);
    }

    TFile *fin = TFile::Open(fv0file);
    TTree *hitsTree = (TTree*)fin->Get("o2sim");

    std::vector<o2::fv0::Hit> fv0hits, *fv0hitsPtr = &fv0hits;
    hitsTree->SetBranchAddress("FV0Hit", &fv0hitsPtr);

	int nEntries = hitsTree->GetEntries();
    for (int ient = 0; ient < nEntries; ient++) {
        hitsTree->GetEntry(ient);

        int nhits = fv0hits.size();
        for (int ihit = 0; ihit < nhits; ihit++) {
            o2::fv0::Hit* hit = &(fv0hits.at(ihit));

            int idet = hit->GetDetectorID();
            double x = hit->GetStartX();
            double y = hit->GetStartY();

            hHitsFV0[idet]->Fill(x, y);
        }
    }

	fin->Close();

	TCanvas *c1 = new TCanvas("c1", "c1");
	c1->Divide(4, 5);

	// Plot hitmaps
	for (int ihist = 0; ihist < (int)nfv0ch/2.; ihist++) {
		c1->cd(ihist+1);
		hHitsFV0[ihist]->Draw("COL");
		gAngle[ihist]->Draw("SAME P");
	}

	c1->SaveAs("hitmaps1.pdf");
    
	TCanvas *c2 = new TCanvas("c2", "c2");
	c2->Divide(4, 5);

	// Plot hitmaps
	int canvas = 1;
	for (int ihist = (int)nfv0ch/2.; ihist < nfv0ch; ihist++) {
		c2->cd(canvas);
		hHitsFV0[ihist]->Draw("COL");
		gAngle[ihist]->Draw("SAME P");
		canvas++;
	}

	c2->SaveAs("hitmaps2.pdf");

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

