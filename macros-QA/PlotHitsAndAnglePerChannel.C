/**
 *	Macro to check if the phi angle in the code alings with the
 *	given channel.
 *
 *	Uses one hits file for FV0 and FT0 from simulation directory
 *	given as argument.
 */

const int nfv0ch = 48;

void PlotHitsAndAnglePerChannel(TString sDir)
{

	TString cmdfv0(Form("ls %s*/o2sim_HitsFV0.root", sDirName.Data()));
	//TString cmdft0(Form("ls %s*/o2sim_HitsFT0.root", sDirName.Data()));
	
    std::vector<TString> fv0Files = dm->GetFileNames(cmdfv0);
    //std::vector<TString> ft0Files = dm->GetFileNames(cmdft0);

    TString tok;
    Ssiz_t from = 0;
    TString output = gSystem->GetFromPipe(cmdfv0);
    while (output.Tokenize(tok, from, "\n")) fv0files.push_back(tok);

    if (fv0files.size()<0) {
        std::cout << "Could not retrieve files, stop macro" << std::endl;
        return;
    }

    TH2D *hHitsFV0[nfv0ch];
    for (int ihist = 0; ihist < nfv0ch; ihist++) {
        hHitsFV0[ihist] = new TH2D(Form("hHitsFV0%d", ihist), Form("Channel %d", ihist), 100, -50., 50, 100, -50., 50.);
    }

    int nfile = fv0Files.size();
    for (int ifile = 0; ifile < nfile; ifile++) {

        TFile *fin = TFile::Open(fv0Files[ifile]);
        TTree *hitsTree = (TTree*)fin->Get("o2sim");

        std::vector<o2::fv0::Hit> fv0hits, *fv0hitsPtr = &fv0hits;
        hitsTree->SetBranchAddress("FV0Hit", &fv0hitsPtr);

        for (UInt_t ient = 0; ient < nEntries; ient++) {
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
    }
    
}
