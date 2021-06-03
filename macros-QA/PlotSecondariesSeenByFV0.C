void PlotSecondariesSeenByFV0(TString sFileName="o2sim_Kine.root")
{
    TFile *fin = TFile::Open(sFileName.Data(), "READ");
    TTree *kineTree = (TTree*)fin->Get("o2sim");

    std::vector<o2::MCTrack> *mctrack = nullptr;
    auto mcbr = kineTree->GetBranch("MCTrack");
    mcbr->SetAddress(&mctrack);

    TH2D *hVrtxXY = new TH2D("hVrtxXY", "hVrtxXY", 301, -20., 20., 301, -20., 20.);
    TH2D *hVrtxXZ = new TH2D("hVrtxXZ", "hVrtxXZ", 301, -80., 80., 301, 314., 334.);

    int nent = kineTree->GetEntries();
    for (int ient = 0; ient < nent; ient++) {
        kineTree->GetEntry(ient);

        std::cout << ient+1 << "/" << nent << std::endl;

        for (auto t : *mctrack) {
           
            int pdg = t.GetPdgCode();
            if (pdg == 50000050) continue; 
            if (!t.leftTrace(o2::detectors::DetID::FV0)) continue;
            double x = t.Vx();
            double y = t.Vy();
            double z = t.Vz();

            if (z > 314. && z < 334.) {
                hVrtxXY->Fill(x, y); 
                hVrtxXZ->Fill(x, z); 
            }
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    hVrtxXY->Draw("COLZ");

    TCanvas *c2 = new TCanvas("c2", "c2");
    hVrtxXZ->Draw("COLZ");

}
