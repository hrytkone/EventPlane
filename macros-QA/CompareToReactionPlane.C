void CompareToReactionPlane(TString sDirName)   
{
    TString digitfile(Form("%s/fv0digits.root", sDirName.Data()));
    TString mcfile(Form("%s/o2sim_Kine.root", sDirName.Data()));

    TFile *fin = TFile::Open(mcfile);
    TTree *kineTree = (TTree*)fin->Get("o2sim");

    std::vector<o2::MCTrack>* mctrack = nullptr;
    auto mcbr = kineTree->GetBranch("MCTrack");
    mcbr->SetAddress(&mctrack);

    o2::dataformats::MCEventHeader* mcheader = nullptr;
    auto hdrbr = kineTree->GetBranch("MCEventHeader.");
    hdrbr->SetAddress(&mcheader);

    int nent = kineTree->GetEntries();
    for (int ient = 0; ient < nent; ient++) {
        
        kineTree->GetEntry(ient);

        int ntrack = 0;
			
		double rp = mcheader->GetRotZ();
		std::cout << rp << std::endl;

    }
}
