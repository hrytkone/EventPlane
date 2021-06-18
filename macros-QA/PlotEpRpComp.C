void PlotEpRpComp(TString file="rpep.root")
{
    TFile *fin = TFile::Open(file.Data(), "READ");
    fin->GetListOfKeys()->Print();

    TH1D *hDiff[8];
    TH1D *hCosDiff[8];
    for (int ih = 0; ih < 8; ih++) {
         hDiff[ih] = (TH1D*)fin->Get(Form("hDiff%i%i", ih, ih+1));
         hCosDiff[ih] = (TH1D*)fin->Get(Form("hCosDiff%i%i", ih, ih+1));
    }

    TString phi[9] = {"-#pi/2", "-3#pi/8", "-#pi/4", "-#pi/8", "0", "#pi/8", "#pi/4", "3#pi/8", "#pi/2"};
    
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 800);
    c1->Divide(4, 2);
    for (int ic = 0; ic < 8; ic++) {
        c1->cd(ic+1);
        //hDiff[ic]->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
        hDiff[ic]->SetTitle(Form("[%s, %s]; #Psi_{RP}-#Psi_{EP}; counts", phi[ic].Data(), phi[ic+1].Data()));
        hDiff[ic]->Draw("HIST");
    }
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 800);
    c2->Divide(4, 2);
    for (int ic = 0; ic < 8; ic++) {
        c2->cd(ic+1);
        hCosDiff[ic]->SetTitle(Form("[%s, %s]; cos(2(#Psi_{RP}-#Psi_{EP})); counts", phi[ic].Data(), phi[ic+1].Data()));
        hCosDiff[ic]->Draw("HIST");
    }

    c1->SaveAs("diff.pdf");
    c2->SaveAs("cosdiff.pdf");
}
