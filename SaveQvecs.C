#include "src/Eventplane.h"

const double bmin = 0.; // To change impact parameter range if data is min bias
const double bmax = 20.;

TFile *fout;

TH2D *hQvecfv0;
TH2D *hQvecft0a;
TH2D *hQvecft0c;

Eventplane *ep;

void InitOutput(TString outfile);
void CloseFiles();

void SaveQvecs(TString kinefile="", TString fv0digitfile="", TString ft0digitfile="", TString outfile="output.root")
{
    ep = new Eventplane(bmin, bmax);
    ep->OpenFiles(kinefile, fv0digitfile, ft0digitfile);
            
    std::vector<TComplex> QvecAfv0 = ep->GetQvecA("FV0");
    for (int i = 0; i < (int)QvecAfv0.size(); i++) hQvecfv0->Fill(QvecAfv0[i].Re(), QvecAfv0[i].Im());            
            
    std::vector<TComplex> QvecAft0a = ep->GetQvecA("FT0A");
    for (int i = 0; i < (int)QvecAft0a.size(); i++) hQvecft0a->Fill(QvecAft0a[i].Re(), QvecAft0a[i].Im());            
            
    std::vector<TComplex> QvecAft0c = ep->GetQvecA("FT0C");
    for (int i = 0; i < (int)QvecAft0c.size(); i++) hQvecft0c->Fill(QvecAft0c[i].Re(), QvecAft0c[i].Im()); 

    fout->Write("", TObject::kOverWrite);
    CloseFiles();
}

//_____________________________________________________________

void InitOutput(TString outfile)
{
    if (!gSystem->AccessPathName(outfile.Data())) {
        fout = TFile::Open(outfile.Data(), "UPDATE");
        hQvecfv0 = (TH1D*)fout->Get("hQvecfv0");
        hQvecft0a = (TH1D*)fout->Get("hQvecft0a");
        hQvecft0c = (TH1D*)fout->Get("hQvecft0c");
    } else {
    	fout = TFile::Open(outfile.Data(), "RECREATE");
    	hQvecfv0 = new TH2D("hQvecfv0", "hQvecfv0", 251, -0.5, 0.5, 251);
    	hQvecft0a = new TH2D("hQvecft0a", "hQvecft0a", 251, -0.5, 0.5, 251);
    	hQvecft0c = new TH2D("hQvecft0c", "hQvecft0c", 251, -0.5, 0.5, 251);    	    	
    }
}

void CloseFiles()
{
	ep->CloseFiles();
	fout->Close();
	delete ep;
}