#include "src/Eventplane.h"
#include "src/Corrections.h"

TFile *fin;

TH2D *hQvecfv0;
TH2D *hQvecft0a;
TH2D *hQvecft0c;

Eventplane *ep;
Corrections *corrFV0;
Corrections *corrFT0A;
Corrections *corrFT0C;

int LoadInput(TString inputfile);
void CloseFiles();

void CalculateCorrections()
{
	int fileopened = LoadInput(inputfile);
	if (!fileopened) return;

    ep = new Eventplane(bmin, bmax);
    corrFV0 = new Corrections("corrections_fv0.txt");
    corrFT0A = new Corrections("corrections_ft0a.txt");
    corrFT0C = new Corrections("corrections_ft0c.txt");	

    corrFV0->SaveCorrections(hQvecfv0);
    corrFT0A->SaveCorrections(hQvecft0a);
    corrFT0C->SaveCorrections(hQvecft0c);
}

//_____________________________________________________________

int LoadInput(TString inputfile)
{
    if (!gSystem->AccessPathName(inputfile.Data())) {
        fin = TFile::Open(inputfile.Data(), "READ");
        hQvecfv0 = (TH1D*)fin->Get("hQvecfv0");
        hQvecft0a = (TH1D*)fin->Get("hQvecft0a");
        hQvecft0c = (TH1D*)fin->Get("hQvecft0c");
    } else {
    	std::cout << "Input file not found! Stopping.." << std::endl;
    	return 0;
    }
    return 1;
}

void CloseFiles()
{
	ep->CloseFiles();
	fout->Close();
}