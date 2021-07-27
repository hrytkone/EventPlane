TH2D *hQvecfv0;
TH2D *hQvecft0a;
TH2D *hQvecft0c;

void CalculateCorrections()
{

    hQvecfv0 = new TH2D("hQvecfv0", "hQvecfv0", 251, -0.5, 0.5, 251);
}

//_____________________________________________________________

void InitOutput(TString outfile)
{
    if (!gSystem->AccessPathName(outfile.Data())) {
        fout = TFile::Open(outfile.Data(), "UPDATE");
    } else {

    }
}
