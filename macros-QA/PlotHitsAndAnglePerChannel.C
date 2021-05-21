/**
 *	Macro to check if the phi angle in the code alings with the
 *	given channel.
 *
 *	Uses one hits file for FV0 and FT0 from simulation directory
 *	given as argument.
 */

void PlotHitsAndAnglePerChannel(TString sDir)
{

	TString cmdfv0(Form("ls %s*/.root", sDirName.Data()));
	TString cmdft0(Form("ls %s*/fv0digits.root", sDirName.Data()));
	std::vector<TString> fv0Files = dm->GetFileNames(cmdfv0);
}
