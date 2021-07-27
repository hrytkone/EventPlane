void RunAnalysis(TString kinefile="", TString clustfile="", TString output="") {
   
    TString lib = "/scratch/project_2003583/focal_test/FOCAL/FOCAL"; 
    gROOT->ProcessLine(Form(".x %s/LoadFOCAL.C(\"%s\")", lib.Data(), lib.Data()));
    //gROOT->ProcessLine(Form(".x /scratch/project_2003583/focal_test/focal_macros/AnalysePythiaData.C(\"%s\", \"%s\", \"%s\")", kinefile.Data(), clustfile.Data(), output.Data()));
    gROOT->ProcessLine(Form(".x /scratch/project_2003583/focal_test/focal_macros/MassAnalysis_matching.C(\"%s\", \"%s\", \"%s\")", kinefile.Data(), clustfile.Data(), output.Data()));
    // THESE NEED TO BE FIXED
    //gROOT->ProcessLine(Form(".x /scratch/project_2003583/focal_test/focal_macros/MatchClusters.C(\"%s\")", input.Data()));
    //gROOT->ProcessLine(Form(".x /scratch/project_2003583/focal_test/focal_macros/TestLeadingCorr.C(\"%s\")", input.Data()));
}
