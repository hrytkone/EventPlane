#include "DataManager.h"

DataManager::DataManager() 
{
}

std::vector<TString> DataManager::GetFileNames(TString cmd)
{
    std::vector<TString> files;

    TString tok;
    Ssiz_t from = 0;
    TString output = gSystem->GetFromPipe(cmd);
    while (output.Tokenize(tok, from, "\n")) files.push_back(tok);

    return files;
}
