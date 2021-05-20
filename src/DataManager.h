#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <iostream>
#include <vector>

#include "TString.h"
#include "TSystem.h"

class DataManager {

public:
    DataManager();
    virtual ~DataManager() {;}

    std::vector<TString> GetFileNames(TString cmd);
};

#endif
