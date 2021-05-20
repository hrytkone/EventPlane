#include "Corrections.h"

Corrections::Corrections(TString filename)
{
    saveFileName = filename;
    for (int i = 0; i < nCorrections; i++) corrections[i] = 0;
}

double *Corrections::CalculateCorrections(TH2D *hQ)
{
    double xmean, ymean, xdev, ydev, aplus, aminus, lambdaplus, lambdaminus;
    GetRecenteringCorrection(hQ, xmean, ymean);
    GetWidthCorrection(hQ, xdev, ydev);
    GetTwistAndRescaleCorrection(hQ, aplus, aminus, lambdaplus, lambdaminus);

    corrections[0] = xmean;
    corrections[1] = ymean;
    corrections[2] = xdev;
    corrections[3] = ydev;
    corrections[4] = aplus;
    corrections[5] = aminus;
    corrections[6] = lambdaplus;
    corrections[7] = lambdaminus;

    return corrections;
}

void Corrections::GetRecenteringCorrection(TH2D *hQ, double &xmean, double &ymean)
{
    xmean = hQ->GetMean(1);
    ymean = hQ->GetMean(2);
}

void Corrections::GetWidthCorrection(TH2D *hQ, double &xdev, double &ydev)
{
    xdev = hQ->GetStdDev(1);
    ydev = hQ->GetStdDev(2);
}

void Corrections::GetTwistAndRescaleCorrection(TH2D *hQ, double &aplus, double &aminus, double &lambdaplus, double &lambdaminus)
{
    double rho = hQ->GetCorrelationFactor();
    double sigmax = hQ->GetStdDev(1);
    double sigmay = hQ->GetStdDev(2);
    
    double b = rho*sigmax*sigmay*TMath::Sqrt(2.0*(sigmax*sigmax + sigmay*sigmay 
                - 2.0*sigmax*sigmay*TMath::Sqrt(1.0 - rho*rho))/((sigmax*sigmax 
                + sigmay*sigmay)*(sigmax*sigmax + sigmay*sigmay) 
                + 4.0*(sigmax*sigmay*rho)*(sigmax*sigmay*rho)));

    aplus = TMath::Sqrt(2.0*sigmax*sigmax - b*b);
    aminus = TMath::Sqrt(2.0*sigmay*sigmay - b*b);

    lambdaplus = b/aplus;
    lambdaminus = b/aminus;
}

void Corrections::SaveCorrections(TH2D *hQvec)
{
    std::cout << "\nCalculating corrections.." << std::endl;
    CalculateCorrections(hQvec);
    outputFile.open(saveFileName.Data());
    for (int i = 0; i < nCorrections; i++) {
        outputFile << corrections[i] << "\n";
    }
    outputFile.close();
    std::cout << "Corrections saved to file " << saveFileName << std::endl;
}

void Corrections::LoadCorrections(TString filename)
{
    inputFile.open(filename);
    std::string line;
    int icorr = 0;
    while (std::getline(inputFile, line)) {
        corrections[icorr] = atof(line.c_str());
        icorr++;
    }
    inputFile.close();
}

void Corrections::DoCorrections(double &qx, double &qy)
{
    // 1) Recentering
    qx -= corrections[0];
    qy -= corrections[1];
    
    // 2) Twist
    qx = (qx - corrections[7]*qy)/(1.0 - corrections[7]*corrections[7]);
    qy = (qy - corrections[6]*qx)/(1.0 - corrections[6]*corrections[6]);

    // 3) Rescale
    if (corrections[4]==0 || corrections[5]==0) {
        std::cout << "Correction warning: aplus or aminus is equal to zero, rescale not applied!" << std::endl;
    } else {
        qx /= corrections[4];
        qy /= corrections[5];
    }
}

void Corrections::Print() 
{
    std::cout << "<Q>      : " << corrections[0] << "  " << corrections[1] << std::endl;
    std::cout << "a+-      : " << corrections[4] << "  " << corrections[5] << std::endl;
    std::cout << "lambda+- : " << corrections[6] << "  " << corrections[7] << std::endl;
}   
