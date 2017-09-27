#include "TFile.h" 
#include "TTree.h" 
#include "TH1F.h" 
#include "TGraphAsymmErrors.h" 

void ComputeEfficiency(std::string inputs,std::string outputs, std::string Name)
{
    TFile* inputFile = TFile::Open(inputs.c_str());
    TFile* outputFile = new TFile(outputs.c_str(),"RECREATE");
    
    TTree* h4 = (TTree*)inputFile->Get("h4");

    TH1F* num = new TH1F("num","",28,1250,4050);
    TH1F* den = new TH1F("den","",28,1250,4050);

    //h4->Draw("HV25-HVAMP >> num","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200)");
    //h4->Draw("HV25-HVAMP >> den","amp_max[MiB2]>200 && amp_max[Rm2]>200");

    h4->Draw("HV10-HVAMP >> num","amp_max[M10]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200)");
    h4->Draw("HV10-HVAMP >> den","amp_max[MiB2]>200 && amp_max[Rm2]>200");
    
    TGraphAsymmErrors* eff = new TGraphAsymmErrors(num,den);

    eff->Draw("AP"); 
     
    outputFile->cd();
    eff->Write(Name.c_str()); 
    outputFile->Close();
}
	
