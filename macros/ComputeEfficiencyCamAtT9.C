#include "TFile.h" 
#include "TTree.h" 
#include "TStyle.h" 
#include "TH1.h" 
#include "TH2.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h" 
#include <iostream>

using namespace std;

void ComputeEfficiencyCamAtT9() {

  gStyle->SetOptStat(0);

  // inputs
  TFile* inFileMib10_1000 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/T9_HVscanMib10mu_mid1000.root");
  TFile* inFileMib10_1200 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/T9_HVscanMib10mu_mid1200.root");

  TFile* inFileMib25_1000 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/T9_HVscanMib25mu_mid1000.root");
  TFile* inFileMib25_1200 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/T9_HVscanMib25mu_mid1200.root");
  TFile* inFileMib25_1300 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/T9_HVscanMib25mu_mid1300.root");

  TTree* h4Mib10_1000 = (TTree*)inFileMib10_1000->Get("h4");
  TTree* h4Mib10_1200 = (TTree*)inFileMib10_1200->Get("h4");

  TTree* h4Mib25_1000 = (TTree*)inFileMib25_1000->Get("h4");
  TTree* h4Mib25_1200 = (TTree*)inFileMib25_1200->Get("h4");
  TTree* h4Mib25_1300 = (TTree*)inFileMib25_1300->Get("h4");

  // pi-ele run
  TFile* inFileRun4594 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/t92016_Ele_Pre_B_4594.root");
  TTree* h4Run4594 = (TTree*)inFileRun4594->Get("h4");

  // comparison with pi run
  TFile* inFileRun4566 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/t92016_MiB_10m_B_4566.root");
  TTree* h4Run4566 = (TTree*)inFileRun4566->Get("h4");

  // histos
  TH1F* numMib10_1000 = new TH1F("numMib10_1000","",1200,1400.,2600.); 
  TH1F* numMib10_1200 = new TH1F("numMib10_1200","",1200,1400.,2600.); 

  TH1F* numMib25_1000 = new TH1F("numMib25_1000","",1200,1400.,2600.); 
  TH1F* numMib25_1200 = new TH1F("numMib25_1200","",1200,1400.,2600.); 
  TH1F* numMib25_1300 = new TH1F("numMib25_1300","",1200,1400.,2600.); 

  TH1F* denMib10_1000 = new TH1F("denMib10_1000","",1200,1400.,2600.); 
  TH1F* denMib10_1200 = new TH1F("denMib10_1200","",1200,1400.,2600.); 

  TH1F* denMib25_1000 = new TH1F("denMib25_1000","",1200,1400.,2600.); 
  TH1F* denMib25_1200 = new TH1F("denMib25_1200","",1200,1400.,2600.); 
  TH1F* denMib25_1300 = new TH1F("denMib25_1300","",1200,1400.,2600.); 

  
  // denominators
  TString commonDen = "amp_max[MiB2]>200 && amp_max[Rm2]>200 && fabs(time_max[MiB2])<150 && fabs(time_max[Rm2])<150 && amp_max[Cerenkov]<80";
  TString commonDen1000 = commonDen + " && HVAMP==1000";
  TString commonDen1200 = commonDen + " && HVAMP==1200";
  TString commonDen1300 = commonDen + " && HVAMP==1300";

  h4Mib10_1000->Project("denMib10_1000","HV10-HVAMP",commonDen1000);
  h4Mib10_1200->Project("denMib10_1200","HV10-HVAMP",commonDen1200);

  h4Mib25_1000->Project("denMib25_1000","HV25-HVAMP",commonDen1000);
  h4Mib25_1200->Project("denMib25_1200","HV25-HVAMP",commonDen1200);
  h4Mib25_1300->Project("denMib25_1300","HV25-HVAMP",commonDen1300);

  // numerators
  TString snumMib10_1000 = commonDen1000 + " && amp_max[M10]>17";
  TString snumMib10_1200 = commonDen1200 + " && amp_max[M10]>17";
  h4Mib10_1000->Project("numMib10_1000","HV10-HVAMP",snumMib10_1000); 
  h4Mib10_1200->Project("numMib10_1200","HV10-HVAMP",snumMib10_1200); 

  TString snumMib25_1000 = commonDen1000 + " && amp_max[M25]>17";
  TString snumMib25_1200 = commonDen1200 + " && amp_max[M25]>17";
  TString snumMib25_1300 = commonDen1300 + " && amp_max[M25]>17";
  h4Mib25_1000->Project("numMib25_1000","HV25-HVAMP",snumMib25_1000); 
  h4Mib25_1200->Project("numMib25_1200","HV25-HVAMP",snumMib25_1200); 
  h4Mib25_1300->Project("numMib25_1300","HV25-HVAMP",snumMib25_1300); 

  // efficiencies
  TGraphAsymmErrors* effMib10_1000 = new TGraphAsymmErrors(numMib10_1000,denMib10_1000);
  TGraphAsymmErrors* effMib10_1200 = new TGraphAsymmErrors(numMib10_1200,denMib10_1200);
  effMib10_1000->SetMarkerStyle(21);
  effMib10_1000->SetMarkerColor(2);
  effMib10_1000->SetLineColor(2);
  effMib10_1200->SetMarkerStyle(20);
  effMib10_1200->SetMarkerColor(1);
  effMib10_1200->SetLineColor(1);

  TGraphAsymmErrors* effMib25_1000 = new TGraphAsymmErrors(numMib25_1000,denMib25_1000);
  TGraphAsymmErrors* effMib25_1200 = new TGraphAsymmErrors(numMib25_1200,denMib25_1200);
  TGraphAsymmErrors* effMib25_1300 = new TGraphAsymmErrors(numMib25_1300,denMib25_1300);
  effMib25_1000->SetMarkerStyle(23);
  effMib25_1000->SetMarkerColor(4);
  effMib25_1000->SetLineColor(4);
  effMib25_1200->SetMarkerStyle(24);
  effMib25_1200->SetMarkerColor(6);
  effMib25_1200->SetLineColor(6);
  effMib25_1300->SetMarkerStyle(25);
  effMib25_1300->SetMarkerColor(7);
  effMib25_1300->SetLineColor(7);

  TH2F* H2 = new TH2F("H2","",1200,1400.,2600.,100,0.,1.);
  H2->GetXaxis()->SetTitle("HV - HVMID [V]");
  H2->GetYaxis()->SetTitle("Efficiency");

  TLegend *leg;
  leg = new TLegend(0.1,0.6,0.4,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(effMib10_1000, "Mib10, MID=1000", "pl");
  leg->AddEntry(effMib10_1200, "Mib10, MID=1200", "pl");
  leg->AddEntry(effMib25_1000, "Mib25, MID=1000", "pl");
  leg->AddEntry(effMib25_1200, "Mib25, MID=1200", "pl");
  leg->AddEntry(effMib25_1300, "Mib25, MID=1300", "pl");

  TCanvas* c1 = new TCanvas("c1","c",1);
  c1->SetGrid();
  H2->Draw();
  effMib10_1000->Draw("Plsame"); 
  effMib10_1200->Draw("Plsame"); 
  effMib25_1000->Draw("Plsame"); 
  effMib25_1200->Draw("Plsame"); 
  effMib25_1300->Draw("Plsame"); 
  leg->Draw();
  c1->SaveAs("effMibAll.png");


  // Comparison ele vs pions
  TString commonDenElePi = "amp_max[MiB2]>200 && amp_max[Rm2]>200 && fabs(time_max[MiB2])<150 && fabs(time_max[Rm2])<150 && HVAMP==1200";

  TString denEle = commonDenElePi + " && amp_max[Cerenkov]>80";
  TString denPi  = commonDenElePi + " && amp_max[Cerenkov]<80";
  TString numEle = denEle + " && amp_max[M10]>17";
  TString numPi  = denPi  + " && amp_max[M10]>17";

  TH1F* HdenEle = new TH1F("HdenEle","",1200,1400.,2600.); 
  TH1F* HdenPi  = new TH1F("HdenPi", "",1200,1400.,2600.); 
  h4Run4594->Project("HdenEle","HV10-HVAMP",denEle);
  h4Run4594->Project("HdenPi", "HV10-HVAMP",denPi);

  TH1F* HnumEle = new TH1F("HnumEle","",1200,1400.,2600.); 
  TH1F* HnumPi  = new TH1F("HnumPi", "",1200,1400.,2600.); 
  h4Run4594->Project("HnumEle","HV10-HVAMP",numEle);
  h4Run4594->Project("HnumPi", "HV10-HVAMP",numPi);

  cout << endl;
  cout << "Run4594: efficiency of electrons = " << HnumEle->Integral()/HdenEle->Integral() << endl;
  cout << "Run4594: efficiency of pions = "     << HnumPi->Integral()/HdenPi->Integral()   << endl;


  // Cross-check on pi run with identical voltage conditions
  TH1F* HdenPi4566 = new TH1F("HdenPi4566", "",1200,1400.,2600.); 
  TH1F* HnumPi4566 = new TH1F("HnumPi4566", "",1200,1400.,2600.); 
  h4Run4566->Project("HdenPi4566", "HV10-HVAMP",denPi);
  h4Run4566->Project("HnumPi4566", "HV10-HVAMP",numPi);
  cout << endl;
  cout << "Run4566: efficiency of pions = " << HnumPi4566->Integral()/HdenPi4566->Integral()   << endl;

}
	
