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

void ComputeEfficiencyCam() {

  gStyle->SetOptStat(0);

  // inputs
  TFile* inFileMib10_700  = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanMib10_MID700.root");
  TFile* inFileMib10_1000 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanMib10_MID1000.root");

  TFile* inFileMib25_700  = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanMib25_MID700.root");
  TFile* inFileMib25_1000 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanMib25_MID1000.root");
  TFile* inFileMib25_1200 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanMib25_MID1200.root");

  TFile* inFileRm5_700  = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanRm5_MID700.root");
  TFile* inFileRm5_1000 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanRm5_MID1000.root");
  TFile* inFileRm5_1200 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanRm5_MID1200.root");

  TFile* inFileRm8_700  = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanRm8_MID700.root");
  TFile* inFileRm8_1000 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanRm8_MID1000.root");
  TFile* inFileRm8_1200 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/btf2016_HVscanRm8_MID1200.root");

  TTree* h4Mib10_700  = (TTree*)inFileMib10_700->Get("h4");
  TTree* h4Mib10_1000 = (TTree*)inFileMib10_1000->Get("h4");

  TTree* h4Mib25_700  = (TTree*)inFileMib25_700->Get("h4");
  TTree* h4Mib25_1000 = (TTree*)inFileMib25_1000->Get("h4");
  TTree* h4Mib25_1200 = (TTree*)inFileMib25_1200->Get("h4");

  TTree* h4Rm5_700  = (TTree*)inFileRm5_700->Get("h4");
  TTree* h4Rm5_1000 = (TTree*)inFileRm5_1000->Get("h4");
  TTree* h4Rm5_1200 = (TTree*)inFileRm5_1200->Get("h4");

  TTree* h4Rm8_700  = (TTree*)inFileRm8_700->Get("h4");
  TTree* h4Rm8_1000 = (TTree*)inFileRm8_1000->Get("h4");
  TTree* h4Rm8_1200 = (TTree*)inFileRm8_1200->Get("h4");

  // histos
  TH1F* numMib10_700  = new TH1F("numMib10_700", "",2600,700.,3300.);
  TH1F* numMib10_1000 = new TH1F("numMib10_1000","",2600,700.,3300.);

  TH1F* numMib25_700  = new TH1F("numMib25_700", "",2600,700.,3300.);
  TH1F* numMib25_1000 = new TH1F("numMib25_1000","",2600,700.,3300.);
  TH1F* numMib25_1200 = new TH1F("numMib25_1200","",2600,700.,3300.);

  TH1F* numRm5_700  = new TH1F("numRm5_700", "",2600,700.,3300.);
  TH1F* numRm5_1000 = new TH1F("numRm5_1000","",2600,700.,3300.);
  TH1F* numRm5_1200 = new TH1F("numRm5_1200","",2600,700.,3300.);

  TH1F* numRm8_700  = new TH1F("numRm8_700", "",2600,700.,3300.);
  TH1F* numRm8_1000 = new TH1F("numRm8_1000","",2600,700.,3300.);
  TH1F* numRm8_1200 = new TH1F("numRm8_1200","",2600,700.,3300.);

  TH1F* denMib10_700  = new TH1F("denMib10_700", "",2600,700.,3300.);
  TH1F* denMib10_1000 = new TH1F("denMib10_1000","",2600,700.,3300.);

  TH1F* denMib25_700  = new TH1F("denMib25_700", "",2600,700.,3300.);
  TH1F* denMib25_1000 = new TH1F("denMib25_1000","",2600,700.,3300.);
  TH1F* denMib25_1200 = new TH1F("denMib25_1200","",2600,700.,3300.);

  TH1F* denRm5_700  = new TH1F("denRm5_700", "",2600,700.,3300.);
  TH1F* denRm5_1000 = new TH1F("denRm5_1000","",2600,700.,3300.);
  TH1F* denRm5_1200 = new TH1F("denRm5_1200","",2600,700.,3300.);

  TH1F* denRm8_700  = new TH1F("denRm8_700", "",2600,700.,3300.);
  TH1F* denRm8_1000 = new TH1F("denRm8_1000","",2600,700.,3300.);
  TH1F* denRm8_1200 = new TH1F("denRm8_1200","",2600,700.,3300.);
  
  // denominators
  TString commonDen = "adc_data>200 && adc_data<700 && n_hitsX>0 && n_hitsY>0 && n_hitsX<3 && n_hitsY<3 && amp_max[MiB2]>200 && amp_max[Rm2]>200 && amp_max[Rm2]<1200 && fabs(time_max[MiB2])<150 && fabs(time_max[Rm2])<150";
  TString commonDen700  = commonDen + " && HVAMP==700";
  TString commonDen1000 = commonDen + " && HVAMP==1000";
  TString commonDen1200 = commonDen + " && HVAMP==1200";

  h4Mib10_700 ->Project("denMib10_700", "HV10-HVAMP",commonDen700);
  h4Mib10_1000->Project("denMib10_1000","HV10-HVAMP",commonDen1000);

  h4Mib25_700 ->Project("denMib25_700", "HV25-HVAMP",commonDen700);
  h4Mib25_1000->Project("denMib25_1000","HV25-HVAMP",commonDen1000);
  h4Mib25_1200->Project("denMib25_1200","HV25-HVAMP",commonDen1200);

  h4Rm5_700 ->Project("denRm5_700", "HV5-HVAMP",commonDen700);    
  h4Rm5_1000->Project("denRm5_1000","HV5-HVAMP",commonDen1000);
  h4Rm5_1200->Project("denRm5_1200","HV5-HVAMP",commonDen1200);

  h4Rm8_700 ->Project("denRm8_700", "HV8-HVAMP",commonDen700);    
  h4Rm8_1000->Project("denRm8_1000","HV8-HVAMP",commonDen1000);
  h4Rm8_1200->Project("denRm8_1200","HV8-HVAMP",commonDen1200);

  // numerators
  TString snumMib10_700  = commonDen700  + " && amp_max[M10]>20";
  TString snumMib10_1000 = commonDen1000 + " && amp_max[M10]>20";
  h4Mib10_700 ->Project("numMib10_700", "HV10-HVAMP",snumMib10_700); 
  h4Mib10_1000->Project("numMib10_1000","HV10-HVAMP",snumMib10_1000); 

  TString snumMib25_700  = commonDen700  + " && amp_max[M25]>17";
  TString snumMib25_1000 = commonDen1000 + " && amp_max[M25]>17";
  TString snumMib25_1200 = commonDen1200 + " && amp_max[M25]>17";
  h4Mib25_700 ->Project("numMib25_700", "HV25-HVAMP",snumMib25_700); 
  h4Mib25_1000->Project("numMib25_1000","HV25-HVAMP",snumMib25_1000); 
  h4Mib25_1200->Project("numMib25_1200","HV25-HVAMP",snumMib25_1200); 

  TString snumRm5_700  = commonDen700  + " && amp_max[M5]>21"; 
  TString snumRm5_1000 = commonDen1000 + " && amp_max[M5]>21";
  TString snumRm5_1200 = commonDen1200 + " && amp_max[M5]>21";
  h4Rm5_700 ->Project("numRm5_700", "HV5-HVAMP",snumRm5_700);    
  h4Rm5_1000->Project("numRm5_1000","HV5-HVAMP",snumRm5_1000); 
  h4Rm5_1200->Project("numRm5_1200","HV5-HVAMP",snumRm5_1200); 

  TString snumRm8_700  = commonDen700  + " && amp_max[M8]>16";
  TString snumRm8_1000 = commonDen1000 + " && amp_max[M8]>16";
  TString snumRm8_1200 = commonDen1200 + " && amp_max[M8]>16";
  h4Rm8_700 ->Project("numRm8_700", "HV8-HVAMP",snumRm8_700);    
  h4Rm8_1000->Project("numRm8_1000","HV8-HVAMP",snumRm8_1000); 
  h4Rm8_1200->Project("numRm8_1200","HV8-HVAMP",snumRm8_1200); 


  // efficiencies
  TGraphAsymmErrors* effMib10_700  = new TGraphAsymmErrors(numMib10_700, denMib10_700);
  TGraphAsymmErrors* effMib10_1000 = new TGraphAsymmErrors(numMib10_1000,denMib10_1000);
  effMib10_700->SetMarkerStyle(20);
  effMib10_700->SetMarkerColor(1);
  effMib10_700->SetLineColor(1);
  effMib10_1000->SetMarkerStyle(21);
  effMib10_1000->SetMarkerColor(2);
  effMib10_1000->SetLineColor(2);

  TGraphAsymmErrors* effMib25_700  = new TGraphAsymmErrors(numMib25_700, denMib25_700);
  TGraphAsymmErrors* effMib25_1000 = new TGraphAsymmErrors(numMib25_1000,denMib25_1000);
  TGraphAsymmErrors* effMib25_1200 = new TGraphAsymmErrors(numMib25_1200,denMib25_1200);
  effMib25_700->SetMarkerStyle(23);
  effMib25_700->SetMarkerColor(4);
  effMib25_700->SetLineColor(4);
  effMib25_1000->SetMarkerStyle(24);
  effMib25_1000->SetMarkerColor(6);
  effMib25_1000->SetLineColor(6);
  effMib25_1200->SetMarkerStyle(25);
  effMib25_1200->SetMarkerColor(7);
  effMib25_1200->SetLineColor(7);

  TGraphAsymmErrors* effRm5_700  = new TGraphAsymmErrors(numRm5_700, denRm5_700);
  TGraphAsymmErrors* effRm5_1000 = new TGraphAsymmErrors(numRm5_1000,denRm5_1000);
  TGraphAsymmErrors* effRm5_1200 = new TGraphAsymmErrors(numRm5_1200,denRm5_1200);
  effRm5_700->SetMarkerStyle(20);
  effRm5_700->SetMarkerColor(1);
  effRm5_700->SetLineColor(1);
  effRm5_1000->SetMarkerStyle(21);
  effRm5_1000->SetMarkerColor(2);
  effRm5_1000->SetLineColor(2);
  effRm5_1200->SetMarkerStyle(22);
  effRm5_1200->SetMarkerColor(3);
  effRm5_1200->SetLineColor(3);

  TGraphAsymmErrors* effRm8_700  = new TGraphAsymmErrors(numRm8_700, denRm8_700);
  TGraphAsymmErrors* effRm8_1000 = new TGraphAsymmErrors(numRm8_1000,denRm8_1000);
  TGraphAsymmErrors* effRm8_1200 = new TGraphAsymmErrors(numRm8_1200,denRm8_1200);
  effRm8_700->SetMarkerStyle(23);
  effRm8_700->SetMarkerColor(4);
  effRm8_700->SetLineColor(4);
  effRm8_1000->SetMarkerStyle(24);
  effRm8_1000->SetMarkerColor(6);
  effRm8_1000->SetLineColor(6);
  effRm8_1200->SetMarkerStyle(25);
  effRm8_1200->SetMarkerColor(7);
  effRm8_1200->SetLineColor(7);


  // Rescaling for pi/ele ratio
  // Eff on electrons = 0.933
  // Eff on pions     = 0.862
  // Multiplicity: eff(n) = 1 - (1-eff1)^n ==> n = ln(1.-effPi-deltaEff) / ln(1-effPi)
  // Correction factor ==> eff(1part) = 1 - pow( (1-effEle), (1/n) );
  float myN = log(1.-0.862-(0.933-0.862)) / log(1-0.862);

  TGraphAsymmErrors* effMib25Corr_1200 = new TGraphAsymmErrors(numMib25_1200,denMib25_1200);
  effMib25Corr_1200->SetMarkerStyle(20);
  effMib25Corr_1200->SetMarkerColor(1);
  effMib25Corr_1200->SetLineColor(1);
  int N_effMib25Corr_1200 = effMib25Corr_1200->GetN();
  for (int ii = 0; ii<N_effMib25Corr_1200; ii++ ) {
    float thisValue = effMib25Corr_1200->GetY()[ii];
    float theCorr   = 1 - pow( (1-thisValue), (1/myN) );  
    if (effMib25Corr_1200->GetY()[ii]>0) effMib25Corr_1200->GetY()[ii] *= theCorr/thisValue; 
  }

  TGraphAsymmErrors* effRm8Corr_1200 = new TGraphAsymmErrors(numRm8_1200,denRm8_1200);
  effRm8Corr_1200->SetMarkerStyle(20);
  effRm8Corr_1200->SetMarkerColor(1);
  effRm8Corr_1200->SetLineColor(1);
  int N_effRm8Corr_1200 = effRm8Corr_1200->GetN();
  for (int ii = 0; ii<N_effRm8Corr_1200; ii++ ) {
    float thisValue = effRm8Corr_1200->GetY()[ii];
    float theCorr   = 1 - pow( (1-thisValue), (1/myN) );  
    if (effRm8Corr_1200->GetY()[ii]>0) effRm8Corr_1200->GetY()[ii] *= theCorr/thisValue; 
  }

  TGraphAsymmErrors* effRm5Corr_1200 = new TGraphAsymmErrors(numRm5_1200,denRm5_1200);
  effRm5Corr_1200->SetMarkerStyle(20);
  effRm5Corr_1200->SetMarkerColor(1);
  effRm5Corr_1200->SetLineColor(1);
  int N_effRm5Corr_1200 = effRm5Corr_1200->GetN();
  for (int ii = 0; ii<N_effRm5Corr_1200; ii++ ) {
    float thisValue = effRm5Corr_1200->GetY()[ii];
    float theCorr   = 1 - pow( (1-thisValue), (1/myN) );  
    if (effRm5Corr_1200->GetY()[ii]>0) effRm5Corr_1200->GetY()[ii] *= theCorr/thisValue; 
  }
  
  TFile myOutputFile("myOutputFileCameroniBtf.root", "RECREATE");
  myOutputFile.cd();
  effMib25_1200->Write("effMib25_1200");
  effMib25Corr_1200->Write("effMib25Corr_1200");
  effRm8_1200->Write("effRm8_1200");
  effRm8Corr_1200->Write("effRm8Corr_1200");
  effRm5_1200->Write("effRm5_1200");
  effRm5Corr_1200->Write("effRm5Corr_1200");
  myOutputFile.Close();

  TH2F* H2 = new TH2F("H2","",2600,700.,3300.,100,0.,1.);
  H2->GetXaxis()->SetTitle("HV - HVMID [V]");
  H2->GetYaxis()->SetTitle("Efficiency");

  TLegend *legA;
  legA = new TLegend(0.1,0.6,0.4,0.9);
  legA->SetFillStyle(0);
  legA->SetBorderSize(0);
  legA->SetTextSize(0.05);
  legA->SetFillColor(0);
  legA->AddEntry(effMib10_700,  "90x2+40x1, 10#mu (Mib10). HVMID=700", "pl");
  legA->AddEntry(effMib10_1000, "90x2+40x1, 10#mu (Mib10). HVMID=1000", "pl");

  TLegend *legB;
  legB = new TLegend(0.1,0.6,0.4,0.9);
  legB->SetFillStyle(0);
  legB->SetBorderSize(0);
  legB->SetTextSize(0.05);
  legB->SetFillColor(0);
  legB->AddEntry(effMib25_700,  "90x2+40x1, 25#mu (Mib25). HVMID=700", "pl");
  legB->AddEntry(effMib25_1000, "90x2+40x1, 25#mu (Mib25). HVMID=1000", "pl");
  legB->AddEntry(effMib25_1200, "90x2+40x1, 25#mu (Mib25). HVMID=1200", "pl");

  TLegend *legC;
  legC = new TLegend(0.1,0.6,0.4,0.9);
  legC->SetFillStyle(0);
  legC->SetBorderSize(0);
  legC->SetTextSize(0.05);
  legC->SetFillColor(0);
  legC->AddEntry(effMib10_700,  "Mib10, MID=700", "pl");
  legC->AddEntry(effMib10_1000, "Mib10, MID=1000", "pl");
  legC->AddEntry(effMib25_700,  "Mib25, MID=700", "pl");
  legC->AddEntry(effMib25_1000, "Mib25, MID=1000", "pl");
  legC->AddEntry(effMib25_1200, "Mib25, MID=1200", "pl");

  TLegend *legA2;
  legA2 = new TLegend(0.1,0.6,0.4,0.9);
  legA2->SetFillStyle(0);
  legA2->SetBorderSize(0);
  legA2->SetTextSize(0.05);
  legA2->SetFillColor(0);
  legA2->AddEntry(effRm8_700,  "90x1+40x1, 8#mu (Rm8). HVMID=700", "pl");
  legA2->AddEntry(effRm8_1000, "90x1+40x1, 8#mu (Rm8). HVMID=1000", "pl");
  legA2->AddEntry(effRm8_1200, "90x1+40x1, 8#mu (Rm8). HVMID=1200", "pl");

  TLegend *legB2;
  legB2 = new TLegend(0.1,0.6,0.4,0.9);
  legB2->SetFillStyle(0);
  legB2->SetBorderSize(0);
  legB2->SetTextSize(0.05);
  legB2->SetFillColor(0);
  legB2->AddEntry(effRm5_700,  "90x1+40x1, 5#mu (Rm5). HVMID=700", "pl");
  legB2->AddEntry(effRm5_1000, "90x1+40x1, 5#mu (Rm5). HVMID=1000", "pl");
  legB2->AddEntry(effRm5_1200, "90x1+40x1, 5#mu (Rm5). HVMID=1200", "pl");

  TLegend *legC2;
  legC2 = new TLegend(0.6,0.6,0.9,0.9);
  legC2->SetFillStyle(0);
  legC2->SetBorderSize(0);
  legC2->SetTextSize(0.05);
  legC2->SetFillColor(0);
  legC2->AddEntry(effRm5_700,  "Rm5, MID=700",  "pl");
  legC2->AddEntry(effRm5_1000, "Rm5, MID=1000", "pl");
  legC2->AddEntry(effRm5_1200, "Rm5, MID=1200", "pl");
  legC2->AddEntry(effRm8_700,  "Rm8, MID=700",  "pl");
  legC2->AddEntry(effRm8_1000, "Rm8, MID=1000", "pl");
  legC2->AddEntry(effRm8_1200, "Rm8, MID=1200", "pl");

  TCanvas* c1 = new TCanvas("c1","c",1);
  c1->SetGrid();
  H2->Draw();
  effMib10_700->Draw("Plsame"); 
  effMib10_1000->Draw("Plsame"); 
  legA->Draw();
  c1->SaveAs("effMib10.png");

  TCanvas* c2 = new TCanvas("c2","c",1);
  c2->SetGrid();
  H2->Draw();
  effMib25_700->Draw("Plsame"); 
  effMib25_1000->Draw("Plsame"); 
  effMib25_1200->Draw("Plsame"); 
  legB->Draw();
  c2->SaveAs("effMib25.png");

  TCanvas* cc2 = new TCanvas("cc2","c",1);
  cc2->SetGrid();
  H2->Draw();
  effMib25_1200->Draw("Plsame"); 
  effMib25Corr_1200->Draw("Plsame"); 
  cc2->SaveAs("effMib25_prePostMultipCorr.png");

  TCanvas* c3 = new TCanvas("c3","c",1);
  c3->SetGrid();
  H2->Draw();
  effMib10_700->Draw("Plsame"); 
  effMib10_1000->Draw("Plsame"); 
  effMib25_700->Draw("Plsame"); 
  effMib25_1000->Draw("Plsame"); 
  effMib25_1200->Draw("Plsame"); 
  legC->Draw();
  c3->SaveAs("effMibAll.png");

  TCanvas* c2a = new TCanvas("c2a","c",1);
  c2a->SetGrid();
  H2->Draw();
  effRm8_700->Draw("Plsame"); 
  effRm8_1000->Draw("Plsame"); 
  effRm8_1200->Draw("Plsame"); 
  legA2->Draw();
  c2a->SaveAs("effRm8.png");

  TCanvas* cc2a = new TCanvas("cc2a","c",1);
  cc2a->SetGrid();
  H2->Draw();
  effRm8_1200->Draw("Plsame"); 
  effRm8Corr_1200->Draw("Plsame"); 
  cc2a->SaveAs("effRm8_prePostMultipCorr.png");

  TCanvas* c2b = new TCanvas("c2b","c",1);
  c2b->SetGrid();
  H2->Draw();
  effRm5_700->Draw("Plsame"); 
  effRm5_1000->Draw("Plsame"); 
  effRm5_1200->Draw("Plsame"); 
  legB2->Draw();
  c2b->SaveAs("effRm5.png");

  TCanvas* cc2b = new TCanvas("cc2b","c",1);
  cc2b->SetGrid();
  H2->Draw();
  effRm5_1200->Draw("Plsame"); 
  effRm5Corr_1200->Draw("Plsame"); 
  cc2b->SaveAs("effRm5_prePostMultipCorr.png");

  TCanvas* c3b = new TCanvas("c3b","c",1);
  c3b->SetGrid();
  H2->Draw();
  effRm5_700->Draw("Plsame"); 
  effRm5_1000->Draw("Plsame"); 
  effRm5_1200->Draw("Plsame"); 
  effRm8_700->Draw("Plsame"); 
  effRm8_1000->Draw("Plsame"); 
  effRm8_1200->Draw("Plsame"); 
  legC2->Draw();
  c3b->SaveAs("effRmAll.png");


  // Spectra comparison
  TH1F *ampRm8  = new TH1F("ampRm8", "ampRm8", 50,0.,4000);
  TH1F *ampRm5  = new TH1F("ampRm5", "ampRm5", 50,0.,4000);
  TH1F *ampMi10 = new TH1F("ampMi10","ampMi10",50,0.,4000);
  TH1F *ampMi25 = new TH1F("ampMi25","ampMi25",50,0.,4000);
  ampRm8  -> SetLineColor(1);
  ampRm5  -> SetLineColor(2);
  ampMi10 -> SetLineColor(3);
  ampMi25 -> SetLineColor(4);
  ampRm8  -> SetLineWidth(2);
  ampRm5  -> SetLineWidth(2);
  ampMi10 -> SetLineWidth(2);
  ampMi25 -> SetLineWidth(2);

  h4Rm8_1200  ->Project("ampRm8", "amp_max[M8]", commonDen1200+" && HV8==2800");
  h4Rm5_1200  ->Project("ampRm5", "amp_max[M5]", commonDen1200+" && HV5==2800");
  h4Mib25_1200->Project("ampMi25","amp_max[M25]",commonDen1200+" && HV25==3600");

  TLegend *legAmp;
  legAmp = new TLegend(0.5,0.6,0.9,0.9);
  legAmp->SetFillStyle(0);
  legAmp->SetBorderSize(0);
  legAmp->SetTextSize(0.05);
  legAmp->SetFillColor(0);
  legAmp->AddEntry(ampRm8,  "Rm8",   "pl");
  legAmp->AddEntry(ampRm5,  "Rm5",   "pl");
  legAmp->AddEntry(ampMi25, "MiB25", "pl");

  TCanvas* c4 = new TCanvas("c4","c",1);
  c4->SetGrid();
  ampRm8->SetTitle("MID=1200");
  ampRm5->SetTitle("MID=1200");
  ampMi25->SetTitle("MID=1200");
  ampRm8->DrawNormalized(); 
  ampRm5->DrawNormalized("same"); 
  ampMi25->DrawNormalized("same"); 
  legAmp->Draw();
  c4->SetLogy();
  c4->SaveAs("amplitude.png");

}
	
