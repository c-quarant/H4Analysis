#include "TFile.h" 
#include "TTree.h" 
#include "TStyle.h" 
#include "TH1.h" 
#include "TH2.h" 
#include "TF1.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h" 
#include "TGraphAsymmErrors.h" 
#include <iostream>

using namespace std;

void ComputeEfficiencyH4() {

  gStyle->SetOptStat(0);

  // inputs
  TFile* inFileSee  = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/H4_scan5_SEE.root"); 
  TFile* inFileMcp2 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/H4_scan6_MCP2.root");
  TFile* inFileMcp5 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/H4_scan7_MCP5.root"); 
  TFile* inFileMcp4 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/H4_scan8_MCP4.root"); 
  TFile* inFilePMTMcp = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/H4_HVscan2_ZS2.root"); 
  TFile* inFileMgO = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/MgO.root"); 

  //
  TTree* h4See  = (TTree*)inFileSee->Get("h4");
  TTree* h4binpMcp2 = (TTree*)inFileMcp2->Get("h4");
  TTree* h4binpMcp5 = (TTree*)inFileMcp5->Get("h4");
  TTree* h4binpMcp4 = (TTree*)inFileMcp4->Get("h4");
  TTree* h4PMTMcp = (TTree*)inFilePMTMcp->Get("h4");
  TTree* h4MgO = (TTree*)inFileMgO->Get("h4");
  
  // input absorber scan     
  TFile* inFileAbsorb = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/MCP_MergedSamples/H4_absScan2_ZS2.root"); 
  TTree* h4Absorb = (TTree*)inFileAbsorb->Get("h4");


  // ----------------------------------------

  TH1F* numSee  = new TH1F("numSee", "",3500,0.,3500.);
  TH1F* numMcp2 = new TH1F("numMcp2","",3500,0.,3500.);
  TH1F* numMcp5 = new TH1F("numMcp5","",3500,0.,3500.);
  TH1F* numMcp4 = new TH1F("numMcp4","",3500,0.,3500.);
  TH1F* numPMTMcp = new TH1F("numPMTMcp","",3500,0.,3500.);
  TH1F* numMgO_PadA = new TH1F("numMgO_PadA","",3500,0.,3500.);
  TH1F* numMgO_PadB = new TH1F("numMgO_PadB","",3500,0.,3500.);
  TH1F* numMgO_PadC = new TH1F("numMgO_PadC","",3500,0.,3500.);
  TH1F* numMgO_PadD = new TH1F("numMgO_PadD","",3500,0.,3500.);

  TH1F* denSee  = new TH1F("denSee", "",3500,0.,3500.);
  TH1F* denMcp2 = new TH1F("denMcp2","",3500,0.,3500.);
  TH1F* denMcp5 = new TH1F("denMcp5","",3500,0.,3500.);
  TH1F* denMcp4 = new TH1F("denMcp4","",3500,0.,3500.);
  TH1F* denPMTMcp = new TH1F("denPMTMcp","",3500,0.,3500.);
  TH1F* denMgO_PadA = new TH1F("denMgO_PadA","",3500,0.,3500.);
  TH1F* denMgO_PadB = new TH1F("denMgO_PadB","",3500,0.,3500.);
  TH1F* denMgO_PadC = new TH1F("denMgO_PadC","",3500,0.,3500.);
  TH1F* denMgO_PadD = new TH1F("denMgO_PadD","",3500,0.,3500.);
  
  // Nominal denominator: selection based on ref-mcp amplitudes and hodo
  TString commonDen = "amp_max[MiB2]>200 && fabs(time_max[MiB2])<150 && X[1]>-800 && Y[1]>-800";
  TString sdenSee    = commonDen + " && X[0]>=-9 && X[0]<=-1 && Y[0]>=-1 && Y[0]<=5 && !( amp_max[SEE]>25 && amp_max[SEE]<55 && (charge_sig[SEE]/charge_tot[SEE])<0.05)";
  TString sdenZS1    = commonDen + " && X[0]>=-9 && X[0]<=-1 && Y[0]>=-2 && Y[0]<=4 && !( amp_max[ZS1]>25 && amp_max[ZS1]<55 && (charge_sig[ZS1]/charge_tot[ZS1])<0.05)";
  TString sdenZS1_Abs    = commonDen + " && X[0]>=-9 && X[0]<=-2 && Y[0]>=-2 && Y[0]<=4 && !( amp_max[ZS1]>25 && amp_max[ZS1]<55 && (charge_sig[ZS1]/charge_tot[ZS1])<0.05)";
  TString sdenZS2    = commonDen + " && X[0]>=-9 && X[0]<=-1 && Y[0]>=-2 && Y[0]<=4 && !( amp_max[ZS2]>25 && amp_max[ZS2]<55 && (charge_sig[ZS2]/charge_tot[ZS2])<0.05)";
  TString sdenZS2_Abs    = commonDen + " && X[0]>=-9 && X[0]<=-1 && Y[0]>=-2 && Y[0]<=4 && !( amp_max[ZS2]>25 && amp_max[ZS2]<55 && (charge_sig[ZS2]/charge_tot[ZS2])<0.05)";
  TString sdenMiB3   = commonDen + " && X[0]>=-9 && X[0]<=-1 && Y[0]>=-2 && Y[0]<=4 && !( amp_max[MiB3]>25 && amp_max[MiB3]<55 && (charge_sig[MiB3]/charge_tot[MiB3])<0.05)";
  TString sdenMiB3_Abs   = commonDen + " && X[0]>=-9 && X[0]<=-1 && Y[0]>=-2 && Y[0]<=4 && !( amp_max[MiB3]>25 && amp_max[MiB3]<55 && (charge_sig[MiB3]/charge_tot[MiB3])<0.05)";
  TString sdenMgO   = "amp_max[MiB2]>100. && fabs(X)<50 && fabs(Y)<50 && run<3506 && (amp_max[PadA] > 0. || amp_max[PadB] > 0. || amp_max[PadC] > 0. || amp_max[PadD] > 0.)";
  TString sdenMgO_PadA   = "fabs(X)<50 && fabs(Y)<50 && run<3506 && X>-4 && X<0 && Y>-10 && Y<-6 && amp_max[PadA]>0."; 
  TString sdenMgO_PadB   = "fabs(X)<50 && fabs(Y)<50 && run<3506 && X>-7 && X<-3 && Y>12 && Y<16 && amp_max[PadB]>0."; 
  TString sdenMgO_PadC   = "fabs(X)<50 && fabs(Y)<50 && run<3506 && X>12 && X<16 && Y>11 && Y<15 && amp_max[PadC]>0."; 
  TString sdenMgO_PadD   = "fabs(X)<50 && fabs(Y)<50 && run<3506 && X>13 && X<17 && Y>-7 && Y<-3 && amp_max[PadD]>0."; 
  h4See ->Project("denSee", "HVSEE", sdenSee);
  h4binpMcp2->Project("denMcp2","HVZS2", sdenZS2);
  h4binpMcp5->Project("denMcp5","HVMiB3",sdenMiB3);
  h4binpMcp4->Project("denMcp4","HVZS1", sdenZS1);
  h4PMTMcp->Project("denPMTMcp","HVZS2", sdenZS2);
  h4MgO->Project("denMgO_PadA","HV2", sdenMgO_PadA);
  h4MgO->Project("denMgO_PadB","HV2", sdenMgO_PadB);
  h4MgO->Project("denMgO_PadC","HV2", sdenMgO_PadC);
  h4MgO->Project("denMgO_PadD","HV2", sdenMgO_PadD);

  // Nominal numerators
  TString snumSee  = sdenSee  + " && amp_max[SEE]>23";
  TString snumZS1  = sdenZS1  + " && amp_max[ZS1]>23";
  TString snumZS1_Abs  = sdenZS1_Abs  + " && amp_max[ZS1]>23";
  TString snumZS2  = sdenZS2  + " && amp_max[ZS2]>21.";
  TString snumZS2_Abs  = sdenZS2_Abs  + " && amp_max[ZS2]>21.";
  TString snumMiB3 = sdenMiB3 + " && amp_max[MiB3]>14.5";
  TString snumMiB3_Abs = sdenMiB3_Abs + " && amp_max[MiB3]>14.5";
  TString snumMgO = "amp_max[MiB2]>100. && fabs(X)<50 && fabs(Y)<50 && run<3506 && (amp_max[PadA] > 18.6 || amp_max[PadB] > 17.4 || amp_max[PadC] > 16.4 || amp_max[PadD] > 17.9)";
  TString snumMgO_PadA = "fabs(X)<50 && fabs(Y)<50 && run<3506 && X>-4 && X<0 && Y>-10 && Y<-6 && amp_max[PadA]>19.3";
  TString snumMgO_PadB = "fabs(X)<50 && fabs(Y)<50 && run<3506 && X>-7 && X<-3 && Y>12 && Y<16 && amp_max[PadB]>18.";
  TString snumMgO_PadC = "fabs(X)<50 && fabs(Y)<50 && run<3506 && X>12 && X<16 && Y>11 && Y<15 && amp_max[PadC]>18."; 
  TString snumMgO_PadD = "fabs(X)<50 && fabs(Y)<50 && run<3506 && X>13 && X<17 && Y>-7 && Y<-3 && amp_max[PadD]>18."; 
  h4See ->Project("numSee", "HVSEE", snumSee);
  h4binpMcp2->Project("numMcp2","HVZS2", snumZS2);
  h4binpMcp5->Project("numMcp5","HVMiB3",snumMiB3);
  h4binpMcp4->Project("numMcp4","HVZS1", snumZS1);
  h4PMTMcp->Project("numPMTMcp","HVZS2", snumZS2);
  h4MgO->Project("numMgO_PadA","HV2", snumMgO_PadA);
  h4MgO->Project("numMgO_PadB","HV2", snumMgO_PadB);
  h4MgO->Project("numMgO_PadC","HV2", snumMgO_PadC);
  h4MgO->Project("numMgO_PadD","HV2", snumMgO_PadD);
  
  cout << endl;
  cout << "-------------------------------------" << endl;
  cout << "check chiara: nominal " << endl;
  cout << "denSee = " << sdenSee << endl;
  cout << "numSee = " << snumSee  << endl;
  cout << endl;
  cout << "denZS1 = " << sdenZS1 << endl;
  cout << "numZS1 = " << snumZS1  << endl;
  cout << endl;
  cout << "denZS2 = " << sdenZS2 << endl;
  cout << "numZS2 = " << snumZS2  << endl;
  cout << endl;
  cout << "denMiB3 = " << sdenMiB3 << endl;
  cout << "numMiB3 = " << snumMiB3  << endl;
  cout << endl;
  cout << "denMgO = " << sdenMgO << endl;
  cout << "numMgO = " << snumMgO  << endl;
  cout << "-------------------------------------" << endl;
  cout << endl;
  cout << endl;

  // efficiencies
  TGraphAsymmErrors* effSee  = new TGraphAsymmErrors(numSee, denSee);
  TGraphAsymmErrors* effMcp2 = new TGraphAsymmErrors(numMcp2,denMcp2);
  TGraphAsymmErrors* effMcp5 = new TGraphAsymmErrors(numMcp5,denMcp5);
  TGraphAsymmErrors* effMcp4 = new TGraphAsymmErrors(numMcp4,denMcp4);
  TGraphAsymmErrors* effPMTMcp = new TGraphAsymmErrors(numPMTMcp,denPMTMcp);
  TGraphAsymmErrors* effMgO_PadA = new TGraphAsymmErrors(numMgO_PadA,denMgO_PadA);
  TGraphAsymmErrors* effMgO_PadB = new TGraphAsymmErrors(numMgO_PadB,denMgO_PadB);
  TGraphAsymmErrors* effMgO_PadC = new TGraphAsymmErrors(numMgO_PadC,denMgO_PadC);
  TGraphAsymmErrors* effMgO_PadD = new TGraphAsymmErrors(numMgO_PadD,denMgO_PadD);
  effSee->SetMarkerStyle(20);
  effSee->SetMarkerColor(6);
  effSee->SetLineColor(6);
  effMcp2->SetMarkerStyle(21);
  effMcp2->SetMarkerColor(2);
  effMcp2->SetLineColor(2);
  effMcp5->SetMarkerStyle(22);
  effMcp5->SetMarkerColor(3);
  effMcp5->SetLineColor(3);
  effMcp4->SetMarkerStyle(23);
  effMcp4->SetMarkerColor(4);
  effMcp4->SetLineColor(4);
  effPMTMcp->SetMarkerStyle(24);
  effPMTMcp->SetMarkerColor(5);
  effPMTMcp->SetLineColor(5);
  /*effMgO->SetMarkerStyle(25);
  effMgO->SetMarkerColor(6);
  effMgO->SetLineColor(6);*/

  // Rescaling for voltage partitor
  TGraphAsymmErrors* effMcp2Corr = new TGraphAsymmErrors(numMcp2,denMcp2);
  effMcp2Corr->SetMarkerStyle(20);
  effMcp2Corr->SetMarkerColor(2);
  effMcp2Corr->SetLineColor(2);
  int N_effMcp2Corr = effMcp2Corr->GetN();
  for (int ii = 0; ii<N_effMcp2Corr; ii++ ) {
    float thisValue = effMcp2Corr->GetX()[ii];
    float theCorr   = 10./11.;
    if (effMcp2Corr->GetX()[ii]>0) effMcp2Corr->GetX()[ii] *= theCorr;
  }

  TGraphAsymmErrors* effPMTMcpCorr = new TGraphAsymmErrors(numPMTMcp,denPMTMcp);
  effPMTMcpCorr->SetMarkerStyle(24);
  effPMTMcpCorr->SetMarkerColor(5);
  effPMTMcpCorr->SetLineColor(5);
  int N_effPMTMcpCorr = effPMTMcpCorr->GetN();
  for (int ii = 0; ii<N_effPMTMcpCorr; ii++ ) {
    float thisValue = effPMTMcpCorr->GetX()[ii];
    float theCorr   = 10./11.;
    if (effPMTMcpCorr->GetX()[ii]>0) effPMTMcpCorr->GetX()[ii] *= theCorr;
  }

  TGraphAsymmErrors* effMcp5Corr = new TGraphAsymmErrors(numMcp5,denMcp5);
  effMcp5Corr->SetMarkerStyle(20);
  effMcp5Corr->SetMarkerColor(3);
  effMcp5Corr->SetLineColor(3);
  int N_effMcp5Corr = effMcp5Corr->GetN();
  for (int ii = 0; ii<N_effMcp5Corr; ii++ ) {
    float thisValue = effMcp5Corr->GetX()[ii];
    float theCorr   = 10./11.;
    cout << "pre: " << effMcp5Corr->GetX()[ii] << endl;
    if (effMcp5Corr->GetX()[ii]>0) effMcp5Corr->GetX()[ii] *= theCorr;
    cout << "post: " << effMcp5Corr->GetX()[ii] << endl;
  }

  TLegend *leg;
  leg = new TLegend(0.1,0.5,0.4,0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(effSee,  "See",  "pl");
  leg->AddEntry(effMcp2, "ZS2",  "pl");
  leg->AddEntry(effPMTMcp, "PMT-ZS2",  "pl");
  leg->AddEntry(effMcp5, "MiB3", "pl");
  leg->AddEntry(effMcp4, "ZS1",  "pl");

  TCanvas* c1c = new TCanvas("c1c","c",1);
  c1c->SetGrid();
  TH2F* H2 = new TH2F("H2","",3500,0.,3500.,100,0.,1.);
  H2->GetXaxis()->SetTitle("Bias voltage [V]");
  H2->GetYaxis()->SetTitle("Efficiency");
  H2->Draw();
  effSee->Draw("Plsame"); 
  effMcp2->Draw("Plsame"); 
  effMcp5->Draw("Plsame"); 
  effMcp4->Draw("Plsame");
  effPMTMcp->Draw("Plsame");
  //effMgO->Draw("Plsame");
  leg->Draw();
  c1c->SaveAs("effH4.png");
  c1c->SaveAs("effH4.pdf");

  TCanvas* c1d = new TCanvas("c1d","c",1);
  c1d->SetGrid();
  H2->Draw();
  effMcp2->Draw("Plsame"); 
  effMcp5->Draw("Plsame"); 
  effMcp2Corr->Draw("Plsame"); 
  effMcp5Corr->Draw("Plsame"); 
  effPMTMcpCorr->Draw("Plsame"); 
  //effMgO->Draw("Plsame"); 
  c1d->SaveAs("effH4corrHV.png");


  // ----------------------------------------------------------------------------------

  // Efficiency vs absorber length
  float numZS2[5], numZS1[5], numMiB3[5];
  float denZS2[5], denZS1[5], denMiB3[5];
  
  // denominator
  for (int ii=0; ii<5; ii++) {
    TH1F* denAbsZS2 = new TH1F("denAbsZS2", "",6,-0.5,5.5);
    if (ii==0) h4Absorb->Project("denAbsZS2","0", sdenZS2_Abs+" && run==657 ");
    if (ii==1) h4Absorb  ->Project("denAbsZS2","1", sdenZS2_Abs+" && run==650");
    if (ii==2) h4Absorb  ->Project("denAbsZS2","2", sdenZS2_Abs+" && run==649");
    if (ii==3) h4Absorb  ->Project("denAbsZS2","3", sdenZS2_Abs+" && run==641");
    if (ii==4) h4Absorb  ->Project("denAbsZS2","5", sdenZS2_Abs+" && run==640");
    denZS2[ii] = denAbsZS2->Integral();
    delete denAbsZS2;

    TH1F* numAbsZS2 = new TH1F("numAbsZS2", "",6,-0.5,5.5);
    if (ii==0) h4Absorb->Project("numAbsZS2","0", snumZS2_Abs+" && run==657");
    if (ii==1) h4Absorb  ->Project("numAbsZS2","1", snumZS2_Abs+" && run==650");
    if (ii==2) h4Absorb  ->Project("numAbsZS2","2", snumZS2_Abs+" && run==649");
    if (ii==3) h4Absorb  ->Project("numAbsZS2","3", snumZS2_Abs+" && run==641");
    if (ii==4) h4Absorb  ->Project("numAbsZS2","5", snumZS2_Abs+" && run==640");
    numZS2[ii] = numAbsZS2->Integral();
    delete numAbsZS2;

    TH1F* denAbsZS1 = new TH1F("denAbsZS1", "",6,-0.5,5.5);
    if (ii==0) h4Absorb->Project("denAbsZS1","0", sdenZS1_Abs+" && run==657 ");
    if (ii==1) h4Absorb  ->Project("denAbsZS1","1", sdenZS1_Abs+" && run==650");
    if (ii==2) h4Absorb  ->Project("denAbsZS1","2", sdenZS1_Abs+" && run==649");
    if (ii==3) h4Absorb  ->Project("denAbsZS1","3", sdenZS1_Abs+" && run==641");
    if (ii==4) h4Absorb  ->Project("denAbsZS1","5", sdenZS1_Abs+" && run==640");
    denZS1[ii] = denAbsZS1->Integral();
    delete denAbsZS1;

    TH1F* numAbsZS1 = new TH1F("numAbsZS1", "",6,-0.5,5.5);
    if (ii==0) h4Absorb->Project("numAbsZS1","0", snumZS1_Abs+" && run==657");
    if (ii==1) h4Absorb  ->Project("numAbsZS1","1", snumZS1_Abs+" && run==650");
    if (ii==2) h4Absorb  ->Project("numAbsZS1","2", snumZS1_Abs+" && run==649");
    if (ii==3) h4Absorb  ->Project("numAbsZS1","3", snumZS1_Abs+" && run==641");
    if (ii==4) h4Absorb  ->Project("numAbsZS1","5", snumZS1_Abs+" && run==640");
    numZS1[ii] = numAbsZS1->Integral();
    delete numAbsZS1;
      
    TH1F* denAbsMiB3 = new TH1F("denAbsMiB3", "",6,-0.5,5.5);
    if (ii==0) h4Absorb->Project("denAbsMiB3","0", sdenMiB3_Abs+" && run==657 ");
    if (ii==1) h4Absorb  ->Project("denAbsMiB3","1", sdenMiB3_Abs+" && run==650");
    if (ii==2) h4Absorb  ->Project("denAbsMiB3","2", sdenMiB3_Abs+" && run==649");
    if (ii==3) h4Absorb  ->Project("denAbsMiB3","3", sdenMiB3_Abs+" && run==641");
    if (ii==4) h4Absorb  ->Project("denAbsMiB3","5", sdenMiB3_Abs+" && run==640");
    denMiB3[ii] = denAbsMiB3->Integral();
    delete denAbsMiB3;

    TH1F* numAbsMiB3 = new TH1F("numAbsMiB3", "",6,-0.5,5.5);
    if (ii==0) h4Absorb->Project("numAbsMiB3","0", snumMiB3_Abs+" && run==657");
    if (ii==1) h4Absorb  ->Project("numAbsMiB3","1", snumMiB3_Abs+" && run==650");
    if (ii==2) h4Absorb  ->Project("numAbsMiB3","2", snumMiB3_Abs+" && run==649");
    if (ii==3) h4Absorb  ->Project("numAbsMiB3","3", snumMiB3_Abs+" && run==641");
    if (ii==4) h4Absorb  ->Project("numAbsMiB3","5", snumMiB3_Abs+" && run==640");
    numMiB3[ii] = numAbsMiB3->Integral();
    delete numAbsMiB3;

  }

  // efficiencies
  float xaxis[5];
  float effAbsZS2[5], effAbsZS1[5], effAbsMiB3[5];
  for (int ii=0; ii<5; ii++) {
    effAbsZS2[ii]  = numZS2[ii]/denZS2[ii];
    effAbsZS1[ii]  = numZS2[ii]/denZS1[ii];
    effAbsMiB3[ii] = numMiB3[ii]/denMiB3[ii];
    xaxis[ii] = (float)ii;
    if(ii == 4) xaxis[ii] = 5.;
  }
  TGraph* effAbsZS2G  = new TGraph(5,xaxis,effAbsZS2);
  TGraph* effAbsZS1G = new TGraph(5,xaxis,effAbsZS1);
  TGraph* effAbsMiB3G = new TGraph(5,xaxis,effAbsMiB3);
  effAbsZS2G->SetMarkerStyle(20);
  effAbsZS2G->SetMarkerColor(1);
  effAbsZS2G->SetLineColor(1);
  effAbsZS1G->SetMarkerStyle(21);
  effAbsZS1G->SetMarkerColor(2);
  effAbsZS1G->SetLineColor(2);
  effAbsMiB3G->SetMarkerStyle(22);
  effAbsMiB3G->SetMarkerColor(3);
  effAbsMiB3G->SetLineColor(3);

  TLegend *leg2;
  leg2 = new TLegend(0.7,0.1,0.9,0.3);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(effAbsMiB3G, "40x2", "pl");
  leg2->AddEntry(effAbsZS2G,  "40x3",  "pl");
  leg2->AddEntry(effAbsZS1G,  "40x3 PMT",  "pl");

  TCanvas* c2c = new TCanvas("c2c","c",1);
  c2c->SetGrid();
  TH2F* H2b = new TH2F("H2b","",5,0.,5.1,100,0.5,1.);
  H2b->GetXaxis()->SetTitle("Number of X_{0}");
  H2b->GetYaxis()->SetTitle("Efficiency");
  H2b->Draw();
  effAbsZS2G->Draw("Plsame");
  effAbsZS1G->Draw("Plsame"); 
  effAbsMiB3G->Draw("Plsame"); 
  leg2->Draw();
  c2c->SaveAs("effH4vsAbs.png");
  c2c->SaveAs("effH4vsAbs.pdf");

  TFile myOutputFile("myOutputFileH4.root", "RECREATE");
  myOutputFile.cd();
  effSee ->Write("effSee");
  effMcp2->Write("effZS2");
  effMcp5->Write("effMib3");
  effMcp4->Write("effZS1");
  effMcp2Corr->Write("effZS2Corr");
  effMcp5Corr->Write("effMib3Corr");
  effPMTMcpCorr->Write("effZS2PMTCorr");
  effMgO_PadA->Write("effMgO_PadA");
  effMgO_PadB->Write("effMgO_PadB");
  effMgO_PadC->Write("effMgO_PadC");
  effMgO_PadD->Write("effMgO_PadD");
  effAbsZS2G->Write("effAbsZS2G");
  effAbsZS1G->Write("effAbsZS1G");
  effAbsMiB3G->Write("effAbsMiB3G");
  myOutputFile.Close();

}
	
